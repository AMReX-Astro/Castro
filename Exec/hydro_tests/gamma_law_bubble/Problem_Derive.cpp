#include <AMReX_REAL.H>

#include <Derive.H>
#include <Castro.H>
#include <prob_util.H>

using namespace amrex;

using RealVector = amrex::Gpu::ManagedVector<amrex::Real>;

void ca_derpi(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
              const FArrayBox& datfab, const Geometry& geomdata,
              Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

  // derive the dynamic pressure

  const auto dx = geomdata.CellSizeArray();
  const auto problo = geomdata.ProbLoArray();

  auto const dat = datfab.array();
  auto const der = derfab.array();

  // first make a 1D initial model for the entire domain
  int npts_1d = (2.0_rt * problem::center[AMREX_SPACEDIM-1] + 1.e-8_rt) /
    dx[AMREX_SPACEDIM-1];

  RealVector pressure_rv(npts_1d, 0);
  RealVector density_rv(npts_1d, 0);
  RealVector temp_rv(npts_1d, 0);
  RealVector eint_rv(npts_1d, 0);

  Real* const pressure = pressure_rv.dataPtr();
  Real* const density = density_rv.dataPtr();
  Real* const temp = temp_rv.dataPtr();
  Real* const eint = eint_rv.dataPtr();

  gamma_law_initial_model(pressure, density, temp, eint, npts_1d, dx);

  // Compute pressure from the EOS

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
  {

    Real rhoInv = 1.0_rt / dat(i,j,k,URHO);
    Real T = dat(i,j,k,UTEMP);

    eos_t eos_state;

    eos_state.rho = dat(i,j,k,URHO);
    eos_state.T = dat(i,j,k,UTEMP);
    eos_state.e = dat(i,j,k,UEINT) * rhoInv;
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = dat(i,j,k,UFS+n) * rhoInv;
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = dat(i,j,k,UFX+n) * rhoInv;
    }
#endif

    if (eos_state.e <= 0.0_rt) {
      eos(eos_input_rt, eos_state);
      der(i,j,k,0) = eos_state.p;

    } else {
      eos(eos_input_re, eos_state);
      der(i,j,k,0) = eos_state.p;
    }

    der(i,j,k,0) -= pressure[AMREX_D_PICK(i,j,k)];

  });

}


void ca_derpioverp0(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                    const FArrayBox& datfab, const Geometry& geomdata,
                    Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

  const auto dx = geomdata.CellSizeArray();
  const auto problo = geomdata.ProbLoArray();

  auto const dat = datfab.array();
  auto const der = derfab.array();

  // first make a 1D initial model for the entire domain
  int npts_1d = (2.0_rt * problem::center[AMREX_SPACEDIM-1] + 1.e-8_rt) /
    dx[AMREX_SPACEDIM-1];

  RealVector pressure_rv(npts_1d, 0);
  RealVector density_rv(npts_1d, 0);
  RealVector temp_rv(npts_1d, 0);
  RealVector eint_rv(npts_1d, 0);

  Real* const pressure = pressure_rv.dataPtr();
  Real* const density = density_rv.dataPtr();
  Real* const temp = temp_rv.dataPtr();
  Real* const eint = eint_rv.dataPtr();

  gamma_law_initial_model(pressure, density, temp, eint, npts_1d, dx);

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
  {

    Real rhoInv = 1.0_rt/dat(i,j,k,URHO);

    Real e = dat(i,j,k,UEINT) * rhoInv;
    Real T = dat(i,j,k,UTEMP);

    eos_t eos_state;

    eos_state.rho = dat(i,j,k,URHO);
    eos_state.T = T;
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = dat(i,j,k,UFS+n) * rhoInv;
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = dat(i,j,k,UFX+n) * rhoInv;
    }
#endif
    eos_state.e = e;

    if (e <= 0.0_rt) {
      eos(eos_input_rt, eos_state);
      der(i,j,k,0) = eos_state.p;

    } else {
      eos(eos_input_re, eos_state);
      der(i,j,k,0) = eos_state.p;

    }

    der(i,j,k,0) = (der(i,j,k,0) - pressure[AMREX_D_PICK(i,j,k)]) / pressure[AMREX_D_PICK(i,j,k)];

  });

}


void ca_derrhopert(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                    const FArrayBox& datfab, const Geometry& geomdata,
                    Real /*time*/, const int* /*bcrec*/, int /*level*/)

{

  const auto dx = geomdata.CellSizeArray();
  const auto problo = geomdata.ProbLoArray();

  auto const dat = datfab.array();
  auto const der = derfab.array();

  // first make a 1D initial model for the entire domain
  int npts_1d = (2.0_rt * problem::center[AMREX_SPACEDIM-1] + 1.e-8_rt) /
    dx[AMREX_SPACEDIM-1];

  RealVector density_rv(npts_1d, 0);
  Real* const density = density_rv.dataPtr();

  density[0]  = problem::dens_base;

  // only initialize the first species
  Real xn[NumSpec];
  xn[0] = 1.0_rt;

  // compute the pressure scale height (for an isothermal, ideal-gas
  // atmosphere)
  Real H = problem::pres_base / problem::dens_base / std::abs(gravity::const_grav);

  for (int j = 0; j < npts_1d; j++) {

     if (problem::do_isentropic) {
       Real z = static_cast<Real>(j) * dx[AMREX_SPACEDIM-1];
       density[j] = problem::dens_base *
         std::pow((gravity::const_grav * problem::dens_base * (eos_gamma - 1.0_rt) * z/
                   (eos_gamma * problem::pres_base) + 1.0_rt), 1.0_rt/(eos_gamma - 1.0_rt));
     } else {
       Real z = (static_cast<Real>(j) + 0.5_rt) * dx[AMREX_SPACEDIM-1];
       density[j] = problem::dens_base * std::exp(-z/H);
     }

  }

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
  {
    der(i,j,k,0) = dat(i,j,k,URHO) - density[AMREX_D_PICK(i,j,k)];
  });

}


void ca_dertpert(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                 const FArrayBox& datfab, const Geometry& geomdata,
                 Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

  const auto dx = geomdata.CellSizeArray();
  const auto problo = geomdata.ProbLoArray();

  auto const dat = datfab.array();
  auto const der = derfab.array();

  // first make a 1D initial model for the entire domain
  int npts_1d = (2.0_rt * problem::center[AMREX_SPACEDIM-1] + 1.e-8_rt) /
    dx[AMREX_SPACEDIM-1];

  RealVector pressure_rv(npts_1d, 0);
  RealVector density_rv(npts_1d, 0);
  RealVector temp_rv(npts_1d, 0);
  RealVector eint_rv(npts_1d, 0);

  Real* const pressure = pressure_rv.dataPtr();
  Real* const density = density_rv.dataPtr();
  Real* const temp = temp_rv.dataPtr();
  Real* const eint = eint_rv.dataPtr();

  gamma_law_initial_model(pressure, density, temp, eint, npts_1d, dx);

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
  {
    der(i,j,k,0) = dat(i,j,k,UTEMP) - temp[AMREX_D_PICK(i,j,k)];
  });

}
