#include "AMReX_REAL.H"

#include "Derive.H"
#include "Castro.H"
#include "Castro_F.H"
#include "prob_parameters.H"

using namespace amrex;


void ca_derpi(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
              const FArrayBox& datfab, const Geometry& geomdata,
              Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

  // derive the dynamic pressure

  const auto dx = geomdata.CellSizeArray();
  const auto problo = geomdata.ProbLoArray();

  GpuArray<Real, 3> center;
  ca_get_center(center.begin());

  auto const dat = datfab.array();
  auto const der = derfab.array();

  // first make a 1D initial model for the entire domain
  int npts_1d = (2.0_rt*center[AMREX_SPACEDIM] + 1.e-8_rt) /
    dx[AMREX_SPACEDIM];

  Array1D<Real, 0, npts_1d-1> pressure;
  Array1D<Real, 0, npts_1d-1> density;
  Array1D<Real, 0, npts_1d-1> temp;
  Array1D<Real, 0, npts_1d-1> eint;

  pressure[0] = pres_base;
  density[0]  = dens_base;

  // only initialize the first species
  Real xn[NumSpec];
  xn[0] = 1.0_rt;

  // compute the pressure scale height (for an isothermal, ideal-gas
  // atmosphere)
  Real H = pres_base / dens_base / std::abs(Gravity::const_grav);

  for (int j = 0; j < npts_1d; j++) {

    // initial guess
     temp[j] = T_guess;

     if (do_isentropic) {
       Real z = static_cast<Real>(j) * dx[AMREX_SPACEDIM];
       density[j] = dens_base * 
         (Gravity::const_grav * dens_base * (gamma_const - 1.0_rt) * z/
          std::pow((gamma_const*pres_base) + 1.0_rt), 1.0_rt/(gamma_const - 1.0_rt));
     } else {
       Real z = (static_cast<Real>(j)+HALF) * dx[AMREX_SPACEDIM];
       density[j] = dens_base * std::exp(-z/H);
     }

     if (j > 0) {
        pressure[j] = pressure[j-1] -
          dx[AMREX_SPACEDIM] * 0.5_rt *
          (density[j] + density[j-1]) * std::abs(Gravity::const_grav);
     }

     eos_t eos_state;

     eos_state.rho = density[j];
     eos_state.T = temp[j];
     for (int n = 0; n < NumSpec; n++) {
       eos_state.xn[n] = xn[n];
     }
     eos_state.p = pressure[j];

     eos(eos_input_rp, eos_state);

     eint[j] = eos_state.e;
     temp[j] = eos_state.T;

  }

  // Compute pressure from the EOS

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    Real rhoInv = 1.0_rt / dat[i,j,k,URHO];
    Real T = dat(i,j,k,UTEMP);

    eos_t eos_state;

    eos_state.rho = dat(i,j,k,URHO);
    eos_state.T = dat(i,j,k,UTEMP);
    eos_state.e = dat(i,j,k,UEINT) * rhoInv;
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = dat(i,j,k,UFS+n) / dat(i,j,k,URHO);
    }
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = dat(i,j,k,UFX+n) / dat(i,j,k,URHO);
    }

    if (eos_state.e <= 0.0_rt) {
      eos(eos_input_rt, eos_state);
      der(i,j,k,0) = eos_state.p;

    } else {
      eos(eos_input_re, eos_state);
      der(i,j,k,0) = eos_state.p;
    }

#if AMREX_SPACEDIM == 1
    p(i,j,k,0) -= pressure(i);
#elif AMREX_SPACEDIM == 2
    p(i,j,k,0) -= pressure(j);
#else
    p(i,j,k,0) -= pressure(k);
#endif

  });

}


void ca_derpioverp0(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                    const FArrayBox& datfab, const Geometry& geomdata,
                    Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

  const auto dx = geomdata.CellSizeArray();
  const auto problo = geomdata.ProbLoArray();

  GpuArray<Real, 3> center;
  ca_get_center(center.begin());

  auto const dat = datfab.array();
  auto const der = derfab.array();

  // first make a 1D initial model for the entire domain
  int npts_1d = (2.0_rt*center[AMREX_SPACEDIM] + 1.e-8_rt) /
    dx[AMREX_SPACEDIM];

  Array1D<Real, 0, npts_1d-1> pressure;
  Array1D<Real, 0, npts_1d-1> density;
  Array1D<Real, 0, npts_1d-1> temp;
  Array1D<Real, 0, npts_1d-1> eint;

  pressure[0] = pres_base;
  density[0]  = dens_base;

  // only initialize the first species
  Real xn[NumSpec];
  xn[0] = 1.0_rt;

  // compute the pressure scale height (for an isothermal, ideal-gas
  // atmosphere)
  Real H = pres_base / dens_base / std::abs(const_grav);

  for (int j = 0; j < npts_1d; j++) {

    // initial guess
     temp[j] = T_guess;

     if (do_isentropic) {
       Real z = static_cast<Real>(j) * dx[AMREX_SPACEDIM];
       density[j] = dens_base * 
         (const_grav * dens_base * (gamma_const - 1.0_rt) * z/
          std::pow((gamma_const*pres_base) + 1.0_rt), 1.0_rt/(gamma_const - 1.0_rt));
     } else {
       Real z = (static_cast<Real>(j)+HALF) * dx[AMREX_SPACEDIM];
       density[j] = dens_base * std::exp(-z/H);
     }

     if (j > 0) {
        pressure[j] = pressure[j-1] -
          dx[AMREX_SPACEDIM] * 0.5_rt *
          (density[j] + density[j-1]) * std::abs(const_grav);
     }

     eos_t eos_state;

     eos_state.rho = density[j];
     eos_state.T = temp[j];
     for (int n = 0; n < NumSpec; n++) {
       eos_state.xn[n] = xn[n];
     }
     eos_state.p = pressure[j];

     eos(eos_input_rp, eos_state);

     eint[j] = eos_state.e;
     temp[j] = eos_state.T;

  }


  // Compute pressure from the EOS

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    Real rhoInv = 1.0_rt/dat[i,j,k,URHO];

    Real e = dat(i,j,k,UEINT) * rhoInv;
    Real T = dat(i,j,k,UTEMP);

    eos_t eos_state;

    eos_state.rho = dat(i,j,k,URHO);
    eos_state.T = T;
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = dat(i,j,k,UFS+n) * rhoInv;
    }
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = dat(i,j,k,UFX+n) * rhoInv;
    }
    eos_state.e = e;

    if (e <= 0.0_rt) {
      eos(eos_input_rt, eos_state);
      der(i,j,k,0) = eos_state.p;

    } else {
      eos(eos_input_re, eos_state);
      der(i,j,k,0) = eos_state.p;

    }

#if AMREX_SPACEDIM == 1
    der(i,j,k,0) = (der(i,j,k,0) - pressure(i)) / pressure(i);
#elif AMREX_SPACEDIM == 2
    der(i,j,k,0) = (der(i,j,k,0) - pressure(j)) / pressure(j);
#else
    der(i,j,k,0) = (der(i,j,k,0) - pressure(k)) / pressure(k);
#endif

  });

}


void ca_derrhopert(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                    const FArrayBox& datfab, const Geometry& geomdata,
                    Real /*time*/, const int* /*bcrec*/, int /*level*/)

{

  const auto dx = geomdata.CellSizeArray();
  const auto problo = geomdata.ProbLoArray();

  GpuArray<Real, 3> center;
  ca_get_center(center.begin());

  auto const dat = datfab.array();
  auto const der = derfab.array();

  // first make a 1D initial model for the entire domain
  int npts_1d = (2.0_rt*center[AMREX_SPACEDIM] + 1.e-8_rt) /
    dx[AMREX_SPACEDIM];

  Array1D<Real, 0, npts_1d-1> density;

  density[0]  = dens_base;

  // only initialize the first species
  Real xn[NumSpec];
  xn[0] = 1.0_rt;

  // compute the pressure scale height (for an isothermal, ideal-gas
  // atmosphere)
  Real H = pres_base / dens_base / std::abs(const_grav);

  for (int j = 0; j < npts_1d; j++) {

     if (do_isentropic) {
       Real z = static_cast<Real>(j) * dx[AMREX_SPACEDIM];
       density[j] = dens_base * 
         (const_grav * dens_base * (gamma_const - 1.0_rt) * z/
          std::pow((gamma_const*pres_base) + 1.0_rt), 1.0_rt/(gamma_const - 1.0_rt));
     } else {
       Real z = (static_cast<Real>(j)+HALF) * dx[AMREX_SPACEDIM];
       density[j] = dens_base * std::exp(-z/H);
     }

  }

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

#if AMREX_SPACEDIM == 1
    der(i,j,k,0) = dat(i,j,k,URHO) - density(i);
#elif AMREX_SPACEDIM == 2
    der(i,j,k,0) = dat(i,j,k,URHO) - density(j);
#else
    der(i,j,k,0) = dat(i,j,k,URHO) - density(k);
#endif
  });

}


void ca_dertpert(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                 const FArrayBox& datfab, const Geometry& geomdata,
                 Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

  const auto dx = geomdata.CellSizeArray();
  const auto problo = geomdata.ProbLoArray();

  GpuArray<Real, 3> center;
  ca_get_center(center.begin());

  auto const dat = datfab.array();
  auto const der = derfab.array();

  // first make a 1D initial model for the entire domain
  int npts_1d = (2.0_rt*center[AMREX_SPACEDIM] + 1.e-8_rt) /
    dx[AMREX_SPACEDIM];

  Array1D<Real, 0, npts_1d-1> pressure;
  Array1D<Real, 0, npts_1d-1> density;
  Array1D<Real, 0, npts_1d-1> temp;
  Array1D<Real, 0, npts_1d-1> eint;

  pressure[0] = pres_base;
  density[0]  = dens_base;

  // only initialize the first species
  Real xn[NumSpec];
  xn[0] = 1.0_rt;

  // compute the pressure scale height (for an isothermal, ideal-gas
  // atmosphere)
  Real H = pres_base / dens_base / std::abs(const_grav);

  for (int j = 0; j < npts_1d; j++) {

    // initial guess
     temp[j] = T_guess;

     if (do_isentropic) {
       Real z = static_cast<Real>(j) * dx[AMREX_SPACEDIM];
       density[j] = dens_base * 
         (const_grav * dens_base * (gamma_const - 1.0_rt) * z/
          std::pow((gamma_const*pres_base) + 1.0_rt), 1.0_rt/(gamma_const - 1.0_rt));
     } else {
       Real z = (static_cast<Real>(j)+HALF) * dx[AMREX_SPACEDIM];
       density[j] = dens_base * std::exp(-z/H);
     }

     if (j > 0) {
        pressure[j] = pressure[j-1] -
          dx[AMREX_SPACEDIM] * 0.5_rt *
          (density[j] + density[j-1]) * std::abs(const_grav);
     }

     eos_t eos_state;

     eos_state.rho = density[j];
     eos_state.T = temp[j];
     for (int n = 0; n < NumSpec; n++) {
       eos_state.xn[n] = xn[n];
     }
     eos_state.p = pressure[j];

     eos(eos_input_rp, eos_state);

     eint[j] = eos_state.e;
     temp[j] = eos_state.T;

  }


  // Compute pressure from the EOS

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {
#if AMREX_SPACEDIM == 1
    der(i,j,k,0) = dat(i,j,k,UTEMP) - temp(i);
#elif AMREX_SPACEDIM == 2
    der(i,j,k,0) = dat(i,j,k,UTEMP) - temp(j);
#else
    der(i,j,k,0) = dat(i,j,k,UTEMP) - temp(k);
#endif
  });

}
