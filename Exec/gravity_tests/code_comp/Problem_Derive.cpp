#include <AMReX_REAL.H>

#include <Derive.H>
#include <Problem_Derive_F.H>
#include <Castro.H>
#include <Castro_F.H>
#include <prob_parameters.H>

using namespace amrex;

void ca_derextheating(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                      const FArrayBox& datfab, const Geometry& geomdata,
                      Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

  const auto dx = geomdata.CellSizeArray();
  const auto problo = geomdata.ProbLoArray();

  auto const dat = datfab.array();
  auto const der = derfab.array();

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {
    Real height = problo[AMREX_SPACEDIM-1] + 
      (static_cast<Real>(AMREX_D_PICK(i,j,k)) + 0.5_rt) * dx[AMREX_SPACEDIM-1] ;

    if (height < 1.125e0_rt * 4.e8_rt) {

      Real fheat = std::sin(8.e0_rt * M_PI * (height/ 4.e8_rt - 1.0_rt));
      der(i,j,k,0) = problem::heating_factor * fheat;

    } else {
      der(i,j,k,0) = 0.0_rt;
    }

  });
}


void ca_deradinvariant(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                       const FArrayBox& datfab, const Geometry& geomdata,
                       Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

  auto const dat = datfab.array();
  auto const der = derfab.array();

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    Real rhoInv = 1.0_rt / dat(i,j,k,URHO);

    eos_t eos_state;

    eos_state.e = dat(i,j,k,UEINT) * rhoInv;
    eos_state.T = dat(i,j,k,UTEMP);
    eos_state.rho = dat(i,j,k,URHO);
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n]  = dat(i,j,k,UFS+n) * rhoInv;
    }
#if NAUX > 0
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = dat(i,j,k,UFX+n) * rhoInv;
    }
#endif

    eos(eos_input_re, eos_state);

    der(i,j,k,0) = eos_state.p / std::pow(dat(i,j,k,URHO), 5.0e0_rt / 3.0e0_rt);

  });
}


void ca_derenthalpyfluxy(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                         const FArrayBox& datfab, const Geometry& geomdata,
                         Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

  auto const dat = datfab.array();
  auto const der = derfab.array();

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    Real rhoInv = 1.0_rt / dat(i,j,k,URHO);

    eos_t eos_state;

    eos_state.e = dat(i,j,k,UEINT) * rhoInv;
    eos_state.T = dat(i,j,k,UTEMP);
    eos_state.rho = dat(i,j,k,URHO);
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n]  = dat(i,j,k,UFS+n) * rhoInv;
    }
#if NAUX > 0
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = dat(i,j,k,UFX+n) * rhoInv;
    }
#endif
    eos(eos_input_re, eos_state);

    der(i,j,k,0) = eos_state.h * dat(i,j,k,UMY) / dat(i,j,k,URHO);
  });

}

void ca_derkinengfluxy(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                       const FArrayBox& datfab, const Geometry& geomdata,
                       Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

  auto const dat = datfab.array();
  auto const der = derfab.array();

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    der(i,j,k,0) = 0.5_rt / dat(i,j,k,0) * 
      (dat(i,j,k,1) * dat(i,j,k,1) +
       dat(i,j,k,2) * dat(i,j,k,2) +
       dat(i,j,k,3) * dat(i,j,k,3)) * dat(i,j,k,2);

  });

}
