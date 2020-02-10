#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

using namespace amrex;

AMREX_GPU_HOST_DEVICE void
Castro::add_hybrid_advection_source_c(const int i, const int j, const int k,
                                      const Real dt, GpuArray<Real, 3> const dx,
                                      GpuArray<Real, 3> const center,
                                      Array4<Real> const update,
                                      Array4<Real const> const qx,
                                      Array4<Real const> const qy,
                                      Array4<Real const> const qz) {

  auto loc = position_c(i,j,k, true, true, true);
#if 0
  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    loc[idir] -= center[idir];
  }

  Real R = std::sqrt(loc[0]*loc[0] + loc[1]*loc[1]);

  update(i,j,k,UMR) += - ( (loc[0] / R) * (qx(i+1,j,k,GDPRES) - qx(i,j,k,GDPRES)) / dx[0] +
                           (loc[1] / R) * (qy(i,j+1,k,GDPRES) - qy(i,j,k,GDPRES)) / dx[1] );

#endif
}
