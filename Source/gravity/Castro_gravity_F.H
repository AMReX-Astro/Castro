#ifndef CASTRO_GRAVITY_F_H
#define CASTRO_GRAVITY_F_H
#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void pm_compute_delta_mass
    (const int* lo, const int* hi, amrex::Real* pm_delta_mass,
     const BL_FORT_FAB_ARG_3D(state),
     const BL_FORT_FAB_ARG_3D(stateout),
     const BL_FORT_FAB_ARG_3D(volume),
     const amrex::Real* problo, const amrex::Real* dx,
     const amrex::Real time, const amrex::Real dt);


  void pm_fix_solution
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(stateold),
     BL_FORT_FAB_ARG_3D(statenew),
     const amrex::Real* problo, const amrex::Real* dx,
     const amrex::Real time, const amrex::Real dt);


  void set_pointmass(const amrex::Real* pointmass);

#ifdef __cplusplus
}
#endif

#endif
