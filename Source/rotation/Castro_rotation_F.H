#ifndef CASTRO_ROTATION_F_H
#define CASTRO_ROTATION_F_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void set_rot_period(const amrex::Real* period);

  void get_omega(amrex::Real* omega);

  void ca_fill_rotational_potential
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(phi),
     const amrex::Real* dx, const amrex::Real time);

  void ca_fill_rotational_acceleration
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(rot),
     const BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, const amrex::Real time);

  void ca_fill_rotational_psi
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(psi),
     const amrex::Real* dx, const amrex::Real time);

  void rotational_acceleration
     (const Real* r, const Real* v,
      const int coriolis, Real* Sr);

#ifdef __cplusplus
}
#endif

#endif
