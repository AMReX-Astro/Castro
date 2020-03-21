#ifndef _Castro_rotation_F_H_
#define _Castro_rotation_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void set_rot_period(const amrex::Real* period);

  void get_omega(const amrex::Real time, amrex::Real* omega);

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

  void ca_rsrc
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const BL_FORT_FAB_ARG_3D(phi),
     const BL_FORT_FAB_ARG_3D(rot),
     const BL_FORT_FAB_ARG_3D(s_old),
     BL_FORT_FAB_ARG_3D(source),
     const BL_FORT_FAB_ARG_3D(vol),
     const amrex::Real* dx, const amrex::Real dt, const amrex::Real time);

  void ca_corrrsrc
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const BL_FORT_FAB_ARG_3D(phi_old),
     const BL_FORT_FAB_ARG_3D(phi_new),
     const BL_FORT_FAB_ARG_3D(rot_old),
     const BL_FORT_FAB_ARG_3D(rot_new),
     const BL_FORT_FAB_ARG_3D(S_old),
     const BL_FORT_FAB_ARG_3D(S_new),
     BL_FORT_FAB_ARG_3D(source),
     const BL_FORT_FAB_ARG_3D(xflux),
     const BL_FORT_FAB_ARG_3D(yflux),
     const BL_FORT_FAB_ARG_3D(zflux),
     const amrex::Real* dx, const amrex::Real dt, const amrex::Real time,
     const BL_FORT_FAB_ARG_3D(volume));

#ifdef __cplusplus
}
#endif

#endif
