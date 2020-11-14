#ifndef GRAVITY_F_H
#define GRAVITY_F_H
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void ca_integrate_grav
    (const amrex::Real* avgmass, const amrex::Real* avgden,
     amrex::Real* radial_grav,
     const amrex::Real* max_radius,
     const amrex::Real* dr, const int* numpts_1d);

  void ca_prescribe_phi
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(p),
     const amrex::Real* dx);

  void ca_prescribe_grav
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(g),
     const amrex::Real* dx);

#ifdef GR_GRAV
  void ca_integrate_gr_grav
    (const amrex::Real* avgden, const amrex::Real* avgmass, 
     const amrex::Real* avgpres, amrex::Real* radial_grav, 
     const amrex::Real* dr, const int* numpts_1d);

  void ca_compute_avgpres
    (const int* lo, const int* hi,
     const amrex::Real* dx, amrex::Real dr,
     const BL_FORT_FAB_ARG_3D(S), 
     amrex::Real* avgpres, 
     const amrex::Real* problo, int numpts_1d, 
     int drdxfac, int level); 
#endif

#ifdef __cplusplus
}
#endif

#endif
