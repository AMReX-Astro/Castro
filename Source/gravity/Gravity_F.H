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

#ifdef GR_GRAV
  void ca_integrate_gr_grav
    (const amrex::Real* avgden, const amrex::Real* avgmass, 
     const amrex::Real* avgpres, amrex::Real* radial_grav, 
     const amrex::Real* dr, const int* numpts_1d);
#endif

#ifdef __cplusplus
}
#endif

#endif
