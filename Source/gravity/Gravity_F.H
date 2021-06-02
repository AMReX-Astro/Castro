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
#ifdef GR_GRAV
     amrex::Real* radial_pres,
#endif
     const amrex::Real* max_radius,
     const amrex::Real* dr, const int* numpts_1d);

#ifdef __cplusplus
}
#endif

#endif
