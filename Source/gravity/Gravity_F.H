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

  void ca_integrate_phi
    (const amrex::Real* avgden, amrex::Real* radial_grav,
     amrex::Real* radial_phi, 
     const amrex::Real* dr, const int* numpts_1d);

  void ca_put_radial_grav
    (const int* lo, const int* hi, 
     const amrex::Real* dx, amrex::Real dr,
     BL_FORT_FAB_ARG_3D(grav_vector),
     const amrex::Real* radial_grav, 
     const amrex::Real* problo, 
     int n1d, int level); 

  void ca_put_radial_phi
    (const int* lo, const int* hi, 
     const int* domlo, const int* domhi, 
     const amrex::Real* dx, amrex::Real dr,
     BL_FORT_FAB_ARG_3D(phi),
     const amrex::Real* radial_phi,
     const amrex::Real* problo, 
     int numpts_1d, int fill_interior);

  void ca_put_multipole_phi
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const amrex::Real* dx,
     BL_FORT_FAB_ARG_3D(phi),
     const int lnum,
     const amrex::Real* qL0, const amrex::Real* qLC, const amrex::Real* qLS,
     const amrex::Real* qU0, const amrex::Real* qUC, const amrex::Real* qUS,
     const int npts, const int boundary_only); 

  void ca_compute_multipole_moments
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const amrex::Real* dx,
     const BL_FORT_FAB_ARG_3D(rho),
     const BL_FORT_FAB_ARG_3D(vol),
     const int lnum,
     amrex::Real* qL0, amrex::Real* qLC, amrex::Real* qLS,
     amrex::Real* qU0, amrex::Real* qUC, amrex::Real* qUS,
     const int npts, const int boundary_only); 

  void pm_add_to_grav
    (const int* lo, const int* hi,
     const amrex::Real point_mass, BL_FORT_FAB_ARG_3D(phi), BL_FORT_FAB_ARG_3D(grav_vector),
     const amrex::Real* problo, const amrex::Real* dx);

#ifdef __cplusplus
}
#endif

#endif
