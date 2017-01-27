#ifndef _Gravity_F_H_
#define _Gravity_F_H_
#include <AMReX_BLFort.H>

extern "C"
{

  void get_grav_const(amrex::Real* Gconst);

  void ca_test_residual
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(rhs),
     D_DECL(const BL_FORT_FAB_ARG(xgrad),
            const BL_FORT_FAB_ARG(ygrad),
            const BL_FORT_FAB_ARG(zgrad)), 
            const amrex::Real* dx, const amrex::Real* problo, 
            const int* coord_type);

  void ca_apply_metric
    (const int* lo, const int* hi,
     D_DECL(const int* xlo, const int *ylo, const int *zlo),
     D_DECL(const int* xhi, const int *yhi, const int *zhi),
     BL_FORT_FAB_ARG(cc),
     D_DECL(const BL_FORT_FAB_ARG(xedge),
            const BL_FORT_FAB_ARG(yedge),
            const BL_FORT_FAB_ARG(zedge)), 
            const amrex::Real* dx, const int* coord_type);

  void ca_weight_cc
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(cc),
     const amrex::Real* dx, const int* coord_type);

  void ca_unweight_cc
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(cc),
     const amrex::Real* dx, const int* coord_type);

  void ca_unweight_edges
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(edge),
     const amrex::Real* dx, const int* coord_type, const int* idir);

  void ca_compute_radial_mass
    (const int lo[], const int hi[], 
     const amrex::Real* dx, const amrex::Real* dr,
     const BL_FORT_FAB_ARG(S), 
     const amrex::Real* avgmass, const amrex::Real* avgvol, 
     const amrex::Real* problo, const int* numpts_1d, 
     const int* drdxfac, const int* level); 

  void ca_compute_avgpres
    (const int lo[], const int hi[], 
     const amrex::Real* dx, const amrex::Real* dr,
     const BL_FORT_FAB_ARG(S), 
     const amrex::Real* avgpres, 
     const amrex::Real* problo, const int* numpts_1d, 
     const int* drdxfac, const int* level); 

  void ca_integrate_grav
    (const amrex::Real* avgmass, const amrex::Real* avgden,
     const amrex::Real* radial_grav,
     const amrex::Real* max_radius,
     const amrex::Real* dr, const int* numpts_1d);

  void ca_prescribe_phi
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(p),
     const amrex::Real* dx);

  void ca_prescribe_grav
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(g),
     const amrex::Real* dx);

#ifdef GR_GRAV
  void ca_integrate_gr_grav
    (const amrex::Real* avgden, const amrex::Real* avgmass, 
     const amrex::Real* avgpres, const amrex::Real* radial_grav, 
     const amrex::Real* dr, const int* numpts_1d);
#endif

  void ca_integrate_phi
    (const amrex::Real* avgden, const amrex::Real* radial_grav,
     const amrex::Real* radial_phi, 
     const amrex::Real* dr, const int* numpts_1d);

  void ca_put_radial_grav
    (const int lo[], const int hi[], 
     const amrex::Real* dx, const amrex::Real* dr,
     const BL_FORT_FAB_ARG(grav_vector),
     const amrex::Real* radial_grav, 
     const amrex::Real* problo, 
     const int* n1d, const int* level); 

  void ca_put_radial_phi
    (const int lo[], const int hi[], 
     const int domlo[], const int domhi[], 
     const amrex::Real* dx, const amrex::Real* dr,
     const BL_FORT_FAB_ARG(phi),
     const amrex::Real* radial_phi,
     const amrex::Real* problo, 
     const int* numpts_1d, const int* fill_interior);

  void init_multipole_gravity
    (const int* lnum, const int* lo_bc, const int* hi_bc);

  void ca_put_multipole_phi
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const amrex::Real* dx,
     const BL_FORT_FAB_ARG_3D(phi),
     const int* lnum,
     amrex::Real* qL0, amrex::Real* qLC, amrex::Real* qLS,
     amrex::Real* qU0, amrex::Real* qUC, amrex::Real* qUS,
     const int* npts, const int* boundary_only); 

  void ca_compute_multipole_moments
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const amrex::Real* dx,
     const BL_FORT_FAB_ARG_3D(rho),
     const BL_FORT_FAB_ARG_3D(vol),
     const int* lnum,
     amrex::Real* qL0, amrex::Real* qLC, amrex::Real* qLS,
     amrex::Real* qU0, amrex::Real* qUC, amrex::Real* qUS,
     const int* npts, const int* boundary_only); 

  void ca_compute_direct_sum_bc
    (const int* lo, const int* hi, const amrex::Real* dx,
     const int* symmetry_type, const int* lo_bc, const int* hi_bc,
     const amrex::Real* rho, const int* r_lo, const int* r_hi,
     const amrex::Real* vol, const int* v_lo, const int* v_hi,
     const amrex::Real* problo, const amrex::Real* probhi,
     amrex::Real* bcXYLo, amrex::Real* bcXYHi,
     amrex::Real* bcXZLo, amrex::Real* bcXZHi,
     amrex::Real* bcYZLo, amrex::Real* bcYZHi,
     const int* bclo, const int* bchi, const amrex::Real* bcdx);

  void ca_put_direct_sum_bc
    (const int* lo, const int* hi, 
     amrex::Real* phi, const int* p_lo, const int* p_hi,
     const amrex::Real* bcXYLo, const amrex::Real* bcXYHi,
     const amrex::Real* bcXZLo, const amrex::Real* bcXZHi,
     const amrex::Real* bcYZLo, const amrex::Real* bcYZHi,
     const int* bclo, const int* bchi);

#ifdef POINTMASS
  void pm_add_to_grav
    (const amrex::Real* point_mass, BL_FORT_FAB_ARG_3D(phi), BL_FORT_FAB_ARG_3D(grav_vector),
     const amrex::Real* problo, const amrex::Real* dx, const int* lo, const int* hi);
#endif

}

#endif
