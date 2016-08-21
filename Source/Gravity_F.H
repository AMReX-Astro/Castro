#ifndef _Gravity_F_H_
#define _Gravity_F_H_
#include <BLFort.H>

extern "C"
{

  void get_grav_const(Real* Gconst);

  void ca_test_residual
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(rhs),
     D_DECL(const BL_FORT_FAB_ARG(xgrad),
            const BL_FORT_FAB_ARG(ygrad),
            const BL_FORT_FAB_ARG(zgrad)), 
            const Real* dx, const Real* problo, 
            const int* coord_type);

  void ca_apply_metric
    (const int* lo, const int* hi,
     D_DECL(const int* xlo, const int *ylo, const int *zlo),
     D_DECL(const int* xhi, const int *yhi, const int *zhi),
     BL_FORT_FAB_ARG(cc),
     D_DECL(const BL_FORT_FAB_ARG(xedge),
            const BL_FORT_FAB_ARG(yedge),
            const BL_FORT_FAB_ARG(zedge)), 
            const Real* dx, const int* coord_type);

  void ca_weight_cc
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(cc),
     const Real* dx, const int* coord_type);

  void ca_unweight_cc
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(cc),
     const Real* dx, const int* coord_type);

  void ca_unweight_edges
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(edge),
     const Real* dx, const int* coord_type, const int* idir);

  void ca_compute_radial_mass
    (const int lo[], const int hi[], 
     const Real* dx, const Real* dr,
     const BL_FORT_FAB_ARG(S), 
     const Real* avgmass, const Real* avgvol, 
     const Real* problo, const int* numpts_1d, 
     const int* drdxfac, const int* level); 

  void ca_compute_avgpres
    (const int lo[], const int hi[], 
     const Real* dx, const Real* dr,
     const BL_FORT_FAB_ARG(S), 
     const Real* avgpres, 
     const Real* problo, const int* numpts_1d, 
     const int* drdxfac, const int* level); 

  void ca_integrate_grav
    (const Real* avgmass, const Real* avgden,
     const Real* radial_grav,
     const Real* max_radius,
     const Real* dr, const int* numpts_1d);

  void ca_prescribe_phi
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(p),
     const Real* dx);

  void ca_prescribe_grav
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(g),
     const Real* dx);

#ifdef GR_GRAV
  void ca_integrate_gr_grav
    (const Real* avgden, const Real* avgmass, 
     const Real* avgpres, const Real* radial_grav, 
     const Real* dr, const int* numpts_1d);
#endif

  void ca_integrate_phi
    (const Real* avgden, const Real* radial_grav,
     const Real* radial_phi, 
     const Real* dr, const int* numpts_1d);

  void ca_put_radial_grav
    (const int lo[], const int hi[], 
     const Real* dx, const Real* dr,
     const BL_FORT_FAB_ARG(grav_vector),
     const Real* radial_grav, 
     const Real* problo, 
     const int* n1d, const int* level); 

  void ca_put_radial_phi
    (const int lo[], const int hi[], 
     const int domlo[], const int domhi[], 
     const Real* dx, const Real* dr,
     const BL_FORT_FAB_ARG(phi),
     const Real* radial_phi,
     const Real* problo, 
     const int* numpts_1d, const int* fill_interior);

  void init_multipole_gravity
    (const int* lnum, const int* lo_bc, const int* hi_bc);

  void ca_put_multipole_phi
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const Real* dx,
     const BL_FORT_FAB_ARG_3D(phi),
     const int* lnum,
     Real* qL0, Real* qLC, Real* qLS,
     Real* qU0, Real* qUC, Real* qUS,
     const int* npts, const int* boundary_only); 

  void ca_compute_multipole_moments
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const Real* dx,
     const BL_FORT_FAB_ARG_3D(rho),
     const BL_FORT_FAB_ARG_3D(vol),
     const int* lnum,
     Real* qL0, Real* qLC, Real* qLS,
     Real* qU0, Real* qUC, Real* qUS,
     const int* npts, const int* boundary_only); 

  void ca_compute_direct_sum_bc
    (const int lo[], const int hi[],
     const int domlo[], const int domhi[],
     const int* symmetry_type, const int lo_bc[], const int hi_bc[],
     const Real* dx, const BL_FORT_FAB_ARG(rho),
     const Real* problo, const Real* probhi,
     Real* bcXYLo, Real* bcXYHi,
     Real* bcXZLo, Real* bcXZHi,
     Real* bcYZLo, Real* bcYZHi);

  void ca_put_direct_sum_bc
    (const int lo[], const int hi[], 
     const int domlo[], const int domhi[],
     const BL_FORT_FAB_ARG(rho),
     Real* bcXYLo, Real* bcXYHi,
     Real* bcXZLo, Real* bcXZHi,
     Real* bcYZLo, Real* bcYZHi);

#ifdef POINTMASS
  void pm_add_to_grav
    (const Real* point_mass, BL_FORT_FAB_ARG_3D(phi), BL_FORT_FAB_ARG_3D(grav_vector),
     const Real* problo, const Real* dx, const int* lo, const int* hi);
#endif

}

#endif
