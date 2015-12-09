#ifndef _Gravity_F_H_
#define _Gravity_F_H_
#include <BLFort.H>

BL_FORT_PROC_DECL(GET_GRAV_CONST,get_grav_const)(Real* Gconst);

BL_FORT_PROC_DECL(CA_EDGE_INTERP,ca_edge_interp)
    (const int* flo, const int* fhi,
     const int* nc, const int* refRatio, const int* dir,
     BL_FORT_FAB_ARG(fine));
    
BL_FORT_PROC_DECL(CA_PC_EDGE_INTERP,ca_pc_edge_interp)
    (const int* lo, const int* hi,
     const int* nc, const int* refRatio, const int* dir,
     const BL_FORT_FAB_ARG(crse),
     BL_FORT_FAB_ARG(fine));

BL_FORT_PROC_DECL(CA_AVERAGE_EC,ca_average_ec)
    (const BL_FORT_FAB_ARG(grad_fine),
           BL_FORT_FAB_ARG(grad_crse),
     const int lo[], const int hi[], const int rat[], const int& idir);

BL_FORT_PROC_DECL(CA_TEST_RESIDUAL,ca_test_residual)
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(rhs),
     D_DECL(const BL_FORT_FAB_ARG(xgrad),
            const BL_FORT_FAB_ARG(ygrad),
            const BL_FORT_FAB_ARG(zgrad)), 
            const Real* dx, const Real* problo, 
            const int* coord_type);

BL_FORT_PROC_DECL(CA_APPLY_METRIC,ca_apply_metric)
    (const int* lo, const int* hi,
     D_DECL(const int* xlo, const int *ylo, const int *zlo),
     D_DECL(const int* xhi, const int *yhi, const int *zhi),
     BL_FORT_FAB_ARG(cc),
     D_DECL(const BL_FORT_FAB_ARG(xedge),
            const BL_FORT_FAB_ARG(yedge),
            const BL_FORT_FAB_ARG(zedge)), 
            const Real* dx, const int* coord_type);

BL_FORT_PROC_DECL(CA_WEIGHT_CC,ca_weight_cc)
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(cc),
     const Real* dx, const int* coord_type);

BL_FORT_PROC_DECL(CA_UNWEIGHT_CC,ca_unweight_cc)
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(cc),
     const Real* dx, const int* coord_type);

BL_FORT_PROC_DECL(CA_UNWEIGHT_EDGES,ca_unweight_edges)
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(edge),
     const Real* dx, const int* coord_type, const int* idir);

BL_FORT_PROC_DECL(CA_COMPUTE_1D_GRAV,ca_compute_1d_grav)
    (const BL_FORT_FAB_ARG(den), const int& lo, const int& hi,
     const Real* grav, const Real* phi, const Real* dx, const Real* problo);

#ifdef GR_GRAV
BL_FORT_PROC_DECL(CA_COMPUTE_1D_GR_GRAV,ca_compute_1d_gr_grav)
    (const BL_FORT_FAB_ARG(state),
     const Real* grav, const Real* dx,
     const Real* problo);
#endif

BL_FORT_PROC_DECL(CA_COMPUTE_RADIAL_MASS,ca_compute_radial_mass)
    (const int lo[], const int hi[], 
     const Real* dx, const Real* dr,
     const BL_FORT_FAB_ARG(S), 
     const Real* avgmass, const Real* avgvol, 
     const Real* problo, const int* numpts_1d, 
     const int* drdxfac, const int* level); 

BL_FORT_PROC_DECL(CA_COMPUTE_AVGPRES,ca_compute_avgpres)
    (const int lo[], const int hi[], 
     const Real* dx, const Real* dr,
     const BL_FORT_FAB_ARG(S), 
     const Real* avgpres, 
     const Real* problo, const int* numpts_1d, 
     const int* drdxfac, const int* level); 

BL_FORT_PROC_DECL(CA_INTEGRATE_GRAV,ca_integrate_grav)
    (const Real* avgmass, const Real* avgden,
     const Real* radial_grav,
     const Real* max_radius,
     const Real* dr, const int* numpts_1d);

BL_FORT_PROC_DECL(CA_PRESCRIBE_GRAV,ca_prescribe_grav)
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(S),
     const Real* dx);

#ifdef GR_GRAV
BL_FORT_PROC_DECL(CA_INTEGRATE_GR_GRAV,ca_integrate_gr_grav)
    (const Real* avgden, const Real* avgmass, 
     const Real* avgpres, const Real* radial_grav, 
     const Real* dr, const int* numpts_1d);
#endif

BL_FORT_PROC_DECL(CA_INTEGRATE_PHI,ca_integrate_phi)
    (const Real* avgden, const Real* radial_grav,
     const Real* radial_phi, 
     const Real* dr, const int* numpts_1d);

BL_FORT_PROC_DECL(CA_PUT_RADIAL_GRAV,ca_put_radial_grav)
    (const int lo[], const int hi[], 
     const Real* dx, const Real* dr,
     const BL_FORT_FAB_ARG(grav_vector),
     const Real* radial_grav, 
     const Real* problo, 
     const int* n1d, const int* level); 

BL_FORT_PROC_DECL(CA_PUT_RADIAL_PHI,ca_put_radial_phi)
    (const int lo[], const int hi[], 
     const int domlo[], const int domhi[], 
     const Real* dx, const Real* dr,
     const BL_FORT_FAB_ARG(phi),
     const Real* radial_phi,
     const Real* problo, 
     const int* numpts_1d, const int* fill_interior);

BL_FORT_PROC_DECL(CA_PUT_MULTIPOLE_PHI,ca_put_multipole_phi)
    (const int lo[], const int hi[], 
     const int domlo[], const int domhi[], 
     const Real* dx,
     const BL_FORT_FAB_ARG(phi),
     const int* lnum,
     Real* qL0, Real* qLC, Real* qLS,
     Real* qU0, Real* qUC, Real* qUS,
     const int* npts, const int* boundary_only); 

BL_FORT_PROC_DECL(CA_COMPUTE_MULTIPOLE_MOMENTS,ca_compute_multipole_moments)
    (const int lo[], const int hi[], 
     const int domlo[], const int domhi[], 
     const int* symmetry_type, const int lo_bc[], const int hi_bc[],
     const Real* dx,
     const BL_FORT_FAB_ARG(rho),
     const int* lnum,
     Real* qL0, Real* qLC, Real* qLS,
     Real* qU0, Real* qUC, Real* qUS,
     const int* npts, const int* boundary_only); 

BL_FORT_PROC_DECL(CA_COMPUTE_DIRECT_SUM_BC,ca_compute_direct_sum_bc)
   (const int lo[], const int hi[],
    const int domlo[], const int domhi[],
    const int* symmetry_type, const int lo_bc[], const int hi_bc[],
    const Real* dx, const BL_FORT_FAB_ARG(rho),
    const Real* problo, const Real* probhi,
    Real* bcXYLo, Real* bcXYHi,
    Real* bcXZLo, Real* bcXZHi,
    Real* bcYZLo, Real* bcYZHi);

BL_FORT_PROC_DECL(CA_PUT_DIRECT_SUM_BC,ca_put_direct_sum_bc)
   (const int lo[], const int hi[], 
    const int domlo[], const int domhi[],
    const BL_FORT_FAB_ARG(rho),
    Real* bcXYLo, Real* bcXYHi,
    Real* bcXZLo, Real* bcXZHi,
    Real* bcYZLo, Real* bcYZHi);

#ifdef POINTMASS
BL_FORT_PROC_DECL(PM_ADD_TO_GRAV,pm_add_to_grav)
    (const Real* point_mass, BL_FORT_FAB_ARG(grav_vector),
     const Real* problo, const Real* dx, const int* lo, const int* hi);
#endif

#endif
