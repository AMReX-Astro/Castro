#ifndef _Derive_F_H_
#define _Derive_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef DIFFUSION
  void dercond
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

  void derdiffcoeff
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

  void derdiffterm
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const BL_FORT_FAB_ARG_3D(coeff_x),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(coeff_y),
#endif
#if AMREX_SPACEDIM == 3
     const BL_FORT_FAB_ARG_3D(coeff_z),
#endif
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);
#endif

#ifdef MHD
  void ca_dermagcenx
    (BL_FORT_FAB_ARG(der),const int* nvar,
     const BL_FORT_FAB_ARG(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
     const int* level, const int* grid_no);

  void ca_dermagceny
    (BL_FORT_FAB_ARG(der), const int* nvar,
     const BL_FORT_FAB_ARG(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
     const int* level, const int* grid_no);

  void ca_dermagcenz
    (BL_FORT_FAB_ARG(der), const int* nvar,
     const BL_FORT_FAB_ARG(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
     const int* level, const int* grid_no);

  void ca_derex
    (BL_FORT_FAB_ARG(der), const int* nvar,
     const BL_FORT_FAB_ARG(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
     const int* level, const int* grid_no);

  void ca_derey
    (BL_FORT_FAB_ARG(der), const int* nvar,
     const BL_FORT_FAB_ARG(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
     const int* level, const int* grid_no);

  void ca_derez
    (BL_FORT_FAB_ARG(der), const int* nvar,
     const BL_FORT_FAB_ARG(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
     const int* level, const int* grid_no);

  void ca_derdivb
    (BL_FORT_FAB_ARG(der), const int* nvar,
     const BL_FORT_FAB_ARG(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
     const int* level, const int* grid_no);

#endif



#ifdef __cplusplus
}
#endif
  
/* problem-specific includes */
#include <Problem_Derive_F.H>

#endif
