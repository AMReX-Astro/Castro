#ifndef _Derive_F_H_
#define _Derive_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void derpres
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

  void dereint1
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

  void dereint2
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

  void derlogden
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

  void deruplusc
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

  void deruminusc
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

  void dersoundspeed
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

  void dergamma1
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

  void dermachnumber
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

  void derentropy
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

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

#ifdef REACTIONS
  void derenuctimescale
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);
#endif

  void dervel
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

  void dermagvel
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

  void dermaggrav
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

  void derradialvel
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

  void dermagmom
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

  void derangmomx
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

  void derangmomy
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

  void derangmomz
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

  void derspec
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

  void derabar
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

  void dermagvort
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

  void derdivu
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

  void derstate
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta);

#ifdef __cplusplus
}
#endif
  
/* problem-specific includes */
#include <Problem_Derive_F.H>

#endif
