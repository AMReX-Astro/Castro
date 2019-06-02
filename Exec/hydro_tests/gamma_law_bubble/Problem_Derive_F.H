#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void derpi
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* dx, const amrex::Real time);

  void derpioverp0
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* dx, const amrex::Real time);

  void derrhopert
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* dx, const amrex::Real time);

  void dertpert
    (BL_FORT_FAB_ARG_3D(der), const int nvar,
     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* dx, const amrex::Real time);

#ifdef __cplusplus
}
#endif
