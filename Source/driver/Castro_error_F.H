#ifndef _Castro_error_F_H_
#define _Castro_error_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void ca_denerror
    (const int* lo, const int* hi,
     char* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(den), const int ncomp,
     const amrex::Real* dx, const amrex::Real* problo,
     const char tagval, const char clearval,
     const amrex::Real time, const int level);

  void ca_velerror
    (const int* lo, const int* hi,
     char* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(vel), const int ncomp,
     const amrex::Real* dx, const amrex::Real* problo,
     const char tagval, const char clearval,
     const amrex::Real time, const int level);

  void ca_temperror
    (const int* lo, const int* hi,
     char* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(temp), const int ncomp,
     const amrex::Real* dx, const amrex::Real* problo,
     const char tagval, const char clearval,
     const amrex::Real time, const int level);

  void ca_presserror
    (const int* lo, const int* hi,
     char* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(press), const int ncomp,
     const amrex::Real* dx, const amrex::Real* problo,
     const char tagval, const char clearval,
     const amrex::Real time, const int level);

  void ca_nucerror
    (const int* lo, const int* hi,
     char* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(state), const int ncomp,
     const amrex::Real* dx, const amrex::Real* problo,
     const char tagval, const char clearval,
     const amrex::Real time, const int level);

  void ca_enucerror
    (const int* lo, const int* hi,
     char* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(state), const int ncomp,
     const amrex::Real* dx, const amrex::Real* problo,
     const char tagval, const char clearval,
     const amrex::Real time, const int level);

#ifdef RADIATION
  void ca_raderror
    (const int* lo, const int* hi,
     char* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(rad), const int ncomp,
     const amrex::Real* dx, const amrex::Real* problo,
     const char tagval, const char clearval,
     const amrex::Real time, const int level);
#endif

#ifdef AMREX_DIMENSION_AGNOSTIC
  void set_problem_tags
    (const int* lo, const int* hi,
     char* tag, const int* tag_lo, const int* tag_hi,
     const BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, const amrex::Real* problo,
     const char tagval, const char clearval,
     const amrex::Real time, const int level);
#else
  void set_problem_tags
    (const int* lo, const int* hi,
     char* tag, const int* tag_lo, const int* tag_hi,
     const BL_FORT_FAB_ARG(state),
     const amrex::Real* dx, const amrex::Real* problo,
     const char tagval, const char clearval,
     const amrex::Real time, const int level);
#endif

#ifdef __cplusplus
}
#endif

#endif
