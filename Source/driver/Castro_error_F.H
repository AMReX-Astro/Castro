#ifndef CASTRO_ERROR_F_H
#define CASTRO_ERROR_F_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  AMREX_GPU_DEVICE
  void ca_denerror
    (int i, int j, int k,
     int8_t* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(den), const int ncomp,
     const amrex::Real* dx, const amrex::Real* problo,
     const int8_t tagval, const int8_t clearval,
     const amrex::Real time, const int level);

  AMREX_GPU_DEVICE
  void ca_velerror
    (int i, int j, int k,
     int8_t* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(vel), const int ncomp,
     const amrex::Real* dx, const amrex::Real* problo,
     const int8_t tagval, const int8_t clearval,
     const amrex::Real time, const int level);

  AMREX_GPU_DEVICE
  void ca_temperror
    (int i, int j, int k,
     int8_t* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(temp), const int ncomp,
     const amrex::Real* dx, const amrex::Real* problo,
     const int8_t tagval, const int8_t clearval,
     const amrex::Real time, const int level);

  AMREX_GPU_DEVICE
  void ca_presserror
    (int i, int j, int k,
     int8_t* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(press), const int ncomp,
     const amrex::Real* dx, const amrex::Real* problo,
     const int8_t tagval, const int8_t clearval,
     const amrex::Real time, const int level);

  AMREX_GPU_DEVICE
  void ca_nucerror
    (int i, int j, int k,
     int8_t* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(state), const int ncomp,
     const amrex::Real* dx, const amrex::Real* problo,
     const int8_t tagval, const int8_t clearval,
     const amrex::Real time, const int level);

  AMREX_GPU_DEVICE
  void ca_enucerror
    (int i, int j, int k,
     int8_t* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(state), const int ncomp,
     const amrex::Real* dx, const amrex::Real* problo,
     const int8_t tagval, const int8_t clearval,
     const amrex::Real time, const int level);

#ifdef RADIATION
  AMREX_GPU_DEVICE
  void ca_raderror
    (int i, int j, int k,
     int8_t* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(rad), const int ncomp,
     const amrex::Real* dx, const amrex::Real* problo,
     const int8_t tagval, const int8_t clearval,
     const amrex::Real time, const int level);
#endif

  void set_problem_tags
    (const int* lo, const int* hi,
     int8_t* tag, const int* tag_lo, const int* tag_hi,
     const BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, const amrex::Real* problo,
     const int8_t tagval, const int8_t clearval,
     const amrex::Real time, const int level);

#ifdef __cplusplus
}
#endif

#endif
