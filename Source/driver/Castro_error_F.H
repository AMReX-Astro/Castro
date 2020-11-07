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

  void get_denerr_params
    (amrex::Real* denerr, int* max_denerr_lev,
     amrex::Real* dengrad, int* max_dengrad_lev,
     amrex::Real* dengrad_rel, int* max_dengrad_rel_lev);

  void get_temperr_params
    (amrex::Real* temperr, int* max_temperr_lev,
     amrex::Real* tempgrad, int* max_tempgrad_lev,
     amrex::Real* tempgrad_rel, int* max_tempgrad_rel_lev);

  void get_velerr_params
    (amrex::Real* velerr, int* max_velerr_lev,
     amrex::Real* velgrad, int* max_velgrad_lev,
     amrex::Real* velgrad_rel, int* max_velgrad_rel_lev);

  void get_presserr_params
    (amrex::Real* presserr, int* max_presserr_lev,
     amrex::Real* pressgrad, int* max_pressgrad_lev,
     amrex::Real* pressgrad_rel, int* max_pressgrad_rel_lev);

  void get_raderr_params
    (amrex::Real* raderr, int* max_raderr_lev,
     amrex::Real* radgrad, int* max_radgrad_lev,
     amrex::Real* radgrad_rel, int* max_radgrad_rel_lev);

  void get_dxnuc_params
    (amrex::Real* dxnuc_min, amrex::Real* dxnuc_max, int* max_dxnuc_lev);

  void get_enuc_params
    (amrex::Real* enucerr, int* max_enucerr_lev);

#ifdef __cplusplus
}
#endif

#endif
