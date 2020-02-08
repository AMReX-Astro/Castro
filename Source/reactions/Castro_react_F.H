#ifndef _Castro_react_F_H_
#define _Castro_react_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void ca_react_state_simplified_sdc
    (const int* lo, const int* hi,
     const amrex::Real* uold, const int* uo_lo, const int* uo_hi,
     amrex::Real* unew, const int* un_lo, const int* un_hi,
     const amrex::Real* asrc, const int* as_lo, const int* as_hi,
     amrex::Real* reactions, const int* r_lo, const int* r_hi,
     const int* mask, const int* m_lo, const int* m_hi,
     const amrex::Real time, const amrex::Real dt_react, const int sdc_iter,
     amrex::Real* burn_failure);

  void ca_react_state
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(reactions),
     BL_FORT_FAB_ARG_3D(weights),
     const BL_FORT_IFAB_ARG_3D(mask),
     const amrex::Real time, const amrex::Real dt_react, const int strang_half,
     amrex::Real* burn_failure);


#ifdef __cplusplus
}
#endif

#endif
