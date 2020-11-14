#ifndef CASTRO_REACT_F_H
#define CASTRO_REACT_F_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void ca_react_state
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(reactions),
     const amrex::Real time, const amrex::Real dt_react,
     amrex::Real* burn_failure);


#ifdef __cplusplus
}
#endif

#endif
