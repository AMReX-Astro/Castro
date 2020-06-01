#ifndef CASTRO_BC_EXT_FILL_ND_F_H
#define CASTRO_BC_EXT_FILL_ND_F_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void ext_fill
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real time, const int* bc);

  void ext_denfill
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real time, const int* bc);

#ifdef __cplusplus
}
#endif

#endif
