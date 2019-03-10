#ifndef _Castro_generic_fill_F_H
#define _Castro_generic_fill_F_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void prepare_nvar_bc(const int* bc);

  void prepare_bc(const int* bc);

  void generic_single_fill
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo);

  void generic_multi_fill
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo);

#ifdef __cplusplus
}
#endif

#endif
