#ifndef _Castro_bc_ext_fill_nd_H_
#define _Castro_bc_ext_fill_nd_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void ca_ext_fill
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

    void ca_ext_denfill
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

#ifdef GRAVITY
    void ca_ext_gravxfill
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

    void ca_ext_gravyfill
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

    void ca_ext_gravzfill
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);
#endif

#ifdef __cplusplus
}
#endif

#endif
