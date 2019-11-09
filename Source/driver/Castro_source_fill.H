#ifndef _Castro_source_fill_H_
#define _Castro_source_fill_H_

#include <AMReX_BLFort.H>
#include <Castro.H>

#ifdef AMREX_USE_CUDA
#include <cuda_runtime_api.h>
#include <AMReX_Arena.H>
#endif

#ifdef __cplusplus
extern "C"
{
#endif

    void ca_source_single_fill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

    void ca_source_multi_fill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

#ifdef __cplusplus
}
#endif

#endif
