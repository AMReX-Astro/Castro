#ifndef _Castro_generic_fill_H_
#define _Castro_generic_fill_H_

#include <AMReX_BLFort.H>
#include <Castro.H>

#ifdef AMREX_USE_CUDA
#include <cuda_runtime_api.h>
#endif

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef AMREX_USE_CUDA
    static void set_bc_launch_config(const int* lo, const int* hi, const int* domlo, const int* domhi)
    {
        amrex::Gpu::Device::setNumThreadsMin(Castro::minBCThreads(0), Castro::minBCThreads(1), Castro::minBCThreads(2));
    }

    static void clean_bc_launch_config()
    {
        amrex::Gpu::Device::setNumThreadsMin(1, 1, 1);
    }

    // Return a pointer to bc valid for use in Fortran. For the CPU this is a no-op.

    static int* prepare_bc(const int* bc, const int nvar, const int* lo, const int* hi, const int* domlo, const int* domhi)
    {
        set_bc_launch_config(lo, hi, domlo, domhi);

        int* bc_f;
        AMREX_GPU_SAFE_CALL(cudaMalloc((void**) &bc_f, AMREX_SPACEDIM * 2 * nvar * sizeof(int)));
        AMREX_GPU_SAFE_CALL(cudaMemcpy(bc_f, bc, AMREX_SPACEDIM * 2 * nvar * sizeof(int), cudaMemcpyHostToDevice));
        return bc_f;
    }

    static void clean_bc(int* bc_f)
    {
        AMREX_GPU_SAFE_CALL(cudaFree((void*) bc_f));

        clean_bc_launch_config();
    }
#endif

    void ca_generic_single_fill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

    void ca_generic_multi_fill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

#ifdef __cplusplus
}
#endif

#endif
