#ifndef _Castro_generic_fill_H_
#define _Castro_generic_fill_H_

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

#ifdef AMREX_USE_CUDA
    static void set_bc_launch_config()
    {
        amrex::Gpu::Device::setNumThreadsMin(Castro::minBCThreads(0), Castro::minBCThreads(1), Castro::minBCThreads(2));
    }

    static void clean_bc_launch_config()
    {
        amrex::Gpu::Device::setNumThreadsMin(1, 1, 1);
    }

    // Return a pointer to bc valid for use in Fortran. For the CPU this is a no-op.

    static int* prepare_bc(const int* bc, const int nvar)
    {
        int* bc_f = (int*) amrex::The_Arena()->alloc(AMREX_SPACEDIM * 2 * nvar * sizeof(int));
        AMREX_GPU_SAFE_CALL(cudaMemcpyAsync(bc_f, bc, AMREX_SPACEDIM * 2 * nvar * sizeof(int), cudaMemcpyHostToDevice, amrex::Gpu::Device::cudaStream()));
        return bc_f;
    }

    static void clean_bc(int* bc_f)
    {
        amrex::Gpu::Device::streamSynchronize();
        amrex::The_Arena()->free(bc_f);
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
