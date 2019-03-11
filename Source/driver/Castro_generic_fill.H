#ifndef _Castro_generic_fill_H_
#define _Castro_generic_fill_H_

#include <AMReX_BLFort.H>

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
        // Ensure that our threadblock size is such that it is
        // evenly divisible by the number of zones in the box,
        // and is at least one larger than the number of ghost zones.
        // This ensures that the corners plus one interior zone
        // are all on the same threadblock.

        int left[3] = {0, 0, 0};
        int rght[3] = {0, 0, 0};
        int sz[3] = {0, 0, 0};

        int ng[3] = {0, 0, 0};

        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            left[n] = domlo[n] - lo[n];
            rght[n] = hi[n] - domhi[n];
            sz[n] = hi[n] - lo[n] + 1;
            ng[n] = std::max(0, std::max(left[n], rght[n]));
        }

        int numThreadsMin[3] = {ng[0] + 1, ng[1] + 1, ng[2] + 1};

        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            while (sz[n] % numThreadsMin[n] != 0) {
                ++numThreadsMin[n];
            }
        }

        if (std::min({numThreadsMin[0], numThreadsMin[1], numThreadsMin[2]}) < 1) {
            amrex::Error("Minimum number of CUDA threads must be positive.");
        }

        amrex::Gpu::Device::setNumThreadsMin(numThreadsMin[0], numThreadsMin[1], numThreadsMin[2]);

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
