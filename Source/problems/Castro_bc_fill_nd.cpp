
#include <AMReX_BLFort.H>
#include <Castro_bc_fill_nd.H>
#include <Castro_bc_fill_nd_F.H>

using namespace amrex;

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef AMREX_USE_CUDA
    void set_bc_launch_config(const int* lo, const int* hi, const int* domlo, const int* domhi)
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

        Gpu::Device::setNumThreadsMin(numThreadsMin[0], numThreadsMin[1], numThreadsMin[2]);

    }

    void clean_bc_launch_config()
    {
        Gpu::Device::setNumThreadsMin(1, 1, 1);
    }
#endif

    // Note that these are called with dimension agnostic macros like
    // AMREX_ZFILL and AMREX_INT_ANYD already, so we should expect that
    // everything has three entries, not AMREX_SPACEDIM entries.
    // We still choose to use the macros anyway below, for compatibility
    // with the GPU pragma script.

    void ca_hypfill(Real* adv, const int* adv_lo, const int* adv_hi,
                    const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                    const Real* time, const int* bc)
    {
        int lo[3];
        int hi[3];

        for (int i = 0; i < 3; ++i) {
            lo[i] = adv_lo[i];
            hi[i] = adv_hi[i];
        }

        prepare_nvar_bc(bc);

#ifdef AMREX_USE_CUDA
        set_bc_launch_config(lo, hi, domlo, domhi);
#endif

#pragma gpu
        hypfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                adv, AMREX_INT_ANYD(adv_lo), AMREX_INT_ANYD(adv_hi),
                AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
                AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time);

#ifdef AMREX_USE_CUDA
        clean_bc_launch_config();
#endif
    }

    void ca_denfill(Real* adv, const int* adv_lo, const int* adv_hi,
                    const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                    const Real* time, const int* bc)
    {
        int lo[3];
        int hi[3];

        for (int i = 0; i < 3; ++i) {
            lo[i] = adv_lo[i];
            hi[i] = adv_hi[i];
        }

        prepare_bc(bc);

#ifdef AMREX_USE_CUDA
        set_bc_launch_config(lo, hi, domlo, domhi);
#endif

#pragma gpu
        denfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                adv, AMREX_INT_ANYD(adv_lo), AMREX_INT_ANYD(adv_hi),
                AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
                AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time);

#ifdef AMREX_USE_CUDA
        clean_bc_launch_config();
#endif
    }

#ifdef GRAVITY
    void ca_phigravfill(Real* phi, const int* phi_lo, const int* phi_hi,
                        const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                        const Real* time, const int* bc)
    {
        int lo[3];
        int hi[3];

        for (int i = 0; i < 3; ++i) {
            lo[i] = phi_lo[i];
            hi[i] = phi_hi[i];
        }

        prepare_bc(bc);

        phigravfill(lo, hi, phi, phi_lo, phi_hi, domlo, domhi, dx, xlo, *time);
    }

    void ca_gravxfill(Real* grav, const int* grav_lo, const int* grav_hi,
                      const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                      const Real* time, const int* bc)
    {
        int lo[3];
        int hi[3];

        for (int i = 0; i < 3; ++i) {
            lo[i] = grav_lo[i];
            hi[i] = grav_hi[i];
        }

        prepare_bc(bc);

        gravxfill(lo, hi, grav, grav_lo, grav_hi, domlo, domhi, dx, xlo, *time);
    }

    void ca_gravyfill(Real* grav, const int* grav_lo, const int* grav_hi,
                      const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                      const Real* time, const int* bc)
    {
        int lo[3];
        int hi[3];

        for (int i = 0; i < 3; ++i) {
            lo[i] = grav_lo[i];
            hi[i] = grav_hi[i];
        }

        prepare_bc(bc);

        gravyfill(lo, hi, grav, grav_lo, grav_hi, domlo, domhi, dx, xlo, *time);
    }

    void ca_gravzfill(Real* grav, const int* grav_lo, const int* grav_hi,
                      const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                      const Real* time, const int* bc)
    {
        int lo[3];
        int hi[3];

        for (int i = 0; i < 3; ++i) {
            lo[i] = grav_lo[i];
            hi[i] = grav_hi[i];
        }

        prepare_bc(bc);

        gravzfill(lo, hi, grav, grav_lo, grav_hi, domlo, domhi, dx, xlo, *time);
    }
#endif

#ifdef ROTATION
    void ca_phirotfill(Real* phi, const int* phi_lo, const int* phi_hi,
                       const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                       const Real* time, const int* bc)
    {
        int lo[3];
        int hi[3];

        for (int i = 0; i < 3; ++i) {
            lo[i] = phi_lo[i];
            hi[i] = phi_hi[i];
        }

        prepare_bc(bc);

        phirotfill(lo, hi, phi, phi_lo, phi_hi, domlo, domhi, dx, xlo, *time);
    }

    void ca_rotxfill(Real* rot, const int* rot_lo, const int* rot_hi,
                     const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                     const Real* time, const int* bc)
    {
        int lo[3];
        int hi[3];

        for (int i = 0; i < 3; ++i) {
            lo[i] = rot_lo[i];
            hi[i] = rot_hi[i];
        }

        rotxfill(lo, hi, rot, rot_lo, rot_hi, domlo, domhi, dx, xlo, *time);
    }

    void ca_rotyfill(Real* rot, const int* rot_lo, const int* rot_hi,
                     const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                     const Real* time, const int* bc)
    {
        int lo[3];
        int hi[3];

        for (int i = 0; i < 3; ++i) {
            lo[i] = rot_lo[i];
            hi[i] = rot_hi[i];
        }

        prepare_bc(bc);

        rotyfill(lo, hi, rot, rot_lo, rot_hi, domlo, domhi, dx, xlo, *time);
    }

    void ca_rotzfill(Real* rot, const int* rot_lo, const int* rot_hi,
                     const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                     const Real* time, const int* bc)
    {
        int lo[3];
        int hi[3];

        for (int i = 0; i < 3; ++i) {
            lo[i] = rot_lo[i];
            hi[i] = rot_hi[i];
        }

        prepare_bc(bc);

        rotzfill(lo, hi, rot, rot_lo, rot_hi, domlo, domhi, dx, xlo, *time);
    }
#endif

#ifdef REACTIONS
    void ca_reactfill(Real* react, const int* react_lo, const int* react_hi,
                      const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                      const Real* time, const int* bc)
    {
        int lo[3];
        int hi[3];

        for (int i = 0; i < 3; ++i) {
            lo[i] = react_lo[i];
            hi[i] = react_hi[i];
        }

        prepare_bc(bc);

        reactfill(lo, hi, react, react_lo, react_hi, domlo, domhi, dx, xlo, *time);
    }
#endif

#ifdef RADIATION
    void ca_radfill(Real* rad, const int* rad_lo, const int* rad_hi,
                    const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                    const Real* time, const int* bc)
    {
        int lo[3];
        int hi[3];

        for (int i = 0; i < 3; ++i) {
            lo[i] = rad_lo[i];
            hi[i] = rad_hi[i];
        }

        prepare_bc(bc);

        radfill(lo, hi, rad, rad_lo, rad_hi, domlo, domhi, dx, xlo, *time);
    }
#endif

#ifdef __cplusplus
}
#endif
