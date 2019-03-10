
#include <AMReX_BLFort.H>
#include <Castro.H>
#include <Castro_generic_fill.H>
#include <Castro_generic_fill_F.H>

using namespace amrex;

#ifdef __cplusplus
extern "C"
{
#endif

    // Note that these are called with dimension agnostic macros like
    // AMREX_ZFILL and AMREX_INT_ANYD already, so we should expect that
    // everything has three entries, not AMREX_SPACEDIM entries.
    // We still choose to use the macros anyway below, for compatibility
    // with the GPU pragma script.

    void ca_generic_single_fill(Real* adv, const int* adv_lo, const int* adv_hi,
                                const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                                const Real* time, const int* bc)
    {
        int lo[3];
        int hi[3];

        for (int i = 0; i < 3; ++i) {
            lo[i] = adv_lo[i];
            hi[i] = adv_hi[i];
        }

        int* bc_f = prepare_bc(bc, 1, lo, hi, domlo, domhi);

#pragma gpu
        generic_single_fill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                            adv, AMREX_INT_ANYD(adv_lo), AMREX_INT_ANYD(adv_hi),
                            AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
                            AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), bc_f);

        clean_bc(bc_f);
    }

    void ca_generic_multi_fill(Real* adv, const int* adv_lo, const int* adv_hi,
                               const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                               const Real* time, const int* bc)
    {
        int lo[3];
        int hi[3];

        for (int i = 0; i < 3; ++i) {
            lo[i] = adv_lo[i];
            hi[i] = adv_hi[i];
        }

        int* bc_f = prepare_bc(bc, Castro::numState(), lo, hi, domlo, domhi);

#pragma gpu
        generic_multi_fill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                           adv, AMREX_INT_ANYD(adv_lo), AMREX_INT_ANYD(adv_hi),
                           AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
                           AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), bc_f);

        clean_bc(bc_f);
    }

#ifdef __cplusplus
}
#endif
