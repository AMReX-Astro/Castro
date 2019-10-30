#include <AMReX_BLFort.H>
#include <Castro.H>
#include <Castro_bc_fill_nd.H>
#include <Castro_bc_ext_fill_nd_F.H>

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

    void ca_ext_fill(const int* lo, const int* hi,
                     Real* adv, const int* adv_lo, const int* adv_hi,
                     const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                     const Real* time, const int* bc)
    {

#pragma gpu
        ext_fill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                 adv, AMREX_INT_ANYD(adv_lo), AMREX_INT_ANYD(adv_hi),
                 AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
                 AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc);

    }

    void ca_ext_denfill(const int* lo, const int* hi,
                        Real* adv, const int* adv_lo, const int* adv_hi,
                        const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                        const Real* time, const int* bc)
    {

#pragma gpu
      ext_denfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                  adv, AMREX_INT_ANYD(adv_lo), AMREX_INT_ANYD(adv_hi),
                  AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
                  AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc);

    }


#ifdef GRAVITY
    void ca_ext_gravxfill(const int* lo, const int* hi,
                          Real* adv, const int* adv_lo, const int* adv_hi,
                          const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                          const Real* time, const int* bc)
    {

#pragma gpu
      ext_gravxfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                    adv, AMREX_INT_ANYD(adv_lo), AMREX_INT_ANYD(adv_hi),
                    AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
                    AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc);

    }

    void ca_ext_gravyfill(const int* lo, const int* hi,
                          Real* adv, const int* adv_lo, const int* adv_hi,
                          const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                          const Real* time, const int* bc)
    {

#pragma gpu
      ext_gravyfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                    adv, AMREX_INT_ANYD(adv_lo), AMREX_INT_ANYD(adv_hi),
                    AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
                    AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc);

    }

    void ca_ext_gravzfill(const int* lo, const int* hi,
                          Real* adv, const int* adv_lo, const int* adv_hi,
                          const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                          const Real* time, const int* bc)
    {

#pragma gpu
      ext_gravzfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                    adv, AMREX_INT_ANYD(adv_lo), AMREX_INT_ANYD(adv_hi),
                    AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
                    AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc);

    }
#endif


#ifdef __cplusplus
}
#endif
