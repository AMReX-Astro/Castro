
#include <AMReX_BLFort.H>
#include <Castro.H>
#include <Castro_bc_fill_nd.H>
#include <Castro_bc_fill_nd_F.H>
#include <Castro_bc_ext_fill_nd.H>
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

  void ca_hypfill(Real* adv, const int* adv_lo, const int* adv_hi,
                  const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                  const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = adv_lo[i];
      hi[i] = adv_hi[i];
    }

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bc, NUM_STATE);
    set_bc_launch_config();
#else
    const int* bc_f = bc;
#endif

    IntVect ilo(D_DECL(lo[0], lo[1], lo[2]));
    IntVect ihi(D_DECL(hi[0], hi[1], hi[2]));

    Box bx(ilo, ihi);

#pragma gpu box(bx)
    hypfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
            adv, AMREX_INT_ANYD(adv_lo), AMREX_INT_ANYD(adv_hi),
            AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
            AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc_launch_config();
#endif

    ca_ext_fill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, dx, xlo, time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc(bc_f);
#endif
  }


  void ca_denfill(Real* adv, const int* adv_lo, const int* adv_hi,
                  const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                  const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = adv_lo[i];
      hi[i] = adv_hi[i];
    }

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bc, 1);
    set_bc_launch_config();
#else
    const int* bc_f = bc;
#endif

    IntVect ilo(D_DECL(lo[0], lo[1], lo[2]));
    IntVect ihi(D_DECL(hi[0], hi[1], hi[2]));

    Box bx(ilo, ihi);

#pragma gpu box(bx)
    denfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
            adv, AMREX_INT_ANYD(adv_lo), AMREX_INT_ANYD(adv_hi),
            AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
            AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc_launch_config();
#endif

    ca_ext_denfill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, dx, xlo, time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc(bc_f);
#endif
  }

#ifdef __cplusplus
}
#endif
