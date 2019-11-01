
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
    int* bc_f = prepare_bc(bc, Castro::numState());
    set_bc_launch_config();
#else
    const int* bc_f = bc;
#endif

#pragma gpu
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

#pragma gpu
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

#ifdef MHD
  void ca_face_fillx(Real* var, const int* var_lo, const int* var_hi,
                     const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                     const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = var_lo[i];
      hi[i] = var_hi[i];
    }

    const int* bc_f = bc;


    face_fillx(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
               var, AMREX_INT_ANYD(var_lo), AMREX_INT_ANYD(var_hi),
               AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
               AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

  }

  void ca_face_filly(Real* var, const int* var_lo, const int* var_hi,
                     const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                     const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = var_lo[i];
      hi[i] = var_hi[i];
    }

    const int* bc_f = bc;


    face_filly(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
               var, AMREX_INT_ANYD(var_lo), AMREX_INT_ANYD(var_hi),
               AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
               AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

  }

  void ca_face_fillz(Real* var, const int* var_lo, const int* var_hi,
                     const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                     const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = var_lo[i];
      hi[i] = var_hi[i];
    }

    const int* bc_f = bc;


    face_fillz(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
               var, AMREX_INT_ANYD(var_lo), AMREX_INT_ANYD(var_hi),
               AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
               AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

  }
#endif  

#ifdef GRAVITY
  void ca_phigravfill(Real* phi, const int* phi_lo, const int* phi_hi,
                      const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                      const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = phi_lo[i];
      hi[i] = phi_hi[i];
    }

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bc, 1);
    set_bc_launch_config();
#else
    const int* bc_f = bc;
#endif

#pragma gpu
    phigravfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                phi, AMREX_INT_ANYD(phi_lo), AMREX_INT_ANYD(phi_hi),
                AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
                AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc_launch_config();
    clean_bc(bc_f);
#endif
  }

  void ca_gravxfill(Real* grav, const int* grav_lo, const int* grav_hi,
                    const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                    const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = grav_lo[i];
      hi[i] = grav_hi[i];
    }

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bc, 1);
    set_bc_launch_config();
#else
    const int* bc_f = bc;
#endif

#pragma gpu
    gravxfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
              grav, AMREX_INT_ANYD(grav_lo), AMREX_INT_ANYD(grav_hi),
              AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
              AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc_launch_config();
#endif

    ca_ext_gravxfill(lo, hi, grav, grav_lo, grav_hi, domlo, domhi, dx, xlo, time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc(bc_f);
#endif

  }

  void ca_gravyfill(Real* grav, const int* grav_lo, const int* grav_hi,
                    const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                    const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = grav_lo[i];
      hi[i] = grav_hi[i];
    }

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bc, 1);
    set_bc_launch_config();
#else
    const int* bc_f = bc;
#endif

#pragma gpu
    gravyfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
              grav, AMREX_INT_ANYD(grav_lo), AMREX_INT_ANYD(grav_hi),
              AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
              AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc_launch_config();
#endif

    ca_ext_gravyfill(lo, hi, grav, grav_lo, grav_hi, domlo, domhi, dx, xlo, time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc(bc_f);
#endif

  }

  void ca_gravzfill(Real* grav, const int* grav_lo, const int* grav_hi,
                    const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                    const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = grav_lo[i];
      hi[i] = grav_hi[i];
    }

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bc, 1);
    set_bc_launch_config();
#else
    const int* bc_f = bc;
#endif

#pragma gpu
    gravzfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
              grav, AMREX_INT_ANYD(grav_lo), AMREX_INT_ANYD(grav_hi),
              AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
              AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc_launch_config();
#endif

    ca_ext_gravzfill(lo, hi, grav, grav_lo, grav_hi, domlo, domhi, dx, xlo, time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc(bc_f);
#endif

  }
#endif

#ifdef ROTATION
  void ca_phirotfill(Real* phi, const int* phi_lo, const int* phi_hi,
                     const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                     const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = phi_lo[i];
      hi[i] = phi_hi[i];
    }

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bc, 1);
    set_bc_launch_config();
#else
    const int* bc_f = bc;
#endif

#pragma gpu
    phirotfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
               phi, AMREX_INT_ANYD(phi_lo), AMREX_INT_ANYD(phi_hi),
               AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
               AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc_launch_config();
    clean_bc(bc_f);
#endif
  }

  void ca_rotxfill(Real* rot, const int* rot_lo, const int* rot_hi,
                   const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                   const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = rot_lo[i];
      hi[i] = rot_hi[i];
    }

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bc, 1);
    set_bc_launch_config();
#else
    const int* bc_f = bc;
#endif

#pragma gpu
    rotxfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
             rot, AMREX_INT_ANYD(rot_lo), AMREX_INT_ANYD(rot_hi),
             AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
             AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc_launch_config();
    clean_bc(bc_f);
#endif
  }

  void ca_rotyfill(Real* rot, const int* rot_lo, const int* rot_hi,
                   const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                   const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = rot_lo[i];
      hi[i] = rot_hi[i];
    }

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bc, 1);
    set_bc_launch_config();
#else
    const int* bc_f = bc;
#endif

#pragma gpu
    rotyfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
             rot, AMREX_INT_ANYD(rot_lo), AMREX_INT_ANYD(rot_hi),
             AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
             AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc_launch_config();
    clean_bc(bc_f);
#endif
  }

  void ca_rotzfill(Real* rot, const int* rot_lo, const int* rot_hi,
                   const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                   const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = rot_lo[i];
      hi[i] = rot_hi[i];
    }

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bc, 1);
    set_bc_launch_config();
#else
    const int* bc_f = bc;
#endif

#pragma gpu
    rotzfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
             rot, AMREX_INT_ANYD(rot_lo), AMREX_INT_ANYD(rot_hi),
             AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
             AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc_launch_config();
    clean_bc(bc_f);
#endif
  }
#endif

#ifdef REACTIONS
  void ca_reactfill(Real* react, const int* react_lo, const int* react_hi,
                    const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                    const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = react_lo[i];
      hi[i] = react_hi[i];
    }

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bc, 1);
    set_bc_launch_config();
#else
    const int* bc_f = bc;
#endif

#pragma gpu
    reactfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
              react, AMREX_INT_ANYD(react_lo), AMREX_INT_ANYD(react_hi),
              AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
              AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc_launch_config();
    clean_bc(bc_f);
#endif
  }
#endif

#ifdef RADIATION
  void ca_radfill(Real* rad, const int* rad_lo, const int* rad_hi,
                  const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                  const Real* time, const int* bc)
  {
    int lo[3] = {0};
    int hi[3] = {0};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      lo[i] = rad_lo[i];
      hi[i] = rad_hi[i];
    }

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bc, 1);
    set_bc_launch_config();
#else
    const int* bc_f = bc;
#endif

#pragma gpu
    radfill(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
            rad, AMREX_INT_ANYD(rad_lo), AMREX_INT_ANYD(rad_hi),
            AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
            AMREX_REAL_ANYD(dx), AMREX_REAL_ANYD(xlo), *time, bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc_launch_config();
    clean_bc(bc_f);
#endif
  }
#endif

#ifdef __cplusplus
}
#endif
