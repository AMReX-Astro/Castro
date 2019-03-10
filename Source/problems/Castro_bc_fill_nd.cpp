
#include <AMReX_BLFort.H>
#include <Castro_bc_fill_nd.H>
#include <Castro_bc_fill_nd_F.H>

using namespace amrex;

#ifdef __cplusplus
extern "C"
{
#endif

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

        hypfill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, dx, xlo, time, bc);
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

        denfill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, dx, xlo, time, bc);
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

        phigravfill(lo, hi, phi, phi_lo, phi_hi, domlo, domhi, dx, xlo, time, bc);
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

        gravxfill(lo, hi, grav, grav_lo, grav_hi, domlo, domhi, dx, xlo, time, bc);
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

        gravyfill(lo, hi, grav, grav_lo, grav_hi, domlo, domhi, dx, xlo, time, bc);
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

        gravzfill(lo, hi, grav, grav_lo, grav_hi, domlo, domhi, dx, xlo, time, bc);
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

        phirotfill(lo, hi, phi, phi_lo, phi_hi, domlo, domhi, dx, xlo, time, bc);
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

        rotxfill(lo, hi, rot, rot_lo, rot_hi, domlo, domhi, dx, xlo, time, bc);
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

        rotyfill(lo, hi, rot, rot_lo, rot_hi, domlo, domhi, dx, xlo, time, bc);
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

        rotzfill(lo, hi, rot, rot_lo, rot_hi, domlo, domhi, dx, xlo, time, bc);
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

        reactfill(lo, hi, react, react_lo, react_hi, domlo, domhi, dx, xlo, time, bc);
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

        radfill(lo, hi, rad, rad_lo, rad_hi, domlo, domhi, dx, xlo, time, bc);
    }
#endif

#ifdef __cplusplus
}
#endif
