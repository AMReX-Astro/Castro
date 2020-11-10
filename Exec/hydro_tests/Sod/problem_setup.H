#ifndef problem_setup_H
#define problem_setup_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_INLINE
void problem_initialize ()
{
    const Geometry& dgeom = DefaultGeometry();

    const int coord_type = dgeom.Coord();
    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    problem::split[0] = problem::frac * (problo[0] + probhi[0]);
#if AMREX_SPACEDIM >= 2
    problem::split[1] = problem::frac * (problo[1] + probhi[1]);
#else
    problem::split[1] = 0.0;
#endif
#if AMREX_SPACEDIM == 3
    problem::split[2] = problem::frac * (problo[2] + probhi[2]);
#else
    problem::split[2] = 0.0;
#endif

    // compute the internal energy (erg/cc) for the left and right state
    Real xn[NumSpec] = {0.0};
    xn[0] = 1.0_rt;

    eos_t eos_state;

    if (problem::use_Tinit) {

        eos_state.rho = problem::rho_l;
        eos_state.T = problem::T_l;
        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = xn[n];
        }

        eos(eos_input_rt, eos_state);

        problem::rhoe_l = problem::rho_l * eos_state.e;
        problem::p_l = eos_state.p;

        eos_state.rho = problem::rho_r;
        eos_state.T = problem::T_r;
        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = xn[n];
        }

        eos(eos_input_rt, eos_state);

        problem::rhoe_r = problem::rho_r * eos_state.e;
        problem::p_r = eos_state.p;

    } else {

        eos_state.rho = problem::rho_l;
        eos_state.p = problem::p_l;
        eos_state.T = 100000.0_rt;  // initial guess
        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = xn[n];
        }

        eos(eos_input_rp, eos_state);

        problem::rhoe_l = problem::rho_l * eos_state.e;
        problem::T_l = eos_state.T;

        eos_state.rho = problem::rho_r;
        eos_state.p = problem::p_r;
        eos_state.T = 100000.0_rt;  // initial guess
        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = xn[n];
        }

        eos(eos_input_rp, eos_state);

        problem::rhoe_r = problem::rho_r * eos_state.e;
        problem::T_r = eos_state.T;

    }

}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{
    int coord_type = geomdata.Coord();

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

    Real y = 0.0;
#if AMREX_SPACEDIM >= 2
    y = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
#endif

    Real z = 0.0;
#if AMREX_SPACEDIM == 3
    z = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
#endif

    if (problem::idir == 1) {
        if (x <= problem::split[0]) {
            state(i,j,k,URHO) = problem::rho_l;
            state(i,j,k,UMX) = problem::rho_l * problem::u_l;
            state(i,j,k,UMY) = 0.0_rt;
            state(i,j,k,UMZ) = 0.0_rt;
            state(i,j,k,UEDEN) = problem::rhoe_l +
                0.5 * problem::rho_l * problem::u_l * problem::u_l;
            state(i,j,k,UEINT) = problem::rhoe_l;
            state(i,j,k,UTEMP) = problem::T_l;
        } else {
            state(i,j,k,URHO) = problem::rho_r;
            state(i,j,k,UMX) = problem::rho_r * problem::u_r;
            state(i,j,k,UMY) = 0.0_rt;
            state(i,j,k,UMZ) = 0.0_rt;
            state(i,j,k,UEDEN) = problem::rhoe_r +
                0.5 * problem::rho_r * problem::u_r * problem::u_r;
            state(i,j,k,UEINT) = problem::rhoe_r;
            state(i,j,k,UTEMP) = problem::T_r;
        }

    } else if (problem::idir == 2) {
        if (y <= problem::split[1]) {
            state(i,j,k,URHO) = problem::rho_l;
            state(i,j,k,UMX) = 0.0_rt;
            state(i,j,k,UMY) = problem::rho_l * problem::u_l;
            state(i,j,k,UMZ) = 0.0_rt;
            state(i,j,k,UEDEN) = problem::rhoe_l +
                0.5 * problem::rho_l * problem::u_l * problem::u_l;
            state(i,j,k,UEINT) = problem::rhoe_l;
            state(i,j,k,UTEMP) = problem::T_l;
        } else {
            state(i,j,k,URHO) = problem::rho_r;
            state(i,j,k,UMX) = 0.0_rt;
            state(i,j,k,UMY) = problem::rho_r * problem::u_r;
            state(i,j,k,UMZ) = 0.0_rt;
            state(i,j,k,UEDEN) = problem::rhoe_r +
                0.5 * problem::rho_r * problem::u_r * problem::u_r;
            state(i,j,k,UEINT) = problem::rhoe_r;
            state(i,j,k,UTEMP) = problem::T_r;
        }

    } else if (problem::idir == 3) {
        if (z <= problem::split[2]) {
            state(i,j,k,URHO) = problem::rho_l;
            state(i,j,k,UMX) = 0.0_rt;
            state(i,j,k,UMY) = 0.0_rt;
            state(i,j,k,UMZ) = problem::rho_l * problem::u_l;
            state(i,j,k,UEDEN) = problem::rhoe_l +
                0.5 * problem::rho_l * problem::u_l * problem::u_l;
            state(i,j,k,UEINT) = problem::rhoe_l;
            state(i,j,k,UTEMP) = problem::T_l;
        } else {
            state(i,j,k,URHO) = problem::rho_r;
            state(i,j,k,UMX) = 0.0_rt;
            state(i,j,k,UMY) = 0.0_rt;
            state(i,j,k,UMZ) = problem::rho_r * problem::u_r;
            state(i,j,k,UEDEN) = problem::rhoe_r +
                0.5 * problem::rho_r * problem::u_r * problem::u_r;
            state(i,j,k,UEINT) = problem::rhoe_r;
            state(i,j,k,UTEMP) = problem::T_r;
        }

    } else {
#ifndef AMREX_USE_GPU
        amrex::Error("invalid idir");
#endif
    }

    state(i,j,k,UFS) = state(i,j,k,URHO);
    for (int n = 1; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = 0.0_rt;
    }

}

#endif
