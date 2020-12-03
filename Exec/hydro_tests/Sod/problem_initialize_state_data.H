#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>

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
