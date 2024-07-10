#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

    amrex::ignore_unused(geomdata);

    state(i,j,k,URHO) = problem::rho_0;
    state(i,j,k,UMX) = problem::rho_0 * problem::u_x;
    state(i,j,k,UMY) = problem::rho_0 * problem::u_y;
    state(i,j,k,UMZ) = 0.0_rt;

    state(i,j,k,UEDEN) = problem::rhoe_0 +
        0.5_rt * problem::rho_0 * (problem::u_x * problem::u_x +
                                   problem::u_y * problem::u_y);
    state(i,j,k,UEINT) = problem::rhoe_0;
    state(i,j,k,UTEMP) = problem::T_0;

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = 0.0_rt;
    }
    state(i,j,k,UFS) = state(i,j,k,URHO);
}
#endif
