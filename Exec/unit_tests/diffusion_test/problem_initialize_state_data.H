#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>
#include <prob_util.H>
#include <conductivity.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();
    int coord_type = geomdata.Coord();

    Real r[3] = {0.0_rt};

    r[0] = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt) - problem::center[0];

#if AMREX_SPACEDIM >= 2
    r[1] = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt) - problem::center[1];
#endif

#if AMREX_SPACEDIM == 3
    r[2] = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt) - problem::center[2];
#endif

    // set the composition

    Real X[NumSpec] = {0.0_rt};
    X[0] = 1.0_rt;

    state(i,j,k,URHO) = problem::rho0;

    Real temp = analytic(r, 0.0_rt, coord_type);

    state(i,j,k,UTEMP) = temp;

    // compute the internal energy and temperature

    eos_t eos_state;
    eos_state.T = temp;
    eos_state.rho = state(i,j,k,URHO);
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = X[n];
    }

    eos(eos_input_rt, eos_state);

    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

    state(i,j,k,UEDEN) = problem::rho0 * eos_state.e;
    state(i,j,k,UEINT) = problem::rho0 * eos_state.e;

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = problem::rho0 * X[n];
    }

}

#endif
