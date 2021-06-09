#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k, Array4<Real> const& state, const GeometryData& geomdata)
{
    int coord_type = geomdata.Coord();

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();
    const Real* probhi = geomdata.ProbHi();

    // set the composition

    Real xn[NumSpec] = {0.0_rt};
    xn[0] = 1.0_rt;

    Real xc = 0.5_rt * (probhi[0] - problo[0]);
    Real yc = 0.5_rt * (probhi[1] - problo[1]);

    Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

    Real y = 0.0;
#if AMREX_SPACEDIM >= 2
    y = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
#endif

    Real z = 0.0;
#if AMREX_SPACEDIM == 3
    z = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
#endif

    state(i,j,k,URHO) = 1.0_rt;

    state(i,j,k,UMX) = -problem::vel_amp * (y - yc);
    state(i,j,k,UMY) =  problem::vel_amp * (x - xc);
    state(i,j,k,UMZ) = 0.0_rt;

    // compute the internal energy and temperature
    eos_t eos_state;

    eos_state.p = 1.0_rt;
    eos_state.rho = state(i,j,k,URHO);
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = xn[n];
    }

    eos(eos_input_rp, eos_state);

    state(i,j,k,UEDEN) = state(i,j,k,URHO) * eos_state.e +
        0.5_rt * (state(i,j,k,UMX) * state(i,j,k,UMX) +
                  state(i,j,k,UMY) * state(i,j,k,UMY) +
                  state(i,j,k,UMZ) * state(i,j,k,UMZ)) / state(i,j,k,URHO);

    state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state.e;

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * xn[n];
    }
}
#endif
