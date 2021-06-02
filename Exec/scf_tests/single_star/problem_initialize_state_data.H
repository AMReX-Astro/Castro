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

    Real xx = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

    Real yy = 0.0;
#if AMREX_SPACEDIM >= 2
    yy = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
#endif

    Real zz = 0.0;
#if AMREX_SPACEDIM == 3
    zz = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
#endif

    // Establish the initial guess: a uniform density sphere

    if (std::sqrt(xx * xx + yy * yy + zz * zz) < castro::scf_equatorial_radius) {
        state(i,j,k,URHO) = castro::scf_maximum_density;
    } else {
        state(i,j,k,URHO) = problem::ambient_dens;
    }

    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

    state(i,j,k,UTEMP) = 1.0e7_rt;
    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) / static_cast<Real>(NumSpec);
    }

    eos_t eos_state;
    eos_state.rho = state(i,j,k,URHO);
    eos_state.T = state(i,j,k,UTEMP);
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state(i,j,k,UFS+n) / state(i,j,k,URHO);
    }

    eos(eos_input_rt, eos_state);

    state(i,j,k,UEINT) = eos_state.e * state(i,j,k,URHO);
    state(i,j,k,UEDEN) = state(i,j,k,UEINT) +
        0.5_rt * (state(i,j,k,UMX) * state(i,j,k,UMX) +
                  state(i,j,k,UMY) * state(i,j,k,UMY) +
                  state(i,j,k,UMZ) * state(i,j,k,UMZ)) / state(i,j,k,URHO);

}
#endif
