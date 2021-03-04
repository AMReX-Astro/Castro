#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <Castro_util.H>
#include <eos.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k, Array4<Real> const& state, const GeometryData& geomdata)
{
    GpuArray<Real, 3> loc;
    position(i, j, k, geomdata, loc);

    // Uniform density, negligible pressure.

    eos_t zone_state;

    zone_state.rho = 1.0e0_rt;
    zone_state.p   = 1.0e-6_rt;
    for (int n = 0; n < NumSpec; ++n) {
        zone_state.xn[n] = 1.0_rt / static_cast<Real>(NumSpec);
    }

    eos(eos_input_rp, zone_state);

    // Radial inflow with |v| = 1.

    Real r = std::sqrt(loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2]);

    Real vel[3] = {-loc[0] / r, -loc[1] / r, -loc[2] / r};

    state(i,j,k,URHO)  = zone_state.rho;
    state(i,j,k,UTEMP) = zone_state.T;
    state(i,j,k,UEINT) = zone_state.e * zone_state.rho;
    for (int n = 0; n < NumSpec; ++n) {
        state(i,j,k,UFS+n) = zone_state.xn[n] * zone_state.rho;
    }

    state(i,j,k,UMX) = state(i,j,k,URHO) * vel[0];
    state(i,j,k,UMY) = state(i,j,k,URHO) * vel[1];
    state(i,j,k,UMZ) = state(i,j,k,URHO) * vel[2];

    state(i,j,k,UEDEN) = state(i,j,k,UEINT) + (state(i,j,k,UMX) * state(i,j,k,UMX) +
                                               state(i,j,k,UMY) * state(i,j,k,UMY) +
                                               state(i,j,k,UMZ) * state(i,j,k,UMZ)) / (2.0_rt * state(i,j,k,URHO));
}

#endif
