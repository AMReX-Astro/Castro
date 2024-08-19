#ifndef problem_bc_fill_H
#define problem_bc_fill_H

#include <Castro_util.H>
#include <eos.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_bc_fill(int i, int j, int k,
                     Array4<Real> const& state,
                     Real time,
                     const Array1D<BCRec, 0, NUM_STATE-1>& bcs,
                     const GeometryData& geomdata)
{
    const Real pres_init = 1.0e-6_rt;
    const Real rho_init = 1.0e0_rt;

    const int* domhi = geomdata.Domain().hiVect();

    // Overwrite the outer boundary conditions

    if (!(i > domhi[0] || j > domhi[1] || k > domhi[2])) {
        return;
    }

    GpuArray<Real, 3> loc;
    position(i, j, k, geomdata, loc);

    Real r = std::sqrt(loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2]);

    eos_t zone_state;
    zone_state.rho = rho_init * std::pow(1.0_rt + time / r, AMREX_SPACEDIM - 1);
    zone_state.p   = pres_init * std::pow(zone_state.rho / rho_init, 1.0_rt + eos_rp::eos_gamma);
    for (int n = 0; n < NumSpec; ++n) {
        zone_state.xn[n] = 1.0_rt / static_cast<Real>(NumSpec);
    }

    eos(eos_input_rp, zone_state);

    // Radial inflow with |v| = 1.

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
