#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <model_parser.H>
#include <ambient.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{
    GpuArray<Real, 3> loc;
    position(i, j, k, geomdata, loc);

    Real pos[3];
    for (int n = 0; n < 3; ++n) {
        pos[n] = loc[n] - problem::center[n];
    }

    Real dist = std::sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);

    const Real* dx = geomdata.CellSize();

    eos_t zone_state;

    if (dist < problem::wd_radius) {
        int nsub = 1;
        zone_state.rho = interpolate_3d(pos, dx, model::idens, nsub, 0);
        zone_state.T   = interpolate_3d(pos, dx, model::itemp, nsub, 0);
        for (int n = 0; n < NumSpec; ++n) {
            zone_state.xn[n] = interpolate_3d(pos, dx, model::ispec + n, nsub, 0);
        }

        eos(eos_input_rt, zone_state);
    }
    else {
        zone_state.rho = ambient::ambient_state[URHO];
        zone_state.T   = ambient::ambient_state[UTEMP];
        zone_state.e   = ambient::ambient_state[UEINT] / ambient::ambient_state[URHO];
        for (int n = 0; n < NumSpec; ++n) {
            zone_state.xn[n] = ambient::ambient_state[UFS+n] / ambient::ambient_state[URHO];
        }
    }

    state(i,j,k,URHO)  = zone_state.rho;
    state(i,j,k,UMX)   = zone_state.rho * problem::wd_speed;
    state(i,j,k,UMY)   = 0.0_rt;
    state(i,j,k,UMZ)   = 0.0_rt;
    state(i,j,k,UTEMP) = zone_state.T;
    state(i,j,k,UEINT) = zone_state.e * zone_state.rho;
    state(i,j,k,UEDEN) = zone_state.e * zone_state.rho +
                         0.5_rt * state(i,j,k,URHO) * problem::wd_speed * problem::wd_speed;
    for (int n = 0; n < NumSpec; ++n) {
        state(i,j,k,UFS+n) = zone_state.rho * zone_state.xn[n];
    }
}

#endif
