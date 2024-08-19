#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <AMReX_Geometry.H>
#include <prob_parameters.H>
#include <eos.H>
#include <network.H>
#include <model_parser.H>
#include <initial_model.H>
#include <ambient.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt) - problem::center[0];

    Real z = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);

    // lateral distance and height

    Real r;
    Real height;

    // we want the flame front to be parallel to the y-axis
    r = x;
    height = z;

    // blending factor

    Real f;

    if (r < problem::x_half_max) {
        f = 1.0_rt;

    } else if (r > problem::x_half_max + problem::x_half_width) {
        f = 0.0_rt;

    } else {
        f = -(r - problem::x_half_max) / problem::x_half_width + 1.0_rt;
    }

    state(i,j,k,URHO) = f * interpolate(height, model::idens, 1) +
             (1.0_rt - f) * interpolate(height, model::idens, 0);

    state(i,j,k,UTEMP) = f * interpolate(height, model::itemp, 1) +
              (1.0_rt - f) * interpolate(height, model::itemp, 0);

    Real temppres = f * interpolate(height, model::ipres, 1) +
         (1.0_rt - f) * interpolate(height, model::ipres, 0);

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = f * interpolate(height, model::ispec+n, 1) +
                  (1.0_rt - f) * interpolate(height, model::ispec+n, 0);
    }

    eos_t eos_state;
    eos_state.rho = state(i,j,k,URHO);
    eos_state.T = state(i,j,k,UTEMP);
    eos_state.p = temppres;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state(i,j,k,UFS+n);
    }

    eos(eos_input_rp, eos_state);

    state(i,j,k,UTEMP) = eos_state.T;
    state(i,j,k,UEINT) = eos_state.rho * eos_state.e;
    state(i,j,k,UEDEN) = state(i,j,k,UEINT);

    // Initial velocities = 0

    state(i,j,k,UMX) = 0.e0_rt;
    state(i,j,k,UMY) = 0.e0_rt;
    state(i,j,k,UMZ) = 0.e0_rt;

    // convert to partial densities

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * state(i,j,k,UFS+n);
    }
}

#endif
