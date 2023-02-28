#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>
#include <model_parser.H>
#include <ambient.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

#if AMREX_SPACEDIM == 3
    Real zcen = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt) - problem::center[2];
#else
    Real zcen = 0.0_rt;
#endif

    Real ycen = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt) - problem::center[1];

    Real xcen = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt) - problem::center[0];

#if AMREX_SPACEDIM == 2
    Real dist = std::sqrt(xcen * xcen + ycen * ycen);
#else
    Real dist = std::sqrt(xcen * xcen + ycen * ycen + zcen * zcen);
#endif

    state(i,j,k,URHO) = interpolate(dist, model::idens);
    state(i,j,k,UTEMP) = interpolate(dist, model::itemp);

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = interpolate(dist, model::ispec+n);
    }

    eos_t eos_state;
    eos_state.rho = state(i,j,k,URHO);
    eos_state.T = state(i,j,k,UTEMP);
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state(i,j,k,UFS+n);
    }

    eos(eos_input_rt, eos_state);

    state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state.e;
    state(i,j,k,UEDEN) = state(i,j,k,URHO) * eos_state.e;

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * state(i,j,k,UFS+n);
    }

    // Initial velocities = 0

    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

}
#endif
