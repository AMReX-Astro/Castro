#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>
#include <network.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();
    const Real* probhi = geomdata.ProbHi();

    Real width = problem::w_T * (probhi[0] - problo[0]);
    Real c_T = problo[0] + problem::center_T * (probhi[0] - problo[0]);

    Real xcen = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

    state(i,j,k,URHO) = problem::dens;

    Real sigma = 1.0_rt / (1.0_rt + std::exp(-(c_T - xcen)/ width));
    Real freq = problem::frequency * 2.0 * 3.1415926535897932384;

    state(i,j,k,UTEMP) = problem::T_l + (problem::T_r - problem::T_l) * (1.0_rt - sigma);

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * (problem::xn_env[n] * std::sin(freq * i) + problem::xn_core[n] * (1.0 - std::sin(freq * i)));
    }

    eos_t eos_state;
    eos_state.rho = state(i,j,k,URHO);
    eos_state.T = state(i,j,k,UTEMP);
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = problem::xn_core[n] + problem::xn_env[n];
    }

    // normalize eos_state.xn
    Real sum_xn = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
        sum_xn += eos_state.xn[n]
    }

    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] /= sum_xn;
    }

    eos(eos_input_rt, eos_state);

    state(i,j,k,UMX) = state(i,j,k,URHO) * (problem::vel - 2 * problem::vel * (1.0_rt - sigma));
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;
    state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state.e;
    state(i,j,k,UEDEN) = state(i,j,k,UEINT) + 0.5_rt * state(i,j,k,URHO) * problem::vel * problem::vel;

}

#endif
