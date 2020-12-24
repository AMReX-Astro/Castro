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

    Real L = probhi[0] - problo[0];
    Real x_int = problo[0] + problem::pert_frac * L;

    Real pert_width = problem::pert_delta * L;

    Real xx = problo[0] + dx[0]* (static_cast<Real>(i) + 0.5_rt);

    // blend the fuel and ash state, keeping the pressure constant

    eos_t eos_state;
    eos_state.rho = (problem::rho_ash - problem::rho_fuel) * 0.5_rt *
        (1.0_rt - std::tanh((xx - x_int) / pert_width)) + problem::rho_fuel;
    eos_state.T = (problem::T_ash - problem::T_fuel) * 0.5_rt *
        (1.0_rt - std::tanh((xx - x_int) / pert_width)) + problem::T_fuel;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = (problem::xn_ash[n] - problem::xn_fuel[n]) * 0.5_rt *
            (1.0_rt - std::tanh((xx - x_int) / pert_width)) + problem::xn_fuel[n];
    }
    eos_state.p = problem::p_fuel;

    eos(eos_input_tp, eos_state);

    state(i,j,k,URHO ) = eos_state.rho;

    state(i,j,k,UMX) = problem::mass_flux;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

    state(i,j,k,UEDEN) = eos_state.rho * eos_state.e +
        0.5_rt * (state(i,j,k,UMX) * state(i,j,k,UMX) +
                  state(i,j,k,UMY) * state(i,j,k,UMY) +
                  state(i,j,k,UMZ) * state(i,j,k,UMZ)) / state(i,j,k,URHO);
    state(i,j,k,UEINT) = eos_state.rho * eos_state.e;
    state(i,j,k,UTEMP) = eos_state.T;
    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = eos_state.rho * eos_state.xn[n];
    }
}
#endif
