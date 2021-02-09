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

    int coord_type = geomdata.Coord();

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();
    const Real* probhi = geomdata.ProbHi();

    Real L = probhi[0] - problo[0];
    Real x_int = problo[0] + problem::pert_frac * L;

    Real pert_width = problem::pert_delta * L;

    // fuel state

    Real xn_fuel[NumSpec];
    for (int n = 0; n < NumSpec; n++) {
        xn_fuel[n] = 0.0_rt;
    }
    xn_fuel[problem::ifuel] = 1.0_rt;

    eos_t eos_state;
    eos_state.rho = problem::rho_fuel;
    eos_state.T = problem::T_fuel;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = xn_fuel[n];
    }

    eos(eos_input_rt, eos_state);

    Real e_fuel = eos_state.e;
    Real p_fuel = eos_state.p;

    // compute the ash state

    Real rho_ash = problem::rho_fuel / (1.0_rt + 0.6_rt * specific_q_burn / e_fuel);
    Real e_ash = e_fuel - p_fuel*(1.0_rt/rho_ash - 1.0_rt / problem::rho_fuel) + specific_q_burn;
    Real xn_ash[NumSpec];
    for (int n = 0; n < NumSpec; n++) {
        xn_ash[n] = 0.0_rt;
    }
    xn_ash[problem::iash] = 1.0_rt;

    eos_state.rho = rho_ash;
    eos_state.e = e_ash;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = xn_ash[n];
    }

    eos(eos_input_re, eos_state);

    Real T_ash = eos_state.T;

    Real xx = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

    if (xx <= x_int) {

        // ash

        state(i,j,k,URHO ) = rho_ash;
        state(i,j,k,UEDEN) = rho_ash * e_ash;
        state(i,j,k,UEINT) = rho_ash * e_ash;
        state(i,j,k,UTEMP) = T_ash;
        for (int n = 0; n < NumSpec; n++) {
            state(i,j,k,UFS+n) = rho_ash * xn_ash[n];
        }

    } else if (xx > x_int && xx < x_int + pert_width) {

        // linearly interpolate

        Real f = (xx - x_int) / pert_width;
        eos_state.e = (1.0_rt - f) * e_ash + f * e_fuel;
        eos_state.rho = (1.0_rt - f) * rho_ash + f * problem::rho_fuel;
        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = (1.0_rt - f) * xn_ash[n] + f * xn_fuel[n];
        }

        eos(eos_input_re, eos_state);

        state(i,j,k,URHO ) = eos_state.rho;
        state(i,j,k,UEDEN) = eos_state.rho * eos_state.e;
        state(i,j,k,UEINT) = eos_state.rho * eos_state.e;
        state(i,j,k,UTEMP) = eos_state.T;
        for (int n = 0; n < NumSpec; n++) {
            state(i,j,k,UFS+n) = eos_state.rho * eos_state.xn[n];
        }

    } else {

        // fuel

        state(i,j,k,URHO ) = problem::rho_fuel;
        state(i,j,k,UEDEN) = problem::rho_fuel * e_fuel;
        state(i,j,k,UEINT) = problem::rho_fuel * e_fuel;
        state(i,j,k,UTEMP) = problem::T_fuel;
        for (int n = 0; n < NumSpec; n++) {
            state(i,j,k,UFS+n) = problem::rho_fuel * xn_fuel[n];
        }

    }

    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

}

#endif
