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

    Real sigma = 1.0_rt / (1.0_rt + std::exp(-(c_T - xcen) / width));

    // define the WD-accretion interface around the center of the domain
    Real delta = problem::delta;
    Real center = problo[0] + (probhi[0] - problo[0]) / 2; // center of domain
    Real x_start = center - delta / 2.0_rt;
    Real x_end = center + delta / 2.0_rt;

    // linear interpolation factor
    Real interp_factor = (xcen - x_start) / (x_end - x_start);

    // clamp interp_factor between 0 and 1
    interp_factor = std::min(std::max(interp_factor, 0.0_rt), 1.0_rt);

    Real temperature = problem::T_l + (problem::T_r - problem::T_l) * interp_factor;

    if (xcen < problem::hotspot_size) {
        temperature += problem::T_ign;
    }

    // create a smooth transition from xn_env (left side of domain)
    // to xn_core (right side of domain)
    state(i,j,k,UTEMP) = temperature;

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * ((1.0_rt - interp_factor) * problem::xn_env[n] +
                                              interp_factor * problem::xn_core[n]);
    }

    eos_t eos_state;
    eos_state.rho = state(i,j,k,URHO);
    eos_state.T = state(i,j,k,UTEMP);
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state(i,j,k,UFS+n) / state(i,j,k,URHO);
    }

    // normalize eos_state.xn
    Real sum_xn = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
        sum_xn += eos_state.xn[n];
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
