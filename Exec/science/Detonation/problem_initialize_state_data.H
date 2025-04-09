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

    Real sigma_arg = -(c_T - xcen) / width;
    Real sigma;
    // need to avoid FP overflow for sigma_arg >= 709
    // 1/(1 + exp(100)) ~= 3e-44, which is much smaller than machine epsilon
    if (sigma_arg < 100) {
        sigma = 1.0_rt / (1.0_rt + std::exp(sigma_arg));
    } else {
        sigma = 0.0_rt;
    }

    state(i,j,k,UTEMP) = problem::T_l + (problem::T_r - problem::T_l) * (1.0_rt - sigma);

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * problem::xn[n];
    }

    eos_t eos_state;
    eos_state.rho = state(i,j,k,URHO);
    eos_state.T = state(i,j,k,UTEMP);
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = problem::xn[n];
    }

#ifdef AUX_THERMO
    // set the aux quantities -- we need to do this if we are using the NSE network
    set_aux_comp_from_X(eos_state);

    for (int n = 0; n < NumAux; n++) {
        state(i,j,k,UFX+n) = state(i,j,k,URHO) * eos_state.aux[n];
    }
#endif

#ifdef NSE_NET
    // set initial conditions for NSE_NET
    state(i,j,k,UMUP) = problem::mu_p;
    state(i,j,k,UMUN) = problem::mu_n;

    burn_t burn_state;
    burn_state.T = state(i,j,k,UTEMP);
    burn_state.rho = state(i,j,k,URHO);
    for (int n = 0; n < NumSpec; ++n) {
      burn_state.xn[n] = problem::xn[n];
    }
    // initialize y_e.
    burn_state.y_e = 0.5_rt;
    burn_state.mu_p = problem::mu_p;
    burn_state.mu_n = problem::mu_n;

    if (burn_state.T > 2.5e9) {
        auto nse_state = get_actual_nse_state(burn_state, 1.0e-10_rt, true);

      // for (int n = 0; n < NumSpec; ++n) {
      //   state(i,j,k,UFS+n) = state(i,j,k,URHO) * nse_state.xn[n];
      // }

      state(i,j,k,UMUP) = nse_state.mu_p;
      state(i,j,k,UMUN) = nse_state.mu_n;
    }
#endif

    eos(eos_input_rt, eos_state);

    state(i,j,k,UMX) = state(i,j,k,URHO) * (problem::vel - 2 * problem::vel * (1.0_rt - sigma));
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;
    state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state.e;
    state(i,j,k,UEDEN) = state(i,j,k,UEINT) + 0.5_rt * state(i,j,k,URHO) * problem::vel * problem::vel;

}

#endif
