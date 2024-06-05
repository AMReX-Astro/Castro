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

    // semi-major axis
    Real a = problem::a_T * (probhi[0] - problo[0]);

    // semi-minor axis
    Real b = a * std::sqrt(1.0 - problem::ecc_T * problem::ecc_T);

    // compute distance from the center
    Real xx = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt) - problem::center[0];
    Real yy = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt) - problem::center[1];

    Real r = std::sqrt((xx / a) * (xx / a) +
                       (yy / b) * (yy / b));

    state(i,j,k,URHO) = problem::dens;

    Real sigma = 1.0_rt / (1.0_rt + std::exp(-(1.0_rt - r)/ problem::w_T));

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

    eos(eos_input_rt, eos_state);

    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;
    state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state.e;
    state(i,j,k,UEDEN) = state(i,j,k,UEINT);

}

#endif
