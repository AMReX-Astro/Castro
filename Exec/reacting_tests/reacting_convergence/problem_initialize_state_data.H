#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>
#include <network.H>
#include <network_properties.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();
    const Real* probhi = geomdata.ProbHi();

    Real xx = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

    Real yy = 0.0;
#if AMREX_SPACEDIM >= 2
    yy = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
#endif

    Real zz = 0.0;
#if AMREX_SPACEDIM == 3
    zz = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
#endif

    Real dist = 0;
    dist = std::sqrt((problem::center[0] - xx) * (problem::center[0] - xx) +
                     (problem::center[1] - yy) * (problem::center[1] - yy) +
                     (problem::center[2] - zz) * (problem::center[2] - zz));

    Real p;
    if (dist <= problem::center[0]) {
        p = problem::p0 * (1.0_rt + problem::dp_fact *
                           std::exp(-(dist/problem::L_pert)*(dist/problem::L_pert)) *
                           std::pow(std::cos(M_PI*(dist/(probhi[0]-problo[0]))), 6));
    } else {
        p = problem::p0;
    }

    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

    // we are isentropic, so find rho
    eos_t eos_state;
    eos_state.p =  p;
    eos_state.T = 1.e4_rt;  // initial guess
    eos_state.rho = problem::rho0;  // initial guess
    eos_state.s = problem::s0;

    Real xn[NumSpec] = {network_rp::small_x};
    xn[0] = 1.0_rt - (NumSpec-1) * network_rp::small_x;

    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = xn[n];
    }

#ifdef AUX_THERMO
    set_aux_comp_from_X(eos_state);
#endif

    eos(eos_input_ps, eos_state);

    state(i,j,k,URHO) = eos_state.rho;

    state(i,j,k,UEDEN) = eos_state.rho * eos_state.e;
    state(i,j,k,UEINT) = eos_state.rho * eos_state.e;

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * xn[n];
    }

#ifdef AUX_THERMO
    for (int n = 0; n < NumAux; n++) {
        state(i,j,k,UFX+n) = state(i,j,k,URHO) * eos_state.aux[n];
    }
#endif

    state(i,j,k,UTEMP) = eos_state.T;
}

#endif
