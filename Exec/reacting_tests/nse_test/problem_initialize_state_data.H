#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>
#include <network.H>
#if defined(NSE_TABLE)
#include <nse_table_type.H>
#include <nse_table.H>
#elif defined(NSE_NET)
#include <nse_solver.H>
#endif

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();
    const Real* probhi = geomdata.ProbHi();

    Real xx = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt) - problem::center[0];

    Real yy = 0.0;
#if AMREX_SPACEDIM >= 2
    yy = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt) - problem::center[1];
#endif

    Real zz = 0.0;
#if AMREX_SPACEDIM == 3
    zz = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt) - problem::center[2];
#endif

    constexpr Real ye0 = 0.5_rt;
#if defined(NSE_TABLE)
    constexpr Real dye = -0.05_rt;
#elif defined(NSE_NET)
    constexpr Real dye = -0.02_rt;
#endif
    Real dist = std::sqrt(xx * xx + yy * yy + zz * zz);

    Real T;
    Real ye;

    if (dist <= problem::center[0]) {
        T = problem::T0 *
            (1.0_rt + problem::dT_fact *
             std::exp(-(dist/problem::L_pert) * (dist/problem::L_pert)) *
             std::pow(std::cos(M_PI*(dist/(probhi[0]-problo[0]))), 6));
        ye = ye0 *
             (1.0_rt + dye * std::exp(-(dist/problem::L_pert) * (dist/problem::L_pert) ) *
              std::pow(std::cos(M_PI*(dist/(probhi[0]-problo[0]))), 6));
    } else {
        T = problem::T0;
        ye = ye0;
    }

    state(i,j,k,UMX) = problem::rho0 * problem::u0;
    state(i,j,k,UMY) = problem::rho0 * problem::v0;
    state(i,j,k,UMZ) = problem::rho0 * problem::w0;

    burn_t burn_state;
    burn_state.T = T;
    burn_state.rho = problem::rho0;
    burn_state.y_e = ye;

#if defined(NSE_TABLE)
    nse_table_t nse_state;
    nse_state.T = T;
    nse_state.rho = problem::rho0;
    nse_state.Ye = ye;
    nse_interp(nse_state);

    for (int n = 0; n < NumSpec; ++n) {
        burn_state.xn[n] = nse_state.X[n];
    }
#elif defined(NSE_NET)
    Real eps = 1.e-10_rt;
    bool input_ye_is_valid = true;

    // manually set the tolerance. Set a really large number to assure we're in NSE
    // since here we only test nse_burn() in this unit test.
    ase_tol = problem::nse_tol;

    // we need an initial guess for the mass fractions -- make them all the same
    for (int n = 0; n < NumSpec; ++n) {
        burn_state.xn[n] = 1.0_rt / static_cast<Real>(NumSpec);
    }

    // we also need a guess for the initial chemical potential
    burn_state.mu_p = -3.0_rt;
    burn_state.mu_n = -12.0_rt;

    auto nse_state = get_actual_nse_state(burn_state, eps, input_ye_is_valid);

    for (int n = 0; n < NumSpec; ++n) {
        burn_state.xn[n] = nse_state.xn[n];
    }
#endif

    // since the species are interpolated, normalize them

    Real sumX = 0.0_rt;
    for (auto X : burn_state.xn) {
        sumX += X;
    }
    for (auto & X : burn_state.xn) {
        X /= sumX;
    }

#ifdef NSE_NET
    state(i,j,k,UMUP) = burn_state.mu_p;
    state(i,j,k,UMUN) = burn_state.mu_n;
#endif

#ifdef AUX_THERMO
    burn_state.aux[AuxZero::iye] = ye;
    burn_state.aux[AuxZero::iabar] = nse_state.abar;
    burn_state.aux[AuxZero::ibea] = nse_state.bea;
#endif

    eos(eos_input_rt, burn_state);

    state(i,j,k,URHO) = burn_state.rho;

    state(i,j,k,UEDEN) = burn_state.rho * burn_state.e +
        0.5_rt * burn_state.rho * (problem::u0 * problem::u0 +
                                   problem::v0 * problem::v0 +
                                   problem::w0 * problem::w0);

    state(i,j,k,UEINT) = burn_state.rho * burn_state.e;

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * burn_state.xn[n];
    }

#ifdef AUX_THERMO
    for (int n = 0; n < NumAux; n++) {
        state(i,j,k,UFX+n) = state(i,j,k,URHO) * burn_state.aux[n];
    }
#endif

    state(i,j,k,UTEMP) = burn_state.T;

}

#endif
