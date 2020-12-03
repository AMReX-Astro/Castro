#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>
#include <network.H>
#include <actual_network.H>

AMREX_INLINE
void problem_initialize ()
{
    const Geometry& dgeom = DefaultGeometry();

    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        problem::center[n] = 0.5_rt * (problo[n] + probhi[n]);
    }

    Real xn[NumSpec] = {small_x};
    xn[0] = 1.0_rt - (NumSpec-1)*small_x;

    eos_t eos_state;

    eos_state.rho = problem::rho0;
    eos_state.T = problem::T0;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = xn[n];
    }

#ifdef NSE_THERMO
    Real aux[NumAux] = {0.0_rt};

    for (int n = 0; n < NumSpec; n++) {
        aux[AuxZero::iye] += xn[n] * zion[n] * aion_inv[n];
        aux[AuxZero::iabar] += xn[n] * aion_inv[n];
        aux[AuxZero::ibea] += xn[n] * aprox19::bion(n+1) * aion_inv[n];
    }
    aux[AuxZero::iabar] = 1.0_rt/aux[AuxZero::iabar];

    for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = aux[n];
    }
#endif

    eos(eos_input_rt, eos_state);

    problem::p0 = eos_state.p;
    problem::s0 = eos_state.s;
}

#endif
