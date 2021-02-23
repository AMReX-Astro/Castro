#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_INLINE
void problem_initialize ()
{
    const Geometry& dgeom = DefaultGeometry();

    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        problem::center[n] = 0.5_rt * (problo[n] + probhi[n]);
    }

    Real xn[NumSpec] = {0.0};

    xn[0] = 1.0_rt;

    eos_t eos_state;

    eos_state.rho = problem::rho0;
    eos_state.T = problem::T0;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = xn[n];
    }

    eos(eos_input_rt, eos_state);

    problem::p0 = eos_state.p;
    problem::s0 = eos_state.s;
}

#endif
