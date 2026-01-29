#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_INLINE
void problem_initialize ()
{

    eos_t eos_state;

    // set some defaults we will need later in the BCs

    eos_state.rho = problem::rho0;
    eos_state.T = problem::T0;
    for (auto & X : eos_state.xn) {
        X = 0.0_rt;
    }
    eos_state.xn[0]  = 1.0_rt;

    eos(eos_input_rt, eos_state);

    problem::eint0 = problem::rho0 * eos_state.e;
    problem::etot0 = problem::eint0 +
        0.5_rt * problem::rho0 * (problem::v0 * problem::v0);

    eos_state.rho = problem::rho1;
    eos_state.T = problem::T1;

    eos(eos_input_rt, eos_state);

    problem::eint1 = problem::rho1 * eos_state.e;
    problem::etot1 = problem::eint1 +
        0.5_rt * problem::rho1 * (problem::v1 * problem::v1);
}
#endif
