#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_INLINE
void problem_initialize ()
{

    eos_t eos_state;

    eos_state.rho = problem::rho_0;
    eos_state.T = problem::T_0;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n]  = 0.0_rt;
    }
    eos_state.xn[0] = 1.0_rt;

    eos(eos_input_rt, eos_state);

    problem::rhoe_0 = problem::rho_0 * eos_state.e;
}
#endif

