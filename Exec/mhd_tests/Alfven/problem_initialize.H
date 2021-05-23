#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_INLINE
void problem_initialize ()
{

    // give value to B
    problem::B_x  = problem::B_0 / std::sqrt(2.0_rt);
    problem::B_y  = problem::B_0 / std::sqrt(2.0_rt);
    problem::B_z  = 0.0_rt;

    // compute the internal energy (erg/cc)

    eos_t eos_state;
    eos_state.rho = problem::rho_0;
    eos_state.p = problem::p_0;
    eos_state.T = 100.0_rt;  // initial guess
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = 0.0_rt;
    }
    eos_state.xn[0] = 1.0_rt;

    eos(eos_input_rp, eos_state);

    problem::rhoe_0 = problem::rho_0 * eos_state.e;
    problem::T_0 = eos_state.T;
}
#endif
