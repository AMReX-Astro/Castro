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

    // get T_0 corresponding to rhoe_0 and rho_0 through the EOS
    Real e_0 = problem::rhoe_0 / problem::rho_0;

    eos_t eos_state;
    eos_state.rho = problem::rho_0;
    eos_state.e = e_0;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = 0.0_rt;
    }
    eos_state.xn[0] = 1.0_rt;

    eos(eos_input_re, eos_state);

    problem::T_0 = eos_state.T;

    // the 'center' variables are the location of the middle of the
    // domain -- this is where we put the interface

    problem::center[0] = 0.5_rt * (problo[0] + probhi[0]);
}
#endif
