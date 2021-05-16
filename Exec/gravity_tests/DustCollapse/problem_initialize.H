#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_INLINE
void problem_initialize ()
{
    problem::center[0] = problem::center_x;

#if AMREX_SPACEDIM >= 2
    problem::center[1] = problem::center_y;
#else
    problem::center[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
    problem::center[2] = problem::center_z;
#else
    problem::center[2] = 0.0_rt;
#endif

    const Geometry& dgeom = DefaultGeometry();
    const Real* problo = dgeom.ProbLo();

    if (problo[0] != 0.0_rt) {
        amrex::Error("ERROR: xmin should be 0!");
    }

#if AMREX_SPACEDIM >= 2
    if (problo[1] != 0.0_rt) {
        amrex::Error("ERROR: ymin should be 0!");
    }
#endif

#if AMREX_SPACEDIM == 3
    if (problo[2] != 0.0_rt) {
        amrex::Error("ERROR: zmin should be 0!");
    }
#endif

    // set the composition to be uniform
    for (int n = 0; n < NumSpec; ++n) {
        problem::X_0[n] = 0.0_rt;
    }
    problem::X_0[0] = 1.0_rt;

    // get the ambient temperature and sphere temperature, T_0

    eos_t eos_state;

    eos_state.rho = problem::rho_0;
    eos_state.p   = problem::p_0;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = problem::X_0[n];
    }
    eos_state.T = castro::small_temp; // Initial guess for the EOS

    eos(eos_input_rp, eos_state);

    problem::T_0 = eos_state.T;

    eos_state.rho = problem::rho_ambient;

    eos(eos_input_rp, eos_state);

    problem::T_ambient = eos_state.T;
}
#endif
