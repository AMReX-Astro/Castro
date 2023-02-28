#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_INLINE
void problem_initialize ()
{
    const Geometry& dgeom = DefaultGeometry();

    const int coord_type = dgeom.Coord();
    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    problem::split[0] = problem::frac * (problo[0] + probhi[0]);
#if AMREX_SPACEDIM >= 2
    problem::split[1] = problem::frac * (problo[1] + probhi[1]);
#else
    problem::split[1] = 0.0;
#endif
#if AMREX_SPACEDIM == 3
    problem::split[2] = problem::frac * (problo[2] + probhi[2]);
#else
    problem::split[2] = 0.0;
#endif

    // compute the internal energy (erg/cc) for the left and right state

    Real xn[NumSpec] = {0.0};
    xn[0] = 1.0_rt;

    eos_t eos_state;
    eos_state.rho = problem::rho_l;
    eos_state.p = problem::p_l;
    eos_state.T = 100000.0_rt;  // initial guess
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = xn[n];
    }

    eos(eos_input_rp, eos_state);

    problem::rhoe_l = problem::rho_l * eos_state.e;
    problem::T_l = eos_state.T;

    eos_state.rho = problem::rho_r;
    eos_state.p = problem::p_r;
    eos_state.T = 100000.0_rt;  // initial guess

    eos(eos_input_rp, eos_state);

    problem::rhoe_r = problem::rho_r * eos_state.e;
    problem::T_r = eos_state.T;

    // we overwrite the values of the B components from the probin file here

    problem::B_x_l = 4.0_rt / (2.0 * std::sqrt(M_PI));
    problem::B_y_l = 3.6_rt / (2.0 * std::sqrt(M_PI));
    problem::B_z_l = 2.0_rt / (2.0 * std::sqrt(M_PI));
    problem::B_x_r = 4.0_rt / (2.0 * std::sqrt(M_PI));
    problem::B_y_r = 4.0_rt / (2.0 * std::sqrt(M_PI));
    problem::B_z_r = 2.0_rt / (2.0 * std::sqrt(M_PI));

}
#endif
