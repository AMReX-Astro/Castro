// This sets up the Gresho vortex problem as described in
// Miczek, Roeple, and Edelmann 2015
//
// By choosing the reference pressure, p0, we can specify the
// Mach number

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

    // characteristic scales
    problem::x_r = probhi[0] - problo[0];
    problem::q_r = 0.4_rt * M_PI * problem::x_r / problem::t_r;

    // pressure peaks at r = 1/5 from the center
    // where p = p0 + 25/2 r^2, so peak pressure is p = p0 + 1/2
    //
    // Mach number is |u_phi| / sqrt(gamma p / rho), so we get solve
    // for p0.  Note that the peak velocity is q,
    //
    // p0 = rho q^2 / (gamma M^2) - 1/2

    problem::p0 = problem::rho0 * std::pow(problem::q_r, 2.0) /
        (eos_rp::eos_gamma * std::pow(problem::M0, 2.0)) - 0.5_rt;

}

#endif
