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

}

#endif
