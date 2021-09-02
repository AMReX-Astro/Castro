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

    problem::L_x = probhi[0] - problo[0];

}
#endif
