/* Implementations of functions in Problem.H go here */

#include <Castro.H>

using namespace amrex;

Real Castro::initial_moment_of_inertia = 0.0;

Real
Castro::compute_moment_of_inertia()
{

    // Compute the moment of inertia for the mass
    // distribution (mass * r**2). This gives us a
    // quality measure of how well the torus holds
    // together, by comparing the moment of inertia
    // at t == 0 to the current moment of inertia.

    int finest_level = parent->finestLevel();
    Real time = state[State_Type].curTime();
    Real dt = parent->dtLevel(0);

    Real moment = 0.0;

    bool local_flag = true;

    for (int lev = 0; lev <= finest_level; ++lev) {

        // Get the Castro level
        Castro& ca_lev = getLevel(lev);

        // Add up the moment of inertia on this level
        int idir = -1; // So we do r**2, rather than any particular direction
        moment += ca_lev.locSquaredSum("density", time, idir, local_flag);

    }

    // Reduce over all ranks.
    ParallelDescriptor::ReduceRealSum(moment);

    return moment;

}

void
Castro::problem_post_init() {

    if (level != 0) return;

    initial_moment_of_inertia = compute_moment_of_inertia();

    amrex::Print() << std::endl << "  Moment of inertia = " << initial_moment_of_inertia << std::endl;

}

void
Castro::problem_post_restart() {

    if (level != 0) return;

    initial_moment_of_inertia = compute_moment_of_inertia();

    amrex::Print() << std::endl << "  Moment of inertia = " << initial_moment_of_inertia << std::endl;

}

void
Castro::problem_post_timestep() {

    if (level != 0) return;

    Real moment_of_inertia = compute_moment_of_inertia();

    amrex::Print() << std::endl << "  Moment of inertia = " << moment_of_inertia << std::endl;
    amrex::Print() << std::endl << "  Moment of inertia / initial moment of inertia = " << moment_of_inertia / initial_moment_of_inertia << std::endl;

}
