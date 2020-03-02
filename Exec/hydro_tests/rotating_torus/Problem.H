
// Preprocessor macros we need.

#ifndef DO_PROBLEM_POST_INIT
#define DO_PROBLEM_POST_INIT
#endif

#ifndef DO_PROBLEM_POST_RESTART
#define DO_PROBLEM_POST_RESTART
#endif

#ifndef DO_PROBLEM_POST_TIMESTEP
#define DO_PROBLEM_POST_TIMESTEP
#endif

// Post initialization analytics.

void problem_post_init();

// End of timestep analytics.

void problem_post_timestep();

// Post restart analytics.

void problem_post_restart();

// Sum up the moment of inertia on the grid.

amrex::Real compute_moment_of_inertia();

// Save the initial moment of inertia for comparison.

static amrex::Real initial_moment_of_inertia;
