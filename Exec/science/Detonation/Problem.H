
// Preprocessor directives for adding custom hooks.

#ifndef DO_PROBLEM_POST_INIT
#define DO_PROBLEM_POST_INIT
#endif

#ifndef DO_PROBLEM_POST_RESTART
#define DO_PROBLEM_POST_RESTART
#endif

#ifndef DO_PROBLEM_POST_TIMESTEP
#define DO_PROBLEM_POST_TIMESTEP
#endif

// Input parameters

static int use_stopping_criterion;
static amrex::Real ts_te_stopping_criterion;
static amrex::Real T_stopping_criterion;

// Problem post-initialization routine.

void problem_post_init();

// Problem post-restart routine.

void problem_post_restart();

// End of timestep routine.

void problem_post_timestep();
