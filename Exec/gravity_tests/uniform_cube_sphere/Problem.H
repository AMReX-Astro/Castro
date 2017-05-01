// Preprocessor directive for allowing us to do a post-initialization update.

#ifndef DO_PROBLEM_POST_INIT
#define DO_PROBLEM_POST_INIT
#endif

static int problem;
static amrex::Real density;
static amrex::Real diameter;

void problem_post_init();
