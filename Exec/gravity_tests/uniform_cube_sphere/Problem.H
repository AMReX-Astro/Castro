// Preprocessor directive for allowing us to do a post-initialization update.

#ifndef DO_PROBLEM_POST_INIT
#define DO_PROBLEM_POST_INIT
#endif

static int problem;
static Real density;
static Real diameter;

void problem_post_init();
