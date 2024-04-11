#ifndef PROBLEM_H
#define PROBLEM_H

/* problem-specific Castro:: declarations go here */

#ifdef DO_PROBLEM_POST_TIMESTEP
void problem_post_timestep();
#endif

#ifdef DO_PROBLEM_POST_RESTART
void problem_post_restart();
#endif

#ifdef DO_PROBLEM_POST_INIT
void problem_post_init();
#endif

#endif
