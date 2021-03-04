#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_INLINE
void problem_initialize ()
{

    // Force a different pressure choice for problem 5

    if (problem::problem == 5) {
        problem::pressure = 10.0;
    }
}
#endif
