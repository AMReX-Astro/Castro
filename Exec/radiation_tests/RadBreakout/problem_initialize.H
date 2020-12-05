#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>
#include <model_parser.H>
#include <ambient.H>

AMREX_INLINE
void problem_initialize ()
{

    // Read initial model

    read_model_file(problem::model_name);

    if (NumAux != 2) {
        amrex::Error("NumAux in network is not equal to 2");
    }

    if (NumSpec != 1) {
        amrex::Error("NumSpec in network is not equal to 1");
    }

}
#endif
