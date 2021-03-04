#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>
#include <model_parser.H>

AMREX_INLINE
void problem_initialize ()
{

    const Geometry& dgeom = DefaultGeometry();

    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();


    // Read initial model

    read_model_file(problem::model_name);


#if AMREX_SPACEDIM == 1
    problem::center[0] = 0.0_rt;

#elif AMREX_SPACEDIM == 2
    //  assume axisymmetric
    problem::center[0] = 0.0_rt;
    problem::center[1] = 0.5_rt * (problo[1] + probhi[1]);

#else
    problem::center[0] = 0.5_rt * (problo[0] + probhi[0]);
    problem::center[1] = 0.5_rt * (problo[1] + probhi[1]);
    problem::center[2] = 0.5_rt * (problo[2] + probhi[2]);
#endif

}

#endif
