#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>
#include <network.H>
#include <model_parser.H>
#include <ambient.H>

AMREX_INLINE
void problem_initialize ()
{

    const Geometry& dgeom = DefaultGeometry();

    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    // set center variable

#if AMREX_SPACEDIM == 2
    // axisymmetric
    problem::center[0] = 0.0_rt;
    problem::center[1] = 0.5_rt * (problo[1] + probhi[1]);
    problem::center[2] = 0.0_rt;

#elif AMREX_SPACEDIM == 3
    problem::center[0] = 0.5_rt * (problo[0] + probhi[0]);
    problem::center[1] = 0.5_rt * (problo[1] + probhi[1]);
    problem::center[2] = 0.5_rt * (problo[2] + probhi[2]);
#endif


    problem::ihe4 = network_spec_index("helium-4");

    if (problem::ihe4 < 0) {
        amrex::Error("Error: no He4 in the network");
    }

    // Read initial model

    read_model_file(problem::model_name);

    // find the distance (starting from the center) where the He layer begins

    problem::R_He_base = 0.0;

    for (int n = 0; n < model::npts; n++) {
        if (model::profile(0).state(n, model::ispec+problem::ihe4) > 0.5_rt) {
            problem::R_He_base = model::profile(0).r(n);
            break;
        }
    }

    amrex::Print() << "base of He layer found at " << problem::R_He_base << std::endl;
}

#endif
