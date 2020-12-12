#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>
#include <model_parser.H>
#include <ambient.H>

AMREX_INLINE
void problem_initialize ()
{

    const Geometry& dgeom = DefaultGeometry();

    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    if (problem::num_vortices > problem::max_num_vortices) {
        amrex::Error("num_vortices too large, please increase max_num_vortices and the size of xloc_vortices");
    }

    // Read initial model

    read_model_file(problem::model_name);

    // HSE check

    Real max_hse_err = -1.e30_rt;
    Real dr_model = model::profile(0).r(1) - model::profile(0).r(0);

    for (int i = 1; i < model::npts-2; i++) {
        Real dpdr = (model::profile(0).state(i, model::ipres) -
                     model::profile(0).state(i-1, model::ipres)) / dr_model;
        Real rhog = 0.5_rt * (model::profile(0).state(i, model::idens) +
                              model::profile(0).state(i-1, model::idens)) * gravity::const_grav;
        Real hse_err = std::abs(dpdr - rhog);

        // only count the error if we have a non-zero gradient in the next cell

        if (model::profile(0).state(i+1, model::ipres) /=
            model::profile(0).state(i, model::ipres)) {
            max_hse_err = amrex::max(max_hse_err, hse_err);
        }
    }

    std::cout << "maximum hse error = " << max_hse_err << std::endl;

    // set local variable defaults

    for (int d = 0; d < AMREX_SPACEDIM; d++) {
        problem::center[d] = 0.5_rt * (problo[d] + probhi[d]);
    }

    // velocity perturbation stuff

    Real offset = (probhi[0] - problo[0]) / problem::num_vortices;

    for (int i = 0; i < problem::num_vortices; i++) {
        problem::xloc_vortices[i] = (static_cast<Real>(i) + 0.5_rt) * offset + problo[0];
    }

}
#endif
