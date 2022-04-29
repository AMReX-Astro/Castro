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

    if (problem::num_vortices > problem::max_num_vortices) {
        amrex::Error("num_vortices too large, please increase max_num_vortices and the size of xloc_vortices");
    }

    // Read initial model
    read_model_file(problem::model_name);

    if (ParallelDescriptor::IOProcessor()) {
        for (int i = 0; i < model::npts; i++) {
            std::cout << i << " " << model::profile(0).r(i) << " " << model::profile(0).state(i,model::idens) << std::endl;
        }
    }

    // velocity perturbation x-axis vortex location

    Real offset = (probhi[0] - problo[0]) / problem::num_vortices;

    for (int i = 0; i < problem::num_vortices; i++) {
        problem::xloc_vortices[i] = (static_cast<Real>(i) + 0.5_rt) * offset + problo[0];
    }

    problem::ih1 = network_spec_index("hydrogen-1");
    problem::layer_h1_base = 0.0_rt;

    for (int n = 0; n < model::npts; n++) {
        if (model::profile(0).state(n, model::ispec+problem::ih1) > 0.5_rt) {
            problem::layer_h1_base = model::profile(0).r(n);
            break;
        }
    }
}

#endif
