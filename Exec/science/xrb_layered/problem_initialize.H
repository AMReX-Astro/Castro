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

    // Read initial model
    read_model_file(problem::model_name);

    if (ParallelDescriptor::IOProcessor()) {
        for (int i = 0; i < model::npts; i++) {
            std::cout << i << " " << model::profile(0).r(i) << " " << model::profile(0).state(i,model::idens) << std::endl;
        }
    }

    // Set up Castro data logs for this problem

    if (castro::sum_interval > 0 && amrex::ParallelDescriptor::IOProcessor()) {

        Castro::problem_data_logs.resize(1);

        Castro::problem_data_logs[0].reset(new std::fstream);
        Castro::problem_data_logs[0]->open("xrb_diag.out", std::ios::out | std::ios::app);
        if (!Castro::problem_data_logs[0]->good()) {
            amrex::FileOpenFailed("xrb_diag.out");
        }

    }

    // set the ambient state for the upper boundary condition

    ambient::ambient_state[URHO] = model::profile(0).state(model::npts-1, model::idens);
    ambient::ambient_state[UTEMP] = model::profile(0).state(model::npts-1, model::itemp);
    for (int n = 0; n < NumSpec; n++) {
        ambient::ambient_state[UFS+n] =
            ambient::ambient_state[URHO] * model::profile(0).state(model::npts-1, model::ispec+n);
    }

    ambient::ambient_state[UMX] = 0.0_rt;
    ambient::ambient_state[UMY] = 0.0_rt;
    ambient::ambient_state[UMZ] = 0.0_rt;

    // make the ambient state thermodynamically consistent

    eos_t eos_state;
    eos_state.rho = ambient::ambient_state[URHO];
    eos_state.T = ambient::ambient_state[UTEMP];
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = ambient::ambient_state[UFS+n] / eos_state.rho;
    }

    eos(eos_input_rt, eos_state);

    ambient::ambient_state[UEINT] = eos_state.rho * eos_state.e;
    ambient::ambient_state[UEDEN] = eos_state.rho * eos_state.e;

}
#endif
