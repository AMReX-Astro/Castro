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

    problem::ihe4 = network_spec_index("helium-4");

    if (problem::ihe4 < 0) {
        amrex::Error("Error: no He4 in the network");
    }

    // Read initial model
    read_model_file(problem::model_name);

    Real max_hse_error = std::numeric_limits<Real>::lowest();

    for (int i = 0; i < model::npts-2; ++i) {
        Real dpdr = (model::profile(0).state(i+1, model::ipres) -
                     model::profile(0).state(i, model::ipres)) /
            (model::profile(0).r(i+1) - model::profile(0).r(i));

        Real rhog = 0.5_rt * (model::profile(0).state(i+1, model::idens) +
                              model::profile(0).state(i, model::idens)) * gravity::const_grav;

        // we'll check that dpdr in the next zone up is also negative, to ensure we are
        // not hitting the part where we stopped integrating HSE
        Real p_check = model::profile(0).state(i+2, model::ipres) -
                       model::profile(0).state(i+1, model::ipres);

        if (dpdr < 0 && p_check < 0) {
            max_hse_error = amrex::max(std::abs(dpdr - rhog) / std::abs(dpdr),
                                       max_hse_error);
        }
    }

    amrex::Print() << Font::Bold << FGColor::Yellow
                   << "HSE error in initial model = " << max_hse_error
                   << ResetDisplay << std::endl;

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
