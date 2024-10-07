#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>
#include <network.H>

AMREX_INLINE
void problem_initialize ()
{

    // get the species indices

    problem::ihe4 = network_spec_index("helium-4");
    problem::ic12 = network_spec_index("carbon-12");
    problem::in14 = network_spec_index("nitrogen-14");
    problem::io16 = network_spec_index("oxygen-16");

    if (problem::ihe4 < 0 || problem::ic12 < 0 || problem::io16 < 0) {
        amrex::Error("ERROR: species indices not found");
    }

    // make sure that the carbon fraction falls between 0 and 1
    // in the core and the envelope

    if (problem::cfrac_core > 1.0_rt || problem::cfrac_core < 0.0_rt
        || problem::cfrac_env > 1.0_rt || problem::cfrac_env < 0.0_rt) {
        amrex::Error("ERROR: cfrac must fall between 0 and 1");
    }

    // make sure that the nitrogen fraction falls between 0 and 1
    // in the core and the envelope

    if (problem::nfrac_core > 1.0_rt || problem::nfrac_core < 0.0_rt
        || problem::nfrac_env > 1.0_rt || problem::nfrac_env < 0.0_rt) {
        amrex::Error("ERROR: nfrac must fall between 0 and 1");
    }

    // make sure that the oxygen fraction falls between 0 and 1
    // in the core and the envelope

    if (problem::ofrac_core > 1.0_rt || problem::ofrac_core < 0.0_rt
        || problem::ofrac_env > 1.0_rt || problem::ofrac_env < 0.0_rt) {
        amrex::Error("ERROR: ofrac must fall between 0 and 1");
    }

    // make sure that the C/O fraction sums to no more than 1
    // in the core and the envelope

    if (problem::cfrac_core + problem::nfrac_core + problem::ofrac_core > 1.0_rt
        || problem::cfrac_env + problem::nfrac_env + problem::ofrac_env > 1.0_rt) {
        amrex::Error("ERROR: cfrac + nfrac + ofrac cannot exceed 1.");
    }

    // set the default mass fractions

    for (int n = 0; n < NumSpec; n++) {
        problem::xn_core[n] = problem::smallx;
        problem::xn_env[n] = problem::smallx;
    }

    problem::xn_core[problem::ic12] = amrex::max(problem::cfrac_core, problem::smallx);
    problem::xn_core[problem::io16] = amrex::max(problem::ofrac_core, problem::smallx);

    problem::xn_env[problem::ic12] = amrex::max(problem::cfrac_env, problem::smallx);
    problem::xn_env[problem::io16] = amrex::max(problem::ofrac_env, problem::smallx);

    if (problem::in14 >= 0) {
       problem::xn_core[problem::in14] = amrex::max(problem::nfrac_core, problem::smallx);
       problem::xn_core[problem::ihe4] = 1.0_rt - problem::xn_core[problem::ic12]
                                           - problem::xn_core[problem::in14]
                                           - problem::xn_core[problem::io16]
                                           - (NumSpec - 4) * problem::smallx;

       problem::xn_env[problem::in14] = amrex::max(problem::nfrac_env, problem::smallx);
       problem::xn_env[problem::ihe4] = 1.0_rt - problem::xn_env[problem::ic12]
                                           - problem::xn_env[problem::in14]
                                           - problem::xn_env[problem::io16]
                                           - (NumSpec - 4) * problem::smallx;
    }
    else {
       problem::xn_core[problem::ihe4] = 1.0_rt - problem::xn_core[problem::ic12]
                                           - problem::xn_core[problem::io16]
                                           - (NumSpec - 3) * problem::smallx;

       problem::xn_env[problem::ihe4] = 1.0_rt - problem::xn_env[problem::ic12]
                                           - problem::xn_env[problem::io16]
                                           - (NumSpec - 3) * problem::smallx;
    }

    // Set the ambient material

    problem::ambient_dens = problem::dens;
    for (int n = 0; n < NumSpec; n++) {
        problem::ambient_comp[n] = problem::xn_env[n];
    }

    eos_t eos_state;
    eos_state.rho = problem::ambient_dens;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n]  = problem::ambient_comp[n];
    }

    eos_state.T = problem::T_l;

    eos(eos_input_rt, eos_state);

    problem::ambient_e_l = eos_state.e;

    eos_state.T = problem::T_r;

    eos(eos_input_rt, eos_state);

    problem::ambient_e_r = eos_state.e;

}

#endif
