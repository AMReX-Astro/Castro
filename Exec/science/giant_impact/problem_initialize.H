#ifndef problem_initialize_H
#define problem_initialize_H

#include <model_parser.H>
#include <initial_model.H>
#include <fundamental_constants.H>
#include <prob_parameters.H>

const Real M_earth = 5.9722e27_rt;

AMREX_INLINE
void problem_initialize ()
{
    // Set some default values for these quantities;
    // we'll update them soon.

    for (int n = 0; n < 3; ++n) {
        problem::target_initial_center[n] = problem::center[n];
        problem::impactor_initial_center[n] = problem::center[n];

        problem::target_initial_velocity[n] = 0.0_rt;
        problem::impactor_initial_velocity[n] = 0.0_rt;
    }

    // Integrate to reach the desired mass for both target and impactor.

    if (problem::target_mass <= 0.0_rt) {
        amrex::Error("Target mass must be positive");
    }
    problem::target_mass *= M_earth;

    Real target_comp[NumSpec];

    // Arbitrarily put all the composition in the first component

    for (int n = 0; n < NumSpec; ++n) {
        target_comp[n] = 0.0_rt;
    }
    target_comp[0] = 1.0_rt;

    Real target_central_density;

    generate_initial_model(problem::target_mass, problem::target_temperature,
                           target_central_density, problem::target_radius,
                           target_comp, problem::initial_model_dx,
                           0);

    if (problem::impactor_mass <= 0.0_rt) {
        amrex::Error("Impactor mass must be positive");
    }
    problem::impactor_mass *= M_earth;

    amrex::Print() << std::endl;

    amrex::Print() << "Generated initial model for target of mass " << std::setprecision(3) << problem::target_mass / M_earth
                   << " earth masses, central density " << std::setprecision(3) << std::scientific << target_central_density
                   << " g cm**-3, and radius " << std::setprecision(3) << std::scientific << problem::target_radius << " cm."
                   << std::endl << std::endl;

    Real impactor_comp[NumSpec];

    // Arbitrarily put all the composition in the first component

    for (int n = 0; n < NumSpec; ++n) {
        impactor_comp[n] = 0.0_rt;
    }
    impactor_comp[0] = 1.0_rt;

    Real impactor_central_density;

    generate_initial_model(problem::impactor_mass, problem::impactor_temperature,
                           impactor_central_density, problem::impactor_radius,
                           impactor_comp, problem::initial_model_dx,
                           1);

    amrex::Print() << "Generated initial model for impactor of mass " << std::setprecision(3) << problem::impactor_mass / M_earth
                   << " earth masses, central density " << std::setprecision(3) << std::scientific << impactor_central_density
                   << " g cm**-3, and radius " << std::setprecision(3) << std::scientific << problem::impactor_radius << " cm."
                   << std::endl << std::endl;

    // Set ambient properties
    
    for (int n = 0; n < NumSpec; ++n) {
        ambient::ambient_state[UFS+n] = ambient::ambient_state[URHO] * (target_comp[n] + impactor_comp[n]) / 2;
    }

    // We have completed (rho, T, xn) for the ambient state, so we can call the EOS.

    eos_t eos_state;
    eos_state.rho = ambient::ambient_state[URHO];
    eos_state.T   = ambient::ambient_state[UTEMP];
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = ambient::ambient_state[UFS+n] / ambient::ambient_state[URHO];
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; ++n) {
        eos_state.aux[n] = ambient::ambient_state[UFX+n] / ambient::ambient_state[URHO];
    }
#endif

    eos(eos_input_rt, eos_state);

    ambient::ambient_state[UEINT] = ambient::ambient_state[URHO] * eos_state.e;
    ambient::ambient_state[UEDEN] = ambient::ambient_state[UEINT];

    // Set up the initial distances and velocities.
    // This follows the prescription of Meier, Reinhardt and Stadel (2021) except
    // that we have the center of mass at the system located at the origin and use
    // the center of mass frame such that the total momentum is zero.

    Real r_crit = problem::target_radius + problem::impactor_radius;
    Real v_esc_targ = std::sqrt(2.0_rt * C::Gconst * problem::target_mass / r_crit);
    Real v_esc_sys = std::sqrt(2.0_rt * C::Gconst * (problem::target_mass + problem::impactor_mass) / r_crit);

    Real theta = (v_esc_targ * v_esc_targ) / ((problem::impact_velocity * v_esc_sys) *
                                              (problem::impact_velocity * v_esc_sys));
    Real b_max = r_crit * std::sqrt(1.0_rt + theta);

    // Now calculate the effective impact parameter and impact velocity at the current location (in physical units).

    Real v_imp = std::sqrt((problem::impact_velocity * v_esc_sys) *
                           (problem::impact_velocity * v_esc_sys) +
                           v_esc_targ * v_esc_targ);

    Real b_imp = (problem::impact_parameter * b_max  * problem::impact_velocity * v_esc_sys) / (v_imp * r_crit);

    const Real impactor_mass_ratio = problem::target_mass / (problem::target_mass + problem::impactor_mass);

    problem::target_initial_center[0] = -impactor_mass_ratio * b_imp * r_crit;
    problem::target_initial_center[1] =  impactor_mass_ratio * r_crit * std::sqrt(1.0_rt - b_imp * b_imp);
    problem::target_initial_center[2] =  0.0_rt;

    problem::target_initial_velocity[0] =  0.0_rt;
    problem::target_initial_velocity[1] = -impactor_mass_ratio * v_imp;
    problem::target_initial_velocity[2] =  0.0_rt;

    const Real target_mass_ratio = 1.0_rt - impactor_mass_ratio;

    problem::impactor_initial_center[0] =  target_mass_ratio * b_imp * r_crit;
    problem::impactor_initial_center[1] = -target_mass_ratio * r_crit * std::sqrt(1.0_rt - b_imp * b_imp);
    problem::impactor_initial_center[2] =  0.0_rt;

    problem::impactor_initial_velocity[0] =  0.0_rt;
    problem::impactor_initial_velocity[1] =  target_mass_ratio * v_imp;
    problem::impactor_initial_velocity[2] =  0.0_rt;

    // Make sure the domain is big enough to hold both the target and impactor.

    const Real* problo = DefaultGeometry().ProbLo();
    const Real* probhi = DefaultGeometry().ProbHi();

    if (problem::target_initial_center[0] - problem::target_radius <= problo[0] ||
        problem::target_initial_center[1] - problem::target_radius <= problo[1] ||
        problem::target_initial_center[2] - problem::target_radius <= problo[2]) {
        amrex::Error("ERROR: The domain width is too small to include the target.");
    }

    if (problem::impactor_initial_center[0] + problem::impactor_radius >= probhi[0] ||
        problem::impactor_initial_center[1] + problem::impactor_radius >= probhi[1] ||
        problem::impactor_initial_center[2] + problem::impactor_radius >= probhi[2]) {
        amrex::Error("ERROR: The domain width is too small to include the impactor.");
    }

    amrex::Print() << "Target location: (" << std::setprecision(3) << std::scientific
                   << problem::target_initial_center[0] << ", "
                   << problem::target_initial_center[1] << ", "
                   << problem::target_initial_center[2] << ")"
                   << std::endl;

    amrex::Print() << "Target velocity: (" << std::setprecision(3) << std::scientific
                   << problem::target_initial_velocity[0] << ", "
                   << problem::target_initial_velocity[1] << ", "
                   << problem::target_initial_velocity[2] << ")"
                   << std::endl;

    amrex::Print() << "Impactor location: (" << std::setprecision(3) << std::scientific
                   << problem::impactor_initial_center[0] << ", "
                   << problem::impactor_initial_center[1] << ", "
                   << problem::impactor_initial_center[2] << ")"
                   << std::endl;

    amrex::Print() << "Impactor velocity: (" << std::setprecision(3) << std::scientific
                   << problem::impactor_initial_velocity[0] << ", "
                   << problem::impactor_initial_velocity[1] << ", "
                   << problem::impactor_initial_velocity[2] << ")"
                   << std::endl;

    amrex::Print() << std::endl;
}

#endif
