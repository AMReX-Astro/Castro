#include <wdmerger_util.H>
#include <binary.H>
#include <network.H>
#include <ambient.H>
#include <model_parser.H>

// Given total mass of a binary system and the initial separation of
// two point particles, obtain the velocity at this separation
// assuming the point masses fell in from infinity. This will
// be the velocity in the frame where the center of mass is stationary.

void freefall_velocity (Real mass, Real distance, Real& vel)
{
    vel = std::sqrt(2.0_rt * C::Gconst * mass / distance);
}

// Accepts the masses of two stars (in solar masses)
// and the orbital period of a system,
// and returns the semimajor axis of the orbit (in cm),
// as well as the distances a_1 and a_2 from the center of mass.

void kepler_third_law (Real radius_1, Real mass_1, Real radius_2, Real mass_2,
                       Real& period, Real eccentricity, Real phi, Real& a,
                       Real& r_1, Real& r_2, Real& v_1r, Real& v_2r, Real& v_1p, Real& v_2p)
{
    amrex::ignore_unused(radius_1);
    amrex::ignore_unused(radius_2);

    Real M  = mass_1 + mass_2;     // Total mass
    Real mu = mass_1 * mass_2 / M; // Reduced mass

    // First, solve for the orbit in the reduced one-body problem, where
    // an object of mass mu orbits an object with mass M located at r = 0.
    // For this we follow Carroll and Ostlie, Chapter 2, but many texts discuss this.
    // Note that we use the convention that phi measures angle from aphelion,
    // which is opposite to the convention they use.

    if (period > 0.0_rt && a < 0.0_rt)
    {
        a = std::pow(C::Gconst * M * period * period / (4.0_rt * M_PI * M_PI), 1.0_rt / 3.0_rt); // C + O, Equation 2.37
    }
    else if (period < 0.0_rt && a > 0.0_rt)
    {
        period = std::sqrt(a * a * a * 4.0_rt * M_PI * M_PI / (C::Gconst * M));
    }
    else {
        amrex::Error("Error: overspecified Kepler's third law calculation.");
    }

    Real r = a * (1.0_rt - eccentricity * eccentricity) / (1.0_rt - eccentricity * std::cos(phi)); // C + O, Equation 2.3

    // To get the radial and azimuthal velocity, we take the appropriate derivatives of the above.
    // v_r = dr / dt = dr / d(phi) * d(phi) / dt, with d(phi) / dt being derived from
    // C + O, Equation 2.30 for the angular momentum, and the fact that L = mu * r**2 * d(phi) / dt.

    Real v_r   = -2.0_rt * M_PI * a * eccentricity * std::sin(phi) /
                  (period * std::sqrt(1.0_rt - eccentricity * eccentricity));
    Real v_phi =  2.0_rt * M_PI * a * (1.0_rt - eccentricity * std::cos(phi)) /
                  (period * std::sqrt(1.0_rt - eccentricity * eccentricity));

    // Now convert everything back to the binary frame, using C+O, Equation 2.23 and 2.24. This applies
    // to the velocities as well as the positions because the factor in front of r_1 and r_2 is constant.

    r_1  = -(mu / mass_1) * r;
    r_2  =  (mu / mass_2) * r;

    v_1r = -(mu / mass_1) * v_r;
    v_2r =  (mu / mass_2) * v_r;

    v_1p = -(mu / mass_1) * v_phi;
    v_2p =  (mu / mass_2) * v_phi;
}

// Given a WD mass, set its core and envelope composition.

void set_wd_composition (Real mass, Real& envelope_mass, Real core_comp[NumSpec], Real envelope_comp[NumSpec], const std::string& star_type)
{
    int iHe4 = network_spec_index("helium-4");
    int iC12 = network_spec_index("carbon-12");
    int iO16 = network_spec_index("oxygen-16");
    int iNe20 = network_spec_index("neon-20");
    int iMg24 = network_spec_index("magnesium-24");

    for (int n = 0; n < NumSpec; ++n) {
        core_comp[n] = small_x;
        envelope_comp[n] = small_x;
    }

    envelope_mass = 0.0_rt;

    // Here we follow the prescription of Dan et al. 2012.

    if (mass > 0.0_rt && mass < problem::max_he_wd_mass) {

        if (iHe4 < 0) {
            amrex::Error("Must have He4 in the nuclear network.");
        }

        core_comp[iHe4] = 1.0_rt;

        amrex::Print() << "Created a pure He " << star_type << "." << std::endl;

        for (int n = 0; n < NumSpec; ++n) {
            envelope_comp[n] = core_comp[n];
        }

    }
    else if (mass >= problem::max_he_wd_mass && mass < problem::max_hybrid_wd_mass) {

        if (iC12 < 0) {
            amrex::Error("Must have C12 in the nuclear network.");
        }
        if (iO16 < 0) {
            amrex::Error("Must have O16 in the nuclear network.");
        }

        core_comp[iC12] = problem::hybrid_wd_c_frac;
        core_comp[iO16] = problem::hybrid_wd_o_frac;

        envelope_mass = problem::hybrid_wd_he_shell_mass;

        amrex::Print()<< "Creating " << star_type << "with CO core with mass fractions C = "<< problem::hybrid_wd_c_frac << "and O = "
           << problem::hybrid_wd_o_frac <<" and a He shell of solar mass =" << problem::hybrid_wd_he_shell_mass << "." << std::endl;


        if (envelope_mass > 0.0_rt) {
            if (iHe4 < 0) {
                amrex::Error("Must have He4 in the nuclear network.");
            }
            envelope_comp[iHe4] = 1.0_rt;
        }
        else {
            for (int n = 0; n < NumSpec; ++n) {
                envelope_comp[n] = core_comp[n];
            }
        }
    }
    else if (mass >= problem::max_hybrid_wd_mass && mass < problem::max_co_wd_mass) {

        if (iC12 < 0) {
            amrex::Error("Must have C12 in the nuclear network.");
        }
        if (iO16 < 0) {
            amrex::Error("Must have O16 in the nuclear network.");
        }

        core_comp[iC12] = problem::co_wd_c_frac;
        core_comp[iO16] = problem::co_wd_o_frac;

        envelope_mass = problem::co_wd_he_shell_mass;

        amrex::Print()<<"Creating " << star_type << " with CO core with mass fractions C = "<<problem::co_wd_c_frac<<"and O = "
            <<problem::co_wd_o_frac<<" and a He shell of solar mass ="<<problem::co_wd_he_shell_mass<< "." << std::endl;

        if (envelope_mass > 0.0_rt) {
            if (iHe4 < 0) {
                amrex::Error("Must have He4 in the nuclear network.");
            }
            envelope_comp[iHe4] = 1.0_rt;
        }
        else {
            for (int n = 0; n < NumSpec; ++n) {
                envelope_comp[n] = core_comp[n];
            }
        }

    }
    else if (mass > problem::max_co_wd_mass) {

        if (iO16 < 0) {
            amrex::Error("Must have O16 in the nuclear network.");
        }
        if (iNe20 < 0) {
            amrex::Error("Must have Ne20 in the nuclear network.");
        }
        if (iMg24 < 0) {
            amrex::Error("Must have Mg24 in the nuclear network.");
        }

        core_comp[iO16]  = problem::onemg_wd_o_frac;
        core_comp[iNe20] = problem::onemg_wd_ne_frac;
        core_comp[iMg24] = problem::onemg_wd_mg_frac;

        amrex::Print()<<"Creating an ONeMg " << star_type << "." <<std::endl;

        for (int n = 0; n < NumSpec; ++n) {
            envelope_comp[n] = core_comp[n];
        }

    }

    // Normalize compositions so that they sum to one.

    Real core_sum = 0.0_rt;
    Real envelope_sum = 0.0_rt;

    for (int n = 0; n < NumSpec; ++n) {
        core_sum += core_comp[n];
        envelope_sum += envelope_comp[n];
    }

    for (int n = 0; n < NumSpec; ++n) {
        core_comp[n] /= core_sum;
        envelope_comp[n] /= envelope_sum;
    }
}

// This routine checks to see if the primary mass is actually larger
// than the secondary mass, and switches them if not.

void ensure_primary_mass_larger ()
{
    // We want the primary WD to be more massive. If what we're calling
    // the primary is less massive, switch the stars.

    if (problem::mass_P < problem::mass_S) {

        amrex::Print() << "Primary mass is less than secondary mass; switching the stars so that the primary is more massive." << std::endl;

        Real temp_mass = problem::mass_P;
        problem::mass_P = problem::mass_S;
        problem::mass_S = temp_mass;

    }
}

// This routine calls all of the other subroutines at the beginning
// of a simulation to fill in the basic problem data.

void initialize_problem ()
{
    // We have already read in the namelist; do some parameter postprocessing.

    finalize_probdata();

    // Establish binary parameters and create initial models.

    binary_setup();

    // Set small_pres and small_ener.

    set_small();
}



// Postprocess the parameters.

void finalize_probdata ()
{
    // Convert masses from solar masses to grams.

    problem::mass_P *= C::M_solar;
    problem::mass_S *= C::M_solar;

    problem::max_he_wd_mass *= C::M_solar;
    problem::max_hybrid_wd_mass *= C::M_solar;
    problem::max_co_wd_mass *= C::M_solar;

    problem::hybrid_wd_he_shell_mass *= C::M_solar;
    problem::co_wd_he_shell_mass *= C::M_solar;

    if (problem::mass_S < 0.0_rt && problem::central_density_S < 0.0_rt) {
        problem::single_star = 1;
    }

    // Make sure that the primary mass is really the larger mass

    ensure_primary_mass_larger();

    // We enforce that the orbital plane is (z, phi) for two-dimensional problems,
    // with rotation about z = 0 along the radial axis.

    if (AMREX_SPACEDIM == 2) {

        if (DefaultGeometry().Coord() != 1) {
            amrex::Error("We only support cylindrical coordinates in two dimensions. Set coord_type == 1.");
        }

        problem::axis_1 = 2;
        problem::axis_2 = 3;
        castro::rot_axis = 1;

    }

    // Make sure we have a sensible collision impact parameter.

    if (problem::collision_impact_parameter > 1.0) {
        amrex::Error("Impact parameter must be less than one in our specified units.");
    }

    // Safety check: we can't run most problems in one dimension.

    if (AMREX_SPACEDIM == 1 && problem::problem != 0) {
        amrex::Error("Can only run a collision or freefall in 1D. Exiting.");
    }

    // Don't do a collision, free-fall, or TDE in a rotating reference frame.

    if (problem::problem == 0 && castro::do_rotation == 1) {
        amrex::Error("The free-fall/collision problem does not make sense in a rotating reference frame.");
    }

    if (problem::problem == 2 && castro::do_rotation == 1) {
        amrex::Error("The TDE problem does not make sense in a rotating reference frame.");
    }

    // Make sure we have a sensible eccentricity.

    if (problem::orbital_eccentricity >= 1.0) {
        amrex::Error("Orbital eccentricity cannot be larger than one.");
    }

    // Make sure we have a sensible angle. Then convert it to radians.

    if (problem::orbital_angle < 0.0 || problem::orbital_angle > 360.0) {
        amrex::Error("Orbital angle must be between 0 and 360 degrees.");
    }

    problem::orbital_angle *= M_PI / 180.0;

    // If we're doing a relaxation, we need to reset the relaxation_is_done parameter.
    // This will be reset as appropriate from the checkpoint if we're performing a restart.

    if (problem::problem == 1 && problem::relaxation_damping_factor > 0.0_rt) {
        problem::relaxation_is_done = 0;
    }

    // As above, but for radial damping.

    if (problem::problem == 1 && problem::radial_damping_velocity_factor > 0.0_rt) {
        problem::radial_damping_is_done = 0;
    }

    // TDE sanity checks

    if (problem::problem == 2) {

        // We must have a BH point mass defined.

        if (castro::point_mass <= 0.0_rt) {

            amrex::Error("No point mass specified for the TDE problem.");

        }


        // This problem cannot have a secondary WD.

        if (problem::mass_S >= 0.0_rt) {

            amrex::Error("TDE problem cannot have a secondary WD.");

        }

        // Beta parameter must be positive.

        if (problem::tde_beta <= 0.0_rt) {

            amrex::Error("TDE beta must be positive.");

        }

        // Initial distance must be positive.

        if (problem::tde_separation <= 0.0_rt) {

            amrex::Error("TDE separation must be positive.");

        }

    }

}



// Calculate small_pres and small_ener

void set_small ()
{
    eos_t eos_state;

    // Given the inputs of small_dens and small_temp, figure out small_pres.

    eos_state.rho = castro::small_dens;
    eos_state.T   = castro::small_temp;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = ambient::ambient_state[UFS+n] / ambient::ambient_state[URHO];
    }
#ifdef AUX_THERMO
    set_aux_comp_from_X(eos_state);
#endif

    eos(eos_input_rt, eos_state);

    castro::small_pres = eos_state.p;
    castro::small_ener = eos_state.e;
}



// Update the locations of the Roche radii

void update_roche_radii ()
{
    if (problem::mass_P > 0.0_rt && problem::mass_S > 0.0_rt) {
        Real r = std::sqrt(std::pow(problem::com_P[0] - problem::com_S[0], 2) +
                           std::pow(problem::com_P[1] - problem::com_S[1], 2) +
                           std::pow(problem::com_P[2] - problem::com_S[2], 2));

        get_roche_radii(problem::mass_S / problem::mass_P, problem::roche_rad_S, problem::roche_rad_P, r);

        // Beyond a certain point, it doesn't make sense to track the stars separately
        // anymore. We'll set the secondary to a fixed constant and keep it there
        // if its Roche radius becomes smaller than 10% of the primary's. Also, for exactly
        // equal mass systems sometimes it is the primary that disrupts, perhaps
        // just due to numerical noise, so do the same check for the primary.

        if (problem::roche_rad_S < problem::roche_rad_P / 10.0_rt) {
            for (int n = 0; n < 3; ++n) {
                problem::com_S[n] = problem::center[n];
                problem::vel_S[n] = 0.0_rt;
            }
            problem::mass_S = 0.0_rt;
            problem::roche_rad_S = 0.0_rt;
            problem::t_ff_S = 0.0_rt;
        }
        else if (problem::roche_rad_P < problem::roche_rad_S / 10.0_rt) {
            for (int n = 0; n < 3; ++n) {
                problem::com_P[n] = problem::center[n];
                problem::vel_P[n] = 0.0_rt;
            }
            problem::mass_S = 0.0_rt;
            problem::roche_rad_P = 0.0_rt;
            problem::t_ff_P = 0.0_rt;
        }
    }

    // Determine when the radial damping force terminates.

    if (problem::problem == 1 && problem::radial_damping_velocity_factor > 0.0_rt && problem::radial_damping_is_done != 1) {
        // It does not make sense to do damping if there's only one star remaining.

        if (problem::mass_S == 0.0_rt) {
            problem::radial_damping_is_done = 1;
        }

        // Only do radial damping until the stars get sufficiently close. The way we will measure this
        // is when the secondary overflows the Roche lobe. At that point we terminate the damping and
        // let Newtonian gravity do the rest of the work.

        if (problem::roche_rad_S <= problem::radial_damping_roche_factor * problem::radius_S) {
            problem::radial_damping_is_done = 1;
        }

        if (problem::radial_damping_is_done == 1) {
            amrex::Print() << "\n\n  Terminating radial damping force since the secondary is about to overflow its Roche lobe.\n\n";
        }
    }
}



// Set up a binary simulation

void binary_setup ()
{
    // Safety check: ensure that if we have a symmetric lower boundary, that the
    // domain center (and thus the stars) are on that boundary.

    const Real* problo = DefaultGeometry().ProbLo();
    const Real* probhi = DefaultGeometry().ProbHi();

    if (Castro::physbc().lo(0) == amrex::PhysBCType::symmetry && problem::center[0] != problo[0]) {
        amrex::Error("Symmetric lower x-boundary but the center is not on this boundary.");
    }

    if (Castro::physbc().lo(1) == amrex::PhysBCType::symmetry && problem::center[1] != problo[1]) {
        amrex::Error("Symmetric lower y-boundary but the center is not on this boundary.");
    }

    if (Castro::physbc().lo(2) == amrex::PhysBCType::symmetry && problem::center[2] != problo[2]) {
        amrex::Error("Symmetric lower z-boundary but the center is not on this boundary.");
    }

    // Set some default values for these quantities;
    // we'll update them soon.

    for (int n = 0; n < 3; ++n) {
        problem::center_P_initial[n] = problem::center[n];
        problem::center_S_initial[n] = problem::center[n];

        problem::com_P[n] = problem::center[n];
        problem::com_S[n] = problem::center[n];

        problem::vel_P[n] = 0.0_rt;
        problem::vel_S[n] = 0.0_rt;
    }

    // Fill in the model's physical details.
    // If we're integrating to reach a desired mass, set the composition accordingly.
    // If instead we're fixing the central density, then first we'll assume the composition is
    // that of a solar mass WD as a initial guess, and get the corresponding mass.
    // Then we set the composition to match this preliminary mass, and we'll get a final mass later.

    if (problem::central_density_P > 0.0_rt) {
        // Initial guess for the mass
        problem::mass_P = C::M_solar;
    }
    else if (problem::mass_P <= 0.0_rt) {
        amrex::Error("Must specify either a positive primary mass or a positive primary central density.");
    }

    set_wd_composition(problem::mass_P, problem::envelope_mass_P, problem::core_comp_P, problem::envelope_comp_P, "primary");



    if (!problem::single_star) {

        if (problem::central_density_S > 0.0_rt) {
            // Initial guess for the mass
            problem::mass_S = C::M_solar;
        }
        else if (problem::mass_S <= 0.0_rt) {
            amrex::Error("If we are doing a binary calculation, we must specify either a positive secondary mass or a positive secondary central density");
        }

        set_wd_composition(problem::mass_S, problem::envelope_mass_S, problem::core_comp_S, problem::envelope_comp_S, "secondary");

        for (int n = 0; n < NumSpec; ++n) {
            ambient::ambient_state[UFS+n] = ambient::ambient_state[URHO] * (problem::envelope_comp_P[n] + problem::envelope_comp_S[n]) / 2;
        }

    }
    else {

        for (int n = 0; n < NumSpec; ++n) {
            ambient::ambient_state[UFS+n] = ambient::ambient_state[URHO] * problem::envelope_comp_P[n];
        }

    }

    // We have completed (rho, T, xn) for the ambient state, so we can call the EOS.

    eos_t eos_state;
    eos_state.rho = ambient::ambient_state[URHO];
    eos_state.T   = ambient::ambient_state[UTEMP];
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = ambient::ambient_state[UFS+n] / ambient::ambient_state[URHO];
    }
#ifdef AUX_THERMO
    set_aux_comp_from_X(eos_state);
#endif

    eos(eos_input_rt, eos_state);

    ambient::ambient_state[UEINT] = ambient::ambient_state[URHO] * eos_state.e;
    ambient::ambient_state[UEDEN] = ambient::ambient_state[UEINT];



    problem::roche_rad_P = 0.0_rt;
    problem::roche_rad_S = 0.0_rt;

    // Generate primary and secondary WD models.

    establish_hse(problem::mass_P, problem::central_density_P, problem::radius_P,
                  problem::core_comp_P, problem::stellar_temp, problem::initial_model_dx,
                  problem::envelope_mass_P, problem::envelope_comp_P, 0);

    amrex::Print() << std::endl;

    amrex::Print() << "Generated initial model for primary WD of mass " << std::setprecision(3) << problem::mass_P / C::M_solar
                   << " solar masses, central density " << std::setprecision(3) << std::scientific << problem::central_density_P
                   << " g cm**-3, and radius " << std::setprecision(3) << std::scientific << problem::radius_P << " cm."
                   << std::endl << std::endl;

    problem::roche_rad_P = problem::radius_P;

    if (!problem::single_star) {

        establish_hse(problem::mass_S, problem::central_density_S, problem::radius_S,
                      problem::core_comp_S, problem::stellar_temp, problem::initial_model_dx,
                      problem::envelope_mass_S, problem::envelope_comp_S, 1);

        amrex::Print() << "Generated initial model for secondary WD of mass " << std::setprecision(3) << problem::mass_S / C::M_solar
                       << " solar masses, central density " << std::setprecision(3) << std::scientific << problem::central_density_S
                       << " g cm**-3, and radius " << std::setprecision(3) << std::scientific << problem::radius_S << " cm."
                       << std::endl << std::endl;

        problem::roche_rad_S = problem::radius_S;

        // Compute initial Roche radii

        get_roche_radii(problem::mass_S / problem::mass_P, problem::roche_rad_S, problem::roche_rad_P);

        // Set up the stellar distances and velocities according to the problem choice

        if (problem::problem == 0) {

            problem::collision_separation *= problem::radius_S;

            if (problem::collision_velocity < 0.0e0_rt) {

                Real v_ff;
                freefall_velocity(problem::mass_P + problem::mass_S, problem::collision_separation, v_ff);

                problem::vel_P[problem::axis_1-1] =  (problem::mass_P / (problem::mass_S + problem::mass_P)) * v_ff;
                problem::vel_S[problem::axis_1-1] = -(problem::mass_S / (problem::mass_S + problem::mass_P)) * v_ff;

            }
            else {

                problem::vel_P[problem::axis_1-1] =  problem::collision_velocity;
                problem::vel_S[problem::axis_1-1] = -problem::collision_velocity;

            }

            problem::r_P_initial = -(problem::mass_P / (problem::mass_S + problem::mass_P)) * problem::collision_separation;
            problem::r_S_initial =  (problem::mass_S / (problem::mass_S + problem::mass_P)) * problem::collision_separation;

            problem::a = problem::r_S_initial - problem::r_P_initial;

            problem::center_P_initial[problem::axis_1-1] += problem::r_P_initial;
            problem::center_S_initial[problem::axis_1-1] += problem::r_S_initial;

            // We also permit a non-zero impact parameter b in the direction perpendicular
            // to the motion of the stars. This is measured in units of the radius of the
            // primary, so that b > 1 doesn't make any sense as the stars won't collide.
            // Since the secondary's radius is greater than the primary's, measuring in the
            // units of the primary's radius will guarantee contact.

            Real collision_offset = problem::collision_impact_parameter * problem::radius_P;

            problem::center_P_initial[problem::axis_2-1] -= collision_offset;
            problem::center_S_initial[problem::axis_2-1] += collision_offset;

        }
        else if (problem::problem == 1) {

            if (problem::roche_radius_factor < 0.0_rt) {

                // Determine the orbital distance based on the rotational period.

                problem::a = -1.0_rt;

            }
            else {

                // Set the orbital distance, then calculate the rotational period.

                problem::a = problem::roche_radius_factor * (problem::radius_S / problem::roche_rad_S);

                castro::rotational_period = -1.0_rt;

            }

            Real v_P_r, v_S_r, v_P_phi, v_S_phi;

            kepler_third_law(problem::radius_P, problem::mass_P, problem::radius_S, problem::mass_S,
                             castro::rotational_period, problem::orbital_eccentricity, problem::orbital_angle,
                             problem::a, problem::r_P_initial, problem::r_S_initial, v_P_r, v_S_r, v_P_phi, v_S_phi);

            // Make sure the domain is big enough to hold stars in an orbit this size.

            Real length;

            if (Castro::physbc().lo(problem::axis_1-1) == amrex::PhysBCType::symmetry) {

                // In this case we're only modelling the secondary.
                length = problem::r_P_initial + problem::radius_P;

            }
            else {

                length = (problem::r_S_initial - problem::r_P_initial) + problem::radius_P + problem::radius_S;

            }

            if (length > (probhi[problem::axis_1-1] - problo[problem::axis_1-1])) {
                amrex::Error("ERROR: The domain width is too small to include the binary orbit.");
            }

            // Make sure the stars are not touching.

            if (problem::radius_P + problem::radius_S > problem::a) {
                amrex::Error("ERROR: Stars are touching!");
            }

            amrex::Print() << "Generated binary orbit of distance " << std::setprecision(3) << std::scientific
                           << problem::a << " cm = " << problem::a / C::AU << " AU." << std::endl;
            amrex::Print() << "The primary orbits the center of mass at distance " << std::setprecision(3) << std::scientific
                           << problem::r_P_initial << " cm = " << problem::r_P_initial / C::AU << " AU." << std::endl;
            amrex::Print() << "The secondary orbits the center of mass at distance " << std::setprecision(3) << std::scientific
                           << problem::r_S_initial << " cm = " << problem::r_S_initial / C::AU << " AU." << std::endl;
            amrex::Print() << "The initial orbital period is " << std::setprecision(6) << std::fixed
                           << castro::rotational_period << " s." << std::endl;
            amrex::Print() << "The initial orbital speed of the primary is " << std::setprecision(3) << std::scientific
                           << 2.0_rt * M_PI * problem::r_P_initial / castro::rotational_period << " cm/s." << std::endl;
            amrex::Print() << "The initial orbital speed of the secondary is " << std::setprecision(3) << std::scientific
                           << 2.0_rt * M_PI * problem::r_S_initial / castro::rotational_period << " cm/s." << std::endl;
            amrex::Print() << std::endl;

            // Star center positions -- we'll put them in the midplane, with the center of mass at the center of the domain.

            problem::center_P_initial[problem::axis_1-1] += problem::r_P_initial * std::cos(problem::orbital_angle);
            problem::center_P_initial[problem::axis_2-1] += problem::r_P_initial * std::sin(problem::orbital_angle);

            problem::center_S_initial[problem::axis_1-1] += problem::r_S_initial * std::cos(problem::orbital_angle);
            problem::center_S_initial[problem::axis_2-1] += problem::r_S_initial * std::sin(problem::orbital_angle);

            // Star velocities, from Kepler's third law. Note that these are the velocities in the inertial frame.

            problem::vel_P[problem::axis_1-1] = v_P_r   * std::cos(problem::orbital_angle) - v_P_phi * std::sin(problem::orbital_angle);
            problem::vel_P[problem::axis_2-1] = v_P_phi * std::cos(problem::orbital_angle) + v_P_r   * std::sin(problem::orbital_angle);

            problem::vel_S[problem::axis_1-1] = v_S_r   * std::cos(problem::orbital_angle) - v_S_phi * std::sin(problem::orbital_angle);
            problem::vel_S[problem::axis_2-1] = v_S_phi * std::cos(problem::orbital_angle) + v_S_r   * std::sin(problem::orbital_angle);

        }
        else {

            amrex::Error("Error: Unknown problem choice.");

        }

        // Scale the Roche radii by the initial distance.

        problem::roche_rad_P *= problem::a;
        problem::roche_rad_S *= problem::a;

    }
    else {

        if (problem::problem == 2) {

            // The tidal radius is given by (M_BH / M_WD)^(1/3) * R_WD.

            problem::tde_tidal_radius = std::pow(castro::point_mass / problem::mass_P, 1.0_rt / 3.0_rt) *
                                        problem::radius_P;

            // The usual definition for the Schwarzschild radius.

            problem::tde_schwarzschild_radius = 2.0_rt * C::Gconst * castro::point_mass / (C::c_light * C::c_light);

            // The pericenter radius is the distance of closest approach,
            // for a point mass on a parabolic orbit.

            problem::tde_pericenter_radius = problem::tde_tidal_radius / problem::tde_beta;

            // Given the pericenter distance, we can calculate the parameters of
            // the parabolic orbit. A parabolic orbit has E = 0, so KE = PE,
            // or (1/2) M_WD v**2 = G * M_WD * M_BH / r. This simplifies to
            // v = sqrt(2 * G * M_BH / r). The initial distance is set at runtime.

            problem::r_P_initial = problem::tde_separation * problem::tde_tidal_radius;

            Real v_P = std::sqrt(2.0_rt * C::Gconst * castro::point_mass / problem::r_P_initial);

            // Now we need to convert this into angular and radial components. To
            // do this, we need the orbital angle, which comes from the orbit equation,
            // r = r_0 / (1 - eccentricity * cos(phi)), where for a parabolic orbit,
            // r_0 is twice the pericenter distance.

            problem::orbital_eccentricity = 1.0_rt;
            problem::orbital_angle = std::acos(1.0_rt - (2.0_rt * problem::tde_pericenter_radius) / problem::r_P_initial);

            // Now set the x and y components of the position and velocity. The position is
            // straightforward: the orbital angle is the usual angle phi such that x = r cos(phi)
            // and y = r sin(phi). The velocity is a little more involved and depends on the orbit
            // equation. The Cartesian form of the parabolic orbit equation is x = y**2 / (2 * r_0) + r_0 / 2,
            // so dx/dt = (y / r_0) * dy/dt. Given that v**2 = v_x**2 + v_y**2, we have
            // v_x = v / sqrt( 1 + (r_0 / (r * sin(phi))) )**2 , and
            // v_y = v / sqrt( 1 + (r * sin(phi) / r_0)**2 ).

            problem::center_P_initial[problem::axis_1-1] -= problem::r_P_initial * std::cos(problem::orbital_angle);
            problem::center_P_initial[problem::axis_2-1] -= problem::r_P_initial * std::sin(problem::orbital_angle);

            if (problem::tde_initial_velocity == 1) {

                problem::vel_P[problem::axis_1-1] = v_P / std::sqrt(1.0_rt + std::pow(2.0_rt * problem::tde_pericenter_radius /
                                                                                      (problem::r_P_initial * std::sin(problem::orbital_angle)), 2));
                problem::vel_P[problem::axis_2-1] = v_P / std::sqrt(1.0_rt + std::pow(problem::r_P_initial * std::sin(problem::orbital_angle) /
                                                                                      (2.0_rt * problem::tde_pericenter_radius), 2));

            }

        }

    }

    for (int n = 0; n < 3; ++n) {
        problem::com_P[n] = problem::center_P_initial[n];
        problem::com_S[n] = problem::center_S_initial[n];
    }

    // Safety check: make sure the stars are actually inside the computational domain.

    if (!(AMREX_SPACEDIM == 2 && Castro::physbc().lo(1) == amrex::PhysBCType::symmetry)) {

        if ((0.5_rt * (probhi[0] - problo[0]) < problem::radius_P) ||
            (0.5_rt * (probhi[1] - problo[1]) < problem::radius_P) ||
            (0.5_rt * (probhi[2] - problo[2]) < problem::radius_P && AMREX_SPACEDIM == 3)) {
            amrex::Error("Primary does not fit inside the domain.");
        }

        if ((0.5_rt * (probhi[0] - problo[0]) < problem::radius_S) ||
            (0.5_rt * (probhi[1] - problo[1]) < problem::radius_S) ||
            (0.5_rt * (probhi[2] - problo[2]) < problem::radius_S && AMREX_SPACEDIM == 3) ) {
            amrex::Error("Secondary does not fit inside the domain.");
        }

    }
    else {

        if ((probhi[0] - problo[0] < problem::radius_S) ||
            (probhi[1] - problo[1] < 2.0_rt * problem::radius_S)) {
            amrex::Error("Secondary does not fit inside the domain.");
        }

    }
}
