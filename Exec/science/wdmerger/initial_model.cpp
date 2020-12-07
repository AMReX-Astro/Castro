#include <initial_model.H>
#include <castro_params.H>
#include <prob_parameters.H>
#include <model_parser_data.H>
#include <eos.H>
#include <ambient.H>

void establish_hse (model::initial_model_t& model,
                    Real& mass_want, Real& central_density_want,
                    Real envelope_mass, Real& radius,
                    const Real core_comp[NumSpec], const Real envelope_comp[NumSpec])
{
    // Note that if central_density > 0, then this initial model generator will use it in calculating
    // the model. If mass is also provided in this case, we assume it is an estimate used for the purpose of 
    // determining the envelope mass boundary.

    // Check to make sure we've specified at least one of them.

    if (mass_want < 0.0_rt && central_density_want < 0.0_rt) {
        amrex::Error("Error: Must specify either mass or central density in the initial model generator.");
    }

    // If we are specifying the mass, then we don't know what WD central density
    // will give the desired total mass, so we need to do a secant iteration
    // over central density. rho_c_old is the 'old' guess for the central
    // density and rho_c is the current guess.  After two loops, we can
    // start estimating the density required to yield our desired mass.

    // If instead we are specifying the central density, then we only need to do a 
    // single HSE integration.

    const int max_hse_iter = 250;
    int max_mass_iter;

    Real rho_c, rho_c_old, drho_c;
    Real mass, mass_old;
    Real p_want, drho;

    if (central_density_want > 0.0_rt) {

        max_mass_iter = 1;

        rho_c_old = central_density_want;
        rho_c     = central_density_want;

    }
    else {

        max_mass_iter = max_hse_iter;

        rho_c_old = -1.0_rt;
        rho_c     = 1.e7_rt;     // A reasonable starting guess for moderate-mass WDs

    }

    // Check to make sure the initial temperature makes sense.

    if (problem::stellar_temp < castro::small_temp) {
        amrex::Error("Error: WD central temperature is less than small_temp. Aborting.");
    }

    bool mass_converged = false;

    for (int mass_iter = 1; mass_iter <= max_mass_iter; ++mass_iter) {

        bool fluff = false;

        // We start at the center of the WD and integrate outward.  Initialize
        // the central conditions.

        model.state(0, model::itemp) = problem::stellar_temp;
        model.state(0, model::idens) = rho_c;
        for (int n = 0; n < NumSpec; ++n) {
            model.state(0, model::ispec + n) = core_comp[n];
        }

        eos_t eos_state;
        eos_state.rho  = model.state(0, model::idens);
        eos_state.T    = model.state(0, model::itemp);
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = model.state(0, model::ispec + n);
        }

        eos(eos_input_rt, eos_state);

        model.state(0, model::ipres) = eos_state.p;

        model.r(0) = 0.5 * problem::initial_model_dx;

        int icutoff = NPTS_MODEL;

        // Make the initial guess be completely uniform.

        for (int i = 1; i < NPTS_MODEL; ++i) {
            model.state(i, model::idens) = model.state(0, model::idens);
            model.state(i, model::itemp) = model.state(0, model::itemp);
            model.state(i, model::ipres) = model.state(0, model::ipres);
            for (int n = 0; n < NumSpec; ++n) {
                model.state(i, model::ispec + n) = model.state(0, model::ispec + n);
            }
            model.r(i) = model.r(i-1) + problem::initial_model_dx;
        }

        // Keep track of the mass enclosed below the current zone.

        Real rl = 0.0_rt;
        Real rr = rl + problem::initial_model_dx;

        Real M_enclosed = (4.0_rt / 3.0_rt) * M_PI * (std::pow(rr, 3) - std::pow(rl, 3)) * model.state(0, model::idens);
        mass = M_enclosed;

        //-------------------------------------------------------------------------
        // HSE solve
        //-------------------------------------------------------------------------
        for (int i = 1; i < problem::initial_model_npts; ++i) {

            rl += problem::initial_model_dx;
            rr += problem::initial_model_dx;

            // As the initial guess for the density, use the underlying zone.

            model.state(i, model::idens) = model.state(i-1, model::idens);

            if (mass_want > 0.0_rt && M_enclosed >= mass_want - envelope_mass) {
                for (int n = 0; n < NumSpec; ++n) {
                    model.state(i, model::ispec + n) = envelope_comp[n];
                    eos_state.xn[n] = model.state(i, model::ispec + n);
                }
            }
            else {
                for (int n = 0; n < NumSpec; ++n) {
                    model.state(i, model::ispec + n) = core_comp[n];
                    eos_state.xn[n] = model.state(i, model::ispec + n);
                }
            }

            Real g = -C::Gconst * M_enclosed / (std::pow(rl, 2));


            //----------------------------------------------------------------------
            // Iteration loop
            //----------------------------------------------------------------------

            // Start off the Newton loop by assuming that the zone has not converged.

            bool converged_hse = false;

            for (int hse_iter = 1; hse_iter <= max_hse_iter; ++hse_iter) {

                if (fluff) {
                    model.state(i, model::idens) = ambient::ambient_state[URHO];
                    eos_state.rho = ambient::ambient_state[URHO];
                    break;
                }

                // The core is isothermal, so we just need to constrain
                // the density and pressure to agree with the EOS and HSE.

                // We difference HSE about the interface between the current
                // zone and the one just inside.

                Real rho_avg = 0.5_rt * (model.state(i, model::idens) + model.state(i-1, model::idens));
                p_want = model.state(i-1, model::ipres) + problem::initial_model_dx * rho_avg * g;

                eos(eos_input_rt, eos_state);

                drho = (p_want - eos_state.p) / (eos_state.dpdr - 0.5_rt * problem::initial_model_dx * g);

                model.state(i, model::idens) = amrex::max(0.9_rt * model.state(i, model::idens),
                                                     amrex::min(model.state(i, model::idens) + drho,
                                                                1.1_rt * model.state(i, model::idens)));
                eos_state.rho = model.state(i, model::idens);

                if (model.state(i, model::idens) < ambient::ambient_state[URHO]) {
                    icutoff = i;
                    fluff = true;
                }

                if (std::abs(drho) < problem::initial_model_hse_tol * model.state(i, model::idens)) {
                    converged_hse = true;
                    break;
                }

            }

            if (!converged_hse && (!fluff)) {

                std::cout << "Error: zone " <<  i << " did not converge in init_hse()" << std::endl;
                std::cout << model.state(i, model::idens) << " " << model.state(i, model::itemp) << std::endl;
                std::cout << p_want << " " << eos_state.p;
                std::cout << drho << " " << problem::initial_model_hse_tol * model.state(i, model::idens);
                amrex::Error("Error: HSE non-convergence.");

            }

            // Call the EOS to establish the final properties of this zone.

            eos(eos_input_rt, eos_state);

            model.state(i, model::ipres) = eos_state.p;

            // Discretize the mass enclosed as (4 pi / 3) * rho * dr * (rl**2 + rl * rr + rr**2).

            Real dM = (4.0_rt / 3.0_rt) * M_PI * model.state(i, model::idens) * problem::initial_model_dx *
                      (rr * rr + rl * rr + rl * rl);
            M_enclosed += dM;

            if (i <= icutoff) {
                // Also update the final WD mass if we're not in the ambient material.
                mass += dM;
            }

        } // End loop over zones

        radius = model.r(icutoff);

        if (rho_c_old < 0.0_rt) {

            // Not enough iterations yet -- use an arbitrary guess for the next iteration.

            rho_c_old = rho_c;
            rho_c = 0.5_rt * rho_c_old;

        }
        else {

            // Check if we have converged.

            if (std::abs(mass - mass_want) / mass_want < problem::initial_model_mass_tol) {
                mass_converged = true;
                break;
            }

            // Do a secant iteration:
            // M_tot = M(rho_c) + dM/drho |_rho_c x drho + ...

            drho_c = (mass_want - mass) / ((mass  - mass_old) / (rho_c - rho_c_old));

            rho_c_old = rho_c;
            rho_c = amrex::min(1.1e0_rt * rho_c_old, amrex::max((rho_c + drho_c), 0.9e0_rt * rho_c_old));

        }

        mass_old = mass;

    } // End mass constraint loop

    if (!mass_converged && max_mass_iter > 1) {
        amrex::Error("ERROR: WD mass did not converge.");
    }

    central_density_want = model.state(0, model::idens);
    mass_want = mass;

    model::initialized = true;
    model::npts = problem::initial_model_npts;
}
