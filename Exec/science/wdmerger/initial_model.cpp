#include <initial_model.H>

#include <castro_params.H>

#include <eos.H>

using namespace initial_model;

extern "C" {

void establish_hse (model& model,
                    Real rho[initial_model_max_npts],
                    Real T[initial_model_max_npts],
                    Real xn[NumSpec][initial_model_max_npts],
                    Real r[initial_model_max_npts])
{
    // Note that if central_density > 0, then this initial model generator will use it in calculating
    // the model. If mass is also provided in this case, we assume it is an estimate used for the purpose of 
    // determining the envelope mass boundary.

    // Check to make sure we've specified at least one of them.

    if (model.mass < 0.0_rt && model.central_density < 0.0_rt) {
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
    Real mass, mass_old, radius;
    Real p_want, p_last, drho;

    if (model.central_density > 0.0_rt) {

        max_mass_iter = 1;

        rho_c_old = model.central_density;
        rho_c     = model.central_density;

    }
    else {

        max_mass_iter = max_hse_iter;

        rho_c_old = -1.0_rt;
        rho_c     = 1.e7_rt;     // A reasonable starting guess for moderate-mass WDs

    }

    // Check to make sure the initial temperature makes sense.

    if (model.central_temp < castro::small_temp) {
        amrex::Error("Error: WD central temperature is less than small_temp. Aborting.");
    }

    bool mass_converged = false;

    for (int mass_iter = 1; mass_iter <= max_mass_iter; ++mass_iter) {

        bool fluff = false;

        // We start at the center of the WD and integrate outward.  Initialize
        // the central conditions.

        T[0]    = model.central_temp;
        rho[0]  = rho_c;
        for (int n = 0; n < NumSpec; ++n) {
            xn[n][0] = model.core_comp[n];
        }

        eos_t eos_state;
        eos_state.rho  = rho[0];
        eos_state.T    = T[0];
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = xn[n][0];
        }

        eos(eos_input_rt, eos_state);

        p_last = eos_state.p;

        int icutoff;

        // Make the initial guess be completely uniform.

        for (int i = 0; i < initial_model_max_npts; ++i) {
            rho[i] = rho[0];
            T[i]   = T[0];
            for (int n = 0; n < NumSpec; ++n) {
                xn[n][i] = xn[n][0];
            }
        }

        // Keep track of the mass enclosed below the current zone.

        model.M_enclosed[0] = (4.0_rt / 3.0_rt) * M_PI * (std::pow(model.rr[0], 3) - std::pow(model.rl[0], 3)) * rho[0];

        //-------------------------------------------------------------------------
        // HSE solve
        //-------------------------------------------------------------------------
        for (int i = 1; i < model.npts; ++i) {

            // As the initial guess for the density, use the underlying zone.

            rho[i] = rho[i-1];

            if (model.mass > 0.0_rt && model.M_enclosed[i-1] >= model.mass - model.envelope_mass) {
                for (int n = 0; n < NumSpec; ++n) {
                    xn[n][i] = model.envelope_comp[n];
                    eos_state.xn[n] = xn[n][i];
                }
            }
            else{
                for (int n = 0; n < NumSpec; ++n) {
                    xn[n][i] = model.core_comp[n];
                    eos_state.xn[n] = xn[n][i];
                }
            }

            model.g[i] = -C::Gconst * model.M_enclosed[i-1] / (std::pow(model.rl[i], 2));


            //----------------------------------------------------------------------
            // Iteration loop
            //----------------------------------------------------------------------

            // Start off the Newton loop by assuming that the zone has not converged.

            bool converged_hse = false;

            for (int hse_iter = 1; hse_iter <= max_hse_iter; ++hse_iter) {

                if (fluff) {
                    rho[i] = model.min_density;
                    eos_state.rho = model.min_density;
                    break;
                }

                // The core is isothermal, so we just need to constrain
                // the density and pressure to agree with the EOS and HSE.

                // We difference HSE about the interface between the current
                // zone and the one just inside.

                Real rho_avg = 0.5_rt * (rho[i] + rho[i-1]);
                p_want = p_last + model.dx * rho_avg * model.g[i];

                eos(eos_input_rt, eos_state);

                drho = (p_want - eos_state.p) / (eos_state.dpdr - 0.5_rt * model.dx * model.g[i]);

                rho[i] = amrex::max(0.9 * rho[i], min(rho[i] + drho, 1.1 * rho[i]));
                eos_state.rho = rho[i];

                if (rho[i] < model.min_density) {
                    icutoff = i;
                    fluff = true;
                }

                if (std::abs(drho) < model.hse_tol * rho[i]) {
                    converged_hse = true;
                    break;
                }

            }

            if (!converged_hse && (!fluff)) {

                std::cout << "Error: zone " <<  i << " did not converge in init_hse()" << std::endl;
                std::cout << rho[i] << " " << T[i] << std::endl;
                std::cout << p_want << " " << eos_state.p;
                std::cout << drho << " " << model.hse_tol * rho[i];
                amrex::Error("Error: HSE non-convergence.");

            }

            // Call the EOS to establish the final properties of this zone.

            eos(eos_input_rt, eos_state);

            p_last = eos_state.p;

            // Discretize the mass enclose as (4 pi / 3) * rho * dr * (rl**2 + rl * rr + rr**2).

            model.M_enclosed[i] = model.M_enclosed[i-1] +
                                  (4.0_rt / 3.0_rt) * M_PI * rho[i] * model.dx *
                                  (std::pow(model.rr[i], 2) + model.rl[i] * model.rr[i] + std::pow(model.rl[i], 2));

        } // End loop over zones

        mass = model.M_enclosed[icutoff];
        radius = model.r[icutoff];

        if (rho_c_old < 0.0_rt) {

            // Not enough iterations yet -- use an arbitrary guess for the next iteration.

            rho_c_old = rho_c;
            rho_c = 0.5_rt * rho_c_old;

        }
        else {

            // Check if we have converged.

            if (std::abs(mass - model.mass) / model.mass < model.mass_tol) {
                mass_converged = true;
                break;
            }

            // Do a secant iteration:
            // M_tot = M(rho_c) + dM/drho |_rho_c x drho + ...

            drho_c = (model.mass - mass) / ((mass  - mass_old) / (rho_c - rho_c_old));

            rho_c_old = rho_c;
            rho_c = amrex::min(1.1e0_rt * rho_c_old, amrex::max((rho_c + drho_c), 0.9e0_rt * rho_c_old));

        }

        mass_old = mass;

    } // End mass constraint loop

    if (!mass_converged && max_mass_iter > 1) {
        amrex::Error("ERROR: WD mass did not converge.");
    }

    model.central_density = rho[0];
    model.radius = radius;
    model.mass = mass;
}

}
