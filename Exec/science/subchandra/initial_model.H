#ifndef INITIAL_MODEL_H
#define INITIAL_MODEL_H

#include <AMReX_Array.H>

#include <fundamental_constants.H>


namespace subch {

    // TOL_WD_MASS is tolerance used for getting the total WD mass
    // equal to M_tot (defined below).  It can be reasonably small,
    // since there will always be a central density value that can
    // give the desired WD mass on the grid we use

    constexpr Real TOL_WD_MASS = 1.e-3_rt;

    // TOL_HE_MASS is the tolerance used for getting the mass of the
    // He envelope.  This should not be too small, since the values of
    // the He envelope mass we can achieve will depend on our grid
    // spacing.

    constexpr Real TOL_HE_MASS = 2.e-3_rt;

    const int MAX_ITER = 1000;
}

struct model_t {

    amrex::Real M_tot;
    amrex::Real M_He;

    amrex::Real delta;
    amrex::Real T_core;
    amrex::Real T_base;
    amrex::Real low_density_cutoff;
    amrex::Real T_fluff;

    amrex::Real tol_hse;

    bool isothermal_layer;

    amrex::Real xn_core[NumSpec];
    amrex::Real xn_he[NumSpec];

    int ihe4;
    int ic12;
    int in14;
    int io16;

    amrex::Real X_N14;
    amrex::Real X_C12;
    amrex::Real X_O16;
};


// construct the model give a central density (rho_c) and density to
// transition to He (rho_he)

AMREX_INLINE
std::tuple<amrex::Real, amrex::Real, amrex::Real>
build_star(const model_t& model_params, const Real dx,
           const amrex::Real rho_c, const amrex::Real rho_he) {

    bool fluff = false;


    amrex::Array1D<amrex::Real, 0, NPTS_MODEL-1> rl;
    amrex::Array1D<amrex::Real, 0, NPTS_MODEL-1> rr;

    for (int i = 0; i < model::npts; i++) {
        rl(i) = static_cast<Real>(i) * dx;
        rr(i) = (static_cast<Real>(i) + 1.0) * dx;
    }

    amrex::Array1D<amrex::Real, 0, NPTS_MODEL-1> M_enclosed;

    // we start at the center of the WD and integrate outward.
    // Initialize the central conditions.

    eos_t eos_state;
    eos_state.T = model_params.T_core;
    eos_state.rho = rho_c;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = model_params.xn_core[n];
    }

    // (t, rho) -> (p, s)

    eos(eos_input_rt, eos_state);

    // make the initial guess be completely uniform

    for (int i = 0; i < model::npts; ++i) {
        model::profile(0).state(i, model::idens) = eos_state.rho;
        model::profile(0).state(i, model::itemp) = eos_state.T;
        model::profile(0).state(i, model::ipres) = eos_state.p;

        for (int n = 0; n < NumSpec; ++n) {
            model::profile(0).state(i, model::ispec+n) = eos_state.xn[n];
        }
    }


    // keep track of the mass enclosed below the current zone

    M_enclosed(0) = (4.0_rt / 3.0_rt) * M_PI *
        (amrex::Math::powi<3>(rr(0)) - amrex::Math::powi<3>(rl(0))) *
        model::profile(0).state(0, model::idens);

    // this is the zone where we switch composition to He

    int ihe_layer{-1};

    // this will be the zone at the base of the He layer that
    // sets the entropy we will constrain to

    int ihe_entropy{-1};

    // this is the zone where the low density (fluff) begins

    int icutoff{-1};

    // HSE + entropy solve

    bool isentropic;

    Real entropy_base;

    amrex::Real dens_zone;
    amrex::Real temp_zone;
    amrex::Real pres_zone;
    amrex::Real entropy;
    amrex::Real xn[NumSpec];

    for (int i = 1; i < model::npts; ++i) {

        // as the initial guess for the density, use the previous
        // zone

        dens_zone = model::profile(0).state(i-1, model::idens);

        if (dens_zone > rho_he) {
            temp_zone = model_params.T_core;
            for (int n = 0; n < NumSpec; ++n) {
                xn[n] = model_params.xn_core[n];
            }

            isentropic = false;

        } else {

            if (ihe_layer == -1) {
                ihe_layer = i;
            }


            // determine whether we are starting the ramp up.  We
            // will use a tanh profile, centered at
            // (xzn_hse(ihe_layer) + FOUR*delta).  The "+
            // FOUR*delta" enables us to capture the leading edge
            // of the profile.  Since rho_he is computed by
            // considering the integral of He on the grid,
            // shifting the profile by FOUR*delta doesn't affect
            // the overall mass.

            Real test = 0.5_rt * (1.0_rt + std::tanh((model::profile(0).r(i) -
                                                      model::profile(0).r(ihe_layer) -
                                                      4.0_rt * model_params.delta) /
                                                     model_params.delta));

            if (test < 0.999_rt) {

                // small tanh ramp up regime

                for (int n = 0; n < NumSpec; ++n) {
                    xn[n] = model_params.xn_core[n] + test * (model_params.xn_he[n] - model_params.xn_core[n]);
                }

                temp_zone = model_params.T_core + test * (model_params.T_base - model_params.T_core);

                isentropic = false;

            } else {

                if (model_params.isothermal_layer) {
                    // isothermal He layer no matter what
                    temp_zone = model_params.T_base;
                    isentropic = false;
                } else {
                    // fully isentropic

                    if (ihe_entropy == -1) {
                        ihe_entropy = i;
                        temp_zone = model_params.T_base;
                        isentropic = false;
                    } else {
                        temp_zone = model::profile(0).state(i-1, model::itemp);
                        isentropic = true;
                    }
                }

                for (int n = 0; n < NumSpec; ++n) {
                    xn[n] = model_params.xn_he[n];
                }

            }

        }

        amrex::Real g_zone = -C::Gconst * M_enclosed(i-1) / (rl(i) * rl(i));


        // thermodynamic state iteration loop

        // start off the Newton loop by saying that the zone has not converged
        bool converged_hse = false;

        if (! fluff) {

            amrex::Real p_want;
            amrex::Real drho;
            amrex::Real dtemp;

            for (int iter = 0; iter < subch::MAX_ITER; ++iter) {

                if (isentropic) {

                    p_want = model::profile(0).state(i-1, model::ipres) +
                        dx * 0.5_rt * (dens_zone + model::profile(0).state(i-1, model::idens)) * g_zone;

                    // now we have two functions to zero:
                    //   A = p_want - p(rho,T)
                    //   B = entropy_base - s(rho,T)
                    // We use a two dimensional Taylor expansion
                    // and find the deltas for both density and
                    // temperature

                    eos_state.T = temp_zone;
                    eos_state.rho = dens_zone;
                    for (int n = 0; n < NumSpec; ++n) {
                        eos_state.xn[n] = xn[n];
                    }

                    // (t, rho) -> (p, s)
                    eos(eos_input_rt, eos_state);

                    entropy = eos_state.s;
                    pres_zone = eos_state.p;

                    amrex::Real dpT = eos_state.dpdT;
                    amrex::Real dpd = eos_state.dpdr;
                    amrex::Real dsT = eos_state.dsdT;
                    amrex::Real dsd = eos_state.dsdr;

                    amrex::Real A = p_want - pres_zone;
                    amrex::Real B = entropy_base - entropy;

                    amrex::Real dAdT = -dpT;
                    amrex::Real dAdrho = 0.5_rt * dx * g_zone - dpd;
                    amrex::Real dBdT = -dsT;
                    amrex::Real dBdrho = -dsd;

                    dtemp = (B - (dBdrho / dAdrho) * A) /
                        ((dBdrho / dAdrho) * dAdT - dBdT);

                    drho = -(A + dAdT * dtemp) / dAdrho;

                    dens_zone = amrex::Clamp(dens_zone + drho,
                                             0.9_rt * dens_zone, 1.1_rt * dens_zone);

                    temp_zone = amrex::Clamp(temp_zone + dtemp,
                                             0.9_rt * temp_zone, 1.1_rt * temp_zone);

                    // check if the density falls below our minimum
                    // cut-off -- if so, floor it

                    if (dens_zone < model_params.low_density_cutoff) {
                        dens_zone = model_params.low_density_cutoff;
                        temp_zone = model_params.T_fluff;
                        converged_hse = true;
                        fluff = true;
                        break;
                    }

                    if (std::abs(drho) < model_params.tol_hse * dens_zone &&
                        std::abs(dtemp) < model_params.tol_hse * temp_zone) {
                        converged_hse = true;
                        break;
                    }

                } else {

                    // the core is isothermal, so we just need to
                    // constrain the density and pressure to agree
                    // with the EOS and HSE

                    // We difference HSE about the interface
                    // between the current zone and the one just
                    // inside.

                    p_want = model::profile(0).state(i-1, model::ipres) +
                        dx * 0.5_rt * (dens_zone + model::profile(0).state(i-1, model::idens)) * g_zone;

                    eos_state.T = temp_zone;
                    eos_state.rho = dens_zone;
                    for (int n = 0; n < NumSpec; ++n) {
                        eos_state.xn[n] = xn[n];
                    }

                    // (t, rho) -> (p, s)

                    eos(eos_input_rt, eos_state);

                    entropy = eos_state.s;
                    pres_zone = eos_state.p;

                    amrex::Real dpd = eos_state.dpdr;

                    drho = (p_want - pres_zone) / (dpd - 0.5_rt * dx * g_zone);

                    dens_zone = amrex::Clamp(dens_zone + drho,
                                             0.9_rt * dens_zone, 1.1_rt * dens_zone);

                    if (std::abs(drho) < model_params.tol_hse * dens_zone) {
                        converged_hse = true;
                        break;
                    }

                    if (dens_zone < model_params.low_density_cutoff) {

                        icutoff = i;
                        dens_zone = model_params.low_density_cutoff;
                        temp_zone = model_params.T_fluff;
                        converged_hse = true;
                        fluff = true;
                        break;
                    }
                }

                if (temp_zone < model_params.T_fluff && isentropic) {
                    temp_zone = model_params.T_fluff;
                    isentropic = false;
                }

            }  // thermo iteration loop


            if (! converged_hse) {

                std::cout << "Error zone " << i <<  " did not converge in init_1d" << std::endl;
                std::cout << dens_zone << " " << temp_zone << std::endl;
                std::cout << p_want << std::endl;
                std::cout << drho << std::endl;
                amrex::Error("Error: HSE non-convergence");
            }

        } else {
            // fluff region
            dens_zone = model_params.low_density_cutoff;
            temp_zone = model_params.T_fluff;
        }

        // call the EOS one more time for this zone and then go on
        // to the next

        eos_state.T = temp_zone;
        eos_state.rho = dens_zone;
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = xn[n];
        }

        // (t, rho) -> (p, s)

        eos(eos_input_rt, eos_state);

        pres_zone = eos_state.p;

        // determine the entropy that we want to constrain to, if
        // this is the first zone of the He layer

        if (i == ihe_entropy) {
            entropy_base = entropy;
        }

        // update the thermodynamics in this zone

        model::profile(0).state(i, model::idens) = dens_zone;
        model::profile(0).state(i, model::itemp) = temp_zone;
        model::profile(0).state(i, model::ipres) = pres_zone;

        for (int n = 0; n < NumSpec; ++n) {
            model::profile(0).state(i, model::ispec+n) = xn[n];
        }

        M_enclosed(i) = M_enclosed(i-1) +
            (4.0_rt / 3.0_rt) * M_PI * (rr(i) - rl(i)) *
            (amrex::Math::powi<2>(rr(i)) + rl(i) * rr(i) + amrex::Math::powi<2>(rl(i))) * model::profile(0).state(i, model::idens);

    } // end loop over zones

    amrex::Real mass_wd{};
    amrex::Real mass_he{};

    // it might be that we never reach the cutoff density in our
    // domain.  This is especially the case if we do an isothermal
    // model.  Make sure we integrate over everything in that
    // case.
    int max_index = icutoff == -1 ? model::npts : icutoff;

    for (int i = 0; i < max_index; ++i) {

        Real vol{0.0};
        if (i == 0) {
            vol = (4.0_rt / 3.0_rt) * M_PI * (amrex::Math::powi<3>(rr(0)) - amrex::Math::powi<3>(rl(0)));
        } else {
            vol = (4.0_rt / 3.0_rt) * M_PI *
                (rr(i) - rl(i)) * (rr(i) * rr(i) + rl(i) * rr(i) + rl(i) * rl(i));
        }

        // compute the total mass of the He layer and C/O WD

        // note: some He layers can optionally be enriched with C/O, so
        // only count that C and O if we are in the He layer

        amrex::Real layer_X = model::profile(0).state(i, model::ispec+model_params.ihe4);
        if (model_params.X_N14 > 0.0) {
            layer_X += model::profile(0).state(i, model::ispec+model_params.in14);
        }

        if (model_params.X_C12 > 0.0 && layer_X > 10.0 * network_rp::small_x) {
            layer_X += model::profile(0).state(i, model::ispec+model_params.ic12);
        }

        if (model_params.X_O16 > 0.0 && layer_X > 10.0 * network_rp::small_x) {
            layer_X += model::profile(0).state(i, model::ispec+model_params.io16);
        }

        Real core_X{0.0};
        if ((model_params.X_C12 == 0.0 && model_params.X_O16 == 0.0) ||
            (model::profile(0).state(i, model::ispec+model_params.ihe4) <=
             network_rp::small_x)) {
            core_X += model::profile(0).state(i, model::ispec+model_params.ic12) +
                model::profile(0).state(i, model::ispec+model_params.io16);
        }

        mass_he += vol * model::profile(0).state(i, model::idens) * layer_X;
        mass_wd += vol * model::profile(0).state(i, model::idens) * core_X;
    }

    // compute the maximum HSE error

    amrex::Real max_hse_error = -1.e30_rt;

    for (int i = 1; i < model::npts-1; ++i) {
        amrex::Real g_zone = -C::Gconst * M_enclosed(i-1) / (rr(i-1) * rr(i-1));
        amrex::Real dpdr = (model::profile(0).state(i, model::ipres) -
                            model::profile(0).state(i-1, model::ipres)) / dx;
        amrex::Real rhog = 0.5_rt * (model::profile(0).state(i, model::idens) +
                                      model::profile(0).state(i-1, model::idens)) * g_zone;

        if (dpdr != 0.0_rt &&
            model::profile(0).state(i+1, model::idens) >
            model_params.low_density_cutoff) {
            max_hse_error = std::max(max_hse_error, std::abs(dpdr - rhog) / std::abs(dpdr));
        }
    }

    return {mass_wd, mass_he, max_hse_error};

}

// generate an initial model for an arbitrary-mass, isothermal C WD
// with an isentropic He envelope on the surface.

AMREX_INLINE
void
generate_initial_model(const int npts_model, const Real xmin, const Real xmax,
                       const model_t& model_params)
{

    // Create a 1-d uniform grid that is identical to the mesh that we
    // are mapping onto, and then we want to force it into HSE on that
    // mesh.

    // we are doing the full star
    BL_ASSERT(xmin == 0.0);

    model::npts = npts_model;
    model::initialized = true;

    if (npts_model > NPTS_MODEL) {
        amrex::Error("too many zones requested -- increase NPTS_MODEL");
    }

    amrex::Real dx = xmax / static_cast<amrex::Real>(npts_model);

    amrex::Print() << Font::Bold << FGColor::Green;
    amrex::Print() << "generating initial model with " << npts_model << " points and dx = " << dx << " cm" << std::endl;
    amrex::Print() << ResetDisplay;

    // compute the coordinates of the new gridded function

    for (int i = 0; i < npts_model; i++) {
        model::profile(0).r(i) = (static_cast<Real>(i) + 0.5_rt) * dx;
    }


    // We don't know what WD central density will give the desired
    // total mass, so we need to iterate over central density

    // we will do a secant iteration.  rho_c_old is the 'old' guess
    // for the central density and rho_c is the current guess.  After
    // 2 loops, we can start estimating the density required to yield
    // our desired mass

    amrex::Real rho_c_old = -1.0_rt;
    amrex::Real rho_c = 1.e8_rt;  // 1.e8 is a reasonable starting WD central density for a 1 M_sun WD

    // rho_he_old is the old guess for the density to transition to He,
    // where we will be isentropic, and rho_he is the current guess.

    amrex::Real rho_he_old = -1.0_rt;
    amrex::Real rho_he = 0.1_rt * rho_c;

    bool mass_converged = false;

    amrex::Real mass_wd;
    amrex::Real mass_wd_old;
    amrex::Real mass_he;
    amrex::Real mass_he_old;
    amrex::Real max_hse_error;

    for (int iter_mass = 0; iter_mass < subch::MAX_ITER; ++iter_mass) {

        //amrex::Print() << "mass iter = " << iter_mass << " " << rho_c  << std::endl;

        std::tie(mass_wd, mass_he, max_hse_error) = build_star(model_params, dx, rho_c, rho_he);

        if (rho_c_old < 0.0_rt) {
            // not enough iterations yet -- store the old central
            // density and mass and pick a new value

            rho_c_old = rho_c;
            mass_wd_old = mass_wd;

            rho_he_old = rho_he;
            mass_he_old = mass_he;

            rho_c = 1.01_rt * rho_c_old;
            rho_he = 0.95_rt * rho_he_old;

        } else {
            // have we converged

            amrex::Real wd_err = std::abs(mass_wd - model_params.M_tot) / model_params.M_tot;
            amrex::Real he_err = std::abs(mass_he - model_params.M_He) / model_params.M_He;
            if (wd_err < subch::TOL_WD_MASS && he_err < subch::TOL_HE_MASS) {
                mass_converged = true;
                break;
            }

            // do a secant iteration:

            // we'll define 2 expressions:
            //   F = M_wd(rho_c, rho_he) - M_wd_want
            //   G = M_he(rho_c, rho_he) - M_he_want
            //
            // and then do a 2 variable Newton solve

            Real F0 = mass_wd - model_params.M_tot;
            Real G0 = mass_he - model_params.M_He;

            Real dF_drhoc = (mass_wd - mass_wd_old) / (rho_c - rho_c_old);
            Real dF_drhohe = 0.0; // (mass_wd - mass_wd_old) / (rho_he - rho_he_old);

            Real dG_drhoc = 0.0; // (mass_he - mass_he_old) / (rho_c - rho_c_old);
            Real dG_drhohe = (mass_he - mass_he_old) / (rho_he - rho_he_old);

            Real drho_he = -(G0 - dG_drhoc / dF_drhoc * F0) / (dG_drhohe - dG_drhoc * dF_drhohe / dF_drhoc);
            Real drho_c = -(F0 + dF_drhohe * drho_he) / dF_drhoc;

            rho_c_old = rho_c;
            mass_wd_old = mass_wd;

            rho_c = amrex::Clamp(rho_c + drho_c,
                                 0.9_rt * rho_c_old, 1.1_rt * rho_c_old);

            rho_he_old = rho_he;
            mass_he_old = mass_he;

            rho_he = amrex::Clamp(rho_he + drho_he,
                                  0.75_rt * rho_he_old, 1.5_rt * rho_he_old);

        }

    } // end mass constraint loop

    if (! mass_converged) {
        amrex::Error("ERROR: mass did not converge");
    }

    amrex::Print() << Font::Bold << FGColor::Green;
    amrex::Print() << "model generation converged" << std::endl;
    amrex::Print() << ResetDisplay;
    amrex::Print() << "central density of WD: " << rho_c << std::endl;
    amrex::Print() << "density at base of He layer: " << rho_he << std::endl;


    amrex::Print() << Font::Bold << FGColor::Green;
    amrex::Print() << "final masses: ";
    amrex::Print() << ResetDisplay << std::endl;
    amrex::Print() << " mass WD: " << mass_wd / C::M_solar << std::endl;
    amrex::Print() << " mass He: " << mass_he / C::M_solar << std::endl;

    amrex::Print() << Font::Bold << FGColor::Green;
    amrex::Print() << "maximum HSE error = " << max_hse_error << std::endl;
    amrex::Print() << ResetDisplay << std::endl;

}


#endif
