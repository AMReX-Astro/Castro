#ifndef INITIAL_MODEL_H
#define INITIAL_MODEL_H

//  Create a 1-d hydrostatic, atmosphere with an isothermal region
//  (T_star) representing the NS, a hyperbolic tangent rise to a
//  peak temperature (T_hi) representing the base of an accreted
//  layer, an isoentropic profile down to a lower temperature (T_lo),
//  and then isothermal. This can serve as an initial model for a
//  nova or XRB.
//
//  The temperature profile is:
//
//         ^
//  T_hi   +             ^                                     .
//         |            / \                                    .
//         |           /   \                                   .
//         |          /  .  \                                  .
//  T_star +---------+       \                                 .
//         |         .   .    \                                .
//         |                   \                               .
//         |         .   .      \                              .
//  T_lo   +                     +-----------                  .
//         |         .   .                                     .
//         +---------+---+---------------> r                   .
//         |         \  /
//         |       atm_delta
//         |< H_star>|
//
//
//                   ^
//                   |
//                   +-- dens_base
//
//  dens_base is the density at a height H_star -- just below the rise
//  in T up to the peak T_hi.  The composition is "ash" in the lower
//  isothermal region and "fuel" in the isentropic and upper
//  isothermal regions.  In the transition region, we apply the same
//  hyperbolic tangent profile to interpolate the composition.
//
//  The fuel and ash compositions are specified by the fuel?_name,
//  fuel?_frac and ash?_name, ash?_frac parameters (name of the species
//  and mass fraction).  Where ? = 1,2,3.
//
//  The model is placed into HSE by the following differencing:
//
//   (1/dr) [ <P>_i - <P>_{i-1} ] = (1/2) [ <rho>_i + <rho>_{i-1} ] g
//
//  This will be iterated over in tandem with the EOS call,
//  P(i-1) = P_eos(rho(i-1), T(i-1), X(i-1)
//

// this version allows for multiple initial models

namespace fw {
    constexpr Real MAX_ITER = 250;
    constexpr Real TOL = 1.e-10_rt;
}

struct model_t {

    Real xn_base[NumSpec];
    Real xn_star[NumSpec];
    Real xn_perturb[NumSpec];

    Real dens_base;
    Real T_star;
    Real T_hi;
    Real T_lo;

    Real H_star;
    Real atm_delta;

    Real low_density_cutoff;
};


// Evaluate tanh using the exponential form to workaround a PGI bug on Power9

AMREX_INLINE
Real evaluate_tanh(const Real z) {

    Real t;
    if (std::abs(z) <= 4.0_rt) {
        t = (std::exp(z) - std::exp(-z))/(std::exp(z) + std::exp(-z));
    } else if (z < -4.0_rt) {
        t = -1.0_rt;
    } else {
        t = 1.0_rt;
    }

    return t;
}


AMREX_INLINE
void
generate_initial_model(const int npts_model, const Real xmin, const Real xmax,
                       const model_t model_params, const int model_num)
{


    // Create a 1-d uniform grid that is identical to the mesh that we are
    // mapping onto, and then we want to force it into HSE on that mesh.

    // we actually require that the number of points is the same for each
    // model, so we'll just set it each time

    model::npts = npts_model;
    model::initialized = true;

    if (npts_model > NPTS_MODEL) {
        amrex::Error("Error: model has more than NPTS_MODEL points,  Increase MAX_NPTS_MODEL");
    }

    // create the grid -- cell centers

    Real dx = (xmax - xmin) / npts_model;

    for (int i = 0; i < npts_model; i++) {
        model::profile(model_num).r(i) =
            xmin + (static_cast<Real>(i) + 0.5_rt) * dx;
    }


    // find the index of the base height

    int index_base = -1;
    for (int i = 0; i < npts_model; i++) {
        if (model::profile(model_num).r(i) >= xmin + model_params.H_star) {
            index_base = i+1;
            break;
        }
    }

    if (index_base == -1) {
        amrex::Error("ERROR: invalid base_height");
    }


    // put the model onto our new uniform grid

    bool fluff = false;

    // determine the conditions at the base -- this is below the atmosphere

    eos_t eos_state;
    eos_state.T = model_params.T_star;
    eos_state.rho = model_params.dens_base;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = model_params.xn_star[n];
    }

    eos(eos_input_rt, eos_state);

    // store the conditions at the base -- we'll use the entropy later
    // to constrain the isentropic layer

    Real pres_base = eos_state.p;

    // set an initial temperature profile and composition

    for (int i = 0; i < npts_model; i++) {

        Real xc = model::profile(model_num).r(i) -
            (xmin + model_params.H_star) - 1.5_rt * model_params.atm_delta;

        // hyperbolic tangent transition:

        for (int n = 0; n < NumSpec; n++) {
            model::profile(model_num).state(i, model::ispec+n) =
                model_params.xn_star[n] +
                0.5_rt * (model_params.xn_base[n] - model_params.xn_star[n]) *
                (1.0_rt + evaluate_tanh(xc / (0.5_rt * model_params.atm_delta)));
        }

        // force them to sum to 1

        Real sumX = 0.0_rt;
        for (int n = 0; n < NumSpec; n++) {
            sumX += model::profile(model_num).state(i, model::ispec+n);
        }
        for (int n = 0; n < NumSpec; n++) {
            model::profile(model_num).state(i, model::ispec+n) /= sumX;
        }

        // temperature profile -- it is constant below the base

        if (i <= index_base) {
            model::profile(model_num).state(i, model::itemp) = model_params.T_star;
        } else {
            model::profile(model_num).state(i, model::itemp) =
                model_params.T_star +
                0.5_rt * (model_params.T_hi - model_params.T_star) *
                (1.0_rt + evaluate_tanh(xc / (0.5_rt * model_params.atm_delta)));
        }

        // the density and pressure will be determined via HSE,
        // for now, set them to the base conditions

        model::profile(model_num).state(i, model::idens) = model_params.dens_base;
        model::profile(model_num).state(i, model::ipres) = pres_base;

    }

    // make the base thermodynamics consistent for this base point -- that is
    // what we will integrate from!

    eos_state.rho = model::profile(model_num).state(index_base, model::idens);
    eos_state.T = model::profile(model_num).state(index_base, model::itemp);
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = model::profile(model_num).state(index_base, model::ispec+n);
    }

    eos(eos_input_rt, eos_state);

    model::profile(model_num).state(index_base, model::ipres) = eos_state.p;


    //
    // HSE + entropy solve
    //

    // the HSE state will be done putting creating an isentropic state until
    // the temperature goes below T_lo -- then we will do isothermal.
    // also, once the density goes below low_density_cutoff, we stop HSE

    bool isentropic = false;
    bool flipped = false;   // we start out isothermal and then 'flip' to isentropic

    //
    // integrate upward
    //

    Real entropy_base;

    for (int i = index_base+1; i < npts_model; i++) {

       if ((model::profile(model_num).r(i) > xmin + model_params.H_star + 3.0_rt * model_params.atm_delta) && !flipped) {
           if (i == index_base + 1) {
               amrex::Error("atm_delta is too small for the grid spacing (at least one transition zone is required)");
           }
           isentropic = true;
           flipped = true;

           // now we need to know the entropy we are confining ourselves to

           eos_state.rho = model::profile(model_num).state(i-1, model::idens);
           eos_state.T = model::profile(model_num).state(i-1, model::itemp);
           for (int n = 0; n < NumSpec; n++) {
               eos_state.xn[n] = model::profile(model_num).state(i-1, model::ispec+n);
           }

           eos(eos_input_rt, eos_state);

           entropy_base = eos_state.s;

           amrex::Print() << "base density = " << eos_state.rho << " " << eos_state.T << std::endl;
       }

       // we've already set initial guesses for density, temperature, and
       // composition

       Real dens_zone = model::profile(model_num).state(i, model::idens);
       Real temp_zone = model::profile(model_num).state(i, model::itemp);
       Real xn[NumSpec];
       for (int n = 0; n < NumSpec; n++) {
           xn[n] = model::profile(model_num).state(i, model::ispec+n);
       }


       //
       // iteration loop
       //

       // start off the Newton loop by saying that the zone has not converged

       bool converged_hse = false;

       Real pres_zone;
       Real entropy;

       if (!fluff) {

           Real p_want;
           Real drho;
           Real dtemp = 0;

           for (int iter = 0; iter < fw::MAX_ITER; iter++) {

               // get the pressure we want from the HSE equation, just the
               // zone below the current.  Note, we are using an average of
               // the density of the two zones as an approximation of the
               // interface value -- this means that we need to iterate for
               // find the density and pressure that are consistent

               // furthermore, we need to get the entropy that we need,
               // which will come from adjusting the temperature in
               // addition to the density.

               // HSE differencing

               p_want = model::profile(model_num).state(i-1, model::ipres) +
                   dx * 0.5_rt * (dens_zone + model::profile(model_num).state(i-1, model::idens)) * gravity::const_grav;

               if (isentropic) {

                   // now we have two functions to zero:
                   //   A = p_want - p(rho,T)
                   //   B = entropy_base - s(rho,T)
                   // We use a two dimensional Taylor expansion and find the deltas
                   // for both density and temperature

                   // now we know the pressure and the entropy that we want, so we
                   // need to find the temperature and density through a two
                   // dimensional root find

                   // (t, rho) -> (p, s)

                   eos_state.T = temp_zone;
                   eos_state.rho = dens_zone;
                   for (int n = 0; n < NumSpec; n++) {
                       eos_state.xn[n] = xn[n];
                   }

                   eos(eos_input_rt, eos_state);

                   entropy = eos_state.s;
                   pres_zone = eos_state.p;

                   Real dpt = eos_state.dpdT;
                   Real dpd = eos_state.dpdr;
                   Real dst = eos_state.dsdT;
                   Real dsd = eos_state.dsdr;

                   Real A = p_want - pres_zone;
                   Real B = entropy_base - entropy;

                   dtemp = ((dsd / (dpd - 0.5_rt * dx * gravity::const_grav)) * A - B) /
                       (dsd * dpt / (dpd - 0.5_rt * dx * gravity::const_grav) - dst);

                   drho = (A - dpt * dtemp) / (dpd - 0.5_rt * dx * gravity::const_grav);

                   dens_zone = amrex::max(0.9_rt * dens_zone,
                                          amrex::min(dens_zone + drho, 1.1_rt * dens_zone));

                   temp_zone = amrex::max(0.9_rt * temp_zone,
                                          amrex::min(temp_zone + dtemp, 1.1_rt * temp_zone));

                   // check if the density falls below our minimum cut-off --
                   // if so, floor it

                   if (dens_zone < model_params.low_density_cutoff) {

                       dens_zone = model_params.low_density_cutoff;
                       temp_zone = model_params.T_lo;
                       converged_hse = true;
                       fluff = true;
                       break;
                   }

                   // if (A < TOL .and. B < ETOL) then
                   if (std::abs(drho) < fw::TOL * dens_zone && std::abs(dtemp) < fw::TOL * temp_zone) {
                       converged_hse = true;
                       break;
                   }

               } else {

                   // do isothermal

                   if (model::profile(model_num).r(i) > xmin + model_params.H_star + 3.0_rt * model_params.atm_delta) {
                       temp_zone = model_params.T_lo;
                   }

                   // (t, rho) -> (p)

                   eos_state.T = temp_zone;
                   eos_state.rho = dens_zone;
                   for (int n = 0; n < NumSpec; n++) {
                       eos_state.xn[n] = xn[n];
                   }

                   eos(eos_input_rt, eos_state);

                   entropy = eos_state.s;
                   pres_zone = eos_state.p;

                   Real dpd = eos_state.dpdr;

                   drho = (p_want - pres_zone) / (dpd - 0.5_rt * dx * gravity::const_grav);

                   dens_zone = amrex::max(0.9_rt * dens_zone,
                                          amrex::min(dens_zone + drho, 1.1_rt * dens_zone));

                   if (std::abs(drho) < fw::TOL * dens_zone) {
                       converged_hse = true;
                       break;
                   }

                   if (dens_zone < model_params.low_density_cutoff) {
                       dens_zone = model_params.low_density_cutoff;
                       temp_zone = model_params.T_lo;
                       converged_hse = true;
                       fluff = true;
                       break;
                   }

               }

               if (temp_zone < model_params.T_lo) {
                   temp_zone = model_params.T_lo;
                   isentropic = false;
               }

           }

           if (!converged_hse) {
               std::cout << "Error zone " << i << " did not converge in init_1d" << std::endl;
               std::cout << "integrate up" << std::endl;
               std::cout << dens_zone << " " << temp_zone << std::endl;
               std::cout << p_want << " " << entropy_base << " " << entropy << std::endl;
               std::cout << drho << " " << dtemp << std::endl;
               amrex::Error("Error: HSE non-convergence");
          }

       } else {
           // fluff

           dens_zone = model_params.low_density_cutoff;
           temp_zone = model_params.T_lo;
       }


       // call the EOS one more time for this zone and then go on to the next
       // (t, rho) -> (p)

       eos_state.T = temp_zone;
       eos_state.rho = dens_zone;
       for (int n = 0; n < NumSpec; n++) {
           eos_state.xn[n] = xn[n];
       }

       eos(eos_input_rt, eos_state);

       pres_zone = eos_state.p;

       // update the thermodynamics in this zone

       model::profile(model_num).state(i, model::idens) = dens_zone;
       model::profile(model_num).state(i, model::itemp) = temp_zone;
       model::profile(model_num).state(i, model::ipres) = pres_zone;

       // to make this process converge faster, set the density in the
       // next zone to the density in this zone
       // model::profile(model_num).state(i+1, model::idens) = dens_zone;

    }


    //
    // integrate down -- using the temperature profile defined above
    //

    for (int i = index_base-1; i >= 0; --i) {

        // we already set the temperature and composition profiles

        Real temp_zone = model::profile(model_num).state(i, model::itemp);
        Real xn[NumSpec];
        for (int n = 0; n < NumSpec; n++) {
            xn[n] = model::profile(model_num).state(i, model::ispec+n);
        }

        // use our previous initial guess for density

        Real dens_zone = model::profile(model_num).state(i+1, model::idens);

        //
        // iteration loop
        //

        // start off the Newton loop by saying that the zone has not converged

        bool converged_hse = false;

        Real pres_zone;
        Real p_want;
        Real drho;

        for (int iter = 0; iter < fw::MAX_ITER; iter++) {

            // get the pressure we want from the HSE equation, just the
            // zone above the current.  Note, we are using an average of
            // the density of the two zones as an approximation of the
            // interface value -- this means that we need to iterate for
            // find the density and pressure that are consistent

            // HSE differencing

            p_want = model::profile(model_num).state(i+1, model::ipres) -
                dx * 0.5_rt * (dens_zone + model::profile(model_num).state(i+1, model::idens)) * gravity::const_grav;

            // we will take the temperature already defined in gen_model_state
            // so we only need to zero:
            //   A = p_want - p(rho)

            // (t, rho) -> (p)
            eos_state.T = temp_zone;
            eos_state.rho = dens_zone;
            for (int n = 0; n < NumSpec; n++) {
                eos_state.xn[n] = xn[n];
            }

            eos(eos_input_rt, eos_state);

            pres_zone = eos_state.p;

            Real dpd = eos_state.dpdr;

            Real A = p_want - pres_zone;

            drho = A / (dpd + 0.5_rt * dx * gravity::const_grav);

            dens_zone = amrex::max(0.9_rt * dens_zone,
                                   amrex::min(dens_zone + drho, 1.1_rt * dens_zone));

            if (std::abs(drho) < fw::TOL * dens_zone) {
                converged_hse = true;
                break;
            }

        }

        if (!converged_hse) {
            std::cout << "Error zone " << i << " did not converge in init_1d" << std::endl;
            std::cout << "integrate down" << std::endl;
            std::cout << dens_zone << " " << temp_zone << std::endl;
            std::cout << p_want << std::endl;
            std::cout << drho << std::endl;
            amrex::Error("Error: HSE non-convergence");
        }


        // call the EOS one more time for this zone and then go on to the next
        // (t, rho) -> (p)

        eos_state.T = temp_zone;
        eos_state.rho = dens_zone;
        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = xn[n];
        }

        eos(eos_input_rt, eos_state);

        pres_zone = eos_state.p;

        // update the thermodynamics in this zone

        model::profile(model_num).state(i, model::idens) = dens_zone;
        model::profile(model_num).state(i, model::itemp) = temp_zone;
        model::profile(model_num).state(i, model::ipres) = pres_zone;

    }
}
#endif
