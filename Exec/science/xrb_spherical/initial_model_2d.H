#ifndef INITIAL_MODEL_H
#define INITIAL_MODEL_H

#include <AMReX_REAL.H>

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
    constexpr amrex::Real MAX_ITER = 250;
    constexpr amrex::Real TOL = 1.e-10_rt;
}

struct model_t {

    amrex::Real xn_base[NumSpec];
    amrex::Real xn_star[NumSpec];
    amrex::Real xn_perturb[NumSpec];

    amrex::Real dens_base;
    amrex::Real T_star;
    amrex::Real T_hi;
    amrex::Real T_lo;

    amrex::Real H_star;
    amrex::Real atm_delta;

    amrex::Real low_density_cutoff;
};


// Evaluate tanh using the exponential form to workaround a PGI bug on Power9

AMREX_INLINE
amrex::Real evaluate_tanh(const amrex::Real z) {

    amrex::Real t;
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
generate_initial_model(const int npts_model_x, const amrex::Real xmin, const amrex::Real xmax,
                       const int npts_model_y, const amrex::Real ymin, const amrex::Real
                       const model_t& model_params, const int model_num)
{


    // Create a 2-d uniform grid that is identical to the mesh that we are
    // mapping onto, and then we want to force it into HSE on that mesh.

    // we actually require that the number of points is the same for each
    // model, so we'll just set it each time

    model::npts = npts_model_x;
    model::initialized = true;

    if (npts_model_x > NPTS_MODEL) {
        amrex::Error("Error: model has more than NPTS_MODEL points in x-dir,  Increase MAX_NPTS_MODEL");
    }
    if (npts_model_y > NPTS_MODEL) {
        amrex::Error("Error: model has more than NPTS_MODEL points in y-dir,  Increase MAX_NPTS_MODEL");
    }

    // create the grid -- cell centers

    amrex::Real dx = (xmax - xmin) / static_cast<amrex::amrex::Real>(npts_model_x);
    amrex::Real dy = (ymax - ymin) / static_cast<amrex::amrex::Real>(npts_model_y);

    for (int i = 0; i < npts_model_x; i++) {
        model::profile(model_num).x(i) =
            xmin + (static_cast<amrex::Real>(i) + 0.5_rt) * dx;
    }

    for (int j = 0; j < npts_model_y; j++) {
        model::profile(model_num).y(j) =
            ymin + (static_cast<amrex::Real>(j) + 0.5_rt) * dy;
    }

    // find the index of the base height
    // x is assumed to be the r-dir

    int index_base = -1;
    for (int i = 0; i < npts_model_x; i++) {
        if (model::profile(model_num).x(i) >= xmin + model_params.H_star) {
            index_base = i+1;
            break;
        }
    }

    if (index_base == -1) {
        amrex::Error("ERROR: invalid base_height");
    }

    //
    // put the model onto our new uniform grid
    // This is the initial guess conditions.
    //

    // determine the conditions at the base -- this is below the atmosphere
    // this point is assumed to be at the rotating pole, θ=0

    eos_t eos_state;
    eos_state.T = model_params.T_star;
    eos_state.rho = model_params.dens_base;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = model_params.xn_star[n];
    }

    eos(eos_input_rt, eos_state);

    // store the conditions at the base -- we'll use the entropy later
    // to constrain the isentropic layer

    amrex::Real pres_base = eos_state.p;

    // set an initial temperature profile and composition for the 2D grid

    for (int i = 0; i < npts_model_x; i++) {

        amrex::Real xc = model::profile(model_num).x(i) -
            (xmin + model_params.H_star) - 1.5_rt * model_params.atm_delta;

        for (int j = 0; j < npts_model_y; j++) {

            // hyperbolic tangent transition:

            for (int n = 0; n < NumSpec; n++) {
                model::profile(model_num).state(i, j, model::ispec+n) =
                    model_params.xn_star[n] +
                    0.5_rt * (model_params.xn_base[n] - model_params.xn_star[n]) *
                    (1.0_rt + evaluate_tanh(xc / (0.5_rt * model_params.atm_delta)));
            }

            // force them to sum to 1

            amrex::Real sumX = 0.0_rt;
            for (int n = 0; n < NumSpec; n++) {
                sumX += model::profile(model_num).state(i, j, model::ispec+n);
            }
            for (int n = 0; n < NumSpec; n++) {
                model::profile(model_num).state(i, j, model::ispec+n) /= sumX;
            }

            // temperature profile -- it is constant below the base
            // It follows the hyperbolic tangent transition above the base

            if (i <= index_base) {
                model::profile(model_num).state(i, j, model::itemp) = model_params.T_star;
            } else {
                model::profile(model_num).state(i, j, model::itemp) =
                    model_params.T_star +
                    0.5_rt * (model_params.T_hi - model_params.T_star) *
                    (1.0_rt + evaluate_tanh(xc / (0.5_rt * model_params.atm_delta)));
            }

            //
            // We will also have a hyperbolic tangent transition in the theta dir
            //



            // the density and pressure will be determined via HSE,
            // for now, set them to the base conditions

            model::profile(model_num).state(i, j, model::idens) = model_params.dens_base;
            model::profile(model_num).state(i, j, model::ipres) = pres_base;
        }
    }

    // make the base thermodynamics consistent for this base point,
    // this is assumed to be on θ=0 -- that is what we will integrate from!

    eos_state.rho = model::profile(model_num).state(index_base, 0, model::idens);
    eos_state.T = model::profile(model_num).state(index_base, 0, model::itemp);
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = model::profile(model_num).state(index_base, 0, model::ispec+n);
    }

    eos(eos_input_rt, eos_state);

    model::profile(model_num).state(index_base, 0, model::ipres) = eos_state.p;


    //
    // 2D HSE + entropy solve
    //

    bool fluff = false;
    bool isentropic = false;
    bool flipped = false;

    // Start from the base layer at θ=0 and integrate upward
    // We start out isothermal during the hyperbolic transition region
    // then we 'flip' to isentropic
    // the HSE state will be done putting creating an isentropic state until
    // the temperature goes below T_lo -- then we will do isothermal.
    // also, once the density goes below low_density_cutoff, we stop HSE

    amrex::Real entropy_base = -1.0_rt;

    for (int i = index_base+1; i < npts_model_x; i++) {
        if (model::profile(model_num).r(i) > xmin + model_params.H_star + 3.0_rt * model_params.atm_delta &&
            !flipped) {
            // This is to determine when did we get out of the hyperbolic tangent transition phase (isothermal)
            // Then get the appropriate entropy for doing isentropic atmosphere

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

            // Get entropy for isentropic atmosphere

            entropy_base = eos_state.s;

            amrex::Print() << "base density = " << eos_state.rho
                           << " and base Temperature " << eos_state.T << std::endl;
        }

        HSE_integrate(i, 0);
    }


    // Start from the base layer at θ=0 and integrate upward and integrate to the right

    // Now we go through the intersection regions, and integrate upward.
    // Since the integrated result should be the same regardless of the direction,
    // we do an average from both direction

    // Start from the base layer at θ=0 and integrate downward

    // Now go through the intersection regions and integrate downward.

}

AMREX_INLINE
void HSE_integration(const int i, const int j,
                     const amrex::Real dx, const amrex::Real dy,
                     const int idir,
                     const bool isentropic, bool& fluff,
                     const amrex::Real entropy_base,
                     const int index_base,
                     const model_t& model_params, const int model_num) {
    //
    // Goal is to update density, temperature, and pressure of the current zone.
    //

    // Get current state, this is the initial guess

    amrex::Real dens_zone = model::profile(model_num).state(i, j, model::idens);
    amrex::Real temp_zone = model::profile(model_num).state(i, j, model::itemp);
    amrex::Real xn[NumSpec];
    for (int n = 0; n < NumSpec; n++) {
        xn[n] = model::profile(model_num).state(i, model::ispec+n);
    }

    if (fluff) {
        // If in fluff region,
        // i.e. density of the zone is less than low_density_cutoff
        // We stop HSE in this zone and set to a constant density and Temp.

        dens_zone = model_params.low_density_cutoff;
        temp_zone = model_params.T_lo;

        // Call eos to get pressure: (T, rho) -> (p)

       eos_state.T = temp_zone;
       eos_state.rho = dens_zone;
       for (int n = 0; n < NumSpec; n++) {
           eos_state.xn[n] = xn[n];
       }

       eos(eos_input_rt, eos_state);

       pres_zone = eos_state.p;

       // update the thermodynamics in this zone

       model::profile(model_num).state(i, j, model::idens) = dens_zone;
       model::profile(model_num).state(i, j, model::itemp) = temp_zone;
       model::profile(model_num).state(i, j, model::ipres) = pres_zone;
       return;
    }

    //
    // If we are not in fluff region, then we have to do Newton solve
    //

    // start off the Newton loop by saying that the zone has not converged

    bool converged_hse = false;

    amrex::Real pres_zone;
    amrex::Real entropy_zone;
    amrex::Real p_want;
    amrex::Real s_want;
    amrex::Real drho;
    amrex::Real dtemp = 0;

    int i0;
    int j0;

    if (idir == 0) {
        // Integrate upward
        i0 = -1;
        j0 = 0;
    } else if (idir == 1) {
        // Integrate downward
        iX = +1;
        jY = 0;
    } else if (idir == 2) {
        // integrate to the right
        iX = 0;
        jY = -1;
    } else {
        amrex::Error("Invalid integration direction!")
    }

    for (int iter = 0; iter < fw::MAX_ITER; iter++) {
        // Get the pressure from HSE eqn.
        // If integrating in radial (x-) dir, then we have:
        //
        // P_{i,j} = P_{i∓1, j} ± 0.5 (ρ_{i∓1, j} + ρ_{i, j}) dx (g + 0.5 Ω^2 r_{i∓1/2} [1 - cos(2 θ_j)])
        //
        // Note, we are using an average of the density of the two zones as
        // an approximation of the interface value.
        // Also, the reference state depends on the direction of integration.
        // e.g. if integrating upward, reference state is below:   {i-1, j}
        //      if integrating downward, reference state is above: {i+1, j}
        //
        // If integrating in theta (y-) dir, then we have:
        //
        // P_{i,j} = P_{i, j-1} + 0.25 (ρ_{i, j-1} + ρ_{i, j}) dy Ω^2 r_i sin(2 θ_{j-1/2})
        //
        // We assumed integration in the right dir,
        // hence reference state is always from left {i, j-1}

        // First just try integrate either in vertical or horizontal direction

        amrex::Real r = model::profile(model_num).x(i) + static_cast<amrex::Real>(i0) * 0.5_rt * dx;
        amrex::Real twoTheta = 2.0_rt * (model::profile(model_num).y(j) +
                                         static_cast<amrex::Real>(j0) * 0.5_rt * dy);

        if (idir == 0 || idir == 1) {

            // integrating in vertical direction

            p_want = model::profile(model_num).state(i+i0, j, model::ipres) -
                static_cast<amrex::Real>(i0) * 0.5_rt * dx *
                (dens_zone + model::profile(model_num).state(i+i0, j, model::idens)) *
                (gravity::const_grav + 0.5_rt * Omega * Omega * r * (1.0_rt - std::cos(twoTheta)));

        } else if (idir == 2) {

            // integrating in to the right

            p_want = model::profile(model_num).state(i, j+j0, model::ipres) +
                0.25_rt * (dens_zone + model::profile(model_num).state(i, j+j0, model::idens)) *
                Omega * Omega * r * r * dy * std::sin(twoTheta);
        }

        //
        // We can either have isentropic or isothermal case.
        //

        if (isentropic) {
            // isentropic case
            // We have two unknowns: rho and T
            // with two constraint eqns: entropy and pressure constraints

            eos_state.T = temp_zone;
            eos_state.rho = dens_zone;
            for (int n = 0; n < NumSpec; n++) {
                eos_state.xn[n] = xn[n];
            }

            eos(eos_input_rt, eos_state);

            entropy_zone = eos_state.s;
            pres_zone = eos_state.p;

            Real dpt = eos_state.dpdT;
            Real dpd = eos_state.dpdr;
            Real dst = eos_state.dsdT;
            Real dsd = eos_state.dsdr;

            Real A = p_want - pres_zone;
            Real B = entropy_base - entropy_zone;

            dtemp = ((dsd / (dpd - 0.5_rt * dx * gravity::const_grav)) * A - B) /
                (dsd * dpt / (dpd - 0.5_rt * dx * gravity::const_grav) - dst);

            drho = (A - dpt * dtemp) / (dpd - 0.5_rt * dx * gravity::const_grav);

            dens_zone = amrex::Clamp(dens_zone + drho , 0.9_rt * dens_zone, 1.1_rt * dens_zone);
            temp_zone = amrex::Clamp(temp_zone + dtemp, 0.9_rt * temp_zone, 1.1_rt * temp_zone);

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
            // isothermal case
            // Since T is fixed, the only unknown is density
            // with one constraint eqn: pressure constraint

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

            entropy_zone = eos_state.s;
            pres_zone = eos_state.p;

            Real dpd = eos_state.dpdr;

            drho = (p_want - pres_zone) / (dpd - 0.5_rt * dx * gravity::const_grav);

            dens_zone = amrex::Clamp(dens_zone + drho , 0.9_rt * dens_zone, 1.1_rt * dens_zone);

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

        } // End of isothermal case

    } // End of Newton solve

    if (!converged_hse) {
        std::cout << "Error zone " << i << " did not converge in init_1d" << std::endl;
        std::cout << "integrate up" << std::endl;
        std::cout << dens_zone << " " << temp_zone << std::endl;
        std::cout << p_want << " " << entropy_base << " " << entropy_zone << std::endl;
        std::cout << drho << " " << dtemp << std::endl;
        amrex::Error("Error: HSE non-convergence");
    }

    //
    // Once rho and T are determined, call eos to get pressure
    //

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
#endif
