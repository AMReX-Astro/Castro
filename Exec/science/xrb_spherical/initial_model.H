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
void set_initial_condition(const int index_base, const int index_isentropic,
                           const int npts_model_x, const amrex::Real xmin,
                           const int npts_model_y, const amrex::Real ymin,
                           const model_t& model_params, const int model_num) {

    // Set the initial condition of the mesh

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

            // vertical hyperbolic tangent transition

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
            } else if (i >= index_base && i <= index_isentropic) {
                model::profile(model_num).state(i, j, model::itemp) =
                    model_params.T_star +
                    0.5_rt * (model_params.T_hi - model_params.T_star) *
                    (1.0_rt + evaluate_tanh(xc / (0.5_rt * model_params.atm_delta)));
            } else {
                model::profile(model_num).state(i, j, model::itemp) =
                    model_params.T_hi - std::pow(xc, 3);
            }

            // Floor temperature by T_lo

            if (model::profile(model_num).state(i, j, model::itemp) < model_params.T_lo) {
                model::profile(model_num).state(i, j, model::itemp) = model_params.T_lo;
            }

            //
            // We will also have a hyperbolic tangent transition in the theta dir
            // but for now just assume theres is no hot/cold model.
            //



            // the density and pressure will be determined via HSE,
            // for now, set them to the base conditions as initial guess

            model::profile(model_num).state(i, j, model::idens) = model_params.dens_base;
            model::profile(model_num).state(i, j, model::ipres) = pres_base;
        }
    }

    // make the base thermodynamics consistent for this base point,
    // this is assumed to be on θ=0 -- that is what we will integrate from!
    // here 0-index is actually not θ=0, due to extra ghost cell,
    // we might need to determine the actual index.

    eos_state.rho = model::profile(model_num).state(index_base, 0, model::idens);
    eos_state.T = model::profile(model_num).state(index_base, 0, model::itemp);
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = model::profile(model_num).state(index_base, 0, model::ispec+n);
    }

    eos(eos_input_rt, eos_state);

    model::profile(model_num).state(index_base, 0, model::ipres) = eos_state.p;
}


AMREX_INLINE
void HSE_integrate(const int i, const int j,
                   const amrex::Real dx, const amrex::Real dy,
                   const int idir, const int index_isentropic,
                   bool& isentropic, bool& fluff,
                   amrex::Real& entropy_base,
                   const model_t& model_params, const int model_num) {
    //
    // Goal is to update density, temperature, and pressure of the current zone based on integration
    // Also changes the `fluff` condition of the zone if density drops below
    // `low_density_cutoff`.
    //

    // Comment out the following to get pure isothermal profile

    // if (i == index_isentropic + 1) {
    //     // First time reaching isentropic atmosphere from siothermal transition region
    //     // set isentropic to true
    //     isentropic = true;
    // }

    // if (i > index_isentropic && !isentropic) {
    //     // If higher than the atmosphere region and its isothermal,
    //     // then we have reached the low temperature region, floor the temperature
    //     model::profile(model_num).state(i, j, model::itemp) = model_params.T_lo;
    // }

    // Get current state, this is the initial guess

    amrex::Real dens_zone = model::profile(model_num).state(i, j, model::idens);
    amrex::Real temp_zone = model::profile(model_num).state(i, j, model::itemp);
    amrex::Real xn[NumSpec];
    for (int n = 0; n < NumSpec; n++) {
        xn[n] = model::profile(model_num).state(i, j, model::ispec+n);
    }
    amrex::Real pres_zone;

    if (fluff) {
        // If in fluff region,
        // i.e. density of the zone is less than low_density_cutoff
        // We stop HSE in this zone and set to a constant density and temperature.

        dens_zone = model_params.low_density_cutoff;
        temp_zone = model_params.T_lo;

        // Call eos to get pressure: (T, rho) -> (p)

        eos_t eos_state;
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

    amrex::Real entropy_zone;
    amrex::Real p_want;
    amrex::Real drho;
    amrex::Real dtemp;

    amrex::Real Omega2 = 0.0_rt;
    if (castro::rotational_period > 0.0_rt && castro::rotation_include_centrifugal) {
        Omega2 = std::pow(2.0_rt * M_PI / castro::rotational_period, 2);
    }

    // Set index offset based on integration direction
    int i0;
    int j0;
    if (idir == 0) {
        // Integrate upward
        i0 = -1;
        j0 = 0;
    } else if (idir == 1) {
        // Integrate downward
        i0 = +1;
        j0 = 0;
    } else if (idir == 2) {
        // integrate to the right
        i0 = 0;
        j0 = -1;
    } else {
        amrex::Error("Invalid integration direction!");
    }

    // Do Newton solve

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
        // We assumed integration is always in the right dir,
        // hence reference state is always from left {i, j-1}

        eos_t eos_state;
        eos_state.T = temp_zone;
        eos_state.rho = dens_zone;
        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = xn[n];
        }

        eos(eos_input_rt, eos_state);

        entropy_zone = eos_state.s;
        pres_zone = eos_state.p;

        // Get HSE pressure and the density derivative based on integration direction

        amrex::Real r = model::profile(model_num).x(i) + static_cast<amrex::Real>(i0) * 0.5_rt * dx;
        amrex::Real twoTheta = 2.0_rt * model::profile(model_num).y(j) + static_cast<amrex::Real>(j0) * dy;

        amrex::Real dHSEp_drho;

        if (idir == 0 || idir == 1) {
            // integrating in vertical direction

            dHSEp_drho = -static_cast<amrex::Real>(i0) * 0.5_rt * dx *
                (gravity::const_grav + 0.5_rt * Omega2 * r * (1.0_rt - std::cos(twoTheta)));

            p_want = model::profile(model_num).state(i+i0, j, model::ipres) +
                (dens_zone + model::profile(model_num).state(i+i0, j, model::idens)) * dHSEp_drho;

        } else if (idir == 2) {
            // integrating in to the right

            dHSEp_drho = -static_cast<amrex::Real>(j0) * 0.25_rt * dy *
                Omega2 * r * r * std::sin(twoTheta);

            p_want = model::profile(model_num).state(i, j+j0, model::ipres) +
                (dens_zone + model::profile(model_num).state(i, j+j0, model::idens)) * dHSEp_drho;
        }

        amrex::Real A = p_want - pres_zone;
        amrex::Real dAdT = -eos_state.dpdT;
        amrex::Real dAdrho = dHSEp_drho - eos_state.dpdr;

        // Find dtemp if dealing with isentropic atmosphere

        dtemp = 0.0_rt;
        if (isentropic) {
            if (entropy_base <= 0.0_rt) {
                amrex::Error("isentropic atmosphere with negative entropy!");
            }

            amrex::Real B = entropy_base - entropy_zone;
            amrex::Real dBdT = -eos_state.dsdT;
            amrex::Real dBdrho = -eos_state.dsdr;

            dtemp = (B - (dBdrho / dAdrho) * A) /
                ((dBdrho / dAdrho) * dAdT - dBdT);
        }

        drho = -(A + dAdT * dtemp) / dAdrho;

        // Update zone

        dens_zone = amrex::Clamp(dens_zone + drho , 0.9_rt * dens_zone, 1.1_rt * dens_zone);
        temp_zone = amrex::Clamp(temp_zone + dtemp, 0.9_rt * temp_zone, 1.1_rt * temp_zone);

        // if (A < TOL .and. B < ETOL) then
        if (std::abs(drho) < fw::TOL * dens_zone && std::abs(dtemp) < fw::TOL * temp_zone) {
            converged_hse = true;
            break;
        }

        // Check if we entered isothermal region when integrating upward
        if (isentropic && temp_zone < model_params.T_lo && idir == 0) {
            temp_zone = model_params.T_lo;
            isentropic = false;
        }

        // check if the density falls below our minimum cut-off --
        // if so, floor it
        // Also detect whether we've entered the fluff region.
        // Note that, its only possible when integrating upward.

        if (dens_zone < model_params.low_density_cutoff && idir == 0) {
            dens_zone = model_params.low_density_cutoff;
            temp_zone = model_params.T_lo;
            converged_hse = true;
            fluff = true;
            break;
        }

    } // End of Newton solve

    if (!converged_hse) {
        std::cout << "Error zone " << i << ", " << j << " did not converge in 2D HSE" << std::endl;
        std::cout << dens_zone << " " << temp_zone << std::endl;
        std::cout << p_want << " " << entropy_base << " " << entropy_zone << std::endl;
        std::cout << drho << " " << dtemp << std::endl;
        amrex::Error("Error: HSE non-convergence");
    }

    //
    // Once rho and T are determined, call eos to get pressure
    //

    eos_t eos_state;
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

    // If reaching the zone prior to the isentropic atmosphere,
    // save the entropy, which will be used for isentropic integration
    if (i == index_isentropic) {
        entropy_base = eos_state.s;
    }
}


AMREX_INLINE
void
set_boundary_condition(const int npts_model_x, const int npts_model_y,
                       const amrex::Real dx, const amrex::Real dy,
                       const int index_base, const int index_isentropic,
                       bool fluff_mask[NPTS_X_MODEL][NPTS_Y_MODEL],
                       const model_t& model_params, const int model_num) {

    // This function sets the boundary condition along θ_min, θ_max,
    // r_min and r=H_star + 3 delta_atm.
    // This defines the boundary condition for the bottom isothermal region.

    /*
      Integrate the grid, and use the net result as boundary condition
      or the initial guess for the smoothing procedure.
     */


    // Integrate horizontally along the base

    bool fluff = false;
    bool isentropic = false;
    amrex::Real entropy_base = -1.0_rt;

    for (int j = 1; j < npts_model_y; j++) {
        HSE_integrate(index_base, j, dx, dy, 2,
                      index_isentropic,
                      isentropic, fluff, entropy_base,
                      model_params, model_num);
        fluff_mask[index_base][j] = fluff;
    }

    // Integrate vertically

    for (int j = 0; j < npts_model_y; ++j) {

        fluff = false;
        isentropic = false;
        entropy_base = -1.0_rt;

        // First integrate downward from base

        for (int i = index_base-1; i >= 0; --i) {
            HSE_integrate(i, j, dx, dy, 1,
                          index_isentropic,
                          isentropic, fluff, entropy_base,
                          model_params, model_num);
            fluff_mask[i][j] = fluff;
        }

        // Then integrate upward from base

        fluff = false;
        isentropic = false;
        entropy_base = -1.0_rt;

        for (int i = index_base+1; i < npts_model_x; i++) {
            HSE_integrate(i, j, dx, dy, 0,
                          index_isentropic,
                          isentropic, fluff, entropy_base,
                          model_params, model_num);
            fluff_mask[i][j] = fluff;
        }
    }


    //
    // integrate along θ_min from r_base (left boundary)
    //
    // bool fluff = false;
    // bool isentropic = false;
    // amrex::Real entropy_base = -1.0_rt;

    // // First integrate downward from base

    // for (int i = index_base-1; i >= 0; --i) {
    //     HSE_integrate(i, 0, dx, dy, 1,
    //                   index_isentropic,
    //                   isentropic, fluff, entropy_base,
    //                   model_params, model_num);
    // }

    // // Then integrate upward from base

    // fluff = false;
    // isentropic = false;
    // entropy_base = -1.0_rt;

    // for (int i = index_base+1; i < npts_model_x; i++) {
    //     HSE_integrate(i, 0, dx, dy, 0,
    //                   index_isentropic,
    //                   isentropic, fluff, entropy_base,
    //                   model_params, model_num);
    // }

    // //
    // // integrate along r_min from θ_min (bottom boundary)
    // //

    // fluff = false;
    // isentropic = false;
    // entropy_base = -1.0_rt;

    // for (int j = 1; j < npts_model_y; j++) {
    //     HSE_integrate(0, j, dx, dy, 2,
    //                   index_isentropic,
    //                   isentropic, fluff, entropy_base,
    //                   model_params, model_num);
    // }

    // // integrate along θ_max from r_min (right boundary)

    // fluff = false;
    // isentropic = false;
    // entropy_base = -1.0_rt;

    // // First integrate downward from base

    // for (int i = index_base-1; i >= 0; --i) {
    //     HSE_integrate(i, npts_model_y-1, dx, dy, 1,
    //                   index_isentropic,
    //                   isentropic, fluff, entropy_base,
    //                   model_params, model_num);
    // }

    // // Then integrate upward from base

    // fluff = false;
    // isentropic = false;
    // entropy_base = -1.0_rt;

    // for (int i = index_base+1; i < npts_model_x; i++) {
    //     HSE_integrate(i, npts_model_y-1, dx, dy, 0,
    //                   index_isentropic,
    //                   isentropic, fluff, entropy_base,
    //                   model_params, model_num);
    // }

    // // for (int i = 1; i < npts_model_x; i++) {
    // //     HSE_integrate(i, npts_model_y-1, dx, dy, 0,
    // //                   index_isentropic,
    // //                   isentropic, fluff, entropy_base,
    // //                   model_params, model_num);
    // // }

    // // Set upper boundary to fluff condition
    // fluff= true;
    // for (int j = 0; j < npts_model_y; ++j) {
    //     fluff_mask[npts_model_x-1][j] = fluff;

    //     HSE_integrate(npts_model_x-1, j, dx, dy, 2,
    //                   index_isentropic,
    //                   isentropic, fluff, entropy_base,
    //                   model_params, model_num);
    // }


    // // integrate along r = H_star + 3 delta_atm (isothermal top boundary)

    // fluff = false;
    // isentropic = false;
    // entropy_base = -1.0_rt;

    // for (int j = 1; j < npts_model_y; j++) {
    //     HSE_integrate(index_isentropic+1, j, dx, dy, 2,
    //                   index_isentropic,
    //                   isentropic, fluff, entropy_base,
    //                   model_params, model_num);
    // }


}

AMREX_INLINE
amrex::Real compute_relative_error(const int npts_model_x, const int npts_model_y,
                                   const amrex::Real dx, const amrex::Real dy,
                                   const bool fluff_mask[NPTS_X_MODEL][NPTS_Y_MODEL],
                                   const model_t& model_params, const int model_num) {

    // Computes the relative error for the smoothing process

    amrex::Real Omega2 = 0.0_rt;
    if (castro::rotational_period > 0.0_rt && castro::rotation_include_centrifugal) {
        Omega2 = std::pow(2.0_rt * M_PI / castro::rotational_period, 2);
    }

    amrex::Real residual_norm = 0.0_rt;;
    amrex::Real f_norm = 0.0_rt;

    for (int i = 1; i < npts_model_x-1; ++i) {
        for (int j = 1; j < npts_model_y-1; ++j) {

            // Ignore current zone if:
            // 1) the surrounding zone is fluff then we're at the top boundary interface
            // 2) the current zone is fluff
            if (fluff_mask[i][j] ||
                fluff_mask[i+1][j] || fluff_mask[i-1][j] ||
                fluff_mask[i][j-1] || fluff_mask[i][j+1]) {
                continue;
            }

            amrex::Real rc = model::profile(model_num).x(i);
            amrex::Real rci = 1.0_rt / rc;
            amrex::Real rc2i = 1.0_rt / (rc * rc);
            amrex::Real rp = rc + 0.5_rt * dx;
            amrex::Real rm = rc - 0.5_rt * dx;
            amrex::Real rp2 = rp * rp;
            amrex::Real rm2 = rm * rm;

            amrex::Real sinThetac = std::sin(model::profile(model_num).y(j));
            amrex::Real sinThetap = std::sin(model::profile(model_num).y(j) + 0.5_rt * dy);
            amrex::Real sinThetam = std::sin(model::profile(model_num).y(j) - 0.5_rt * dy);

            amrex::Real cosThetac = std::cos(model::profile(model_num).y(j));
            amrex::Real cosThetap = std::cos(model::profile(model_num).y(j) + 0.5_rt * dy);
            amrex::Real cosThetam = std::cos(model::profile(model_num).y(j) - 0.5_rt * dy);

            amrex::Real dxi =  1.0_rt / dx;
            amrex::Real dx2i = 1.0_rt / (dx * dx);
            amrex::Real dySinThetaci = 1.0_rt / (dy * sinThetac);
            amrex::Real dy2sinThetaci = 1.0_rt / (dy * dy * sinThetac);

            // amrex::Real sinTwoThetac = std::sin(2.0_rt * model::profile(model_num).y(j));
            // amrex::Real sinTwoThetap = std::sin(2.0_rt * model::profile(model_num).y(j) + dy);
            // amrex::Real sinTwoThetam = std::sin(2.0_rt * model::profile(model_num).y(j) - dy);
            // amrex::Real oneMinusCosTwoThetac = 1.0_rt - std::cos(2.0_rt * model::profile(model_num).y(j));

            amrex::Real psi_rc = Omega2 * rc * sinThetac * sinThetac + gravity::const_grav;
            amrex::Real psi_rp = Omega2 * rp * sinThetac * sinThetac + gravity::const_grav;
            amrex::Real psi_rm = Omega2 * rm * sinThetac * sinThetac + gravity::const_grav;

            amrex::Real psi_tc = Omega2 * rc * sinThetac * cosThetac;
            amrex::Real psi_tp = Omega2 * rc * sinThetap * cosThetap;
            amrex::Real psi_tm = Omega2 * rc * sinThetam * cosThetam;

            // amrex::Real psi_rc = 0.5_rt * Omega2 * rc * oneMinusCosTwoThetac + gravity::const_grav;
            // amrex::Real psi_rp = 0.5_rt * Omega2 * rp * oneMinusCosTwoThetac + gravity::const_grav;
            // amrex::Real psi_rm = 0.5_rt * Omega2 * rm * oneMinusCosTwoThetac + gravity::const_grav;

            // amrex::Real psi_tc = 0.5_rt * Omega2 * rc * sinTwoThetac;
            // amrex::Real psi_tp = 0.5_rt * Omega2 * rc * sinTwoThetap;
            // amrex::Real psi_tm = 0.5_rt * Omega2 * rc * sinTwoThetam;

            // amrex::Real psi_drho = psi_rc * (1.0_rt + (psi_tc * psi_tc) / (psi_rc * psi_rc)) *
            //     (model::profile(model_num).state(i+1, j, model::idens) -
            //      model::profile(model_num).state(i-1, j, model::idens)) / (2.0_rt * dx);

            // amrex::Real psi_drho = 0.5_rt * ((model::profile(model_num).state(i+1, j  , model::idens) -
            //                                   model::profile(model_num).state(i-1, j  , model::idens)) *
            //                                  psi_rc / dx +
            //                                  (model::profile(model_num).state(i  , j+1, model::idens) -
            //                                   model::profile(model_num).state(i  , j-1, model::idens)) *
            //                                  psi_tc / (rc * dy));

            // amrex::Real dHSEp_drho = 2.0_rt * (Omega2 + gravity::const_grav / rc);

            // amrex::Real f = model::profile(model_num).state(i, j, model::idens) * dHSEp_drho + psi_drho;

            amrex::Real rhoip = 0.5_rt * (model::profile(model_num).state(i+1, j  , model::idens) +
                                          model::profile(model_num).state(i  , j  , model::idens));
            amrex::Real rhoim = 0.5_rt * (model::profile(model_num).state(i  , j  , model::idens) +
                                          model::profile(model_num).state(i-1, j  , model::idens));
            amrex::Real rhojp = 0.5_rt * (model::profile(model_num).state(i  , j+1, model::idens) +
                                          model::profile(model_num).state(i  , j  , model::idens));
            amrex::Real rhojm = 0.5_rt * (model::profile(model_num).state(i  , j  , model::idens) +
                                          model::profile(model_num).state(i  , j-1, model::idens));

            amrex::Real f = (rp2 * rhoip * psi_rp - rm2 * rhoim * psi_rm) * rc2i * dxi +
                            (sinThetap * rhojp * psi_tp - sinThetam * rhojm * psi_tm) * rci * dySinThetaci;

            amrex::Real d2p = rc2i *
                (dx2i * (rp2 * model::profile(model_num).state(i+1, j  , model::ipres) +
                         rm2 * model::profile(model_num).state(i-1, j  , model::ipres) -
                 (rp2 + rm2) * model::profile(model_num).state(i  , j  , model::ipres)) +
  dy2sinThetaci * (sinThetap * model::profile(model_num).state(i  , j+1, model::ipres) +
                   sinThetam * model::profile(model_num).state(i  , j-1, model::ipres) -
     (sinThetap + sinThetam) * model::profile(model_num).state(i  , j  , model::ipres)));

            amrex::Real residual = f - d2p;

            // Do L1 norm
            // if (std::abs(residual) > std::abs(f)) {
            //     std::cout << "f in " << i << " " << j << " is " << f << std::endl;
            //     std::cout << "d2p in " << i << " " << j << " is " << d2p << std::endl;
            // }
            f_norm += std::abs(f);
            residual_norm += std::abs(residual);
        }
    }

    // Return relative error
    return residual_norm / f_norm;
}


AMREX_INLINE
void HSE_smooth(const int i, const int j,
                const amrex::Real dx, const amrex::Real dy,
                bool fluff_mask[NPTS_X_MODEL][NPTS_Y_MODEL],
                const model_t& model_params, const int model_num) {
    // Smoothing procedure for HSE

    // Get current state, this is the initial guess

    amrex::Real dens_zone = model::profile(model_num).state(i, j, model::idens);
    amrex::Real temp_zone = model::profile(model_num).state(i, j, model::itemp);
    amrex::Real xn[NumSpec];
    for (int n = 0; n < NumSpec; n++) {
        xn[n] = model::profile(model_num).state(i, j, model::ispec+n);
    }
    amrex::Real pres_zone;

    // Ignore current zone if:
    // 1) the surrounding zone is fluff then we're at the top boundary interface
    // 2) the current zone is fluff
    if (fluff_mask[i][j] ||
        fluff_mask[i+1][j] || fluff_mask[i-1][j] ||
        fluff_mask[i][j-1] || fluff_mask[i][j+1]) {
        return;
    }

    // if (fluff_mask[i][j]) {
        // If in fluff region,
        // i.e. density of the zone is less than low_density_cutoff
        // We stop HSE in this zone and set to a constant density and temperature.

        // dens_zone = model_params.low_density_cutoff;
        // temp_zone = model_params.T_lo;

        // // Call eos to get pressure: (T, rho) -> (p)

        // eos_t eos_state;
        // eos_state.T = temp_zone;
        // eos_state.rho = dens_zone;
        // for (int n = 0; n < NumSpec; n++) {
        //     eos_state.xn[n] = xn[n];
        // }

        // eos(eos_input_rt, eos_state);

        // pres_zone = eos_state.p;

        // // update the thermodynamics in this zone

        // model::profile(model_num).state(i, j, model::idens) = dens_zone;
        // model::profile(model_num).state(i, j, model::itemp) = temp_zone;
        // model::profile(model_num).state(i, j, model::ipres) = pres_zone;
        // return;
    // }

    //
    // If we are not in fluff region, then we have to do Newton solve
    //

    // start off the Newton loop by saying that the zone has not converged

    bool converged_hse = false;

    amrex::Real entropy_zone;
    amrex::Real p_want;
    amrex::Real drho;
    amrex::Real dtemp;

    // precompute some terms for evaluating HSE pressure.
    amrex::Real Omega2 = 0.0_rt;
    if (castro::rotational_period > 0.0_rt && castro::rotation_include_centrifugal) {
        Omega2 = std::pow(2.0_rt * M_PI / castro::rotational_period, 2);
    }


    amrex::Real rc = model::profile(model_num).x(i);
    amrex::Real rc2 = rc * rc;
    amrex::Real rci = 1.0_rt / rc;
    amrex::Real rc2i = 1.0_rt / (rc * rc);
    amrex::Real rp = rc + 0.5_rt * dx;
    amrex::Real rm = rc - 0.5_rt * dx;
    amrex::Real rp2 = rp * rp;
    amrex::Real rm2 = rm * rm;

    amrex::Real sinThetac = std::sin(model::profile(model_num).y(j));
    amrex::Real sinThetap = std::sin(model::profile(model_num).y(j) + 0.5_rt * dy);
    amrex::Real sinThetam = std::sin(model::profile(model_num).y(j) - 0.5_rt * dy);

    amrex::Real cosThetac = std::cos(model::profile(model_num).y(j));
    amrex::Real cosThetap = std::cos(model::profile(model_num).y(j) + 0.5_rt * dy);
    amrex::Real cosThetam = std::cos(model::profile(model_num).y(j) - 0.5_rt * dy);

    amrex::Real dx2 = dx * dx;
    amrex::Real dxi =  1.0_rt / dx;
    amrex::Real dx2i = 1.0_rt / (dx * dx);
    amrex::Real dy2sinThetac = dy * dy * sinThetac;
    amrex::Real dySinThetaci = 1.0_rt / (dy * sinThetac);
    amrex::Real dy2sinThetaci = 1.0_rt / (dy * dy * sinThetac);

    // amrex::Real sinTwoThetac = std::sin(2.0_rt * model::profile(model_num).y(j));
    // amrex::Real sinTwoThetap = std::sin(2.0_rt * model::profile(model_num).y(j) + dy);
    // amrex::Real sinTwoThetam = std::sin(2.0_rt * model::profile(model_num).y(j) - dy);
    // amrex::Real oneMinusCosTwoThetac = 1.0_rt - std::cos(2.0_rt * model::profile(model_num).y(j));

    amrex::Real psi_rc = Omega2 * rc * sinThetac * sinThetac + gravity::const_grav;
    amrex::Real psi_rp = Omega2 * rp * sinThetac * sinThetac + gravity::const_grav;
    amrex::Real psi_rm = Omega2 * rm * sinThetac * sinThetac + gravity::const_grav;

    amrex::Real psi_tc = Omega2 * rc * sinThetac * cosThetac;
    amrex::Real psi_tp = Omega2 * rc * sinThetap * cosThetap;
    amrex::Real psi_tm = Omega2 * rc * sinThetam * cosThetam;

    // amrex::Real psi_rc = 0.5_rt * Omega2 * rc * oneMinusCosTwoThetac + gravity::const_grav;
    // amrex::Real psi_rp = 0.5_rt * Omega2 * rp * oneMinusCosTwoThetac + gravity::const_grav;
    // amrex::Real psi_rm = 0.5_rt * Omega2 * rm * oneMinusCosTwoThetac + gravity::const_grav;

    // amrex::Real psi_tc = 0.5_rt * Omega2 * rc * sinTwoThetac;
    // amrex::Real psi_tp = 0.5_rt * Omega2 * rc * sinTwoThetap;
    // amrex::Real psi_tm = 0.5_rt * Omega2 * rc * sinTwoThetam;

    // amrex::Real fac = psi_rc * (1.0_rt + (psi_tc * psi_tc) / (psi_rc * psi_rc)) *
    //     (model::profile(model_num).state(i+1, j, model::idens) -
    //      model::profile(model_num).state(i-1, j, model::idens)) / (2.0_rt * dx);

    // amrex::Real dHSEp_drho = 2.0_rt * (Omega2 + gravity::const_grav / rc);

    // Do Newton solve

    for (int iter = 0; iter < fw::MAX_ITER; iter++) {

        eos_t eos_state;
        eos_state.T = temp_zone;
        eos_state.rho = dens_zone;
        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = xn[n];
        }

        eos(eos_input_rt, eos_state);

        entropy_zone = eos_state.s;
        pres_zone = eos_state.p;

        // Get HSE pressure and the density derivative

        // amrex::Real f = dens_zone * dHSEp_drho + fac;

        amrex::Real rhoip = 0.5_rt * (model::profile(model_num).state(i+1, j  , model::idens) +
                                      model::profile(model_num).state(i  , j  , model::idens));
        amrex::Real rhoim = 0.5_rt * (model::profile(model_num).state(i  , j  , model::idens) +
                                      model::profile(model_num).state(i-1, j  , model::idens));
        amrex::Real rhojp = 0.5_rt * (model::profile(model_num).state(i  , j+1, model::idens) +
                                      model::profile(model_num).state(i  , j  , model::idens));
        amrex::Real rhojm = 0.5_rt * (model::profile(model_num).state(i  , j  , model::idens) +
                                      model::profile(model_num).state(i  , j-1, model::idens));

        amrex::Real f = (rp2 * rhoip * psi_rp - rm2 * rhoim * psi_rm) * rc2i * dxi +
            (sinThetap * rhojp * psi_tp - sinThetam * rhojm * psi_tm) * rci * dySinThetaci;

        amrex::Real df_drho = 0.5_rt * ((rp2 * psi_rp - rm2 * psi_rm) * rc2i * dxi +
                                        (sinThetap * psi_tp - sinThetam * psi_tm) * rci * dySinThetaci);
        amrex::Real dHSEp_drho = - rc2 * dx2 * dy2sinThetac * df_drho /
            (dy2sinThetac * (rp2 + rm2) + dx2 * (sinThetap + sinThetam));

        p_want = (
                  dy2sinThetac * (rp2 * model::profile(model_num).state(i+1, j  , model::ipres)  +
                                  rm2 * model::profile(model_num).state(i-1, j  , model::ipres)) +
                     dx2 * (sinThetap * model::profile(model_num).state(i  , j+1, model::ipres)  +
                            sinThetam * model::profile(model_num).state(i  , j-1, model::ipres)) -
                  rc2 * dx2 * dy2sinThetac * f
                  ) / ( dy2sinThetac * (rp2 + rm2) + dx2 * (sinThetap + sinThetam) );

        amrex::Real A = p_want - pres_zone;
        amrex::Real dAdT = -eos_state.dpdT;
        amrex::Real dAdrho = dHSEp_drho - eos_state.dpdr;

        // Find dtemp if dealing with isentropic atmosphere

        dtemp = 0.0_rt;
        drho = -(A + dAdT * dtemp) / dAdrho;

        // Update zone

        dens_zone = amrex::Clamp(dens_zone + drho , 0.9_rt * dens_zone, 1.1_rt * dens_zone);
        temp_zone = amrex::Clamp(temp_zone + dtemp, 0.9_rt * temp_zone, 1.1_rt * temp_zone);

        // if (A < TOL .and. B < ETOL) then
        if (std::abs(drho) < fw::TOL * dens_zone && std::abs(dtemp) < fw::TOL * temp_zone) {
            converged_hse = true;
            break;
        }

        // check if the density falls below our minimum cut-off --
        // if so, floor it
        // Also detect whether we've entered the fluff region.
        // Note that, its only possible when integrating upward.

        if (dens_zone < model_params.low_density_cutoff) {
            dens_zone = model_params.low_density_cutoff;
            temp_zone = model_params.T_lo;
            converged_hse = true;
            fluff_mask[i][j] = true;
            break;
        }

    } // End of Newton solve

    if (!converged_hse) {
        std::cout << "Error zone " << i << ", " << j << " did not converge in 2D HSE" << std::endl;
        std::cout << dens_zone << " " << temp_zone << std::endl;
        std::cout << p_want << " " << " " << entropy_zone << std::endl;
        std::cout << drho << " " << dtemp << std::endl;
        amrex::Error("Error: HSE non-convergence");
    }

    //
    // Once rho and T are determined, call eos to get pressure
    //

    eos_t eos_state;
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
}


AMREX_INLINE
void
generate_initial_model(const int npts_model_x, const amrex::Real xmin, const amrex::Real xmax,
                       const int npts_model_y, const amrex::Real ymin, const amrex::Real ymax,
                       const model_t& model_params, const int model_num)
{


    // Create a 2-d uniform grid that is identical to the mesh that we are
    // mapping onto, and then we want to force it into HSE on that mesh.

    model::npts_x = npts_model_x;
    model::npts_y = npts_model_y;
    model::initialized = true;

    if (npts_model_x > NPTS_X_MODEL) {
        amrex::Error("Error: model has more than NPTS_MODEL points in x-dir,  Increase MAX_NPTS_X_MODEL");
    }
    if (npts_model_y > NPTS_Y_MODEL) {
        amrex::Error("Error: model has more than NPTS_MODEL points in y-dir,  Increase MAX_NPTS_Y_MODEL");
    }

    // create the grid -- cell centers

    amrex::Real dx = (xmax - xmin) / static_cast<amrex::Real>(npts_model_x);
    amrex::Real dy = (ymax - ymin) / static_cast<amrex::Real>(npts_model_y);

    for (int i = 0; i < npts_model_x; i++) {
        model::profile(model_num).x(i) =
            xmin + (static_cast<amrex::Real>(i) + 0.5_rt) * dx;
    }

    for (int j = 0; j < npts_model_y; j++) {
        model::profile(model_num).y(j) =
            ymin + (static_cast<amrex::Real>(j) + 0.5_rt) * dy;
    }

    // find the index of the base height

    int index_base = -1;
    for (int i = 0; i < npts_model_x; i++) {
        if (model::profile(model_num).x(i) >= xmin + model_params.H_star) {
            index_base = i;
            break;
        }
    }

    if (index_base == -1) {
        amrex::Error("ERROR: invalid base_height");
    }

    // Find the index of isentropic atmosphere,
    // this is to find the entropy of the last zone of the
    // isothermal layer prior to the isentropic atmosphere

    int index_isentropic = -1;
    for (int i = 0; i < npts_model_x; i++) {
        if (model::profile(model_num).x(i) > xmin + model_params.H_star + 3.0_rt * model_params.atm_delta) {
            index_isentropic = i - 1;
            if (i == index_base + 1) {
                amrex::Error("atm_delta is too small for the grid spacing (at least one transition zone is required)");
            }
            break;
        }
    }

    //
    // put the model onto our new uniform grid
    // This is the initial guess conditions.
    //

    set_initial_condition(index_base, index_isentropic,
                          npts_model_x, xmin,
                          npts_model_y, ymin,
                          model_params, model_num);

    //
    // set boundary conditions
    //

    bool fluff_mask[NPTS_X_MODEL][NPTS_Y_MODEL] = {false};

    set_boundary_condition(npts_model_x, npts_model_y, dx, dy,
                           index_base, index_isentropic, fluff_mask,
                           model_params, model_num);

    amrex::Real relative_error = compute_relative_error(npts_model_x, npts_model_y,
                                                        dx, dy, fluff_mask,
                                                        model_params, model_num);
    std::cout << "initial relative error is " << relative_error << std::endl;

    //
    // Multigrid Relaxation + Newton solve
    //

    int smooth_count = 0;
    int converging_counter = 0;

    while (converging_counter < 20) {

        // Do smoothing over the interior grid

        for (int i = 1; i < npts_model_x-1; ++i) {
            for (int j = 1; j < npts_model_y-1; ++j) {

            // Ignore current zone if:
            // 1) the surrounding zone is fluff then we're at the top boundary interface
            // 2) the current zone is fluff
            if (fluff_mask[i][j] ||
                fluff_mask[i+1][j] || fluff_mask[i-1][j] ||
                fluff_mask[i][j-1] || fluff_mask[i][j+1]) {
                continue;
            }
                HSE_smooth(i, j, dx, dy,
                           fluff_mask,
                           model_params, model_num);
            }
        }

        // Compute relative error

        if (smooth_count % 10 == 0) {
            amrex::Real new_relative_error = compute_relative_error(npts_model_x, npts_model_y,
                                                                    dx, dy, fluff_mask,
                                                                    model_params, model_num);
            std::cout << "relative error is " << new_relative_error << std::endl;

            if (relative_error < new_relative_error) {
                ++converging_counter;
            }
            relative_error = new_relative_error;
        }

        ++smooth_count;

        if (smooth_count > 1e5) {
            break;
        }
    }
    std::cout << "Number of relaxation iterations took are " << smooth_count << std::endl;
}
#endif
