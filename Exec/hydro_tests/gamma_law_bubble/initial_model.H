#ifndef INITIAL_MODEL_H
#define INITIAL_MODEL_H

#include <prob_parameters.H>
#include <network.H>
#include <eos.H>

struct model_t {
    Real dens_base = -1;
    Real p_base = -1;
    Real xn[NumSpec] = {0.0};
};



///
/// construct an initial model in HSE.  Note: this does not return
/// anything, but rather updates the model_parser globals with the
/// model information.
///
AMREX_INLINE
void
generate_initial_model(const int npts_model, const Real xmin, const Real xmax,
                       const model_t model_params) {

    model::npts = npts_model;
    model::initialized = true;

    if (npts_model > NPTS_MODEL) {
        amrex::Error("Error: model has more than NPTS_MODEL points,  Increase MAX_NPTS_MODEL");
    }

    // compute the pressure scale height (for an isothermal, ideal-gas
    // atmosphere)

    Real H = model_params.p_base / model_params.dens_base / std::abs(gravity::const_grav);

    // create the grid -- cell centers

    Real dx = (xmax - xmin) / npts_model;

    for (int i = 0; i < npts_model; i++) {
        model::profile(0).r(i) = xmin + (static_cast<Real>(i) + 0.5_rt) * dx;
    }

    for (int i = 0; i < npts_model; i++) {

        // initial guess

        Real temp = T_guess;

        Real density;

        // get the zone-center density

        if (problem::do_isentropic == 1) {

            // We integrate dp/dz = rho g with p = K rho**gamma
            // (isentropic).  Take the zero-point of our integration
            // as the first zone center, where p = p_base.  This
            // means that the next zone center is just dx away, so for
            // evaluating the height here, we just use i * dx, not
            // (i+0.5)*dx.

            Real z = xmin + static_cast<Real>(i) * dx;

            density = model_params.dens_base *
                std::pow(gravity::const_grav * model_params.dens_base *
                 (eos_rp::eos_gamma - 1.0) * z /
                 (eos_rp::eos_gamma * model_params.p_base) + 1.0_rt, 1.0_rt/(eos_rp::eos_gamma - 1.0_rt));

        } else {
            Real z = xmin + (static_cast<Real>(i) + 0.5_rt) * dx;

            density = model_params.dens_base * std::exp(-z/H);
        }

        // for zone 0, the pressure at the center is just the base
        // pressure.  For other zones, we will difference HSE.

        Real pressure;

        if (i == 0) {
            pressure = model_params.p_base;

        } else {
            pressure = model::profile(0).state(i-1, model::ipres) -
                dx * 0.5_rt * (density + model::profile(0).state(i-1, model::idens)) * std::abs(gravity::const_grav);
        }

        eos_t eos_state;
        eos_state.p = pressure;
        eos_state.T = temp;
        eos_state.rho = density;
        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = model_params.xn[n];
        }

        eos(eos_input_rp, eos_state);

        model::profile(0).state(i, model::idens) = density;
        model::profile(0).state(i, model::ipres) = pressure;
        model::profile(0).state(i, model::itemp) = eos_state.T;
        for (int n = 0; n < NumSpec; n++) {
            model::profile(0).state(i, model::ispec+n) = model_params.xn[n];
        }

    }

}

#endif
