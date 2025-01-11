#ifndef INITIAL_MODEL_H
#define INITIAL_MODEL_H

#include <prob_parameters.H>
#include <network.H>
#include <eos.H>

struct model_t {
    amrex::Real dens_base = -1;
    amrex::Real scale_height = -1;
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

    amrex::Real pres_base = model_params.scale_height * model_params.dens_base * std::abs(gravity::const_grav);

    // create the grid -- cell centers

    Real dx = (xmax - xmin) / npts_model;

    for (int i = 0; i < npts_model; i++) {
        model::profile(0).r(i) = xmin + (static_cast<Real>(i) + 0.5_rt) * dx;
    }

    for (int i = 0; i < npts_model; i++) {

        amrex::Real density;
        amrex::Real pressure;

        amrex::Real profile = 1.0 - (eos_rp::eos_gamma - 1.0) / eos_rp::eos_gamma * model::profile(0).r(i) / model_params.scale_height;

        if (profile > 0.0) {
            density = std::max(model_params.dens_base * std::pow(profile, 1.0/(eos_rp::eos_gamma - 1.0)),
                               problem::low_density_cutoff);
        } else {
            density = problem::low_density_cutoff;
        }

        if (i == 0) {
            pressure = pres_base;
        } else if (density <= problem::low_density_cutoff + 1.e-30_rt) {
            pressure = model::profile(0).state(i-1, model::ipres);
        } else {
            pressure = pres_base *
                std::pow(density / model_params.dens_base, eos_rp::eos_gamma);
        }

        // initial guess

        amrex::Real temp = T_guess;

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
