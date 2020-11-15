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
        model::r(i) = xmin + (static_cast<Real>(i) + 0.5_rt) * dx;
    }

    // create the constant if we are isentropic

    eos_t eos_state;
    eos_state.rho = model_params.dens_base;
    eos_state.p = model_params.p_base;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = model_params.xn[n];
    }

    // initial guess

    eos_state.T = 1000.0_rt;

    eos(eos_input_rp, eos_state);

    Real gamma_const = model_params.p_base/(model_params.dens_base * eos_state.e) + 1.0_rt;

    model::state(0, model::idens) = model_params.dens_base;
    model::state(0, model::ipres) = model_params.p_base;
    model::state(0, model::itemp) = eos_state.T;

    for (int n = 0; n < NumSpec; n++) {
        model::state(0, model::ispec+n) = model_params.xn[n];
    }

    model::r(0) = xmin + 0.5_rt * dx;

    // integrate up from the base

    for (int j = 0; j < npts_model; j++) {

        model::r(j) = xmin + (static_cast<Real>(j) + 0.5_rt) * dx;

        if (problem::do_isentropic == 1) {
            model::state(j, model::idens) = model_params.dens_base *
                std::pow(gravity::const_grav * model_params.dens_base * (gamma_const - 1.0) *
                         (model::r(j) - model::r(0)) /
                         (gamma_const * model_params.p_base) + 1.0_rt, 1.0_rt/(gamma_const - 1.e0_rt));
        } else {
            model::state(j, model::idens) = model_params.dens_base *
                std::exp(-(model::r(j) - model::r(0)) / H);
        }

        model::state(j, model::ipres) = model::state(j-1, model::ipres) -
            dx * 0.5_rt * (model::state(j, model::idens) + model::state(j-1, model::idens)) *
            std::abs(gravity::const_grav);


        eos_t eos_state;

        eos_state.rho = model::state(j, model::idens);
        eos_state.p = model::state(j, model::ipres);
        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = model_params.xn[n];
        }

        // initial guess

        eos_state.T = 1000.0_rt;

        eos(eos_input_rp, eos_state);

        model::state(j, model::itemp) = eos_state.T;

        for (int n = 0; n < NumSpec; n++) {
            model::state(j, model::ispec+n) = model_params.xn[n];
        }
    }
}
#endif
