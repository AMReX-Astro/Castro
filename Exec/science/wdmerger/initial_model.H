#ifndef initial_model_H
#define initial_model_H

#include <AMReX.H>

#include <network_properties.H>
#include <eos.H>
#include <interpolate.H>
#include <model_parser_data.H>
#include <prob_parameters.H>

void establish_hse (model::initial_model_t& model,
                    Real& mass, Real& central_density,
                    Real envelope_mass, Real& radius,
                    const Real core_comp[NumSpec], const Real envelope_comp[NumSpec]);

// Takes a one-dimensional stellar model and interpolates it to a point in
// 3D Cartesian grid space. Optionally does a sub-grid-scale interpolation
// if nsub > 1 (set in the probin file).

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void interpolate_3d_from_1d (model::initial_model_t& model, Real radius,
                             const Real* loc, const Real* dx, eos_t& state, int nsub = 1)
{
    state.rho = 0.0_rt;
    state.p   = 0.0_rt;
    state.T   = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        state.xn[n] = 0.0_rt;
    }

    // If the model radius is smaller than the zone size, just use the center of the model.

    Real dist = std::sqrt(loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2]);

    Real max_dx = dx[0];
    if (AMREX_SPACEDIM >= 2) {
        max_dx = amrex::max(max_dx, dx[1]);
    }
    if (AMREX_SPACEDIM == 3) {
        max_dx = amrex::max(max_dx, dx[2]);
    }

    if (radius <= max_dx && dist < max_dx) {

        state.rho = model.state(0, model::idens);
        state.T   = model.state(0, model::itemp);
        for (int n = 0; n < NumSpec; ++n) {
            state.xn[n]  = model.state(0, model::ispec + n);
        }

    }
    else {

        // We perform a sub-grid-scale interpolation, where
        // nsub determines the number of intervals we split the zone into.
        // If nsub = 1, we simply interpolate using the cell-center location.
        // As an example, if nsub = 3, then the sampled locations will be
        // k = 0 --> z = loc(3) - dx(3) / 3   (1/6 of way from left edge of zone)
        // k = 1 --> z = loc(3)               (halfway between left and right edge)
        // k = 2 --> z = loc(3) + dx(3) / 3   (1/6 of way from right edge of zone)

        for (int k = 0; k < nsub; ++k) {
            Real z = loc[2] + (static_cast<Real>(k) + 0.5_rt * (1 - nsub)) * dx[2] / nsub;

            for (int j = 0; j < nsub; ++j) {
                Real y = loc[1] + (static_cast<Real>(j) + 0.5_rt * (1 - nsub)) * dx[1] / nsub;

                for (int i = 0; i < nsub; ++i) {
                    Real x = loc[0] + (static_cast<Real>(i) + 0.5_rt * (1 - nsub)) * dx[0] / nsub;

                    Real dist = std::sqrt(x * x + y * y + z * z);

                    state.rho = state.rho + interpolate::interpolate(dist, problem::initial_model_npts, &model.r(0), &model.state(0, model::idens));
                    state.T   = state.T   + interpolate::interpolate(dist, problem::initial_model_npts, &model.r(0), &model.state(0, model::itemp));

                    for (int n = 0; n < NumSpec; ++n) {
                        state.xn[n] += interpolate::interpolate(dist, problem::initial_model_npts, &model.r(0), &model.state(0, model::ispec + n));
                    }
                }
            }
        }

        // Now normalize by the number of intervals.

        state.rho = state.rho / (nsub * nsub * nsub);
        state.T   = state.T   / (nsub * nsub * nsub);
        for (int n = 0; n < NumSpec; ++n) {
            state.xn[n] = state.xn[n] / (nsub * nsub * nsub);
        }

    }

    // Complete the thermodynamics.

    eos(eos_input_rt, state);
}

#endif
