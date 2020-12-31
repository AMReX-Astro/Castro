#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>
#include <model_parser.H>
#include <ambient.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real coord[3] = {0.0};

    coord[0] = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);
#if AMREX_SPACEDIM >= 2
    coord[1] = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
#endif
#if AMREX_SPACEDIM == 3
    coord[2] = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
#endif

    Real shear_height_r = problem::shear_height * dx[2];

    int int_shear_width_x = static_cast<int>(problem::shear_width_x);
    int int_shear_width_y = static_cast<int>(problem::shear_width_y);
    Real velocity_gradient = problem::shear_amplitude / static_cast<Real>(problem::shear_height);

    int shear_bottom_index = static_cast<int>(problem::shear_height_loc / dx[2]) - problem::shear_height;

    state(i,j,k,URHO) = interpolate(coord[AMREX_SPACEDIM-1], model::idens);
    state(i,j,k,UTEMP) = interpolate(coord[AMREX_SPACEDIM-1], model::itemp);
    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = interpolate(coord[AMREX_SPACEDIM-1], model::ispec+n);
    }

    eos_t eos_state;
    eos_state.rho = state(i,j,k,URHO);
    eos_state.T = state(i,j,k,UTEMP);
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state(i,j,k,UFS+n);
    }

    eos(eos_input_rt, eos_state);

    // switch to conserved quantities

    state(i,j,k,UEDEN) = state(i,j,k,URHO) * eos_state.e;
    state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state.e;

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * state(i,j,k,UFS+n);
    }

    // Initial velocities

    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

#if AMREX_SPACEDIM == 3
    if (problem::shear_vel_field) {
        // First give shear velocity at h >= shear_height_loc - shear_height/2

        if (coord[2] > shear_bottom_index &&
            coord[2] <= shear_bottom_index + problem::shear_height) {
            state(i,j,k,UMX) = state(i,j,k,URHO) *
                velocity_gradient * std::abs(k - shear_bottom_index);

        } else if (k > shear_bottom_index + problem::shear_height) {
           state(i,j,k,UMX) = state(i,j,k,URHO) * problem::shear_amplitude;
        }

        if (int_shear_width_y > 2) {
            if (k > shear_bottom_index &&
                k <= shear_bottom_index + problem::shear_height &&
                j % int_shear_width_y <= 4) {
                state(i,j,k,UMX) = 0.0_rt;
                state(i,j,k,UMY) = 0.0_rt;
                state(i,j,k,UMZ) = 0.0_rt;
            }
        }

        if (int_shear_width_x > 2) {
            if (k > shear_bottom_index &&
                k <= shear_bottom_index + problem::shear_height &&
                i % int_shear_width_x <= 4) {
                state(i,j,k,UMX) = 0.0_rt;
                state(i,j,k,UMY) = 0.0_rt;
                state(i,j,k,UMZ) = 0.0_rt;
            }
        }
    }
#endif

    //  Now add the velocity perturbation (update the kinetic energy too)
    if (problem::apply_vel_field) {

#if AMREX_SPACEDIM == 2
        Real ydist = coord[1] - problem::velpert_height_loc;

        if (y >= shear_height_loc) {
            state(i,j,k,UMX) = state(i,j,k,URHO) * problem::shear_amplitude;
            state(i,j,k,UMY) = 0.0_rt;
        }

        Real upert[3] = {0.0_rt};

        // loop over each vortex

        for (int ivortex = 0; ivortex < problem::num_vortices; ivortex++) {

            Real xdist = coord[0] - problem::xloc_vortices[ivortex];

            Real r = std::sqrt(xdist * xdist + ydist * ydist);

            upert[0] += - (ydist / problem::velpert_scale) * problem::velpert_amplitude *
                std::exp(-r * r / std::pow(2.0_rt * problem::velpert_scale, 2)) *
                std::pow(-1.0_rt, ivortex+1);

            upert[1] += (xdist / problem::velpert_scale) * problem::velpert_amplitude *
                std::exp(-r * r / std::pow(2.0_rt * velpert_scale, 2)) *
                std::pow(-1.0_rt, ivortex+1);
        }

        state(i,j,UMX) += state(i,j,URHO) * upert[0];
        state(i,j,UMY) += state(i,j,URHO) * upert[1];

        state(i,j,UEDEN) += 0.5_rt * (state(i,j,k,UMX) * state(i,j,k,UMX) +
                                      state(i,j,k,UMY) * state(i,j,k,UMY) +
                                      state(i,j,k,UMZ) * state(i,j,k,UMZ)) / state(i,j,k,URHO);

#else
        Real zdist = coord[2] - problem::velpert_height_loc;

        Real upert[3] = {0.0_rt};

        // loop over each vortex

        if (problem::velpert_height_loc > 0.0_rt &&
            problem::velpert_scale > 0.0_rt &&
            problem::velpert_amplitude > 0.0 &&
            problem::num_vortices > 0) {

            for (int ivortex = 0; ivortex < problem::num_vortices; ivortex++) {

                Real xdist = coord[0] - problem::xloc_vortices[ivortex];
                Real r = std::sqrt(xdist * xdist + zdist * zdist);

                upert[0] += - (zdist / problem::velpert_scale) * problem::velpert_amplitude *
                    std::exp(-r * r / std::pow(2.0_rt * problem::velpert_scale, 2)) *
                    std::pow(-1.0_rt, ivortex+1);

                upert[2] += (xdist / problem::velpert_scale) * problem::velpert_amplitude *
                    std::exp(-r * r / std::pow(2.0_rt * problem::velpert_scale, 2)) *
                    std::pow(-1.0_rt, ivortex+1);
            }
        }

        state(i,j,k,UMX) += state(i,j,k,URHO) * upert[0];
        state(i,j,k,UMZ) += state(i,j,k,URHO) * upert[2];

        state(i,j,UEDEN) += 0.5_rt * (state(i,j,k,UMX) * state(i,j,k,UMX) +
                                      state(i,j,k,UMY) * state(i,j,k,UMY) +
                                      state(i,j,k,UMZ) * state(i,j,k,UMZ)) / state(i,j,k,URHO);
#endif
    }
}
#endif
