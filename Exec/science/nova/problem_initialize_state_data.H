#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>
#include <model_parser.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

    Real y = 0.0;
#if AMREX_SPACEDIM >= 2
    y = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
#endif

    Real z = 0.0;
#if AMREX_SPACEDIM == 3
    z = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
#endif

#if AMREX_SPACEDIM == 2
    Real height = y;
#else
    Real height = z;
#endif

    state(i,j,k,URHO) = interpolate(height, model::idens);
    state(i,j,k,UTEMP) = interpolate(height, model::itemp);
    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = amrex::max(interpolate(height, model::ispec+n), network_rp::small_x);
    }

    // normalize

    Real sum = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
        sum += state(i,j,k,UFS+n);
    }

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) /= sum;
    }

    // Here we enforce the EOS by calling rho and T (giving us p)
    eos_t eos_state;
    eos_state.rho = state(i,j,k,URHO);
    eos_state.T = state(i,j,k,UTEMP);
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state(i,j,k,UFS+n);
    }

    eos(eos_input_rt, eos_state);

    state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state.e;
    state(i,j,k,UEDEN) = state(i,j,k,UEINT);

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * state(i,j,k,UFS+n);
    }


    // Setting up the initial velocities = 0
    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;


    //  Now add the velocity perturbation
    if (problem::apply_vel_field) {

        Real zdist = 0.0_rt;
        Real ydist = y - problem::y_h1 - problem::velpert_scale_ratio * problem::velpert_scale;

        Real upert[3] = {0.0_rt};

        // loop over each vortex

        for (int vortex = 0; vortex < problem::num_vortices; vortex++) {

            Real xdist = x - problem::xloc_vortices[vortex];

            Real r = std::sqrt(xdist * xdist + ydist * ydist + zdist * zdist);

            upert[0] += -(ydist / problem::velpert_scale) * problem::velpert_amplitude *
                std::exp(-r * r / (2.0_rt * problem::velpert_scale * problem::velpert_scale)) *
                std::pow(-1.0_rt, vortex+1);

            upert[1] += (xdist / problem::velpert_scale) * problem::velpert_amplitude *
                std::exp(-r * r / (2.0_rt * problem::velpert_scale * problem::velpert_scale)) *
                std::pow(-1.0_rt, vortex+1);

            upert[2] = 0.0_rt;

        }

            state(i,j,k,UMX) = state(i,j,k,URHO) * upert[0];
            state(i,j,k,UMY) = state(i,j,k,URHO) * upert[1];
            state(i,j,k,UMZ) = state(i,j,k,URHO) * upert[2];

            state(i,j,k,UEDEN) = state(i,j,k,UEINT) + 0.5_rt *
                                    (state(i,j,k,UMX) * state(i,j,k,UMX) +
                                    state(i,j,k,UMY) * state(i,j,k,UMY) +
                                    state(i,j,k,UMZ) * state(i,j,k,UMZ)) / state(i,j,k,URHO);

    }

    // Now we add the temperature perturbation
    if (problem::apply_temp_field){

#if AMREX_SPACEDIM == 2

        Real ydist = y - problem::y_h1 - problem::width_ratio*problem::width;
        Real delta = problem::amplitude;
        Real temp_pert = 0.0_rt;

        if (problem::num_vortices % 2 == 0){
            temp_pert = delta * (1.0_rt + std::sin(static_cast<Real>(problem::num_vortices)*M_PI*x/problem::L_x) ) *
                std::exp(-std::pow(ydist / problem::width, 2));
        } else {
            temp_pert = delta * (1.0_rt + std::cos(static_cast<Real>(problem::num_vortices)*M_PI*x/problem::L_x) ) *
                std::exp(-std::pow(ydist / problem::width, 2));
        }

        state(i,j,k,UTEMP) = state(i,j,k, UTEMP)*(1.0_rt + temp_pert);

#elif AMREX_SPACEDIM == 3

        Real zdist = z - problem::y_h1 - problem::width_ratio*problem::width;
        Real delta = problem::amplitude;
        Real temp_pert = 0.0_rt;

        if (problem::num_vortices % 2 == 0){
            temp_pert = delta * (1.0_rt + std::sin(static_cast<Real>(problem::num_vortices_x)*M_PI*x/problem::L_x) *
                                 std::cos(static_cast<Real>(problem::num_vortices_y)*M_PI*y/problem::L_y)) *
                std::exp(-std::pow(zdist / problem::width, 2));
        } else {
            temp_pert = delta * (1.0_rt + std::cos(static_cast<Real>(problem::num_vortices)*M_PI*x/problem::L_x)*
                                    std::sin(static_cast<Real>(problem::num_vortices_y)*M_PI*y/problem::L_y)) *
                std::exp(-std::pow(zdist / problem::width, 2));
        }
        state(i,j,k,UTEMP) = state(i,j,k, UTEMP)*(1.0_rt + temp_pert);
#endif

        // Now we enforce the EOS by calling T and p (giving us rho)
        eos_state.T = state(i,j,k,UTEMP);

        eos(eos_input_tp, eos_state);

        state(i,j,k,URHO) = eos_state.rho;

        for (int n = 0; n < NumSpec; n++) {
            state(i,j,k,UFS+n) = eos_state.rho * eos_state.xn[n];
        }

        state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state.e;
        state(i,j,k,UEDEN) = state(i,j,k,UEINT);

    }

}
#endif
