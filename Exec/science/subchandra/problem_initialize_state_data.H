#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>
#include <network.H>
#include <model_parser.H>
#include <ambient.H>
#include <burn_type.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt) - problem::center[0];

    Real y = 0.0;
#if AMREX_SPACEDIM >= 2
    y = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt) - problem::center[1];
#endif

    Real z = 0.0;
#if AMREX_SPACEDIM == 3
    z = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt) - problem::center[2];
#endif

    Real dist = sqrt(x*x + y*y + z*z);

    state(i,j,k,URHO) = interpolate(dist, model::idens);
    state(i,j,k,UTEMP) = interpolate(dist, model::itemp);
    Real sumX{0.0_rt};
    for (int n = 0; n < NumSpec; ++n) {
        state(i,j,k,UFS+n) = amrex::min(1.0_rt, amrex::max(interpolate(dist, model::ispec+n), network_rp::small_x));
        sumX += state(i,j,k,UFS+n);
    }

    // normalize
    for (int n = 0; n < NumSpec; ++n) {
        state(i,j,k,UFS+n) /= sumX;
    }

    burn_t burn_state;

    burn_state.rho = state(i,j,k,URHO);
    burn_state.T = state(i,j,k,UTEMP);
    for (int n = 0; n < NumSpec; n++) {
        burn_state.xn[n] = state(i,j,k,UFS+n);
    }

#ifdef AUX_THERMO
    // Set AUX variables

    set_aux_comp_from_X(burn_state);
    for (int n = 0; n < NumAux; n++) {
        state(i,j,k,UFX+n) = state(i,j,k,URHO) * burn_state.aux[n];
    }
#endif

    eos(eos_input_rt, burn_state);

    state(i,j,k,UEINT) = state(i,j,k,URHO) * burn_state.e;
    state(i,j,k,UEDEN) = state(i,j,k,URHO) * burn_state.e;

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * state(i,j,k,UFS+n);
    }

#ifdef NSE_NET
    state(i,j,k,UMUP) = problem::mu_p;
    state(i,j,k,UMUN) = problem::mu_n;
#endif

    // initial velocities = 0
    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

    // add a perturbation at the north pole

    Real T0 = state(i,j,k,UTEMP);

    // perturbation is on the vertical-axis

    Real pert_center = problem::R_pert + problem::R_He_base;

#if AMREX_SPACEDIM == 1
    Real r1 = std::sqrt((x - pert_center) * (x - pert_center)) /
        (2.5e6_rt * problem::pert_rad_factor);
#elif AMREX_SPACEDIM == 2
    Real r1 = std::sqrt(x * x + (y - pert_center) * (y - pert_center)) /
        (2.5e6_rt * problem::pert_rad_factor);
#else
    Real r1 = std::sqrt(x * x + y * y + (z - pert_center) * (z - pert_center)) /
        (2.5e6_rt * problem::pert_rad_factor);
#endif

    // convolve the temperature perturbation with the amount of He
    Real X_he = burn_state.xn[problem::ihe4];

    Real Tpert = T0 * (1.0_rt + X_he * problem::pert_temp_factor *
                       (0.150e0_rt * (1.0_rt + std::tanh(2.0_rt - r1))));

    Real dT = Tpert - T0;


#if AMREX_SPACEDIM == 3
    // optional second perturbation

    if (problem::second_pert_angle > 0) {

        Real angle_rad = problem::second_pert_angle * M_PI / 180.0_rt;

        Real second_pert_center = problem::second_R_pert + problem::R_He_base;

        // these are measured with respect to the center

        Real x_pert = second_pert_center * std::sin(angle_rad);
        Real y_pert = 0.0_rt;
        Real z_pert = second_pert_center * std::cos(angle_rad);

        Real dx_pert = x - x_pert;
        Real dy_pert = y - y_pert;
        Real dz_pert = z - z_pert;

        Real r2 = std::sqrt(dx_pert * dx_pert + dy_pert * dy_pert + dz_pert * dz_pert) /
            (2.5e6_rt * problem::second_pert_rad_factor);

        // we are assuming here that the 2 perturbations don't overlap

        Tpert = T0 * (1.0_rt + X_he * problem::second_pert_temp_factor *
                      (0.150e0_rt * (1.0_rt + std::tanh(2.0_rt - r2))));

        dT += (Tpert - T0);

    }
#endif

    burn_state.rho = state(i,j,k,URHO);
    burn_state.T = T0 + dT;

    // we don't need to refill xn, since it still holds unchanged from above

    eos(eos_input_rt, burn_state);

    // the internal energy changed

    state(i,j,k,UEINT) = burn_state.e * state(i,j,k,URHO);
    state(i,j,k,UEDEN) = burn_state.e * state(i,j,k,URHO);

}

#endif
