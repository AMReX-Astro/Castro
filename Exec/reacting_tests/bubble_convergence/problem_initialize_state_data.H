#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>
#include <model_parser.H>
#include <initial_model.H>
#include <network.H>

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

#if AMREX_SPACEDIM == 3
    Real z = 0.0;
    z = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
#endif

#if AMREX_SPACEDIM == 1
    Real height = x;
#elif AMREX_SPACEDIM == 2
    Real height = y;
#else
    Real height = z;
#endif

    state(i,j,k,URHO) = interpolate(height, model::idens);
    state(i,j,k,UTEMP) = interpolate(height, model::itemp);
    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = interpolate(height, model::ispec+n);
    }

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

    // Initial velocities = 0

    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

    if (problem::do_pert) {

        // Now add the perturbation
        Real t0 = state(i,j,k,UTEMP);

        Real rho_old = state(i,j,k,URHO);

#if AMREX_SPACEDIM == 1
        Real r = std::abs(x - problem::center[0]) /  problem::pert_width;
#elif AMREX_SPACEDIM == 2
        Real r = std::sqrt(std::pow(x - problem::center[0], 2) +
                           std::pow(y - problem::center[1], 2)) /  problem::pert_width;
#else
        Real r = std::sqrt(std::pow(x - problem::center[0], 2) +
                           std::pow(y - problem::center[1], 2) +
                           std::pow(z - problem::center[2], 2)) /  problem::pert_width;
#endif

        state(i,j,k,UTEMP) = t0 * (1.0_rt + 0.6_rt * (1.0_rt + std::tanh(4.0_rt - r)));

        // update the temperature in the EOS state -- leave the pressure unchanged

        eos_state.T = state(i,j,k,UTEMP);

        // now get the new density from this T, p

        eos(eos_input_tp, eos_state);

        state(i,j,k,URHO) = eos_state.rho;

        // correct the mass fractions and energy with the new density

        state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state.e;
        state(i,j,k,UEDEN) = state(i,j,k,UEINT);

        for (int n = 0; n < NumSpec; n++) {
            state(i,j,k,UFS+n) = (state(i,j,k,URHO) / rho_old) * state(i,j,k,UFS+n);
        }
    }
}
#endif
