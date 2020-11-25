#ifndef problem_setup_H
#define problem_setup_H

#include <prob_parameters.H>
#include <eos.H>
#include <model_parser.H>
#include <ambient.H>

AMREX_INLINE
void problem_initialize ()
{

    // Read initial model

    read_model_file(problem::model_name);

    // set the ambient state for the upper boundary condition

    ambient::ambient_state[URHO] = model::profile(0).state(model::npts-1, model::idens);
    ambient::ambient_state[UTEMP] = model::profile(0).state(model::npts-1, model::itemp);
    for (int n = 0; n < NumSpec; n++) {
        ambient::ambient_state[UFS+n] = ambient::ambient_state[URHO] * model::profile(0).state(model::npts-1, model::ispec+n);
    }

    ambient::ambient_state[UMX] = 0.0_rt;
    ambient::ambient_state[UMY] = 0.0_rt;
    ambient::ambient_state[UMZ] = 0.0_rt;

    // make the ambient state thermodynamically consistent

    eos_t eos_state;
    eos_state.rho = ambient::ambient_state[URHO];
    eos_state.T = ambient::ambient_state[UTEMP];
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = ambient::ambient_state[UFS+n] / eos_state.rho;
    }

    eos(eos_input_rt, eos_state);

    ambient::ambient_state[UEINT] = eos_state.rho * eos_state.e;
    ambient::ambient_state[UEDEN] = eos_state.rho * eos_state.e;

}


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

    // Now add the perturbation

    Real t0 = state(i,j,k,UTEMP);

    Real x1 = 5.0e7_rt;

#if AMREX_SPACEDIM == 2
    Real y1 = 6.5e7_rt;
    Real z1 = 0.0_rt;
#else
    Real y1 = 5.0e7_rt;
    Real z1 = 6.5e7_rt;
#endif
    Real r1 = std::sqrt((x - x1) * (x - x1) +
                        (y - y1) * (y - y1) +
                        (z - z1) * (z - z1)) /
        (2.5e6_rt * problem::pert_rad_factor);

    Real x2 = 1.2e8_rt;

#if AMREX_SPACEDIM == 2
    Real y2 = 8.5e7_rt;
    Real z2 = 0.0_rt;
#else
    Real y2 = 1.2e8_rt;
    Real z2 = 8.5e7_rt;
#endif

    Real r2 = std::sqrt((x - x2) * (x - x2) +
                        (y - y2) * (y - y2) +
                        (z - z2) * (z - z2)) /
        (2.5e6_rt * problem::pert_rad_factor);

    Real x3 = 2.0e8_rt;

#if AMREX_SPACEDIM == 2
    Real y3 = 7.5e7_rt;
    Real z3 = 0.0_rt;
#else
    Real y3 = 2.0e8_rt;
    Real z3 = 7.5e7_rt;
#endif

    Real r3 = std::sqrt((x - x3) * (x - x3) +
                        (y - y3) * (y - y3) +
                        (z - z3) * (z - z3)) /
        (2.5e6_rt * problem::pert_rad_factor);

    state(i,j,k,UTEMP) = t0 * (1.0_rt + problem::pert_temp_factor *
                               (0.150e0_rt * (1.0_rt + std::tanh(2.0_rt - r1)) +
                                0.300e0_rt * (1.0_rt + std::tanh(2.0_rt - r2)) +
                                0.225e0_rt * (1.0_rt + std::tanh(2.0_rt - r3))));

    // we can reuse the EOS state from above -- the pressure and
    // composition are unchanged so we'll just update the temperature
    // and then ask for the corresponding density

    eos_state.T = state(i,j,k,UTEMP);

    eos(eos_input_tp, eos_state);

    state(i,j,k,URHO) = eos_state.rho;

    // reset the species densities with the updated density
    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = eos_state.rho * eos_state.xn[n];
    }

    state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state.e;
    state(i,j,k,UEDEN) = state(i,j,k,UEINT);
}

#endif
