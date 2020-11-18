#ifndef problem_setup_H
#define problem_setup_H

#include <prob_parameters.H>
#include <eos.H>
#include <network.H>
#include <model_parser.H>
#include <ambient.H>

AMREX_INLINE
void problem_initialize ()
{

    const Geometry& dgeom = DefaultGeometry();

    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    // set center variable

#if AMREX_SPACEDIM == 2
    // axisymmetric
    problem::center[0] = 0.0_rt;
    problem::center[1] = 0.5_rt * (problo[1] + probhi[1]);
    problem::center[2] = 0.0_rt;

#elif AMREX_SPACEDIM == 3
    problem::center[0] = 0.5_rt * (problo[0] + probhi[0]);
    problem::center[1] = 0.5_rt * (problo[1] + probhi[1]);
    problem::center[2] = 0.5_rt * (problo[2] + probhi[2]);
#endif


    problem::ihe4 = network_spec_index("helium-4");

    if (problem::ihe4 < 0) {
        amrex::Error("Error: no He4 in the network");
    }

    // Read initial model

    read_model_file(problem::model_name);

}


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
    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = interpolate(dist, model::ispec+n);
    }

    eos_t eos_state;

    eos_state.rho = state(i,j,k,URHO);
    eos_state.T = state(i,j,k,UTEMP);
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state(i,j,k,UFS+n);
    }

    eos(eos_input_rt, eos_state);

    state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state.e;
    state(i,j,k,UEDEN) = state(i,j,k,URHO) * eos_state.e;

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * state(i,j,k,UFS+n);
    }


    // initial velocities = 0

    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

    // add a perturbation

    Real t0 = state(i,j,k,UTEMP);

    // perturbation is on the vertical-axis

#if AMREX_SPACEDIM == 1
    Real r1 = std::sqrt((x - problem::R_pert) * (x - problem::R_pert)) /
        (2.5e6_rt * problem::pert_rad_factor);
#elif AMREX_SPACEDIM == 2
    Real r1 = std::sqrt(x * x + (y - problem::R_pert) * (y - problem::R_pert)) /
        (2.5e6_rt * problem::pert_rad_factor);
#else
    Real r1 = std::sqrt(x * x + y * y + (z - problem::R_pert) * (z - problem::R_pert)) /
        (2.5e6_rt * problem::pert_rad_factor);
#endif

    // convolve the temperature perturbation with the amount of He
    Real X_he = eos_state.xn[problem::ihe4];

    state(i,j,k,UTEMP) = t0 * (1.0_rt + X_he * problem::pert_temp_factor *
                               (0.150e0_rt * (1.0_rt + std::tanh(2.0_rt - r1))));


    eos_state.rho = state(i,j,k,URHO);
    eos_state.T = state(i,j,k,UTEMP);
    // we don't need to refill xn, since it still holds unchanged from above

    eos(eos_input_rt, eos_state);

    // the internal energy changed

    state(i,j,k,UEINT) = eos_state.e * state(i,j,k,URHO);
    state(i,j,k,UEDEN) = eos_state.e * state(i,j,k,URHO);

}
#endif
