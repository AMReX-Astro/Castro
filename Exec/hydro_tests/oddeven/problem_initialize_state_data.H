#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{
    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real xx = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

    // compute the integer location of the center of the domain
    int icen = static_cast<int>(problem::center[0] / dx[0]);

    int jcen = 0;
#if AMREX_SPACEDIM >= 2
    jcen = static_cast<int>(problem::center[1] / dx[1]);
#endif

    int kcen = 0;
#if AMREX_SPACEDIM == 3
    kcen = static_cast<int>(problem::center[2] / dx[2]);
#endif

    Real dens = problem::dens_ambient;
    if (i == icen && j == jcen && k == kcen) {
        dens = problem::dens_ambient * problem::dens_pert_factor;
    }

    state(i,j,k,URHO) = dens;

    // velocity perturbation

    Real xvel = 0.0;
    if (xx < problem::center[0]) {
        xvel = problem::vel_pert;
    } else if (xx > problem::center[0]) {
        xvel = -problem::vel_pert;
    }

    state(i,j,k,UMX) = dens * xvel;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

    // set the composition
    Real X[NumSpec] = {0.0};
    X[0] = 1.0_rt;

    // compute the internal energy and temperature

    eos_t eos_state;
    eos_state.T = 1.0_rt;  // initial guess
    eos_state.rho = dens;
    eos_state.p = problem::p_ambient;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = X[n];
    }

    eos(eos_input_rp, eos_state);

    state(i,j,k,UTEMP) = eos_state.T;

    state(i,j,k,UEDEN) = dens * eos_state.e +
        0.5_rt * (state(i,j,k,UMX) * state(i,j,k,UMX) / state(i,j,k,URHO) +
                  state(i,j,k,UMY) * state(i,j,k,UMY) / state(i,j,k,URHO) +
                  state(i,j,k,UMZ) * state(i,j,k,UMZ) / state(i,j,k,URHO));

    state(i,j,k,UEINT) = dens * eos_state.e;

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = dens * X[n];
    }
}

#endif
