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

    // Provides the simulation to be run in the x,y,or z direction
    // where length direction is the length side in a square prism

    Real length_cell{};
    if (problem::idir == 1) {
        length_cell = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);
    } else if (problem::idir == 2) {
        length_cell = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
    } else if (problem::idir == 3) {
        length_cell = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
    } else {
#ifndef AMREX_USE_GPU
        amrex::Error("Invalid direction please input idir = [1,3]");
#endif
    }

    if (length_cell < 0.e0_rt) {
        state(i,j,k,URHO) = problem::rho0;

        // set the composition to be all in the first species

        for (int n = 0; n < NumSpec; ++n) {
            state(i,j,k,UFS+n) = 0.0_rt;
        }
        state(i,j,k,UFS) = state(i,j,k,URHO);

        state(i,j,k,UTEMP) = problem::T0;

        if (problem::idir == 1) {
            state(i,j,k,UMX) = problem::rho0 * problem::v0;
            state(i,j,k,UMY) = 0.0_rt;
            state(i,j,k,UMZ) = 0.0_rt;

        } else if (problem::idir == 2) {
            state(i,j,k,UMX) = 0.0_rt;
            state(i,j,k,UMY) = problem::rho0 * problem::v0;
            state(i,j,k,UMZ) = 0.0_rt;

        } else if (problem::idir == 3) {
            state(i,j,k,UMX) = 0.0_rt;
            state(i,j,k,UMY) = 0.0_rt;
            state(i,j,k,UMZ) = problem::rho0 * problem::v0;
        }
    } else {
         state(i,j,k,URHO) = problem::rho1;

         //  set the composition to be all in the first species

         for (int n = 0; n < NumSpec; ++n) {
             state(i,j,k,UFS+n) = 0.0_rt;
         }
         state(i,j,k,UFS) = state(i,j,k,URHO);

         state(i,j,k,UTEMP) = problem::T1;

         if (problem::idir == 1) {
             state(i,j,k,UMX) = problem::rho1 * problem::v1;
             state(i,j,k,UMY) = 0.0_rt;
             state(i,j,k,UMZ) = 0.0_rt;

         } else if (problem::idir == 2) {
             state(i,j,k,UMX) = 0.0_rt;
             state(i,j,k,UMY) = problem::rho1 * problem::v1;
             state(i,j,k,UMZ) = 0.0_rt;

         } else if (problem::idir == 3) {
             state(i,j,k,UMX) = 0.0_rt;
             state(i,j,k,UMY) = 0.0_rt;
             state(i,j,k,UMZ) = problem::rho1 * problem::v1;
         }

    }


#if NAUX_NET > 0
    for (int n = 0; n < NumAux; ++n) {
        state(i,j,k,UFX+n) = 0.0_rt;
    }
    state(i,j,k,UFX) = state(i,j,k,URHO);
#endif

    eos_t eos_state;
    eos_state.rho = state(i,j,k,URHO);
    eos_state.T = state(i,j,k,UTEMP);

    Real rhoInv = 1.0_rt / state(i,j,k,URHO);
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n]  = state(i,j,k,UFS+n) * rhoInv;
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; ++n) {
        eos_state.aux[n] = state(i,j,k,UFX+n) * rhoInv;
    }
#endif

    eos(eos_input_rt, eos_state);

    state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state.e;
    state(i,j,k,UEDEN) = state(i,j,k,UEINT) +
        0.5_rt*(state(i,j,k,UMX) * state(i,j,k,UMX) +
                state(i,j,k,UMY) * state(i,j,k,UMY) +
                state(i,j,k,UMZ) * state(i,j,k,UMZ)) / state(i,j,k,URHO);
}
#endif
