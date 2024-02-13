#ifndef problem_bc_fill_H
#define problem_bc_fill_H

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_bc_fill(int i, int j, int k,
                     Array4<Real> const& state,
                     Real time,
                     const Array1D<BCRec, 0, NUM_STATE-1>& bcs,
                     const GeometryData& geomdata)
{
    // Set the velocity to zero at the top physical boundary.

    const int* domhi = geomdata.Domain().hiVect();

    if (AMREX_SPACEDIM == 1) {

        if (bcs(UMX).hi(0) == amrex::BCType::foextrap && i > domhi[0]) {

            Real vel = amrex::max(state(i,j,k,UMX) / state(i,j,k,URHO), 0.0);
            state(i,j,k,URHO)  = problem::hse_rho_top;
            state(i,j,k,UMX)   = state(i,j,k,URHO) * vel;
            state(i,j,k,UMY)   = 0.0;
            state(i,j,k,UMZ)   = 0.0;
            state(i,j,k,UTEMP) = problem::hse_T_top;
            state(i,j,k,UEINT) = problem::hse_rho_top * problem::hse_eint_top;
            state(i,j,k,UEDEN) = problem::hse_rho_top * problem::hse_eint_top +
                                 0.5_rt * (state(i,j,k,UMX) * state(i,j,k,UMX) +
                                           state(i,j,k,UMY) * state(i,j,k,UMY) +
                                           state(i,j,k,UMZ) * state(i,j,k,UMZ)) / state(i,j,k,URHO);
            for (int n = 0; n < NumSpec; ++n) {
                state(i,j,k,UFS+n) = problem::hse_rho_top * problem::hse_X_top[n];
            }

        }

    }
    else if (AMREX_SPACEDIM == 2) {

        if (bcs(UMY).hi(1) == amrex::BCType::foextrap && j > domhi[1]) {

            Real vel = amrex::max(state(i,j,k,UMY) / state(i,j,k,URHO), 0.0);
            state(i,j,k,URHO)  = problem::hse_rho_top;
            state(i,j,k,UMX)   = 0.0;
            state(i,j,k,UMY)   = state(i,j,k,URHO) * vel;
            state(i,j,k,UMZ)   = 0.0;
            state(i,j,k,UTEMP) = problem::hse_T_top;
            state(i,j,k,UEINT) = problem::hse_rho_top * problem::hse_eint_top;
            state(i,j,k,UEDEN) = problem::hse_rho_top * problem::hse_eint_top +
                                 0.5_rt * (state(i,j,k,UMX) * state(i,j,k,UMX) +
                                           state(i,j,k,UMY) * state(i,j,k,UMY) +
                                           state(i,j,k,UMZ) * state(i,j,k,UMZ)) / state(i,j,k,URHO);
            for (int n = 0; n < NumSpec; ++n) {
                state(i,j,k,UFS+n) = problem::hse_rho_top * problem::hse_X_top[n];
            }

        }

    }
    else {

        if (bcs(UMZ).hi(2) == amrex::BCType::foextrap && k > domhi[2]) {

            Real vel = amrex::max(state(i,j,k,UMZ) / state(i,j,k,URHO), 0.0);
            state(i,j,k,URHO)  = problem::hse_rho_top;
            state(i,j,k,UMX)   = 0.0;
            state(i,j,k,UMY)   = 0.0;
            state(i,j,k,UMZ)   = state(i,j,k,URHO) * vel;
            state(i,j,k,UTEMP) = problem::hse_T_top;
            state(i,j,k,UEINT) = problem::hse_rho_top * problem::hse_eint_top;
            state(i,j,k,UEDEN) = problem::hse_rho_top * problem::hse_eint_top +
                                 0.5_rt * (state(i,j,k,UMX) * state(i,j,k,UMX) +
                                           state(i,j,k,UMY) * state(i,j,k,UMY) +
                                           state(i,j,k,UMZ) * state(i,j,k,UMZ)) / state(i,j,k,URHO);
            for (int n = 0; n < NumSpec; ++n) {
                state(i,j,k,UFS+n) = problem::hse_rho_top * problem::hse_X_top[n];
            }

        }

    }
}

#endif
