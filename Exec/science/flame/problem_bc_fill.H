#ifndef problem_bc_fill_H
#define problem_bc_fill_H

#include <AMReX_Real.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_bc_fill(int i, int j, int k,
                     Array4<Real> const& state,
                     Real time,
                     Array1D<BCRec, 0, NUM_STATE-1> bcs,
                     GeometryData geomdata)
{

    // implement an inflow BC at the right interface

    const int* domlo = geomdata.Domain().loVect();
    const int* domhi = geomdata.Domain().hiVect();

    // Note: in these checks we're only looking at the boundary
    // conditions on the first state component (density).

    if (i > domhi[0] && (bcs(0).hi(0) == EXT_DIR)) {

        state(i,j,k,URHO) = problem::rho_fuel;
        state(i,j,k,UTEMP) = problem::T_fuel;
        state(i,j,k,UEINT) = problem::rho_fuel * problem::e_fuel;
        for (int n = 0; n < NumSpec; ++n) {
            state(i,j,k,UFS+n) = problem::rho_fuel * problem::xn_fuel[n];
        }

        state(i,j,k,UMX) = problem::mass_flux;
        state(i,j,k,UMY) = 0.0_rt;
        state(i,j,k,UMZ) = 0.0_rt;

        state(i,j,k,UEDEN) = state(i,j,k,UEINT) +
            0.5_rt * (state(i,j,k,UMX) * state(i,j,k,UMX) +
                      state(i,j,k,UMY) * state(i,j,k,UMY) +
                      state(i,j,k,UMZ) * state(i,j,k,UMZ)) /
            state(i,j,k,URHO);

    }

}

#endif
