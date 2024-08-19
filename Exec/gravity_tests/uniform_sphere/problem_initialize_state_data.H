#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k, Array4<Real> const& state, const GeometryData& geomdata)
{

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real xx = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt) - problem::center[0];
    Real yy = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt) - problem::center[1];
    Real zz = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt) - problem::center[2];

    // Establish the sphere

    if (std::pow(xx * xx + yy * yy + zz * zz, 0.5_rt) < problem::diameter / 2) {
        state(i,j,k,URHO) = problem::density;
    }
    else {
        state(i,j,k,URHO) = problem::ambient_dens;
    }

    // Establish the thermodynamic quantities. They don't have to be
    // valid because this test will never do a hydro step.

    state(i,j,k,UTEMP) = 1.0_rt;
    state(i,j,k,UEINT) = 1.0_rt;
    state(i,j,k,UEDEN) = 1.0_rt;

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) / static_cast<Real>(NumSpec);
    }
}

#endif
