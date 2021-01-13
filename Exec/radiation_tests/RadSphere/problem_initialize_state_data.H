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

    Real xcen = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

    state(i,j,k,URHO ) = problem::rho_0;
    state(i,j,k,UMX  ) = 0.0_rt;
    state(i,j,k,UEDEN) = problem::rhoe_0;
    state(i,j,k,UEINT) = problem::rhoe_0;

    // set the composition to be all in the first species

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = 0.0_rt;
    }
    state(i,j,k,UFS) = state(i,j,k,URHO);

#if NAUX_NET > 0
    for (int n = 0; n < NumAux; n++) {
        state(i,j,k,UFX+n) = 0.0_rt;
    }
    state(i,j,k,UFX) = state(i,j,k,URHO);
#endif

    state(i,j,k,UTEMP) = problem::T_0;
}
#endif
