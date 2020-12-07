#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

    state(i,j,k,URHO) = 1.0_rt;
    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

    // set the composition to be all in the first species
    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = 0.0_rt;
    }
    state(i,j,k,UFS) = state(i,j,k,URHO);

    // Set temperature and energy to arbitrary, positive values
    // so that the Castro state checkers are OK.
    state(i,j,k,UTEMP) = 1.0e-50_rt;
    state(i,j,k,UEINT) = 1.0e-50_rt;
    state(i,j,k,UEDEN) = 1.0e-50_rt;
}
#endif
