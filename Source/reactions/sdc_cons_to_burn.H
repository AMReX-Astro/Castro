#ifndef SDC_CONS_TO_BURN_H
#define SDC_CONS_TO_BURN_H

#include <burn_type.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
copy_cons_to_burn_type(const int i, const int j, const int k,
                       Array4<const Real> const& state,
                       burn_t& burn_state) {

    burn_state.y[SRHO] = state(i,j,k,URHO);
    burn_state.y[SMX] = state(i,j,k,UMX);
    burn_state.y[SMY] = state(i,j,k,UMY);
    burn_state.y[SMZ] = state(i,j,k,UMZ);
    burn_state.y[SEDEN] = state(i,j,k,UEDEN);
    burn_state.y[SEINT] = state(i,j,k,UEINT);
    for (int n = 0; n < NumSpec; n++) {
        burn_state.y[SFS+n] = state(i,j,k,UFS+n);
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; n++) {
        burn_state.y[SFX+n] = state(i,j,k,UFX+n);
    }
#endif

    burn_state.T = state(i,j,k,UTEMP);
    burn_state.rho = state(i,j,k,URHO);

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
copy_cons_to_burn_type(GpuArray<Real, NUM_STATE> const& state,
                       burn_t& burn_state) {

    burn_state.y[SRHO] = state[URHO];
    burn_state.y[SMX] = state[UMX];
    burn_state.y[SMY] = state[UMY];
    burn_state.y[SMZ] = state[UMZ];
    burn_state.y[SEDEN] = state[UEDEN];
    burn_state.y[SEINT] = state[UEINT];
    for (int n = 0; n < NumSpec; n++) {
        burn_state.y[SFS+n] = state[UFS+n];
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; n++) {
        burn_state.y[SFX+n] = state[UFX+n];
    }
#endif

    burn_state.T = state[UTEMP];
    burn_state.rho = state[URHO];

}



#endif
