#ifndef SDC_CONS_TO_BURN_H
#define SDC_CONS_TO_BURN_H

#include <burn_type.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
copy_cons_to_burn_type(const int i, const int j, const int k,
                       Array4<const Real> const& state,
                       burn_t& burn_state) {

    const auto U_cell = state.cellData(i,j,k);
    burn_state.y[SRHO] = U_cell[URHO];
    burn_state.y[SMX] = U_cell[UMX];
    burn_state.y[SMY] = U_cell[UMY];
    burn_state.y[SMZ] = U_cell[UMZ];
    burn_state.y[SEDEN] = U_cell[UEDEN];
    burn_state.y[SEINT] = U_cell[UEINT];
    amrex::constexpr_for<0, NumSpec>([&] (auto n) {
        burn_state.y[SFS+n] = U_cell[UFS+n];
    });
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; n++) {
        burn_state.y[SFX+n] = state(i,j,k,UFX+n);
    }
#endif

    burn_state.T = U_cell[UTEMP];
    burn_state.rho = burn_state.y[SRHO];

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
