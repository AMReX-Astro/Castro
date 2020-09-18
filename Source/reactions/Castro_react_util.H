#ifndef CASTRO_REACT_UTIL_H
#define CASTRO_REACT_UTIL_H

#include <Castro.H>
#include <burn_type.H>
#include <actual_rhs.H>


AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool
okay_to_burn(const int i, const int j, const int k,
             Array4<const Real> const& state) {

    if (state(i,j,k,UTEMP) < react_T_min || state(i,j,k,UTEMP) > react_T_max ||
        state(i,j,k,URHO) < react_rho_min || state(i,j,k,URHO) > react_rho_max) {
            return false;
    }

    return true;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool
okay_to_burn_type(burn_t const& state) {

    if (state.T < react_T_min || state.T > react_T_max ||
        state.rho < react_rho_min || state.rho > react_rho_max) {
            return false;
    }

    return true;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
single_zone_react_source(const int i, const int j, const int k,
             Array4<const Real> const& state,
             Array4<Real> const& R,
             burn_t burn_state) {

    auto rhoInv = 1.0_rt / state(i,j,k,URHO);

    burn_state.rho = state(i,j,k,URHO);
    burn_state.T = state(i,j,k,UTEMP);
    burn_state.e = state(i,j,k,UEINT) * rhoInv;

    for (int n = 0; n < NumSpec; ++n) {
       burn_state.xn[n] = amrex::max(amrex::min(state(i,j,k,UFS+n) * rhoInv, 1.0_rt), small_x);
    }

#if NAUX_NET > 0
    for (int n = 0; n < NumAux; ++n) {
       burn_state.aux[n] = state(i,j,k,UFX+n) * rhoInv;
    }
#endif

    eos_t eos_state;

    // Ensure that the temperature going in is consistent with the internal energy.
    burn_to_eos(burn_state, eos_state);
    eos(eos_input_re, eos_state);
    eos_to_burn(eos_state, burn_state);

    eos_get_small_temp(&small_temp);
    burn_state.T = amrex::min(MAX_TEMP, amrex::max(burn_state.T, small_temp));

    burn_state.self_heat = false;

    Array1D<Real, 1, neqs> ydot;

    actual_rhs(burn_state, ydot);

    // store the instantaneous R
    for (int n = 0; n < NUM_STATE; ++n) {
        R(i,j,k,n) = 0.0_rt;
    }

    // species rates come back in terms of molar fractions
    for (int n = 0; n < NumSpec; ++n) {
        R(i,j,k,UFS+n) = state(i,j,k,URHO) * aion[n] * ydot(n);
    }

    R(i,j,k,UEDEN) = state(i,j,k,URHO) * ydot(net_ienuc);
    R(i,j,k,UEINT) = state(i,j,k,URHO) * ydot(net_ienuc);
}

#endif