#ifndef CASTRO_REACT_UTIL_H
#define CASTRO_REACT_UTIL_H

#include <Castro.H>
#include <burn_type.H>
#include <actual_network.H> //this is needed to get the definition of NEW_NETWORK_IMPLEMENTATION
#ifdef NEW_NETWORK_IMPLEMENTATION
#include <rhs.H>
#endif
#include <numerical_jacobian.H>


// #include <actual_rhs.H>

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
okay_to_burn(GpuArray<Real, NUM_STATE> const& state) {

    if (state[UTEMP] < react_T_min || state[UTEMP] > react_T_max ||
        state[URHO] < react_rho_min || state[URHO] > react_rho_max) {
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

#endif
