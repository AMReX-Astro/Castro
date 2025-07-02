#ifndef SDC_REACT_UTIL_H
#define SDC_REACT_UTIL_H

#include <burn_type.H>
#include <sdc_cons_to_burn.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
single_zone_react_source(burn_t& burn_state,
                         GpuArray<Real, NUM_STATE>& R) {

    // note: we want this burn_state to come out of here so we can
    // reuse it elsewhere

    // we assume that burn_state.y[] holds the current state, except
    // for temperature guess, and we will initialize the pieces needed
    // to call the RHS

    auto rhoInv = 1.0_rt / burn_state.y[SRHO];

    burn_state.rho = burn_state.y[SRHO];
    burn_state.e = burn_state.y[SEINT] * rhoInv;

    for (int n = 0; n < NumSpec; ++n) {
       burn_state.xn[n] = amrex::max(amrex::min(burn_state.y[SFS+n] * rhoInv, 1.0_rt), small_x);
    }

#if NAUX_NET > 0
    for (int n = 0; n < NumAux; ++n) {
       burn_state.aux[n] = burn_state.y[SFX+n] * rhoInv;
    }
#endif

    // Ensure that the temperature going in is consistent with the internal energy.
    eos(eos_input_re, burn_state);

    // eos_get_small_temp(&small_temp);
    burn_state.T = amrex::Clamp(burn_state.T, small_temp, integrator_rp::MAX_TEMP);

    Array1D<Real, 1, neqs> ydot;

    actual_rhs(burn_state, ydot);

    // store the instantaneous R
    for (int n = 0; n < NUM_STATE; ++n) {
        R[n] = 0.0_rt;
    }

    // species rates come back in terms of molar fractions
    for (int n = 0; n < NumSpec; ++n) {
        R[UFS+n] = burn_state.rho * aion[n] * ydot(1+n);
    }

    R[UEDEN] = burn_state.rho * ydot(net_ienuc);
    R[UEINT] = burn_state.rho * ydot(net_ienuc);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
single_zone_react_source(const int i, const int j, const int k,
                         Array4<const Real> const& state,
                         Array4<Real> const& R) {

    GpuArray<Real, NUM_STATE> R_arr;

    burn_t burn_state;

    copy_cons_to_burn_type(i, j, k, state, burn_state);

    single_zone_react_source(burn_state, R_arr);

    for (int n = 0; n < NUM_STATE; ++n) {
        R(i,j,k,n) = R_arr[n];
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
single_zone_jac(burn_t& burn_state, const Real dt,
                JacNetArray2D& Jac) {

#ifdef SIMPLIFIED_SDC
#ifndef AMREX_USE_GPU
    amrex::Error("we shouldn't be here with the simplified SDC method (USE_SIMPLIFIED_SDC=TRUE)");
#endif
#else

    // Jac has the derivatives with respect to the native
    // network variables, X, e.  Note: the e derivative

    if (integrator_rp::jacobian == 2) {
        // note the numerical Jacobian will be returned in terms of X
        // and will already have the corrections to convert into d/de
        // instead of d/dT

        jac_info_t jac_info;
        jac_info.h = dt;
        numerical_jac(burn_state, jac_info, Jac);

    } else {
        actual_jac(burn_state, Jac);

        // the above will have already multiplied the d/dT terms with 1/c_v

        // The Jacobian from the nets is in terms of dYdot/dY, but we want
        // it was dXdot/dX, so convert here.
        for (int n = 1; n <= NumSpec; n++) {
            for (int m = 1; m <= neqs; m++) {
                // d(X_n)/dq_m
                Jac(n,m) = Jac(n,m) * aion[n-1];
            }
        }

        for (int m = 1; m <= neqs; m++) {
            for (int n = 1; n <= NumSpec; n++) {
                // d RHS_m / dX_n
                Jac(m,n) = Jac(m,n) * aion_inv[n-1];
            }
        }

        // now correct the species derivatives
        // this constructs dy/dX_k |_e = dy/dX_k |_T - e_{X_k} |_T dy/dT / c_v

        // TODO: I think that this can be accomplished via burn_state directly!

        eos_re_extra_t eos_state;
        eos_state.rho = burn_state.rho;
        eos_state.T = burn_state.T;
        eos_state.e = burn_state.e;
        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = burn_state.xn[n];
        }

        eos(eos_input_re, eos_state);
        eos_xderivs_t eos_xderivs = composition_derivatives(eos_state);

        for (int m = 1; m <= neqs; m++) {
            for (int n = 1; n <= NumSpec; n++) {
                Jac(m, n) -= eos_xderivs.dedX[n-1] * Jac(m, net_ienuc);
            }
        }
    }
#endif

    // at this point, our Jacobian should be entirely in terms of X,
    // not Y.

}

#endif
