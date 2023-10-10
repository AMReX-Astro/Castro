#ifndef VODE_RHS_TRUE_SDC_H
#define VODE_RHS_TRUE_SDC_H

#include <network.H>
#include <burn_type.H>
#include <extern_parameters.H>

#ifdef NEW_NETWORK_IMPLEMENTATION
#include <rhs.H>
#else
#include <actual_network.H>
#include <actual_rhs.H>
#endif
#include <Castro_react_util.H>
#include <sdc_cons_to_burn.H>
#include <burner.H>


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
sdc_vode_solve(const Real dt_m,
               GpuArray<Real, NUM_STATE> const& U_old,
               GpuArray<Real, NUM_STATE>& U_new,
               GpuArray<Real, NUM_STATE> const& C,
               const int sdc_iteration) {

    // The purpose of this function is to solve the system the
    // approximate system dU/dt = R + C using the VODE ODE
    // integrator.
    // The solution we get here will then be used as the
    // initial guess to the Newton solve on the real system.

    // We will do the implicit update of only the terms that
    // have reactive sources
    //
    // 1:NumSpec : (rho X)
    // NumSpec+1 : (rho e)


    // Update the density and momenta for this zone -- they don't react

    U_new[URHO] = U_old[URHO] + dt_m * C[URHO];

    for (int n = 0; n < 3; ++n) {
        U_new[UMX+n] = U_old[UMX+n] + dt_m * C[UMX+n];
    }

    // Now only save the subset that participates in the
    // nonlinear solve -- note: we include the old state
    // in f_source

    burn_t burn_state;

    burn_state.success = true;

    // store the state

    copy_cons_to_burn_type(U_old, burn_state);

    // we need an initial T guess for the EOS
    burn_state.T = U_old[UTEMP];

    burn_state.T_fixed = -1.e30_rt;
    burn_state.rho = burn_state.y[SRHO];

    // If we are solving the system as an ODE, then we are
    // solving
    //   dU/dt = R(U) + C
    // so we simply pass in C

    // Note: we need to use the indexing order that is defined in the burn_t
    // so we are compatible with the integrator logic

    burn_state.ydot_a[SRHO] = C[URHO];
    burn_state.ydot_a[SMX] = C[UMX];
    burn_state.ydot_a[SMY] = C[UMY];
    burn_state.ydot_a[SMZ] = C[UMZ];
    burn_state.ydot_a[SEDEN] = C[UEDEN];
    burn_state.ydot_a[SEINT] = C[UEINT];
    for (int n = 0; n < NumSpec; n++) {
        burn_state.ydot_a[SFS+n] = C[UFS+n];
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; n++) {
        burn_state.ydot_a[SFX+n] = C[UFX+n];
    }
#endif

    burn_state.sdc_iter = sdc_iteration;
    burn_state.num_sdc_iters = sdc_order + sdc_extra;

    burner(burn_state, dt_m);

    if (! burn_state.success) {
        amrex::Error("Error: integration failed");
    }

    // update the full U_new -- we already did density and momentum

    U_new[UEDEN] = burn_state.y[SEDEN];
    U_new[UEINT] = burn_state.y[SEINT];
    for (int n = 0; n < NumSpec; n++) {
        U_new[UFS+n] = burn_state.y[SFS+n];
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; n++) {
        U_new[UFX+n] = burn_state.y[SFX+n];
    }
#endif

    // keep our temperature guess
    U_new[UTEMP] = U_old[UTEMP];

}

#endif

