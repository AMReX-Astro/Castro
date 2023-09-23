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
#include <sdc_react_util.H>

#include <vode_dvode.H>

// The f_rhs routine provides the right-hand-side for the DVODE solver.
// This is a generic interface that calls the specific RHS routine in the
// network you're actually using.

template <typename DvodeT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void rhs(const Real time, burn_t& burn_state, DvodeT& vode_state, RArray1D& dUdt, const bool in_jacobian=false)
{

    // note: dUdt is 1-based

    GpuArray<Real, NUM_STATE> U_full;
    GpuArray<Real, NUM_STATE> R_full;

    // evaluate R

    // update the density -- we'll need this for the EOS call in the
    // react update

    U_full[URHO] = burn_state.rho_orig + time * burn_state.ydot_a[SRHO];

    // we are not solving the momentum equations
    // and never need total energy, so we don't worry about filling those

    for (int n = 0; n < NumSpec; n++) {
        U_full[UFS+n] = vode_state.y(1+n);
    }
    U_full[UEINT] = vode_state.y(NumSpec+1);

    // initialize the temperature -- a better value will be found when we do the EOS
    // call in single_zone_react_source

    U_full[UTEMP] = burn_state.T;

    // create a temporary burn_t for this call.  It is just going to
    // get loaded up with U_full

    burn_t burn_state_pass;
    single_zone_react_source(U_full, R_full, burn_state_pass);

    // update our temperature for next time

    burn_state.T = burn_state_pass.T;

    // now construct the RHS -- note this is 1-based

    for (int n = 0; n < NumSpec; ++n) {
        dUdt(n+1) = R_full[UFS+n] + burn_state.ydot_a[SFS+n];
    }
    dUdt(NumSpec+1) = R_full[UEINT] + burn_state.ydot_a[SEINT];

}


template<class MatrixType, typename DvodeT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void jac (const Real time, burn_t& burn_state, DvodeT& vode_state, MatrixType& pd)
{
    // Call the specific network routine to get the Jacobian.

    // note the Jacobian is 1-based at the moment

    GpuArray<Real, NUM_STATE> U_full;
    GpuArray<Real, NUM_STATE> R_full;

    // we are not solving the momentum equations
    // create a full state -- we need this for some interfaces

    U_full[URHO] = burn_state.rho_orig + time * burn_state.ydot_a[URHO];

    for (int n = 0; n < NumSpec; n++) {
        U_full[UFS+n] = vode_state.y(1+n);
    }
    U_full[UEINT] = vode_state.y(NumSpec+1);

    U_full[UTEMP] = burn_state.T;

    burn_t burn_state_pass;

    // this will initialize burn_state_pass and fill the temperature
    // via the EOS

    single_zone_react_source(U_full, R_full, burn_state_pass);

    single_zone_jac(U_full, burn_state_pass, vode_state.H, pd);

}


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
    // 0 : rho
    // 1:NumSpec : species
    // NumSpec+1 : (rho E) or (rho e)

#if (INTEGRATOR == 0)

    // Update the momenta for this zone -- they don't react

    for (int n = 0; n < 3; ++n) {
        U_new[UMX+n] = U_old[UMX+n] + dt_m * C[UMX+n];
    }

    // Now only save the subset that participates in the
    // nonlinear solve -- note: we include the old state
    // in f_source

    // load rpar


    dvode_t<NumSpec+1> dvode_state;

    burn_t burn_state;

    burn_state.rho_orig = U_old[URHO];
    burn_state.rhoe_orig = U_old[UEINT];

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

    burn_state.success = true;

    // temperature will be used as an initial guess in the EOS
    burn_state.T = U_old[UTEMP];


    // store the subset for the nonlinear solve. We only
    // consider (rho e), not (rho E). This is because at
    // present we do not have a method of updating the
    // velocities during the multistep integration.

    // Also note that the dvode_state is 1-based, but we'll
    // access it as 0-based in our implementation of the
    // RHS routine

    for (int n = 0; n < NumSpec; ++n) {
        dvode_state.y(1+n) = U_old[UFS+n];
    }
    dvode_state.y(NumSpec+1) = U_old[UEINT];

    // // set the maximum number of steps allowed
    // dvode_state.MXSTEP = 25000;

    dvode_state.t = 0.0_rt;
    dvode_state.tout = dt_m;
    dvode_state.HMXI = 1.0_rt / ode_max_dt;

    dvode_state.jacobian_type = integrator_rp::jacobian;

    // relative tolerances
    dvode_state.rtol_spec = integrator_rp::rtol_spec;
    dvode_state.rtol_enuc = integrator_rp::rtol_enuc;

    // absolute tolerances
    dvode_state.atol_spec = integrator_rp::atol_spec * U_old[URHO]; // this way, atol is the minimum x
    dvode_state.atol_enuc = integrator_rp::atol_enuc * U_old[UEINT];

    int istate = dvode(burn_state, dvode_state);

    if (istate < 0) {
        Abort("Vode terminated poorly");
    }

    // update the full U_new
    for (int n = 0; n < NumSpec; ++n) {
        U_new[UFS+n] = dvode_state.y(1+n);
    }
    U_new[UEINT] = dvode_state.y(NumSpec+1);

    auto v2 = 0.0_rt;
    for (int n = 0; n < 3; ++n) {
        v2 += U_new[UMX+n] * U_new[UMX+n];
    }
    U_new[UEDEN] = U_new[UEINT] + 0.5_rt * v2 / U_new[URHO];

    // keep our temperature guess
    U_new[UTEMP] = U_old[UTEMP];

#endif
}

#endif

