#ifndef SDC_POWELL_SOLVE_H
#define SDC_POWEKK_SOLVE_H

#include <sdc_react_util.H>
#include <hybrj.H>
#include <hybrj_type.H>

namespace powell {

#ifdef REACTIONS

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void
    f_sdc(amrex::Array1D<amrex::Real, 1, NumSpec+1>& x,  // solution vector
          Array1D<Real, 1, NumSpec+1>& f,  // constraint equations to zero
          const burn_t& burn_state_in,
          int& iflag) {

        // This is used with the hybrid Powell solver to compute f
        // evaluated with x
        //
        // x is {rho X_1, rho X_2, ... rho X_N, rho e}

        GpuArray<Real, NUM_STATE> R_full;

        // we'll use the passed-in burn_t as our working state.  We need to
        // sync it with the input state x

        auto burn_state = burn_state_in;

        // compute the temperature -- this can be removed shortly, since
        // we already do this in single_zone_react_source

        for (int n = 0; n < NumSpec; ++n) {
            burn_state.xn[n] = amrex::Clamp(x(n+1) / burn_state.rho,
                                            network_rp::small_x, 1.0_rt);
        }
#if NAUX_NET > 0
        amrex::Error("error: aux data not currently supported in true SDC");
#endif

        burn_state.e = x(NumSpec+1) / burn_state.rho;

        eos(eos_input_re, burn_state);

        // get the reactive term for this state

        single_zone_react_source(burn_state, R_full);

        // where f is Eq. 36 evaluated with the current guess for U
        // we use ydot_a here to be U^n + dt_m C
        // Note: this is -f, it doesn't matter, since we are zeroing
        // it.

        auto dt_m = burn_state.time;

        for (int n = 1; n <= NumSpec; ++n) {
            f(n) = -burn_state.xn[n+1] +
                dt_m * R_full[UFS-1+n] + burn_state.ydot_a[SFS-1+n];
        }
        f(NumSpec+1) = -burn_state.rho * burn_state.e +
            dt_m * R_full[UEINT] + burn_state.ydot_a[SEINT];

    }


    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void
    jac_sdc(amrex::Array1D<amrex::Real, 1, NumSpec+1>& x,  // solution vector
            JacNetArray2D& Jac,  // constraint equations to zero
            const burn_t& burn_state_in,
            int& iflag) {

        // This is used with the hybrid Powell solver to compute J
        // evaluated with x
        //
        // x is {rho X_1, rho X_2, ... rho X_N, rho e}

        GpuArray<Real, NUM_STATE> R_full;

        // we'll use the passed-in burn_t as our working state.  We need to
        // sync it with the input state x

        auto burn_state = burn_state_in;

        // compute the temperature -- this can be removed shortly, since
        // we already do this in single_zone_react_source

        for (int n = 0; n < NumSpec; ++n) {
            burn_state.xn[n] = amrex::Clamp(x(n+1) / burn_state.rho,
                                            network_rp::small_x, 1.0_rt);
        }
#if NAUX_NET > 0
        amrex::Error("error: aux data not currently supported in true SDC");
#endif

        burn_state.e = x(NumSpec+1) / burn_state.rho;

        eos(eos_input_re, burn_state);

        // Instead of the decomposition into dw/dU and dR/dw
        // written out in the original paper Appendix A, we instead use the
        // form from the simplified-SDC paper (Zingale et al. 2022).  Note:
        // we are not including density anymore.

        auto dt_m = burn_state.time;

        single_zone_jac(burn_state, dt_m, Jac);
    }


    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void
    sdc_powell_solve(const Real dt_m,
                     GpuArray<Real, NUM_STATE> const& U_old,
                     GpuArray<Real, NUM_STATE> & U_new,
                     GpuArray<Real, NUM_STATE> const& C,
                     const int sdc_iteration) {

        // the purpose of this function is to solve the system
        //  U - dt R(U) = U_old + dt C
        // using the hybrid-Powell solver.  This is Eq. 36 in the paper.
        //
        // here, U_new should come in as a guess for the new U for
        // iterations > 0, it will be the solution from the previous
        // iteration initially.
        //
        // upon exit, U_new will be returned with the value that satisfied
        // the nonlinear function

        amrex::ignore_unused(sdc_iteration);

        // we will do the implicit update of only the terms that
        // have reactive sources
        //
        // 1:NumSpec : species
        // NumSpec+1 : (rho e)


        // update the density and momenta for this zone -- they don't react

        U_new[URHO] = U_old[URHO] + dt_m * C[URHO];

        for (int n = 0; n < 3; ++n) {
            U_new[UMX+n] = U_old[UMX+n] + dt_m * C[UMX+n];
        }

        burn_t burn_state;

        copy_cons_to_burn_type(U_new, burn_state);
        burn_state.rho = U_new[URHO];

        // we are solving
        //   f(U) = U - dt R(U) - U_old - dt C = 0
        // we define f_source = U_old + dt C so we are solving
        //   f(U) = U - dt R(U) - f_source = 0
        //
        // we'll store this in the burn_t ydot_a[]

        for (int n = 0; n < NumSpec; ++n) {
            burn_state.ydot_a[SFS + n] = U_old[UFS + n] + dt_m * C[UFS + n];
        }
        burn_state.ydot_a[SEINT] = U_old[UEINT] + dt_m * C[UEINT];

        // we need dt in the RHS, pass this through using the time field
        burn_state.time = dt_m;

        // The Powell solver stores the solution vector in its own type
        hybrj_t<NumSpec+1, JacNetArray2D> hj;
        hj.xtol = 1.e-10;
        hj.mode = 1;

        // locate the initial guess
        for (int n = 0; n < NumSpec; ++n) {
            hj.x(n+1) = U_new[UFS+n];
        }
        hj.x(NumSpec+1) = U_new[UEINT];

        // call the powell solver

        hybrj<NumSpec+1, burn_t>(hj, burn_state, f_sdc, jac_sdc);

        if (hj.info != 1) {
            amrex::Error("hybrid Powell solve failed");
        }

        // update the full U_new

        for (int n = 0; n < NumSpec; ++n) {
            U_new[UFS+n] = hj.x(n+1);
        }

        U_new[UEINT] = hj.x(NumSpec+1);

        // we want to do a conservative update for (rho E), so first
        // figure out the energy generation rate

        Real rho_Sdot = (U_new[UEINT] - U_old[UEINT]) / dt_m - C[UEINT];

        U_new[UEDEN] = U_old[UEDEN] + dt_m * (C[UEDEN] + rho_Sdot);

    }

#endif
}

#endif
