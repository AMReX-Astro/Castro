#ifndef SDC_POWELL_SOLVE_H
#define SDC_POWEKK_SOLVE_H

#include <sdc_react_util.H>
#include <hybrj.H>


namespace powell {

#ifdef REACTIONS

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void
    f_sdc_jac(const Real dt_m,
              burn_t& burn_state,
              Array1D<Real, 1, NumSpec+1>& f,
              JacNetArray2D& Jac) {

        // This is used with the Newton solve and returns f and the Jacobian

        GpuArray<Real, NUM_STATE> R_full;

        // create a burn_t -- this will hold the full state, reconstructed
        // from the current solution U


        // compute the temperature -- this can be removed shortly, since
        // we already do this in single_zone_react_source

        for (int n = 0; n < NumSpec; ++n) {
            burn_state.xn[n] = amrex::Clamp(burn_state.y[SFS+n] / burn_state.y[SRHO], network_rp::small_x, 1.0_rt);
        }
#if NAUX_NET > 0
        amrex::Error("error: aux data not currently supported in true SDC");
#endif

        burn_state.e = burn_state.y[SEINT] / burn_state.y[SRHO];

        eos(eos_input_re, burn_state);

        single_zone_react_source(burn_state, R_full);

        // we are solving J dU = -f
        // where f is Eq. 36 evaluated with the current guess for U
        // note: we store -f

        for (int n = 1; n <= NumSpec; ++n) {
            f(n) = -burn_state.y[SFS-1+n] +
                dt_m * R_full[UFS-1+n] + burn_state.ydot_a[SFS-1+n];
        }
        f(NumSpec+1) = -burn_state.y[SEINT] +
            dt_m * R_full[UEINT] + burn_state.ydot_a[SEINT];

        // get the Jacobian.

        // Instead of the decomposition into dw/dU and dR/dw
        // written out in the original paper Appendix A, we instead use the
        // form from the simplified-SDC paper (Zingale et al. 2022).  Note:
        // we are not including density anymore.

        single_zone_jac(burn_state, dt_m, Jac);

        // Our Jacobian has the form:  J = I - dt dR/dw dwdU
        // (Eq. 38), so now we fix that

        for (int n = 1; n <= NumSpec+1; ++n) {
            for (int m = 1; m <= NumSpec+1; ++m) {
                Real coeff = (m == n) ? 1.0_rt : 0.0_rt;
                Jac(n, m) = coeff - dt_m * Jac(n, m);
            }
        }

    }

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void
    sdc_powell_solve(const Real dt_m,
                     GpuArray<Real, NUM_STATE> const& U_old,
                     GpuArray<Real, NUM_STATE> & U_new,
                     GpuArray<Real, NUM_STATE> const& C,
                     const int sdc_iteration,
                     Real& err_out,
                     int& ierr) {

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

        // for the Jacobian solve, we are solving
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

        // call the powell solver

        TODO HERE

        // update the full U_new

        for (int n = 0; n < NumSpec; ++n) {
            U_new[UFS+n] = burn_state.y[SFS+n];
        }

        U_new[UEINT] = burn_state.y[SEINT];

        // we want to do a conservative update for (rho E), so first
        // figure out the energy generation rate

        Real rho_Sdot = (U_new[UEINT] - U_old[UEINT]) / dt_m - C[UEINT];

        U_new[UEDEN] = U_old[UEDEN] + dt_m * (C[UEDEN] + rho_Sdot);

    }

#endif
}

#endif
