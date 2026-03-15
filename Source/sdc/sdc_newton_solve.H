#ifndef SDC_NEWTON_SOLVE_H
#define SDC_NEWTON_SOLVE_H

#include <sdc_react_util.H>
#include <linpack.H>

// error codes
namespace newton {
    constexpr int NEWTON_SUCCESS = 0;
    constexpr int SINGULAR_MATRIX = -1;
    constexpr int CONVERGENCE_FAILURE = -2;
    constexpr int BAD_MASS_FRACTIONS = -3;

    constexpr amrex::Real species_failure_tolerance = 1.e-2_rt;

    constexpr int MAX_ITER = 200;
    constexpr int MAX_NSUB = 64;

};



#ifdef REACTIONS

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
f_sdc_jac(const Real dt_m,
          burn_t& burn_state,
          Array1D<Real, 1, NumSpec+1>& f,
          JacNetArray2D& Jac) {

    // This is used with the Newton solve and returns f and the
    // Jacobian

    // single_zone_react_source will use burn_state.y[] to fill
    // burn_state.rho, burn_state.e, and burn_state.xn[] and then call
    // the EOS to get the temperature.

    GpuArray<Real, NUM_STATE> R_full;

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

    // Instead of the decomposition into dw/dU and dR/dw written out
    // in the original paper Appendix A, we instead use the form from
    // the simplified-SDC paper (Zingale et al. 2022).  Note: we are
    // not including density anymore.

    // note: single_zone_jac assumes that the burn_state is in sync
    // already with the current solution guess.  In this case, that
    // was done in single_zone_react_source.

    single_zone_jac(burn_state, dt_m, Jac);

    // Our Jacobian has the form:  J = I - dt dR/dw dwdU
    // (Eq. 38), so now we fix that

    Jac.mul(-dt_m);
    Jac.add_identity();

}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
sdc_newton_solve(const Real dt_m,
                 GpuArray<Real, NUM_STATE> const& U_old,
                 GpuArray<Real, NUM_STATE> & U_new,
                 GpuArray<Real, NUM_STATE> const& C,
                 const int sdc_iteration,
                 Real& err_out,
                 int& ierr) {

    // the purpose of this function is to solve the system
    //  U - dt R(U) = U_old + dt C
    // using a Newton solve.  This is Eq. 36 in the paper.
    //
    // here, U_new should come in as a guess for the new U for
    // iterations > 0, it will be the solution from the previous
    // iteration initially.
    //
    // upon exit, U_new will be returned with the value that satisfied
    // the nonlinear function

    amrex::ignore_unused(sdc_iteration);

    JacNetArray2D Jac;

    // we will do the implicit update of only the terms that
    // have reactive sources
    //
    // 1:NumSpec : species
    // NumSpec+1 : (rho e)

    Array1D<Real, 1, NumSpec+1> f;

    ierr = newton::NEWTON_SUCCESS;

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

    // do a simple Newton solve

    // iterative loop
    int iter = 0;

    Real err = 1.e30_rt;
    bool converged = false;

    while (!converged && iter < newton::MAX_ITER) {

        // burn_state.y[] will always contain the current guess for the solution
        // only the species and internal energy are updated though

        // initial guess
        burn_state.T = U_old[UTEMP];

        int info = 0;
        f_sdc_jac(dt_m, burn_state, f, Jac);

        // solve the linear system: Jac dU = -f
#ifdef NEW_NETWORK_IMPLEMENTATION
        info = RHS::dgefa(Jac);
#else
        IArray1D ipvt;
        constexpr bool allow_pivot{true};
        dgefa<NumSpec+1, allow_pivot>(Jac, ipvt, info);
#endif
        if (info != 0) {
            ierr = newton::SINGULAR_MATRIX;
            return;
        }

#ifdef NEW_NETWORK_IMPLEMENTATION
        RHS::dgesl(Jac, f);
#else
        dgesl<NumSpec+1, allow_pivot>(Jac, ipvt, f);
#endif

        // on output, f is the correction to the solution (dU)

        // we want to limit the change by some amount k, such that
        //
        // 1/k < (U + dU) / U < k
        //
        // we'll write this as
        //
        // 1/k < 1 + φ < k
        //
        // where φ = dU/U
        //
        // we want to limit this, so we will introduce a term α, such that
        //
        // 1/k < 1 + α φ < k
        //
        // if φ > 0, then we need α < (k - 1) / φ
        //
        // if φ < 0, then we need α < |1/k - 1| / |φ|

        amrex::Real alpha{1.0_rt};
        constexpr amrex::Real k = 10.0;
        constexpr amrex::Real X_cutoff = 1.e-4_rt;

        for (int n = 0; n < NumSpec; ++n) {
            const amrex::Real rhoX = burn_state.y[SFS+n];
            if (rhoX < burn_state.rho * X_cutoff) {
                continue;
            }
            const amrex::Real drhoX = f(n + 1);
            const amrex::Real phi = drhoX / rhoX;

            if (phi > 0.0_rt) {
                alpha = std::min(alpha, (k - 1.0_rt) / phi);
            } else if (phi < 0.0_rt) {
                alpha = std::min(alpha, (1.0_rt / k - 1.0_rt) / phi);
            }
        }

        // make sure we apply atleast some correction
        alpha = std::max(alpha, 5.e-2_rt);

        for (int n = 0; n < NumSpec; ++n) {
            burn_state.y[SFS+n] += alpha * f(n + 1);

        }
        burn_state.y[SEINT] += alpha * f(NumSpec+1);

        // compute the norm of the weighted error, where the
        // weights are 1/eps_tot

        auto err_sum = 0.0_rt;
        for (int n = 1; n <= NumSpec+1; ++n) {
            Real eps;
            if (n < NumSpec+1) {
                // for species, atol is the mass fraction limit, so we
                // multiply by density to get a partial density limit
                eps = integrator_rp::rtol_spec * std::abs(burn_state.y[SFS-1+n]) +
                      integrator_rp::atol_spec * std::abs(U_new[URHO]);
            } else {
                eps = integrator_rp::rtol_enuc * std::abs(burn_state.y[SEINT]) +
                      integrator_rp::atol_enuc;
            }
            err_sum += (alpha * f(n)) * (alpha * f(n)) / (eps * eps);
        }
        err = std::sqrt(err_sum / static_cast<Real>(NumSpec+1));

        if (err < 1.0_rt) {
            converged = true;
        }
        iter++;
    }

    err_out = err;

    if (!converged) {
        ierr = newton::CONVERGENCE_FAILURE;
        return;
    }

    // update the full U_new

    for (int n = 0; n < NumSpec; ++n) {
        U_new[UFS+n] = burn_state.y[SFS+n];
    }

    U_new[UEINT] = burn_state.y[SEINT];

    // we want to do a conservative update for (rho E), so first figure out the
    // energy generation rate

    Real rho_Sdot = (U_new[UEINT] - U_old[UEINT]) / dt_m - C[UEINT];

    U_new[UEDEN] = U_old[UEDEN] + dt_m * (C[UEDEN] + rho_Sdot);

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
sdc_newton_subdivide(const Real dt_m,
                     GpuArray<Real, NUM_STATE> const& U_old,
                     GpuArray<Real, NUM_STATE>& U_new,
                     GpuArray<Real, NUM_STATE> const& C,
                     const int sdc_iteration,
                     Real& err_out,
                     int& ierr) {

    // This is the driver for solving the nonlinear update for the
    // reating/advecting system using Newton's method. It attempts to
    // do the solution for the full dt_m requested, but if it fails,
    // will subdivide the domain until it converges or reaches our
    // limit on the number of subintervals.

    GpuArray<Real, NUM_STATE> U_begin;

    // subdivide the timestep and do multiple Newtons.

    int nsub = 1;
    ierr = newton::CONVERGENCE_FAILURE;

    while (nsub <= newton::MAX_NSUB && ierr != newton::NEWTON_SUCCESS) {

        // For the first interval, we set the starting state to
        // the old-time solution.

        for (int n = 0; n < NUM_STATE; ++n) {
            U_begin[n] = U_old[n];
        }

        // We come in here with an initial guess for the new solution
        // stored in U_new. That only really makes sense for the case
        // where we have 1 substep. Otherwise, we should just use the
        // previous subinterval solution.

        if (nsub > 1) {
            for (int n = 0; n < NUM_STATE; ++n) {
                U_new[n] = U_old[n];
            }
        }
        Real dt_sub = dt_m / nsub;

        for (int isub = 0; isub < nsub; ++isub) {

            // integrate from U_begin to U_new over a step dt_m / nsub

            // normalize species -- we are working with rho X_k here
            Real sum_rhoX = 0.0_rt;
            for (int n = 0; n < NumSpec; ++n) {
                U_begin[UFS + n] = amrex::max(U_begin[URHO] * network_rp::small_x, U_begin[UFS + n]);
                sum_rhoX += U_begin[UFS + n];
            }
            for (int n = 0; n < NumSpec; ++n) {
                U_begin[UFS + n] *= U_begin[URHO] / sum_rhoX;
            }

            sdc_newton_solve(dt_sub, U_begin, U_new, C, sdc_iteration, err_out, ierr);

            // our solve may have resulted in mass fractions outside
            // of [0, 1] -- reject if this is the case
            for (int n = 0; n < NumSpec; ++n) {
                if (U_new[UFS+n] < -newton::species_failure_tolerance * U_new[URHO] ||
                    U_new[UFS+n] > (1.0_rt + newton::species_failure_tolerance) * U_new[URHO]) {
                    ierr = newton::BAD_MASS_FRACTIONS;
                }
            }

            if (ierr != newton::NEWTON_SUCCESS) {
                // no point in continuing this subdivision if one stage failed
                break;
            }

            for (int n = 0; n < NUM_STATE; ++n) {
                U_begin[n] = U_new[n];
            }
        }
        nsub *= 2;
    }
}
#endif

#endif
