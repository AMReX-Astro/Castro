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

    constexpr Real species_failure_tolerance = 1.e-2_rt;
};

#ifdef REACTIONS

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
f_sdc_jac(const Real dt_m,
          Array1D<Real, 1, NumSpec+1> const& U,
          Array1D<Real, 1, NumSpec+1>& f,
          JacNetArray2D& Jac,
          Array1D<Real, 1, NumSpec+1>& f_source,
          const Real rho,
          GpuArray<Real, 3>& mom_source,
          const Real T_old,
          const Real E_var) {

    // This is used with the Newton solve and returns f and the Jacobian

    GpuArray<Real, NUM_STATE> U_full;
    GpuArray<Real, NUM_STATE> R_full;

    // we are not solving the momentum equations
    // create a full state -- we need this for some interfaces

    U_full[URHO] = rho;

    for (int n = 0; n < NumSpec; ++n) {
        U_full[UFS+n] = U(n+1);
    }

    U_full[UEINT] = U(NumSpec+1);
    U_full[UEDEN] = E_var;

    for (int n = 0; n < 3; ++n) {
        U_full[UMX+n] = mom_source[n];
    }

    // normalize the species
    auto sum_rhoX = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        U_full[UFS+n] = amrex::max(network_rp::small_x, U_full[UFS+n]);
        sum_rhoX += U_full[UFS+n];
    }

    for (int n = 0; n < NumSpec; ++n) {
        U_full[UFS+n] *= U_full[URHO] / sum_rhoX;
    }

    // compute the temperature and species derivatives --
    // maybe this should be done using the burn_state
    // returned by single_zone_react_source, since it is
    // more consistent T from e

    eos_extra_t eos_state;
    eos_state.rho = U_full[URHO];
    eos_state.T = T_old;   // initial guess
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = U_full[UFS+n] / U_full[URHO];
    }
#if NAUX_NET > 0
    amrex::Error("error: aux data not currently supported in true SDC");
#endif

    eos_state.e = U_full[UEINT] / U_full[URHO];

    eos(eos_input_re, eos_state);

    U_full[UTEMP] = eos_state.T;

    // we'll create a burn_state to pass stuff from the RHS to the Jac function

    burn_t burn_state;

    single_zone_react_source(U_full, R_full, burn_state);

    // we are solving J dU = -f
    // where f is Eq. 36 evaluated with the current guess for U

    for (int n = 1; n <= NumSpec; ++n) {
        f(n) = U(n) - dt_m * R_full[UFS-1+n] - f_source(n);
    }
    f(NumSpec+1) = U(NumSpec+1) - dt_m * R_full[UEINT] - f_source(NumSpec+1);

    // get the Jacobian.

    // Instead of the decomposition into dw/dU and dR/dw
    // written out in the original paper Appendix A, we instead use the
    // form from the simplified-SDC paper (Zingale et al. 2022).  Note:
    // we are not including density anymore.

    single_zone_jac(U_full, burn_state, dt_m, Jac);

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
sdc_newton_solve(const Real dt_m,
                 GpuArray<Real, NUM_STATE> const& U_old,
                 GpuArray<Real, NUM_STATE> & U_new,
                 GpuArray<Real, NUM_STATE> const& C,
                 const int sdc_iteration,
                 Real& err_out,
                 int& ierr) {
    // the purpose of this function is to solve the system
    // U - dt R(U) = U_old + dt C using a Newton solve.
    // This is Eq. 36 in the paper.
    //
    // here, U_new should come in as a guess for the new U
    // and will be returned with the value that satisfied the
    // nonlinear function

    JacNetArray2D Jac;

    // we will do the implicit update of only the terms that
    // have reactive sources
    //
    // 1:NumSpec : species
    // NumSpec+1 : (rho e)

    Array1D<Real, 1, NumSpec+1> U_react;
    Array1D<Real, 1, NumSpec+1> f_source;
    GpuArray<Real, 3> mom_source;
    Array1D<Real, 1, NumSpec+1> dU_react;
    Array1D<Real, 1, NumSpec+1> f;
    RArray1D f_rhs;

    const int MAX_ITER = 100;

    ierr = newton::NEWTON_SUCCESS;

    // update the density and momenta for this zone -- they don't react

    U_new[URHO] = U_old[URHO] + dt_m * C[URHO];

    for (int n = 0; n < 3; ++n) {
        U_new[UMX+n] = U_old[UMX+n] + dt_m * C[UMX+n];
    }

    // now only save the subset that participates in the
    // nonlinear solve -- note: we include the old state in
    // f_source

    // for the Jacobian solve, we are solving
    //   f(U) = U - dt R(U) - U_old - dt C = 0
    // we define f_source = U_old + dt C so we are solving
    //   f(U) = U - dt R(U) - f_source = 0

    for (int n = 0; n < NumSpec; ++n) {
        f_source(1 + n) = U_old[UFS + n] + dt_m * C[UFS + n];
    }
    f_source(NumSpec+1) = U_old[UEINT] + dt_m * C[UEINT];

    // set the momenta to be U_new
    for (int n = 0; n < 3; ++n) {
        mom_source[n] = U_new[UMX+n];
    }

    // temperature will be used as an initial guess in the EOS

    Real T_old = U_old[UTEMP];

    // we should be able to do an update for this somehow?

    Real E_var = U_new[UEDEN];

    // we need the density too

    Real rho_new = U_new[URHO];

    // store the subset for the nonlinear solve
    // We use an initial guess if possible

    for (int n = 0; n < NumSpec; ++n) {
        U_react(1+n) = U_new[UFS+n];
    }
    U_react(NumSpec+1) = U_new[UEINT];

#if (INTEGRATOR == 0)

    // do a simple Newton solve

    // iterative loop
    int iter = 0;

    Real err = 1.e30_rt;
    bool converged = false;

    while (!converged && iter < MAX_ITER) {
        int info = 0;
        f_sdc_jac(dt_m, U_react, f, Jac, f_source, rho_new, mom_source, T_old, E_var);

        // solve the linear system: Jac dU_react = -f
#ifdef NEW_NETWORK_IMPLEMENTATION
        RHS::dgefa(Jac);
        info = 0;
#else
        IArray1D ipvt;
        dgefa<NumSpec+1>(Jac, ipvt, info);
#endif
        if (info != 0) {
            ierr = newton::SINGULAR_MATRIX;
            return;
        }

        for (int n = 1; n <= NumSpec+1; ++n) {
            f_rhs(n) = -f(n);
        }

#ifdef NEW_NETWORK_IMPLEMENTATION
        RHS::dgesl(Jac, f_rhs);
#else
        dgesl<NumSpec+1>(Jac, ipvt, f_rhs);
#endif

        for (int n = 1; n <= NumSpec+1; ++n) {
            dU_react(n) = f_rhs(n);
        }

        // how much of dU_react should we apply?
        Real eta = 1.0_rt;
        for (int n = 1; n <= NumSpec+1; ++n) {
            U_react(n) += eta * dU_react(n);
        }

        Array1D<Real, 1, NumSpec+1> eps_tot;

        // for species, atol is the mass fraction limit, so we
        // multiply by density to get a partial density limit

        for (int n = 0; n < NumSpec; ++n) {
            eps_tot(1 + n) = integrator_rp::rtol_spec * std::abs(U_react(1 + n)) +
                integrator_rp::atol_spec * std::abs(U_new[URHO]);
        }
        eps_tot(NumSpec+1) = integrator_rp::rtol_enuc * std::abs(U_react(NumSpec+1)) +
            integrator_rp::atol_enuc;

        // compute the norm of the weighted error, where the
        // weights are 1/eps_tot

        auto err_sum = 0.0_rt;
        for (int n = 1; n <= NumSpec+1; ++n) {
            err_sum += dU_react(n) * dU_react(n) / (eps_tot(n) * eps_tot(n));
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

#endif

    // update the full U_new

    for (int n = 0; n < NumSpec; ++n) {
        U_new[UFS+n] = U_react(1+n);
    }
    auto v2 = 0.0_rt;
    for (int m = 0; m < 3; ++m) {
        v2 += U_new[UMX+m] * U_new[UMX+m];
    }

    U_new[UEINT] = U_react(NumSpec+1);

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
    // This is the driver for solving the nonlinear update for
    // the reating/advecting system using Newton's method. It
    // attempts to do the solution for the full dt_m requested,
    // but if it fails, will subdivide the domain until it
    // converges or reaches our limit on the number of
    // subintervals.

    const int MAX_NSUB = 64;
    GpuArray<Real, NUM_STATE> U_begin;

    // subdivide the timestep and do multiple Newtons. We come
    // in here with an initial guess for the new solution
    // stored in U_new. That only really makes sense for the
    // case where we have 1 substep. Otherwise, we should just
    // use the old time solution.

    int nsub = 1;
    ierr = newton::CONVERGENCE_FAILURE;

    for (int n = 0; n < NUM_STATE; ++n) {
        U_begin[n] = U_old[n];
    }

    while (nsub < MAX_NSUB && ierr != newton::NEWTON_SUCCESS) {
        if (nsub > 1) {
            for (int n = 0; n < NUM_STATE; ++n) {
                U_new[n] = U_old[n];
            }
        }
        Real dt_sub = dt_m / nsub;

        for (int isub = 0; isub < nsub; ++isub) {
            // normalize species
            Real sum_rhoX = 0.0_rt;
            for (int n = 0; n < NumSpec; ++n) {
                U_begin[UFS + n] = amrex::max(network_rp::small_x, U_begin[UFS + n]);
                sum_rhoX += U_begin[UFS + n];
            }
            for (int n = 0; n < NumSpec; ++n) {
                U_begin[UFS + n] *= U_begin[URHO] / sum_rhoX;
            }

            sdc_newton_solve(dt_sub, U_begin, U_new, C, sdc_iteration, err_out, ierr);

            if (ierr != newton::NEWTON_SUCCESS) {
                // no point in continuing this subdivision if one stage failed
                break;
            }

            // our solve may have resulted in mass fractions outside
            // of [0, 1] -- reject if this is the case
            for (int n = 0; n < NumSpec; ++n) {
                if (U_new[UFS+n] < -newton::species_failure_tolerance * U_new[URHO] ||
                    U_new[UFS+n] > (1.0_rt + newton::species_failure_tolerance) * U_new[URHO]) {
                    ierr = newton::BAD_MASS_FRACTIONS;
                }
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
