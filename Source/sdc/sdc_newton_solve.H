#ifndef SDC_NEWTON_SOLVE_H
#define SDC_NEWTON_SOLVE_H

#include <sdc_react_util.H>

// error codes
constexpr int NEWTON_SUCCESS = 0;
constexpr int SINGULAR_MATRIX = -1;
constexpr int CONVERGENCE_FAILURE = -2;


#ifdef REACTIONS

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
f_sdc_jac(const Real dt_m,
          GpuArray<Real, NumSpec+2> const& U,
          GpuArray<Real, NumSpec+2>& f,
          RArray2D& Jac,
          GpuArray<Real, NumSpec+2>& f_source,
          GpuArray<Real, 3>& mom_source,
          const Real T_old,
          const Real E_var) {

    // This is used with the Newton solve and returns f and the Jacobian

    GpuArray<Real, NUM_STATE> U_full;
    GpuArray<Real, NUM_STATE> R_full;
    GpuArray<Real, NumSpec+2> R_react;

    Array2D<Real, 0, NumSpec+1, 0, NumSpec+1> dRdw = {0.0_rt};
    Array2D<Real, 0, NumSpec+1, 0, NumSpec+1> dwdU = {0.0_rt};

    Array2D<Real, 0, NumSpec+1, 0, NumSpec+1> dRdU;

    // we are not solving the momentum equations
    // create a full state -- we need this for some interfaces
    U_full[URHO] = U[0];

    for (int n = 0; n < NumSpec; ++n) {
        U_full[UFS+n] = U[1+n];
    }

    if (sdc_solve_for_rhoe == 1) {
        U_full[UEINT] = U[NumSpec+1];
        U_full[UEDEN] = E_var;
    } else {
        U_full[UEDEN] = U[NumSpec+1];
        U_full[UEINT] = E_var;
    }

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
    for (int n = 0; n < NumAux; ++n) {
        eos_state.aux[n] = U_full[UFX+n] / U_full[URHO];
    }
#endif
    eos_state.e = U_full[UEINT] / U_full[URHO];  // (U_full[UEDEN] - 0.5_rt*sum(U_full(UMX:UMZ))/U_full[URHO])/U_full[URHO]

    eos(eos_input_re, eos_state);

    U_full[UTEMP] = eos_state.T;

    // we'll create a burn_state to pass stuff from the RHS to the Jac function
    burn_t burn_state;

    single_zone_react_source(U_full, R_full, burn_state);

    // store the subset of R used in the Jacobian
    R_react[0] = R_full[URHO];
    for (int n = 0; n < NumSpec; ++n) {
        R_react[n+1] = R_full[UFS+n];
    }
    if (sdc_solve_for_rhoe == 1) {
        R_react[NumSpec+1] = R_full[UEINT];
    } else {
        R_react[NumSpec+1] = R_full[UEDEN];
    }

    for (int n = 0; n < NumSpec+2; ++n) {
        f[n] = U[n] - dt_m * R_react[n] - f_source[n];
    }

    // get dRdw -- this may do a numerical approximation or use the
    // network's analytic Jac
    single_zone_jac(U_full, burn_state, dt_m, dRdw);

    // the density row
    dwdU(iwrho, 0) = 1.0_rt;

    // the X_k rows
    for (int m = 1; m < NumSpec+1; ++m) {
        dwdU(iwfs-1+m, 0) = -U[m] / (U[0]*U[0]);
        dwdU(iwfs-1+m, m) = 1.0_rt / U[0];
    }

    auto eos_xderivs = composition_derivatives(eos_state);

    // now the e row -- this depends on whether we are evolving (rho E) or (rho e)
    auto denom = 1.0_rt / eos_state.rho;
    auto xn_sum = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        xn_sum += eos_state.xn[n] * eos_xderivs.dedX[n];
    }

    if (sdc_solve_for_rhoe == 1) {
        dwdU(iwe, 0) = denom * (xn_sum - eos_state.e);
    } else {
        auto u2_sum = 0.0_rt;
        for (auto n = 0; n < 3; ++n) {
            u2_sum += U_full[UMX+n] * U_full[UMX+n];
        }
        dwdU(iwe, 0) = denom * (xn_sum - eos_state.e -
                               0.5_rt * u2_sum / (eos_state.rho * eos_state.rho));
    }

    for (int m = 0; m < NumSpec; ++m) {
        dwdU(iwe,m+1) = -denom * eos_xderivs.dedX[m];
    }

    dwdU(iwe, NumSpec+1) = denom;

    // construct the Jacobian -- we can get most of the
    // terms from the network itself, but we do not rely on
    // it having derivative wrt density

    // Note: Jac is 1-based, because that is what the linear algebra
    // routines expect but the components, dRdw and dwdU are 0-based!

    for (int n = 1; n <= NumSpec+2; ++n) {
        for (int m = 1; m <= NumSpec+2; ++m) {
            Jac(n, m) = 0.0_rt;
        }
    }

    for (int m = 1; m <= NumSpec+2; ++m) {
        Jac(m, m) = 1.0_rt;
    }

    // auto dRdU = matmul(dRdw, dwdU);
    // dRdU is 0-based

    for (int n = 0; n <= NumSpec+1; ++n) {
        for (int m = 0; m <= NumSpec+1; ++m) {
            dRdU(n,m) = 0.0_rt;
            for (int l = 0; l <= NumSpec+1; ++l) {
                dRdU(n,m) += dRdw(n,l) * dwdU(l,m);
            }
        }
    }

    for (int n = 0; n < NumSpec+2; ++n) {
        for (int m = 0; m < NumSpec+2; ++m) {
            Jac(n+1, m+1) -= dt_m * dRdU(n,m);
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
    //
    // here, U_new should come in as a guess for the new U
    // and will be returned with the value that satisfied the
    // nonlinear function

    RArray2D Jac;

    // we will do the implicit update of only the terms that
    // have reactive sources
    //
    // 0 : rho
    // 1:NumSpec : species
    // NumSpec+1 : (rho E) or (rho e)

    GpuArray<Real, NumSpec+2> U_react;
    GpuArray<Real, NumSpec+2> f_source;
    GpuArray<Real, 3> mom_source;
    GpuArray<Real, NumSpec+2> dU_react;
    GpuArray<Real, NumSpec+2> f;
    RArray1D f_rhs;

    const int MAX_ITER = 100;

    ierr = NEWTON_SUCCESS;

    // the tolerance we are solving to may depend on the
    // iteration
    Real relax_fac = std::pow(sdc_solver_relax_factor, sdc_order - sdc_iteration - 1);
    Real tol_dens = sdc_solver_tol_dens * relax_fac;
    Real tol_spec = sdc_solver_tol_spec * relax_fac;
    Real tol_ener = sdc_solver_tol_ener * relax_fac;

    // update the momenta for this zone -- they don't react
    for (int n = 0; n < 3; ++n) {
        U_new[UMX+n] = U_old[UMX+n] + dt_m * C[UMX+n];
    }

    // now only save the subset that participates in the
    // nonlinear solve -- note: we include the old state in
    // f_source

    // load rpar

    // for the Jacobian solve, we are solving
    //   f(U) = U - dt R(U) - U_old - dt C = 0
    // we define f_source = U_old + dt C so we are solving
    //   f(U) = U - dt R(U) - f_source = 0

    f_source[0] = U_old[URHO] + dt_m * C[URHO];
    for (int n = 0; n < NumSpec; ++n) {
        f_source[1 + n] = U_old[UFS + n] + dt_m * C[UFS + n];
    }
    if (sdc_solve_for_rhoe == 1) {
        f_source[NumSpec+1] = U_old[UEINT] + dt_m * C[UEINT];
    } else {
        f_source[NumSpec+1] = U_old[UEDEN] + dt_m * C[UEDEN];
    }

    // set the momenta to be U_new
    for (int n = 0; n < 3; ++n) {
        mom_source[n] = U_new[UMX+n];
    }

    // temperature will be used as an initial guess in the EOS

    Real T_old = U_old[UTEMP];

    // we should be able to do an update for this somehow?

    Real E_var;

    if (sdc_solve_for_rhoe == 1) {
        E_var = U_new[UEDEN];
    } else {
        E_var = U_new[UEINT];
    }

    // store the subset for the nonlinear solve
    // We use an initial guess if possible
    U_react[0] = U_new[URHO];
    for (int n = 0; n < NumSpec; ++n) {
        U_react[1+n] = U_new[UFS+n];
    }
    if (sdc_solve_for_rhoe == 1) {
        U_react[NumSpec+1] = U_new[UEINT];
    } else {
        U_react[NumSpec+1] = U_new[UEDEN];
    }

#if (INTEGRATOR == 0)

    // do a simple Newton solve

    // iterative loop
    int iter = 0;
    int max_newton_iter = MAX_ITER;

    Real err = 1.e30_rt;
    bool converged = false;

    while (!converged && iter < max_newton_iter) {
        int info = 0;
        f_sdc_jac(dt_m, U_react, f, Jac, f_source, mom_source, T_old, E_var);

        IArray1D ipvt;

        // solve the linear system: Jac dU_react = -f
#ifdef NEW_NETWORK_IMPLEMENTATION
        RHS::dgefa(Jac);
        info = 0;
#else
        dgefa<NumSpec+2>(Jac, ipvt, info);
#endif
        if (info != 0) {
            ierr = SINGULAR_MATRIX;
            return;
        }

        for (int n = 1; n <= NumSpec+2; ++n) {
            f_rhs(n) = -f[n-1];
        }

#ifdef NEW_NETWORK_IMPLEMENTATION
        RHS::dgesl(Jac, f_rhs);
#else
        dgesl<NumSpec+2>(Jac, ipvt, f_rhs);
#endif

        for (int n = 0; n < NumSpec+2; ++n) {
            dU_react[n] = f_rhs(n+1);
        }

        // how much of dU_react should we apply?
        Real eta = 1.0_rt;
        for (int n = 0; n < NumSpec+2; ++n) {
            dU_react[n] *= eta;
            U_react[n] += dU_react[n];
        }

        GpuArray<Real, NumSpec+2> eps_tot;

        eps_tot[0] = tol_dens * std::abs(U_react[0]) + sdc_solver_atol;

        // for species, atol is the mass fraction limit, so we
        // multiply by density to get a partial density limit
        for (int n = 0; n < NumSpec; ++n) {
            eps_tot[1 + n] = tol_spec * std::abs(U_react[1 + n]) + sdc_solver_atol * std::abs(U_react[0]);
        }
        eps_tot[NumSpec+1] = tol_ener * std::abs(U_react[NumSpec+1]) + sdc_solver_atol;

        // compute the norm of the weighted error, where the
        // weights are 1/eps_tot
        auto err_sum = 0.0_rt;
        for (int n = 0; n < NumSpec+2; ++n) {
            err_sum += dU_react[n] * dU_react[n] / (eps_tot[n]* eps_tot[n]);
        }
        err = std::sqrt(err_sum / (NumSpec+2));

        if (err < 1.0_rt) {
            converged = true;
        }
        iter++;
    }

    err_out = err;

    if (!converged) {
        ierr = CONVERGENCE_FAILURE;
        return;
    }

#endif

    // update the full U_new
    // if we updated total energy, then correct internal,
    // or vice versa
    U_new[URHO] = U_react[0];
    for (int n = 0; n < NumSpec; ++n) {
        U_new[UFS+n] = U_react[1+n];
    }
    auto v2 = 0.0_rt;
    for (int m = 0; m < 3; ++m) {
        v2 += U_new[UMX+m] * U_new[UMX+m];
    }

    if (sdc_solve_for_rhoe == 1) {
        U_new[UEINT] = U_react[NumSpec+1];
        U_new[UEDEN] = U_new[UEINT] + 0.5_rt * v2 / U_new[URHO];
    } else {
        U_new[UEDEN] = U_react[NumSpec+1];
        U_new[UEINT] = U_new[UEDEN] - 0.5_rt * v2 / U_new[URHO];
    }
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
    ierr = CONVERGENCE_FAILURE;

    for (int n = 0; n < NUM_STATE; ++n) {
        U_begin[n] = U_old[n];
    }

    while (nsub < MAX_NSUB && ierr != NEWTON_SUCCESS) {
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

            for (int n = 0; n < NUM_STATE; ++n) {
                U_begin[n] = U_new[n];
            }
        }
        nsub *= 2;
    }
}
#endif

#endif
