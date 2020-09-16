#ifndef CASTRO_SDC_UTIL_H
#define CASTRO_SDC_UTIL_H

#include <Castro.H>
#include <Castro_react_util.H>
#include <Castro_sdc_vode.H>

// error codes 
constexpr int NEWTON_SUCCESS = 0;
constexpr int SINGULAR_MATRIX = -1;
constexpr int CONVERGENCE_FAILURE = -2;

// solvers 
constexpr int NEWTON_SOLVE = 1;
constexpr int VODE_SOLVE = 2;
constexpr int HYBRID_SOLVE = 3;


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
normalize_species_sdc(const int i, const int j, const int k,
                  Array4<Real> const& u) {

    GpuArray<Real, NumSpec> xn;

    Real xn_sum = 0.0_rt;

    for (int n = 0; n < NumSpec; ++n) {
        xn[n] = u(i,j,k,UFS+n);
        xn[n] = amrex::max(small_x * u(i,j,k,URHO), amrex::min(u(i,j,k,URHO), xn[n]));
        xn_sum += xn[n];
    }

    for (int n = 0; n < NumSpec; ++n) {
        xn[n] = u(i,j,k,URHO) * (xn[n] / xn_sum);
        u(i,j,k,UFS+n) = xn[n];
    }
}

// #ifdef REACTIONS

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
sdc_solve(const int i, const int j, const int k,
          Array4<const Real> const& U_old, 
          Array4<Real> const& U_new, 
          Array4<const Real> const& C, 
          const Real dt_m,
          const int sdc_iteration) {

    int ierr;

    // for debugging
    GpuArray<Real, NUM_STATE> U_orig;

    for (int n = 0; n < NUM_STATE; ++n) {
        U_orig[n] = U_old(i,j,k,n);
    }

    if (sdc_solver == NEWTON_SOLVE) {
        // we are going to assume we already have a good guess for the
        // solving in U_new and just pass the solve onto the main Newton
        // solve
        auto err_out = sdc_newton_subdivide(dt_m, U_old, U_new, C, sdc_iteration, &ierr);

        // failing?
        if (ierr != NEWTON_SUCCESS) {
            AllPrint() << "Newton convergence failure" << std::endl;
            AllPrint() << "convergence failure, error = " << err_out << std::endl;
            AllPrint() << "density: " << U_old(i,j,k,URHO) << std::endl;
            AllPrint() << "(old) temperature: " << U_old(i,j,k,UTEMP) << std::endl;
            AllPrint() << "mass fractions:" 
            for (int n = 0; n < NumSpec; ++n) {
                AllPrint() << " " << U_old(i,j,k,UFS+n)/U_old(i,j,k,URHO);
            }
            AllPrint() << std::endl;
            Abort("Newton subcycling failed in sdc_solve");
        }

    } else if (sdc_solver == VODE_SOLVE) {
        // use VODE to do the solution
        sdc_vode_solve(dt_m, U_old, U_new, C, sdc_iteration);

    } else if (sdc_solver == HYBRID_SOLVE) {
        // if it is the first iteration, we will use VODE to predict
        // the solution.  Otherwise, we will use Newton.
        if (sdc_iteration == 0) {
            sdc_vode_solve(dt_m, U_old, U_new, C, sdc_iteration);
        }

        // now U_new is the update that VODE predicts, so we will use
        // that as the initial guess to the Newton solve
        sdc_newton_subdivide(dt_m, U_old, U_new, C, sdc_iteration, &ierr);

        // failing?
        if (ierr != NEWTON_SUCCESS) {
            Abort("Newton failure in sdc_solve");
        }
    }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real
sdc_newton_subdivide(const int i, const int j, const int k,
                     Array4<const Real> const& U_old, 
                     Array4<Real> const& U_new, 
                     Array4<const Real> const& C, 
                     const Real dt_m,
                     const int sdc_iteration
                     int& ierr) {
    // This is the driver for solving the nonlinear update for the
    // reacting/advecting system using Newton's method.  It attempts to
    // do the solution for the full dt_m requested, but if it fails,
    // will subdivide the domain until it converges or reaches our
    // limit on the number of subintervals.

    // subdivide the timestep and do multiple Newtons.  We come in here
    // with an initial guess for the new solution stored in U_new.
    // That only really makes sense for the case where we have 1
    // substep.  Otherwise, we should just use the old time solution.
    auto nsub = 1;
    auto ierr = CONVERGENCE_FAILURE;
    const int MAX_NSUB = 64;

    Real err_out = 0.0_rt;

    GpuArray<Real, NUM_STATE> U_begin;

    for (int n = 0; n < NUM_STATE; ++n) {
        U_begin[n] = U_old(i,j,k,n);
    }

    while (nsub < MAX_NSUB && ierr != NEWTON_SUCCESS) {
        if (nsub > 1) {
            for (int n = 0; n < NUM_STATE; ++n) {
                U_new(i,j,k,n) = U_old(i,j,k,n);
            }
        }

        auto dt_sub = dt_m / nsub;

        for (int isub = 0; isub < nsub; ++isub) {

            // normalize species
            for (int n = 0; n < NumSpec; ++n) {
                U_begin[UFS+n] = amrex::max(SMALL_X, U_begin[UFS+n]);
            }

            auto sum_rhoX = 0.0_rt;
            for (int n = 0; n < NumSpec; ++n) {
                sum_rhoX += U_begin[UFS+n];
            }

            for (int n = 0; n < NumSpec; ++n) {
                U_begin[UFS+n] *= U_begin[URHO] / sum_rhoX;
            }
   
            err_out = sdc_newton_solve(dt_sub, U_begin, U_new, C, sdc_iteration, ierr);

            for (int n = 0; n < NUM_STATE; ++n) {
                U_begin[n] = U_new(i,j,k,n);
            }
        }
        nsub *= 2;
    }
    return err_out;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real
sdc_newton_solve(const int i, const int j, const int k,
                 Array4<const Real> const& U_old, 
                 Array4<Real> const& U_new, 
                 Array4<const Real> const& C, 
                 const Real dt_m,
                 const int sdc_iteration
                 int& ierr) {
    // the purpose of this function is to solve the system
    // U - dt R(U) = U_old + dt C using a Newton solve.
    //
    // here, U_new should come in as a guess for the new U and will be
    // returned with the value that satisfies the nonlinear function

    // we will do the implicit update of only the terms that have reactive sources
    //
    //   0               : rho
    //   1:nspec  : species
    //   nspec+1  : (rho E) or (rho e)

    ierr = NEWTON_SUCCESS;

    // the tolerance we are solving to may depend on the iteration
    auto relax_fac = std::pow(sdc_solver_relax_factor, sdc_order - sdc_iteration - 1);
    auto tol_dens = sdc_solver_tol_dens * relax_fac;
    auto tol_spec = sdc_solver_tol_spec * relax_fac;
    auto tol_ener = sdc_solver_tol_ener * relax_fac;

    // update the momenta for this zone -- they don't react
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        U_new(i,j,k,UMX+n) = U_old(i,j,k,UMX+n) + dt_m * C(i,j,k,UMX+n);
    }

    // update the non-reacting species
    for (int n = 0; n < NumSpec; ++n) {
        U_new(i,j,k,UFS+n) = U_old(i,j,k,UFS+n) + dt_m * C(i,j,k,UFS+n);
    }

    // now only save the subset that participates in the nonlinear
    // solve -- note: we include the old state in f_source

    // load rpar

    // for the Jacobian solve, we are solving
    //   f(U) = U - dt R(U) - U_old - dt C = 0
    // we define f_source = U_old + dt C so we are solving
    //   f(U) = U - dt R(U) - f_source = 0

    GpuArray<Real, NumSpec+1> f_source;

    f_source[0] = U_old(i,j,k,URHO) + dt_m * C(i,j,k,URHO);

    for (int n = 0; n < NumSpec; ++n) {
        f_source[1+n] = U_old(i,j,k,UFS+n) + dt_m * C(i,j,k,UFS+n);
    }
    if (sdc_solve_for_rhoe == 1) {
        f_source[NumSpec+1] = U_old(i,j,k,UEINT) + dt_m * C(i,j,k,UEINT);
    } else {
        f_source[NumSpec+1] = U_old(i,j,k,UEDEN) + dt_m * C(i,j,k,UEDEN);
    }

    GpuArray<Real, n_rpar_comps> rpar;

    for (int n = 0; n < NumSpec+1; ++n) {
        rpar[irp_f_source+n] = f_source[n];
    }
    rpar[irp_dt] = dt_m;

    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        rpar[irp_mom+n] = U_new(i,j,k,UMX+n);
    }

    // temperature will be used as an initial guess in the EOS
    rpar[irp_temp] = U_old(i,j,k,UTEMP);

    // we should be able to do an update for this somehow?
    if (sdc_solve_for_rhoe == 1) {
       rpar[irp_evar] = U_new(i,j,k,UEDEN);
    } else {
       rpar[irp_evar] = U_new(i,j,k,UEINT);
    }

    // store the subset for the nonlinear solve
    // We use an initial guess if possible
    GpuArray<Real, NumSpec+2> U_react;

    U_react[0] = U_new(i,j,k,URHO);
    for (int n = 0; n < NumSpec; ++n) {
        U_react[1+n] = U_new(i,j,k,UFS+n);
    }

    if (sdc_solve_for_rhoe == 1) {
        U_react[nspec+1] = U_new(i,j,k,UEINT);
    } else {
        U_react[nspec+1] = U_new(i,j,k,UEDEN);
    }

#if (INTEGRATOR == 0)

    // do a simple Newton solve

    // iterative loop
    auto iter = 0;
    auto max_newton_iter = MAX_ITER;

    auto err = 1.e30_rt;
    auto converged = false;

    while (!converged && iter < max_newton_iter) {

        GpuArray<Real, NumSpec+2> f;
        GpuArray<Real, NumSpec+2> f_rhs;

        f_sdc_jac(NumSpec+2, U_react, f, Jac, NumSpec+2, info, rpar);

        // solve the linear system: Jac dU_react = -f
        dgefa(Jac, ipvt, info);

        if (info != 0) {
            ierr = SINGULAR_MATRIX;
            return err;
        }

        for (int n = 0; n < NumSpec+2; ++n) {
            f_rhs[n] = -f[n];
        }

        dgesl(Jac, ipvt, f_rhs);

        GpuArray<Real, NumSpec+2> dU_react;

        for (int n = 0; n < NumSpec+2; ++n) {
            dU_react[n] = f_rhs[n];
        }

        // how much of dU_react should we apply?
        eta = 1.0_rt;

        for (int n = 0; n < NumSpec+2; ++n) {
            dU_react[n] *= eta;
            U_react[n] += dU_react[n];
        }

        // we still need to normalize here
        // xn(1:nspec) = U_react(1:nspec)/U_react(0)

        // do k = 1, nspec
        //    xn(k) = max(small_x, xn(k))
        // end do
        // xn(:) = xn(:)/sum(xn)

        // U_react(1:nspec) = U_react(0) * xn(1:nspec)

        GpuArray<Real, NumSpec+2> eps_tot;

        eps_tot[0] = tol_dens * std::abs(U_react[0]) + sdc_solver_atol;
        // for species, atol is the mass fraction limit, so we multiply by density to get a partial density limit
        for (int n = 0; n < NumSpec; ++n) {
            eps_tot[1+n] = tol_spec * std::abs(U_react[1+n]) + sdc_solver_atol * std::abs(U_react[0]);
        }

        eps_tot[NumSpec+1] = tol_ener * std::abs(U_react[NumSpec+1]) + sdc_solver_atol;

        // compute the norm of the weighted error, where the weights are 1/eps_tot
        auto sum_dU = 0.0_rt;
        for (int n = 0; n < NumSpec+2; ++n) {
            sum_dU += dU_react[n] * dU_react[n] / (eps_tot[n] * eps_tot[n]);
        }

        err = std::sqrt(dU_react / (NumSpec+2));

        if (err < 1.0_rt) {
            converged = true;
        }

        iter++;
    }

    err_out = err;

    if (!converged) {
       ierr = CONVERGENCE_FAILURE;
       return err_out;
    }

#endif

    // update the full U_new
    // if we updated total energy, then correct internal, or vice versa
    U_new(i,j,k,URHO) = U_react[0];
    for (int n = 0; n < NumSpec; ++n) {
        U_new(i,j,k,UFS+n) = U_react[1+n];
    }

    auto sum_u = 0.0_rt;
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        sum_u += U_new(i,j,k,UMX+n) * U_new(i,j,k,UMX+n);
    }

    if (sdc_solve_for_rhoe == 1) {
        U_new(i,j,k,UEINT) = U_react[nspec+1];
        U_new(i,j,k,UEDEN) = U_new(i,j,k,UEINT) + 0.5_rt * sum_u / U_new(i,j,k,URHO);
    } else {
        U_new(i,j,k,UEDEN) = U_react[nspec+1];
        U_new(i,j,k,UEINT) = U_new(i,j,k,UEDEN) - 0.5_rt * sum_u / U_new(i,j,k,URHO);
    }

    return err_out;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
sdc_vode_solve(const int i, const int j, const int k,
                 Array4<const Real> const& U_old, 
                 Array4<Real> const& U_new, 
                 Array4<const Real> const& C, 
                 const Real dt_m,
                 const int sdc_iteration) {
    // the purpose of this function is to solve the system the
    // approximate system dU/dt = R + C using the VODE ODE integrator.
    // the solution we get here will then be used as the initial guess
    // to the Newton solve on the real system.

    // we will do the implicit update of only the terms that have reactive sources
    //
    //   0               : rho
    //   1:nspec  : species
    //   nspec+1  : (rho E) or (rho e)

#if (INTEGRATOR == 0)

    // the tolerance we are solving to may depend on the iteration
    auto relax_fac = std::pow(sdc_solver_relax_factor, sdc_order - sdc_iteration - 1);
    auto tol_dens = sdc_solver_tol_dens * relax_fac;
    auto tol_spec = sdc_solver_tol_spec * relax_fac;
    auto tol_ener = sdc_solver_tol_ener * relax_fac;

    // update the momenta for this zone -- they don't react
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        U_new(i,j,k,UMX+n) = U_old(i,j,k,UMX+n) + dt_m * C(i,j,k,UMX+n);
    }

    // update the non-reacting species
    for (int n = 0; n < 2; ++n) {
        U_new(i,j,k,UFS-1+n) = U_old(i,j,k,UFS-1+n) + dt_m * C(i,j,k,UFS-1+n);
    }

    // now only save the subset that participates in the nonlinear
    // solve -- note: we include the old state in f_source

    // load rpar

    // if we are solving the system as an ODE, then we
    // are solving
    //    dU/dt = R(U) + C
    // so we simply pass in C

    GpuArray<Real, NumSpec+2> C_react;

    C_react[0] = C[URHO];
    for (int n = 0; n < NumSpec; ++n) {
        C_react[1+n] = C[UFS+n];
    }
    C_react[NumSpec+1] = C[UEINT];

    for (int n = 0; n < NumSpec+2; ++n) {
        dvode_state.rpar[irp_f_source+n] = C_react[n];
    }
    dvode_state.rpar[irp_dt] = dt_m;
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        dvode_state.rpar[irp_mom+n] = U_new(i,j,k,UMX+n);
    }

    // temperature will be used as an initial guess in the EOS
    dvode_state.rpar(irp_temp) = U_old(i,j,k,UTEMP);

    // we are always solving for rhoe with the VODE predict
    dvode_state.rpar(irp_evar) = U_new(i,j,k,UEDEN);

    // store the subset for the nonlinear solve.  We only consider (rho
    // e), not (rho E).  This is because at present we do not have a
    // method of updating the velocities during the multistep
    // integration

    // Also note that the dvode_state is 1-based, but we'll access it
    // as 0-based in our implementation of the RHS routine

    dvode_state.y(1) = U_old(i,j,k,URHO)
    for (int n = 0; n < NumSpec; ++n) {
        dvode_state.y(2+n) = U_old(i,j,k,UFS+n);
    }
    dvode_state.y(NumSpec+2) = U_old(i,j,k,UEINT);

    dvode_state.istate = 1;

    // set the maximum number of steps allowed
    dvode_state.MXSTEP = 25000;

    dvode_state.T = 0.0_rt;
    dvode_state.TOUT = dt_m;

    if (sdc_use_analytic_jac == 1) {
       dvode_state.jacobian = 1;
    } else {
       dvode_state.jacobian = 2;
    }

    if (use_jacobian_caching) {
       dvode_state.JSV = 1;
    } else {
       dvode_state.JSV = -1;
    }

    // relative tolerances
    GpuArray<Real, NumSpec+2> rtol;
    rtol[0] = tol_dens;
    for (int n = 0; n < NumSpec; ++n) {
        rtol[1+n] = tol_spec;
    }
    rtol[NumSpec+1] = tol_ener;

    // absolute tolerances
    GpuArray<Real, NumSpec+2> atol;
    atol[0] = sdc_solver_atol * U_old(i,j,k,URHO);
    for (int n = 0; n < NumSpec; ++n) {
        atol[1+n] = sdc_solver_atol * U_old(i,j,k,URHO);   // this way, atol is the minimum x
    }
    if (sdc_solve_for_rhoe == 1) {
       atol[NumSpec+1] = sdc_solver_atol * U_old(i,j,k,UEINT);
    } else {
       atol[NumSpec+1] = sdc_solver_atol * U_old(i,j,k,UEDEN);
    }

    for (int n = 0; n < NumSpec+2; ++n) {
        dvode_state.atol[n] = atol[n];
        dvode_state.rtol[n] = rtol[n];
    }

    dvode(dvode_state);

    if (dvode_state.istate < 0) {
        AllPrint() << "VODE error, istate = " << dvode_state.istate << std::endl;
        Abort("vode termination poorly");
    }

    // update the full U_new
    U_new(i,j,k,URHO) = dvode_state.y(1);
    for (int n = 0; n < NumSpec; ++n) {
        U_new(i,j,k,UFS+n) = dvode_state.y(2+n);
    }
    U_new(i,j,k,UEINT) = dvode_state.y(nspec+2);

    auto sum_u2 = 0.0_rt;
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        sum_u2 += U_new(i,j,k,UMX+n) * U_new(i,j,k,UMX+n);
    }

    U_new(i,j,k,UEDEN) = U_new(i,j,k,UEINT) + 0.5_rt * sum_u2 / U_new(i,j,k,URHO);

    // keep our temperature guess
    U_new(i,j,k,UTEMP) = dvode_state.rpar(irp_temp);

#endif

}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
sdc_update_o2(const int i, const int j, const int k,
                 Array4<const Real> const& k_m, 
                 Array4<Real> const& k_n, 
                 Array4<const Real> const& A_m, 
                 Array4<const Real> const& R_m_old,
                 Array4<const Real> const& C, 
                 const Real dt_m,
                 const int sdc_iteration, const int m_start) {
    // update k_m to k_n via advection -- this is a second-order accurate update

    // Here, dt_m is the timestep between time-nodes m and m+1

    burn_t burn_state;

    GpuArray<Real, NUM_STATE> U_old;
    GpuArray<Real, NUM_STATE> U_new;
    GpuArray<Real, NUM_STATE> R_full;
    GpuArray<Real, NUM_STATE> C_zone;

    for (int n = 0; n < NUM_STATE; ++n) {
        U_old[n] = k_m(i,j,k,n);
        C_zone[n] = C(i,j,k,n);
    }

    // only burn if we are within the temperature and density
    // limits for burning
    if (!okay_to_burn(U_old)) {
        for (int n = 0; n < NUM_STATE; ++n) {
            R_full[n] = 0.0_rt;
        }
    } else {

        // this is the full state -- this will be updated as we
        // solve the nonlinear system.  We want to start with a
        // good initial guess.  For later iterations, we should
        // begin with the result from the previous iteration.  For
        // the first iteration, let's try to extrapolate forward
        // in time.
        if (sdc_iteration == 0) {
            for (int n = 0; n < NUM_STATE; ++n) {
                U_new[n] = U_old[n] + dt_m * A_m(i,j,k,n) + dt_m * R_m_old(i,j,k,n);
            }
        } else {
            for (int n = 0; n < NUM_STATE; ++n) {
                U_new[n] = k_n(i,j,k,n);
            }
        }

        sdc_solve(dt_m, U_old, U_new, C, sdc_iteration);

        // we solved our system to some tolerance, but let's be sure we are conservative by
        // reevaluating the reactions and { doing the full step update
        single_zone_react_source(i, j, k, U_new, R_full, burn_state);

    }

    for (int n = 0; n < NUM_STATE; ++n) {
        U_new[n] = U_old[n] + dt_m * R_full[n] + dt_m * C_zone[n];

        // copy back to k_n
        k_n(i,j,k,n) = U_new[n];
    }
}


// AMREX_GPU_HOST_DEVICE AMREX_INLINE
// void
// sdc_update_centers_o4(const int i, const int j, const int k,
//                         Array4<const Real> const& U_old, 
//                         Array4<Real> const& U_new, 
//                         Array4<const Real> const& C, 
//                         const Real dt_m,
//                         const int sdc_iteration) {
//     // update U_old to U_new on cell-centers.  This is an implicit
//     // solve because of reactions.  Here U_old corresponds to time node
//     // m and U_new is node m+1.  dt_m is the timestep between m and
//     // m+1

//     // we come in with U_new being a guess for the updated solution
//     if (okay_to_burn(U_old(i,j,k,:))) {
//         sdc_solve(i, j, k, dt_m, U_old, U_new, C, sdc_iteration);
//     } else {
//         // no reactions, so it is a straightforward update
//         for (int n = 0; n < NUM_STATE; ++n) {
//             U_new(i,j,k,n) = U_old(i,j,k,n) + dt_m * C(i,j,k,n);
//         }
//     }
// }

// AMREX_GPU_HOST_DEVICE AMREX_INLINE
// void
// f_sdc_jac(const int neq, 
//           GpuArray<const Real, neq> const& U, 
//           GpuArray<const Real, neq> const& f, 
//           Array2d<Real, 0, ldjac, 0, neq> Jac, 
//           const int ldjac, int iflag, 
//           GpuArray<Real, n_rpar_comps> const& rpar) {
//     // this is used with the Newton solve and returns f and the Jacobian

//     GpuArray<Real, NUM_STATE> U_full;
//     GpuArray<Real, NUM_STATE> R_full;
//     GpuArray<Real, neq> R_react;
//     GpuArray<Real, neq> f_source;
//     burn_t burn_state;
//     eos_t eos_state;
//     eos_xderivs_t eos_xderivs;

//     Array2d<Real, 0, NumSpec+1, 0, NumSpec+1> dRdw;
//     Array2d<Real, 0, NumSpec+1, 0, NumSpec+1> dwdU;

//     // we are not solving the momentum equations
//     // create a full state -- we need this for some interfaces
//     U_full[URHO] = U[0];

//     for (int n = 0; n < NumSpec; ++n) {
//         U_full[UFS+n] = U[1+n];
//     }

//     if (sdc_solve_for_rhoe == 1) {
//         U_full[UEINT] = U[NumSpec];
//         U_full[UEDEN] = rpar[irp_evar];
//     } else {
//         U_full[UEDEN] = U[NumSpec];
//         U_full[UEINT] = rpar[irp_evar];
//     }

//     for (int n = 0; n < AMREX_SPACEDIM; ++n) {
//         U_full[UMX+n] = rpar[irp_mom+n];
//     }

//     // normalize the species
//     auto sum_rhoX = 0.0_rt;
//     for (int n = 0; n < NumSpec; ++n) {
//         U_full[UFS+n] = amrex::max(small_x, U_full[UFS+n]);
//         sum_rhoX += U_full[UFS+n];
//     }

//     for (int n = 0; n < NumSpec; ++n) {
//         U_full[UFS+n] *= U_full[URHO] / sum_rhoX;
//     }

//     // unpack rpar
//     auto dt_m = rpar[irp_dt];
//     for (int n = 0; n < neq; ++n) {
//         f_source[n] = rpar[irp_f_source+n];
//     }

//     // compute the temperature and species derivatives --
//     // maybe this should be done using the burn_state
//     // returned by single_zone_react_source, since it is
//     // more consistent T from e
//     eos_state.rho = U_full[URHO]
//     eos_state.T = rpar(irp_temp)   // initial guess
//     for (int n = 0; n < NumSpec; ++n) {
//         eos_state.xn[n] = U_full[UFS+n]/U_full[URHO];
//     }
//     for (int n = 0; n < NumAux; ++n) {
//         eos_state.aux[n] = U_full[UFX+n]/U_full[URHO];
//     }
//     eos_state.e = U_full[UEINT] / U_full[URHO];  // (U_full[UEDEN] - 0.5_rt*sum(U_full(UMX:UMZ))/U_full[URHO])/U_full[URHO]

//     eos(eos_input_re, eos_state);

//     U_full[UTEMP] = eos_state.T;

//     single_zone_react_source(U_full, R_full, 0, 0, 0, burn_state);

//     // store the subset of R used in the Jacobian
//     R_react[0] = R_full[URHO]
//     for (int n = 0; n < NumSpec; ++n) {
//         R_react[n+1] = R_full[UFS+n];
//     }
//     if (sdc_solve_for_rhoe == 1) {
//         R_react[NumSpec+1] = R_full[UEINT];
//     } else {
//         R_react[NumSpec+1] = R_full[UEDEN];
//     }

//     for (int n = 0; n < neq; ++n) {
//         f[n] = U[n] - dt_m * R_react[n] - f_source[n];
//     }

//     // get dRdw -- this may do a numerical approxiation or use the
//     // network's analytic Jac
//     single_zone_jac(U_full, burn_state, dRdw);

//     // construct dwdU
//     for (int n = 0; n < NumSpec+1; ++n) {
//         for (int m = 0; m < NumSpec+1; ++m) {
//             dwdU(n,m) = 0.0_rt;
//         }
//     }

//     // the density row
//     dwdU(iwrho, 0) = 1.0_rt;

//     // the X_k rows
//     for (int m = 1; m < NumSpec+1; ++m) {
//         dwdU(iwfs-1+m,0) = -U[m] / (U[0]*U[0])
//         dwdU(iwfs-1+m,m) = 1.0_rt / U[0];
//     }

//     composition_derivatives(eos_state, eos_xderivs)

//     // now the T row -- this depends on whether we are evolving (rho E) or (rho e)
//     auto denom = 1.0_rt / (eos_state.rho * eos_state.dedT);
//     auto xn_sum = 0.0_rt;
//     for (int n = 0; n < NumSpec; ++n) {
//         xn_sum += eos_state.xn[n] * eos_xderivs.dedX[n];
//     }

//     if (sdc_solve_for_rhoe == 1) {
//         dwdU(iwT,0) = denom*(xn_sum - eos_state.rho * eos_state.dedr - eos_state.e);
//     } else {
//         auto u2_sum = 0.0_rt;
//         for (auto n = 0; n < AMREX_SPACEDIM; ++n) {
//             u2_sum += U_full[UMX+n] * U_full[UMX+n];
//         }
//         dwdU(iwT,0) = denom * (xn_sum - eos_state.rho * eos_state.dedr - eos_state.e - 
//                                0.5_rt * u2_sum / (eos_state.rho*eos_state.rho));
//     }

//     for (int m = 1; m < NumSpec+1; ++m) {
//         dwdU(iwT,m) = -denom * eos_xderivs.dedX[m];
//     }

//     dwdU(iwT, nspec+1) = denom;

//     // construct the Jacobian -- we can get most of the
//     // terms from the network itself, but we do not rely on
//     // it having derivative wrt density
//     for (int n = 0; n < ldjac; ++n) {
//         for (int m = 0; m < neq; ++m) {
//             Jac(n,m) = 0.0_rt;
//         }
//     }

//     for (int m = 1; m < NumSpec+1; ++m) {
//         Jac(m, m) = 1.0_rt;
//     }

//     auto dRdU = matmul(dRdw, dwdU);

//     for (int n = 0; n < ldjac; ++n) {
//         for (int m = 0; m < neq; ++m) {
//             Jac(n,m) -= dt_m * dRdU(n,m);
//         }
//     }
// }


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
instantaneous_react(const int i, const int j, const int k,
                    Array4<const Real> const& state, 
                    Array4<Real> const& R_source) {

    // convert from cons to prim -- show this be here or in C++-land?
    // or should I do things like we do in burn_state and convert it manually?
    // (in that case, I am not sure if I can assume UTEMP is defined)

    if (okay_to_burn(i, j, k, state)) {
        burn_t burn_state;
        single_zone_react_source(i, j, k, state, R_source, burn_state);
    } else {
        for (int n = 0; n < NUM_STATE; ++n) {
            R_source(i,j,k,n) = 0.0_rt;
        }
    }
}
// #endif


#endif