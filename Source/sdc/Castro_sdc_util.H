#ifndef CASTRO_SDC_UTIL_H
#define CASTRO_SDC_UTIL_H

#include <Castro.H>
#ifdef REACTIONS
#ifdef NEW_NETWORK_IMPLEMENTATION
#include <rhs.H>
#endif
#include <Castro_react_util.H>
#include <sdc_newton_solve.H>
#include <vode_dvode.H>
#ifndef NEW_NETWORK_IMPLEMENTATION
#include <linpack.H>
#endif
#endif

// solvers
constexpr int NEWTON_SOLVE = 1;
constexpr int VODE_SOLVE = 2;
constexpr int HYBRID_SOLVE = 3;

constexpr int ldjac = NumSpec + 2;


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
normalize_species_sdc(const int i, const int j, const int k,
                      Array4<Real> const& u) {

    GpuArray<Real, NumSpec> xn;

    Real xn_sum = 0.0_rt;

    for (int n = 0; n < NumSpec; ++n) {
        xn[n] = u(i,j,k,UFS+n);
        xn[n] = amrex::max(network_rp::small_x * u(i,j,k,URHO), amrex::min(u(i,j,k,URHO), xn[n]));
        xn_sum += xn[n];
    }

    for (int n = 0; n < NumSpec; ++n) {
        xn[n] *= u(i,j,k,URHO) / xn_sum;
        u(i,j,k,UFS+n) = xn[n];
    }
}

#ifdef REACTIONS



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

    // The tolerance we are solving to may depend on the
    // iteration
    auto relax_fac = std::pow(sdc_solver_relax_factor, sdc_order - sdc_iteration - 1);
    auto tol_dens = sdc_solver_tol_dens * relax_fac;
    auto tol_spec = sdc_solver_tol_spec * relax_fac;
    auto tol_ener = sdc_solver_tol_ener * relax_fac;

    // Update the momenta for this zone -- they don't react
    for (int n = 0; n < 3; ++n) {
        U_new[UMX+n] = U_old[UMX+n] + dt_m * C[UMX+n];
    }

    // Now only save the subset that participates in the
    // nonlinear solve -- note: we include the old state
    // in f_source

    // load rpar

    // If we are solving the system as an ODE, then we are
    // solving
    //   dU/dt = R(U) + C
    // so we simply pass in C
    GpuArray<Real, NumSpec+2> C_react;

    C_react[0] = C[URHO];
    for (int n = 0; n < NumSpec; ++n) {
        C_react[1+n] = C[UFS+n];
    }
    C_react[NumSpec+1] = C[UEINT];

    dvode_t<NumSpec+2> dvode_state;

    burn_t burn_state;

    burn_state.success = true;

    for (int n = 0; n < NumSpec+2; n++) {
        burn_state.f_source[n] = C_react[n];
    }

    for (int n = 0; n < 3; n++) {
        burn_state.mom[n] = U_new[UMX+n];
    }

    // we are always solving for rhoe with the VODE predict
    burn_state.e = U_new[UEDEN];

    burn_state.E_var = U_old[UEDEN];

    // temperature will be used as an initial guess in the EOS
    burn_state.T = U_old[UTEMP];


    // store the subset for the nonlinear solve. We only
    // consider (rho e), not (rho E). This is because at
    // present we do not have a method of updating the
    // velocities during the multistep integration.

    // Also note that the dvode_state is 1-based, but we'll
    // access it as 0-based in our implementation of the
    // RHS routine

    dvode_state.y(1) = U_old[URHO];
    for (int n = 0; n < NumSpec; ++n) {
        dvode_state.y(2+n) = U_old[UFS+n];
    }
    dvode_state.y(NumSpec+2) = U_old[UEINT];

    // // set the maximum number of steps allowed
    // dvode_state.MXSTEP = 25000;

    dvode_state.t = 0.0_rt;
    dvode_state.tout = dt_m;
    dvode_state.HMXI = 1.0_rt / ode_max_dt;

    // relative tolerances
    dvode_state.rtol_dens = tol_dens;
    dvode_state.rtol_spec = tol_spec;
    dvode_state.rtol_enuc = tol_ener;

    // absolute tolerances
    dvode_state.atol_dens = sdc_solver_atol * U_old[URHO];
    dvode_state.atol_spec = sdc_solver_atol * U_old[URHO]; // this way, atol is the minimum x
    dvode_state.atol_enuc = sdc_solver_atol * U_old[UEINT];

    int istate = dvode(burn_state, dvode_state);

    if (istate < 0) {
        Abort("Vode terminated poorly");
    }

    // update the full U_new
    U_new[URHO] = dvode_state.y(1);
    for (int n = 0; n < NumSpec; ++n) {
        U_new[UFS+n] = dvode_state.y(2+n);
    }
    U_new[UEINT] = dvode_state.y(NumSpec+2);
    auto v2 = 0.0_rt;
    for (int n = 0; n < 3; ++n) {
        v2 += U_new[UMX+n] * U_new[UMX+n];
    }
    U_new[UEDEN] = U_new[UEINT] + 0.5_rt * v2 / U_new[URHO];

    // keep our temperature guess
    U_new[UTEMP] = U_old[UTEMP];

#endif
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
sdc_solve(const Real dt_m,
          GpuArray<Real, NUM_STATE> const& U_old,
          GpuArray<Real, NUM_STATE>& U_new,
          GpuArray<Real, NUM_STATE> const& C,
          const int sdc_iteration) {

    int ierr;
    Real err_out;

    // for debugging
    GpuArray<Real, NUM_STATE> U_orig;

    for (int n = 0; n < NUM_STATE; ++n) {
        U_orig[n] = U_old[n];
    }

    if (sdc_solver == NEWTON_SOLVE) {
        // We are going to assume we already have a good guess
        // for the solve in U_new and just pass the solve onto
        // the main Newton solve
        sdc_newton_subdivide(dt_m, U_old, U_new, C, sdc_iteration, err_out, ierr);

        // failing?
        if (ierr != NEWTON_SUCCESS) {
            Abort("Newton subcycling failed in sdc_solve");
        }
    } else if (sdc_solver == VODE_SOLVE) {
        // Use VODE to do the solution
        sdc_vode_solve(dt_m, U_old, U_new, C, sdc_iteration);

    } else if (sdc_solver == HYBRID_SOLVE) {
        // If it is the first iteration, we will use VODE to
        // predict the solution. Otherwise, we will use Newton.
        if (sdc_iteration == 0) {
            sdc_vode_solve(dt_m, U_old, U_new, C, sdc_iteration);
        }

        // Now U_new is the update that VODE predicts, so we
        // will use that as the initial guess to the Newton solve
        sdc_newton_subdivide(dt_m, U_old, U_new, C, sdc_iteration, err_out, ierr);

        // Failing?
        if (ierr != NEWTON_SUCCESS) {
            Abort("Newton failure in sdc_solve");
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
sdc_solve(const int i, const int j, const int k,
          const Real dt_m,
          Array4<const Real> const& U_old,
          Array4<Real> const& U_new,
          Array4<const Real> const& C,
          const int sdc_iteration) {
    // wrapper for the zone-by-zone version

    GpuArray<Real, NUM_STATE> U_old_zone;
    GpuArray<Real, NUM_STATE> U_new_zone;
    GpuArray<Real, NUM_STATE> C_zone;

    for (int n = 0; n < NUM_STATE; ++n) {
        U_old_zone[n] = U_old(i,j,k,n);
        U_new_zone[n] = U_new(i,j,k,n);
        C_zone[n] = C(i,j,k,n);
    }

    sdc_solve(dt_m, U_old_zone, U_new_zone, C_zone, sdc_iteration);

    for (int n = 0; n < NUM_STATE; ++n) {
        U_new(i,j,k,n) = U_new_zone[n];
    }
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

    // Only burn if we are within the temperature and density
    // limits for burning
    if (!okay_to_burn(U_old)) {
        for (int n = 0; n < NUM_STATE; ++n) {
            R_full[n] = 0.0_rt;
        }
    } else {

        // This is the full state -- this will be updated as we
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

        sdc_solve(dt_m, U_old, U_new, C_zone, sdc_iteration);

        // we solved our system to some tolerance, but let's be sure we are conservative by
        // reevaluating the reactions and { doing the full step update
        single_zone_react_source(U_new, R_full, burn_state);
    }

    for (int n = 0; n < NUM_STATE; ++n) {
        U_new[n] = U_old[n] + dt_m * R_full[n] + dt_m * C_zone[n];

        // copy back to k_n
        k_n(i,j,k,n) = U_new[n];
    }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
sdc_update_centers_o4(const int i, const int j, const int k,
                      Array4<const Real> const& U_old,
                      Array4<Real> const& U_new,
                      Array4<const Real> const& C,
                      const Real dt_m,
                      const int sdc_iteration) {
    // Update U_old to U_new on cell-centers.  This is an implicit
    // solve because of reactions.  Here U_old corresponds to time node
    // m and U_new is node m+1.  dt_m is the timestep between m and
    // m+1

    // We come in with U_new being a guess for the updated solution
    if (okay_to_burn(i, j, k, U_old)) {
        sdc_solve(i, j, k, dt_m, U_old, U_new, C, sdc_iteration);
    } else {
        // no reactions, so it is a straightforward update
        for (int n = 0; n < NUM_STATE; ++n) {
            U_new(i,j,k,n) = U_old(i,j,k,n) + dt_m * C(i,j,k,n);
        }
    }
}

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
            R_source(i, j, k, n) = 0.0_rt;
        }
    }
}

#endif

#endif
