#ifndef CASTRO_SDC_UTIL_H
#define CASTRO_SDC_UTIL_H

#include <Castro.H>

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
normalize_species(const int i, const int j, const int k,
                  Array4<const Real> const& u) {

    GpuArray<Real, NumSpec> xn;

    Real xn_sum = 0.0_rt;

    for (int n = 0; n < NumSpec; ++n) {
        xn[n] = u(i,j,k,UFS+n);
        xn[n] = amrex::max::(small_x * u(i,j,k,URHO), amrex::min(u(i,j,k,URHO), xn[n]));
        xn_sum += xn[n];
    }

    for (int n = 0; n < NumSpec; ++n) {
        xn[n] = u(i,j,k,URHO) * (xn[n] / xn_sum);
        u(i,j,k,UFS+n) = xn[n];
    }
}

#ifdef REACTIONS

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

    GpuArray<Real, NVAR> U_old;
    GpuArray<Real, NVAR> U_new;
    GpuArray<Real, NVAR> R_full;
    GpuArray<Real, NVAR> C_zone;

    for (int n = 0; n < NVAR; ++n) {
        U_old[n] = k_m(i,j,k,n);
        C_zone[n] = C(i,j,k,n);
    }

    // only burn if we are within the temperature and density
    // limits for burning
    if (!okay_to_burn(U_old)) {
        for (int n = 0; n < NVAR; ++n) {
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
            for (int n = 0; n < NVAR; ++n) {
                U_new[n] = U_old[n] + dt_m * A_m(i,j,k,n) + dt_m * R_m_old(i,j,k,n);
            }
        } else {
            for (int n = 0; n < NVAR; ++n) {
                U_new[n] = k_n(i,j,k,n);
            }
        }

        sdc_solve(dt_m, U_old, U_new, C, sdc_iteration);

        // we solved our system to some tolerance, but let's be sure we are conservative by
        // reevaluating the reactions and { doing the full step update
        single_zone_react_source(i, j, k, U_new, R_full, burn_state);

    }

    for (int n = 0; n < NVAR; ++n) {
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
    // update U_old to U_new on cell-centers.  This is an implicit
    // solve because of reactions.  Here U_old corresponds to time node
    // m and U_new is node m+1.  dt_m is the timestep between m and
    // m+1

    // we come in with U_new being a guess for the updated solution
    if (okay_to_burn(U_old(i,j,k,:))) {
        sdc_solve(i, j, k, dt_m, U_old, U_new, C, sdc_iteration);
    } else {
        // no reactions, so it is a straightforward update
        for (int n = 0; n < NVAR; ++n) {
            U_new(i,j,k,n) = U_old(i,j,k,n) + dt_m * C(i,j,k,n);
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
f_sdc_jac(const int neq, 
          GpuArray<const Real, neq> const& U, 
          GpuArray<const Real, neq> const& f, 
          Array2d<Real, 0, ldjac, 0, neq> Jac, 
          const int ldjac, int iflag, 
          GpuArray<Real, n_rpar_comps> const& rpar) {
    // this is used with the Newton solve and returns f and the Jacobian

    GpuArray<Real, NVAR> U_full;
    GpuArray<Real, NVAR> R_full;
    GpuArray<Real, neq> R_react;
    GpuArray<Real, neq> f_source;
    burn_t burn_state;
    eos_t eos_state;
    eos_xderivs_t eos_xderivs;

    Array2d<Real, 0, NumSpec+1, 0, NumSpec+1> dRdw;
    Array2d<Real, 0, NumSpec+1, 0, NumSpec+1> dwdU;

    // we are not solving the momentum equations
    // create a full state -- we need this for some interfaces
    U_full[URHO] = U[0];

    for (int n = 0; n < NumSpec; ++n) {
        U_full[UFS+n] = U[1+n];
    }

    if (sdc_solve_for_rhoe == 1) {
       U_full[UEINT] = U[NumSpec];
       U_full[UEDEN] = rpar[irp_evar];
    } else {
       U_full[UEDEN] = U[NumSpec];
       U_full[UEINT] = rpar[irp_evar];
    }

    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        U_full[UMX+n] = rpar[irp_mom+n];
    }

    // normalize the species
    auto sum_rhoX = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
       U_full[UFS+n] = amrex::max(small_x, U_full[UFS+n]);
       sum_rhoX += U_full[UFS+n];
    }

    for (int n = 0; n < NumSpec; ++n) {
        U_full[UFS+n] *= U_full[URHO] / sum_rhoX;
    }

    // unpack rpar
    auto dt_m = rpar[irp_dt];
    for (int n = 0; n < neq; ++n) {
        f_source[n] = rpar[irp_f_source+n];
    }

    // compute the temperature and species derivatives --
    // maybe this should be done using the burn_state
    // returned by single_zone_react_source, since it is
    // more consistent T from e
    eos_state.rho = U_full[URHO]
    eos_state.T = rpar(irp_temp)   // initial guess
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = U_full[UFS+n]/U_full[URHO];
    }
    for (int n = 0; n < NumAux; ++n) {
        eos_state.aux[n] = U_full[UFX+n]/U_full[URHO];
    }
    eos_state.e = U_full[UEINT] / U_full[URHO];  // (U_full[UEDEN] - HALF*sum(U_full(UMX:UMZ))/U_full[URHO])/U_full[URHO]

    eos(eos_input_re, eos_state);

    U_full[UTEMP] = eos_state.T;

    single_zone_react_source(U_full, R_full, 0, 0, 0, burn_state);

    // store the subset of R used in the Jacobian
    R_react[0] = R_full[URHO]
    for (int n = 0; n < NumSpec; ++n) {
        R_react[n+1] = R_full[UFS+n];
    }
    if (sdc_solve_for_rhoe == 1) {
       R_react[NumSpec+1] = R_full[UEINT];
    } else {
       R_react[NumSpec+1] = R_full[UEDEN];
    }

    for (int n = 0; n < neq; ++n) {
        f[n] = U[n] - dt_m * R_react[n] - f_source[n];
    }

    // get dRdw -- this may do a numerical approxiation or use the
    // network's analytic Jac
    single_zone_jac(U_full, burn_state, dRdw);

    // construct dwdU
    for (int n = 0; n < NumSpec+1; ++n) {
        for (int m = 0; m < NumSpec+1; ++m) {
            dwdU(n,m) = 0.0_rt;
        }
    }

    // the density row
    dwdU(iwrho, 0) = 1.0_rt;

    // the X_k rows
    for (int m = 1; m < NumSpec+1; ++m) {
       dwdU(iwfs-1+m,0) = -U[m]/(U[0]*U[0])
       dwdU(iwfs-1+m,m) = 1.0_rt / U[0];
    }

    composition_derivatives(eos_state, eos_xderivs)

    // now the T row -- this depends on whether we are evolving (rho E) or (rho e)
    auto denom = 1.0_rt/(eos_state.rho * eos_state.dedT);
    if (sdc_solve_for_rhoe == 1) {
       dwdU(iwT,0) = denom*(sum(eos_state.xn(1:nspec) * eos_xderivs.dedX(1:nspec)) - &
                                       eos_state.rho * eos_state.dedr - eos_state.e)
    } else {
       dwdU(iwT,0) = denom*(sum(eos_state.xn(1:nspec) * eos_xderivs.dedX(1:nspec)) - &
                                       eos_state.rho * eos_state.dedr - eos_state.e - &
                                       HALF*sum(U_full(UMX:UMZ)**2)/eos_state.rho**2)
    }

    for (int m = 1; m < NumSpec+1; ++m) {
       dwdU(iwT,m) = -denom * eos_xderivs.dedX[m];
    }

    dwdU(iwT, nspec+1) = denom;

    // construct the Jacobian -- we can get most of the
    // terms from the network itself, but we do not rely on
    // it having derivative wrt density
    for (int n = 0; n < ldjac; ++n) {
        for (int m = 0; m < neq; ++m) {
            Jac(n,m) = 0.0_rt;
        }
    }

    for (int m = 1; m < NumSpec+1; ++m) {
       Jac(m, m) = 1.0_rt;
    }

    auto dRdU = matmul(dRdw, dwdU);

    for (int n = 0; n < ldjac; ++n) {
        for (int m = 0; m < neq; ++m) {
            Jac(n,m) -= dt_m * dRdU(n,m);
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

    if (okay_to_burn(i,j,k,state)) {
        single_zone_react_source(i, j, k, state, R_source, burn_state);
    } else {
        for (int n = 0; n < NVAR; ++n) {
            R_source(i,j,k,n) = 0.0_rt;
        }
    }
}
#endif


#endif