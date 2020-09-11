#ifndef CASTRO_REACT_UTIL_H
#define CASTRO_REACT_UTIL_H

#include <Castro.H>


AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool
okay_to_burn(const int i, const int j, const int k,
             Array4<const Real> const& state) {

    if (state(i,j,k,UTEMP) < react_T_min || state(i,j,k,UTEMP) > react_T_max ||
        state(i,j,k,URHO) < react_rho_min || state(i,j,k,URHO) > react_rho_max) {
            return false;
    }

    return true;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool
okay_to_burn_type(burn_t const& state) {

    if (state.T < react_T_min || state.T > react_T_max ||
        state.rho < react_rho_min || state.rho > react_rho_max) {
            return false;
    }

    return true;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
single_zone_react_source(const int i, const int j, const int k,
             Array4<const Real> const& state,
             Array4<Real> const& R,
             burn_t burn_state) {

    rhoInv = 1.0_rt / state(i,j,k,URHO);

    burn_state.rho = state(i,j,k,URHO);
    burn_state.T = state(i,j,k,UTEMP);
    burn_state.e = state(i,j,k,UEINT) * rhoInv;

    for (int n = 0; n < NumSpec; ++n) {
       burn_state.xn[n] = amrex::max(amrex::min(state(i,j,k,UFS+n) * rhoInv, 1.0_rt), small_x);
    }

#if NAUX_NET > 0
    for (int n = 0; n < NumAux; ++n) {
       burn_state.aux[n] = state(i,j,k,UFX+n) * rhoInv;
    }
#endif

    eos_t eos_state;

    // Ensure that the temperature going in is consistent with the internal energy.
    burn_to_eos(burn_state, eos_state);
    eos(eos_input_re, eos_state);
    eos_to_burn(eos_state, burn_state);

    eos_get_small_temp(small_temp);
    burn_state.T = amrex::min(max_temp, amrex::max(burn_state.T, small_temp));

    burn_state.self_heat = false;

    GpuArray<Real, neqs> ydot;

    actual_rhs(burn_state, ydot);

    // store the instantaneous R
    for (int n = 0; n < NVAR; ++n) {
        R(i,j,k,n) = 0.0_rt;
    }

    // species rates come back in terms of molar fractions
    for (int n = 0; n < NumSpec; ++n) {
    R(i,j,k,UFS+n) = state(i,j,k,URHO) * aion[n] * ydot[n];
    }

    R(i,j,k,UEDEN) = state(i,j,k,URHO) * ydot[net_ienuc];
    R(i,j,k,UEINT) = state(i,j,k,URHO) * ydot[net_ienuc];
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
single_zone_jac(const int i, const int j, const int k,
             Array4<const Real> const& state,
             burn_t burn_state,
             Array2D<Real, 0, NumSpec+1, 0, NumSpec+1> const& dRdw) {

#ifdef SIMPLIFIED_SDC
#ifndef AMREX_USE_GPU
    amrex::Error("we shouldn't be here with the simplified SDC method (USE_SIMPLIFIED_SDC=TRUE)");
#endif
#else

    GpuArray<Real, neqs> ydot;
    GpuArray<Real, neqs> ydot_pert;
    Array2D<Real, 0, neqs, 0, neqs> jac;

    const auto eps = 1.e-8_rt;

    actual_rhs(burn_state, ydot);

    if (sdc_use_analytic_jac == 0) {
        // note the numerical Jacobian will be returned in terms of X
        numerical_jac(burn_state, jac);
    } else {
        actual_jac(burn_state, jac);

        // The Jacobian from the nets is in terms of dYdot/dY, but we want
        // it was dXdot/dX, so convert here.
        for (int n = 0; n < NumSpec; ++n) {
            for (int m = 0; m < NumSpec; ++m) {
                jac(n,m) = jac(n,m) * aion[n];
                jac(m,n) = jac(m,n) * aion_inv[n];
            }
        }

    }
#endif

    // at this point, our Jacobian should be entirely in terms of X,
    // not Y.  Let's now fix the rhs terms themselves to be in terms of
    // dX/dt and not dY/dt.
    for (int n = 0; n < NumSpec; ++n) {
        ydot[n] *= aion[n];
    }

    // Our jacobian, dR/dw has the form:
    //
    //  /      0                  0                  0                       0          \
    //  | d(rho X1dot)/drho  d(rho X1dot)/dX1   d(rho X1dit)/dX2   ...  d(rho X1dot)/dT |
    //  | d(rho X2dot)/drho  d(rho X2dot)/dX1   d(rho X2dot)/dX2   ...  d(rho X2dot)/dT |
    //  |   ...                                                                         |
    //  \ d(rho Edot)/drho   d(rho Edot)/dX1    d(rho Edot)/dX2    ...  d(rho Edot)/dT  /

    for (int n = 0; n < NumSpec; ++n) {
        for (int m = 0; m < NumSpec; ++m) {
            dRdw(n,m) = 0.0_rt;
        }
    }

    // now perturb density and call the RHS to compute the derivative wrt rho
    // species rates come back in terms of molar fractions
    copy_burn_t(burn_state_pert, burn_state);
    burn_state_pert.rho = burn_state.rho * (1.0_rt + eps);

    actual_rhs(burn_state_pert, ydot_pert);

    // make the rates dX/dt and not dY/dt
    for (int n = 0; n < NumSpec; ++n) {
        ydot_pert[n] *= aion[n];
    }

    // fill the column of dRdw corresponding to the derivative
    // with respect to rho
    for (int m = 0; m < NumSpec; ++m) {
       // d( d(rho X_m)/dt)/drho
       dRdw(m, iwrho) = ydot[m] + 
            state(i,j,k,URHO) * (ydot_pert[m] - ydot[m])/(eps * burn_state.rho)
    }

    // d( d(rho E)/dt)/drho
    dRdw(nspec_evolve+1, iwrho) = ydot[net_ienuc] + 
         state(i,j,k,URHO) * (ydot_pert[net_ienuc] - ydot[net_ienuc])/(eps * burn_state.rho);

    // fill the columns of dRdw corresponding to each derivative
    // with respect to species mass fraction
    for (int n = 0; n < NumSpec; ++n) {
       dRdw(0, iwfs+n) = 0.0_rt;  // density source

       for (int m = 0; m < NumSpec; ++m) {
          // d( d(rho X_m)/dt)/dX_n
          dRdw(m, iwfs+n) = state(i,j,k,URHO) * jac(m, n);
       }

       // d( d(rho E)/dt)/dX_n
       dRdw(NumSpec, iwfs+n) = state(i,j,k,URHO) * jac(net_ienuc, n);

    }

    // now fill the column corresponding to derivatives with respect to
    // temperature -- this column is iwT
    dRdw(0, iwT) = 0.0_rt;

    // d( d(rho X_m)/dt)/dT
    for (int m = 0; m < NumSpec; ++m) {
       dRdw(m, iwT) = state(i,j,k,URHO) * jac(m, net_itemp);
    }

    // d( d(rho E)/dt)/dT
    dRdw(NumSpec, iwT) = state(i,j,k,URHO) * jac(net_ienuc, net_itemp);
}

#endif