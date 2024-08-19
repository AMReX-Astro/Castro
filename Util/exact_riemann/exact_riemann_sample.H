#ifndef RIEMANN_SAMPLE_MODULE_H
#define RIEMANN_SAMPLE_MODULE_H

#include <AMReX_REAL.H>

#include <riemann_constants.H>
#include <exact_riemann_rarefaction.H>

using namespace amrex::literals;

AMREX_INLINE
void
riemann_sample(const amrex::Real rho_l, const amrex::Real u_l, const amrex::Real p_l, const amrex::Real* xn_l,
               const amrex::Real rho_r, const amrex::Real u_r, const amrex::Real p_r, const amrex::Real* xn_r,
               const amrex::Real ustar, const amrex::Real pstar,
               const amrex::Real W_l, const amrex::Real W_r,
               const amrex::Real x, const amrex::Real xjump, const amrex::Real time,
               amrex::Real& rho, amrex::Real& u, amrex::Real& p, amrex::Real* xn) {

    // get the initial sound speeds

    eos_t eos_state;
    eos_state.rho = rho_l;
    eos_state.p = p_l;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn_l[n];
    }
    eos_state.T = castro::T_guess;

    eos(eos_input_rp, eos_state);

    amrex::Real cs_l = std::sqrt(eos_state.gam1 * p_l / rho_l);

    eos_state.rho = rho_r;
    eos_state.p = p_r;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn_r[n];
    }
    eos_state.T = castro::T_guess;

    eos(eos_input_rp, eos_state);

    amrex::Real cs_r = std::sqrt(eos_state.gam1 * p_r / rho_r);


    // find the solution as a function of xi = x/t

    // This follows from the discussion around C&G Eq. 15

    // compute xi = x/t -- this is the similarity variable for the
    // solution

    amrex::Real xi = (x - xjump) / time;

    // check which side of the contact we need to worry about

    amrex::Real chi = std::copysign(1.0_rt, xi - ustar);

    amrex::Real rho_s, u_s, p_s, W_s, cs_s;
    amrex::Real uhat_s, xihat, uhat_star;

    if (chi == -1.0_rt) {
        rho_s = rho_l;
        u_s = u_l;
        p_s = p_l;

        W_s = W_l;

        cs_s = cs_l;

        for (int n = 0; n < NumSpec; ++n) {
            xn[n] = xn_l[n];
        }

        uhat_s = chi * u_s;
        xihat = chi * xi;
        uhat_star = chi * ustar;

    } else if (chi == 1.0_rt) {
        rho_s = rho_r;
        u_s = u_r;
        p_s = p_r;

        W_s = W_r;

        cs_s = cs_r;

        for (int n = 0; n < NumSpec; ++n) {
            xn[n] = xn_r[n];
        }

        uhat_s = chi * u_s;
        xihat = chi * xi;
        uhat_star = chi * ustar;

    } else {
        // we should average in this case
        amrex::Error("Not implemented");

    }


    // are we a shock or rarefaction?

    amrex::Real rhostar, Z_temp;

    if (pstar > p_s) {
        // shock

        rhostar = 1.0_rt / (1.0_rt / rho_s - (pstar - p_s)/ (W_s * W_s));

    } else {
        // rarefaction.  Here we need to integrate the Riemann
        // invariant curves to our pstar to get rhostar and cs_star
        amrex::Real Z_temp, W_temp;
        if (chi == -1.0_rt) {
            rarefaction(pstar, rho_s, u_s, p_s, xn, 1, Z_temp, W_temp, rhostar);
        } else {
            rarefaction(pstar, rho_s, u_s, p_s, xn, 3, Z_temp, W_temp, rhostar);
        }

    }

    // get the soundspeed via the EOS (C&G suggest getting it from
    // the jump conditions)

    eos_state.rho = rhostar;
    eos_state.p = pstar;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn[n];
    }
    eos_state.T = castro::T_guess;

    eos(eos_input_rp, eos_state);

    amrex::Real cs_star = std::sqrt(eos_state.gam1 * pstar / rhostar);

    // now deal with the cases where we are not spanning a rarefaction

    amrex::Real lambdahat_s, lambdahat_star;

    if (pstar <= p_s) {
        lambdahat_s = uhat_s + cs_s;
        lambdahat_star = uhat_star + cs_star;
    } else {
        lambdahat_s = uhat_s + W_s / rho_s;
        lambdahat_star = lambdahat_s;
    }

    if (xihat <= lambdahat_star) {
        p = pstar;
        rho = rhostar;
        u = ustar;

    } else if (xihat > lambdahat_s) {
        p = p_s;
        rho = rho_s;
        u = u_s;

    } else {

        // we are inside the rarefaction.  To find the solution here,
        // we need to integrate up to the point where uhat + c = xihat
        // starting from U = U_s

        // for the 1-rarefaction, chi = -1, so this is -u + c = -xi,
        // or u - c = xi, meaning we integrate to u = xi + c

        // for the 3-rarefaction, chi = 1, so this is u + c = xi,
        // so we integrate to u = xi - c

        // Note that c here is c(rho,p) -- we need to compute that
        // self-consistently as we integrate

        if (chi == -1.0_rt) {
            rarefaction_to_u(rho_s, u_s, p_s, xn, 1, xi, rho, p, u);

        } else if (chi == 1.0_rt) {
            rarefaction_to_u(rho_s, u_s, p_s, xn, 3, xi, rho, p, u);

        }
    }

}
#endif
