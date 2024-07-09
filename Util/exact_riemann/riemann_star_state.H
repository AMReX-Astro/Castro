#ifndef RIEMANN_STAR_STATE_H
#define RIEMANN_STAR_STATE_H

#include <AMReX_REAL.H>

#include <riemann_constants.H>
#include <riemann_shock.H>
#include <riemann_rarefaction.H>

using namespace amrex::literals;

AMREX_INLINE
int
riemann_star_state(const amrex::Real rho_l, const amrex::Real u_l, const amrex::Real p_l, const amrex::Real* xn_l,
                   const amrex::Real rho_r, const amrex::Real u_r, const amrex::Real p_r, const amrex::Real* xn_r,
                   amrex::Real& ustar, amrex::Real& pstar, amrex::Real& W_l, amrex::Real& W_r) {

    const amrex::Real tol = 1.e-10_rt;
    const amrex::Real smallp = 1.e-8_rt;
    const amrex::Real SMALL = 1.e-8_rt;


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

    amrex::Real gammaE_l = p_l / (rho_l * eos_state.e) + 1.0_rt;
    amrex::Real gammaC_l = eos_state.gam1;

    eos_state.rho = rho_r;
    eos_state.p = p_r;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn_r[n];
    }
    eos_state.T = castro::T_guess;

    eos(eos_input_rp, eos_state);

    amrex::Real cs_r = std::sqrt(eos_state.gam1 * p_r / rho_r);

    amrex::Real gammaE_r = p_r / (rho_r * eos_state.e) + 1.0_rt;
    amrex::Real gammaC_r = eos_state.gam1;

    amrex::Real gammaE_bar = 0.5_rt * (gammaE_l + gammaE_r);
    amrex::Real gammaC_bar = 0.5_rt * (gammaC_l + gammaC_r);


    // create an initial guess for pstar using a primitive variable
    // Riemann solver

    // We follow the PVRS solver from Toro (Chapter 9)

    // alternative: two shock solver (see Toro 9.42)

    W_l = rho_l * cs_l;
    W_r = rho_r * cs_r;

    // prevent roundoff errors from giving us a pstar that is unphysical
    // if our input states are the same

    if (W_l == W_r) {
        pstar = 0.5_rt * (p_l + p_r + W_l * (u_l - u_r));
    } else {
        pstar = ((W_r * p_l + W_l * p_r) + W_l * W_r * (u_l - u_r)) / (W_l + W_r);
    }

    pstar = amrex::max(pstar, smallp);


    // find the exact pstar and ustar

    // this procedure follows directly from Colella & Glaz 1985, section 1
    // the basic idea is that we want to find the pstar that satisfies:
    //
    //   ustar_l(pstar) - ustar_r(pstar) = 0
    //
    // where ustar_l connects the left state to the star state across
    // the left wave.  this wave can be a shock or rarefaction, so we
    // need to use either the Riemann invariants (rarefaction) or
    // Rankine-Hugoniot jump conditions (shock).
    //
    // We use a Newton method.  It takes the form:
    //
    //  pstar <-- pstar - Z_l Z_r (ustar_r - ustar_l) / (Z_l + Z_r)
    //
    // where Z_s = | dpstar / dustar_s | is the derivative in the p-u plane
    //
    // and we get ustar via:
    //
    //  ustar_s = u_s +/- (pstar - p_s) / W_s
    //
    // where W_s is the wave speed (this also depends on shock or
    // rarefaction)

    amrex::Real ustar_l, ustar_r;
    bool converged = false;

    int iter = 1;
    while (! converged && iter < castro::riemann_shock_maxiter) {

        // compute Z_l, W_l and Z_r, W_r -- the form of these depend
        // on whether the wave is a shock or a rarefaction

        amrex::Real Z_l, Z_r;

        // left wave

        if (pstar - p_l > SMALL * p_l) {
            // left shock
            shock(pstar, rho_l, u_l, p_l, xn_l, gammaE_bar, gammaC_bar, Z_l, W_l);
        } else {
            // left rarefaction
            amrex::Real rhostar_dummy;
            rarefaction(pstar, rho_l, u_l, p_l, xn_l, 1, Z_l, W_l, rhostar_dummy);
        }

        // right wave

        if (pstar - p_r > SMALL * p_r) {
            // right shock
            shock(pstar, rho_r, u_r, p_r, xn_r, gammaE_bar, gammaC_bar, Z_r, W_r);
        } else {
            // right rarefaction
            amrex::Real rhostar_dummy;
            rarefaction(pstar, rho_r, u_r, p_r, xn_r, 3, Z_r, W_r, rhostar_dummy);
        }

        ustar_l = u_l - (pstar - p_l) / W_l;
        ustar_r = u_r + (pstar - p_r) / W_r;

        amrex::Real pstar_new = pstar - Z_l * Z_r * (ustar_r - ustar_l) / (Z_l + Z_r);

        // estimate the error in the current pstar solution

        amrex::Real error = std::abs(pstar_new - pstar);

        if (error < tol * pstar) {
            converged = true;
        }

        // get ready for the next iteration
        pstar = std::clamp(pstar_new, 0.5 * pstar, 2.0 * pstar);

        iter++;
    }

    if (! converged) {
        amrex::Error("star state did not converge");
    }

    ustar = 0.5_rt * (ustar_l + ustar_r);

    return iter;

}
#endif
