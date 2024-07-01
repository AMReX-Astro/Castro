#ifndef RIEMANN_SHOCK_H
#define RIEMANN_SHOCK_H

#include <AMReX_REAL.H>

#include <eos_type.H>
#include <eos.H>

#include <riemann_constants.H>

AMREX_INLINE
void
W_s_shock(const amrex::Real W_s, const amrex::Real pstar,
          const amrex::Real rho_s, const amrex::Real p_s, const amrex::Real e_s,
          const amrex::Real* xn,
          amrex::Real& rhostar_s, eos_rep_t& eos_state, amrex::Real& f, amrex::Real& fprime) {

    // we need rhostar -- get it from the R-H conditions

    amrex::Real taustar_s = (1.0_rt / rho_s) - (pstar - p_s) / (W_s * W_s);
    rhostar_s = 1.0_rt / taustar_s;

    // get the thermodynamics

    eos_state.rho = rhostar_s;
    eos_state.p = pstar;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn[n];
    }
    eos_state.T = castro::T_guess;

    eos(eos_input_rp, eos_state);

    // compute the correction

    f = W_s * W_s * (eos_state.e - e_s) - 0.5_rt * (pstar * pstar - p_s * p_s);

    // we need de/drho at constant p -- this is not returned by the EOS

    amrex::Real dedrho_p = eos_state.dedr - eos_state.dedT * eos_state.dpdr / eos_state.dpdT;

    fprime = 2.0_rt * W_s * (eos_state.e - e_s) -
        2.0_rt * dedrho_p * (pstar - p_s) * std::pow(rhostar_s, 2) / W_s;

}


AMREX_INLINE
void
newton_shock(amrex::Real& W_s, const amrex::Real pstar,
             const amrex::Real rho_s, const amrex::Real p_s, const amrex::Real e_s,
             const amrex::Real* xn,
             amrex::Real* rhostar_hist, amrex::Real* Ws_hist,
             eos_rep_t& eos_state, bool& converged) {

    // Newton iterations -- we are zeroing the energy R-H jump condition
    // W^2 [e] = 1/2 [p^2]
    //
    // we write f(W) = W^2 (e(pstar, rhostar_s) - e_s) - (1/2)(pstar^2 - p_s)
    //
    // and then compute f'

    converged = false;

    int iter = 1;
    while (! converged && iter < castro::riemann_shock_maxiter) {

        amrex::Real rhostar_s, f, fprime;

        W_s_shock(W_s, pstar, rho_s, p_s, e_s, xn,
                  rhostar_s, eos_state, f, fprime);

        amrex::Real dW = -f / fprime;

        if (std::abs(dW) < castro::riemann_pstar_tol * W_s) {
            converged = true;
        }

        W_s = amrex::min(2.0_rt * W_s, amrex::max(0.5_rt * W_s, W_s + dW));

        // store some history

        rhostar_hist[iter] = rhostar_s;
        Ws_hist[iter] = W_s;

        iter++;
    }
}

AMREX_INLINE
void
shock(const amrex::Real pstar,
      const amrex::Real rho_s, const amrex::Real u_s, const amrex::Real p_s,
      const amrex::Real* xn,
      const amrex::Real gammaE_bar, const amrex::Real gammaC_bar,
      amrex::Real& Z_s, amrex::Real& W_s) {

    amrex::ignore_unused(u_s);

    amrex::Real rhostar_hist[riemann_constants::HISTORY_SIZE], Ws_hist[riemann_constants::HISTORY_SIZE];

    // compute the Z_s function for a shock following C&G Eq. 20 and
    // 23.  Here the "_s" variables are the state (L or R) that we are
    // connecting to the star region through a shock.

    // first we need to compute W_s -- this is iterative because of
    // the nonlinear EOS.  We use the R-H jump conditions + the EOS

    // get the s-state energy, e_s

    eos_rep_t eos_state;
    eos_state.rho = rho_s;
    eos_state.p = p_s;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn[n];
    }
    eos_state.T = castro::T_guess;

    eos(eos_input_rp, eos_state);

    amrex::Real e_s = eos_state.e;

    // to kick things off, we need a guess for W_s.  We'll use the
    // approximation from Colella & Glaz (Eq. 34), which in turn
    // makes an approximation for gammaE_star, using equation 31.

    amrex::Real gammaE_s = p_s / (rho_s * e_s) + 1.0_rt;

    amrex::Real gammaE_star = gammaE_s +
        2.0_rt * (1.0_rt - gammaE_bar / gammaC_bar) * (gammaE_bar - 1.0_rt) *
        (pstar - p_s) / (pstar + p_s);

    // there is a pathological case that if p_s - pstar ~ 0, the root finding
    // just doesn't work.  In this case, we use the ideas from CG, Eq. 35, and
    // take W = sqrt(gamma p rho)

    if (pstar - p_s < riemann_solver::tol_p * p_s) {
        W_s = std::sqrt(eos_state.gam1 * p_s * rho_s);
    } else {
       W_s = std::sqrt((pstar - p_s) *
                       (pstar + 0.5_rt * (gammaE_star - 1.0_rt) * (pstar + p_s)) /
                       (pstar / rho_s - (gammaE_star - 1.0_rt) / (gammaE_s - 1.0_rt) * p_s / rho_s));
    }

    // we need rhostar -- get it from the R-H conditions

    amrex::Real taustar_s = (1.0_rt / rho_s) - (pstar - p_s) / (W_s * W_s);
    amrex::Real rhostar_s;

    if (taustar_s < 0.0_rt) {
        rhostar_s = riemann_constants::small;
        W_s = std::sqrt((pstar - p_s) / (1.0_rt / rho_s - 1.0_rt / rhostar_s));
    }

    // newton

    bool converged;
    newton_shock(W_s, pstar, rho_s, p_s, e_s, xn,
                 rhostar_hist, Ws_hist,
                 eos_state, converged);


    // now did we converge?

    if (! converged) {
        for (int i = 0; i < castro::riemann_shock_maxiter; ++i) {
            std::cout << i << " " << rhostar_hist[i] << " " << Ws_hist[i] << std::endl;
        }

        amrex::Error("shock did not converge");
    }


    // now that we have W_s, we can get rhostar from the R-H conditions
    // (C&G Eq. 12)

    taustar_s = (1.0_rt / rho_s) - (pstar - p_s) / (W_s * W_s);
    rhostar_s = 1.0_rt / taustar_s;

    // next we compute the derivative dW_s/dpstar -- the paper gives
    // dW**2/dpstar (Eq. 23), so we take 1/2W of that

    amrex::Real C = std::sqrt(eos_state.gam1 * pstar * rhostar_s);

    amrex::Real p_e = eos_state.dpdT / eos_state.dedT;
    amrex::Real p_rho = eos_state.dpdr - eos_state.dpdT * eos_state.dedr / eos_state.dedT;

    amrex::Real p_tau = -std::pow(rhostar_s, 2) * p_rho;

    amrex::Real dW2dpstar = (C*C - W_s*W_s) * W_s * W_s /
        ((0.5_rt * (pstar + p_s) * p_e - p_tau) * (pstar - p_s));

    amrex::Real dWdpstar = 0.5_rt * dW2dpstar / W_s;

    // finally compute Z_s

    Z_s = W_s * W_s / (W_s - dWdpstar * (pstar - p_s));

}

#endif
