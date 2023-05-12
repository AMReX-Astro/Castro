#ifndef RIEMANN_SUPPORT_H
#define RIEMANN_SUPPORT_H

#include <eos_type.H>
#include <eos.H>

const Real smallrho = 1.e-5_rt;
const int max_iters = 100;
const Real tol = 1.e-6_rt;


AMREX_INLINE
void
W_s_shock(const Real W_s, const Real pstar, const Real rho_s, const Real p_s, const Real e_s, const Real* xn,
          Real& rhostar_s, eos_t& eos_state, Real& f, Real& fprime) {

    // we need rhostar -- get it from the R-H conditions

    Real taustar_s = (1.0_rt / rho_s) - (pstar - p_s) / (W_s * W_s);
    rhostar_s = 1.0_rt / taustar_s;

    // get the thermodynamics

    eos_state.rho = rhostar_s;
    eos_state.p = pstar;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn[n];
    }
    eos_state.T = problem::initial_temp_guess;

    eos(eos_input_rp, eos_state);

    // compute the correction

    f = W_s * W_s * (eos_state.e - e_s) - 0.5_rt * (pstar * pstar - p_s * p_s);

    // we need de/drho at constant p -- this is not returned by the EOS

    Real dedrho_p = eos_state.dedr - eos_state.dedT * eos_state.dpdr / eos_state.dpdT;

    fprime = 2.0_rt * W_s * (eos_state.e - e_s) - 2.0_rt * dedrho_p * (pstar - p_s) * std::pow(rhostar_s, 2) / W_s;

}


AMREX_INLINE
void
newton_shock(Real& W_s, const Real pstar, const Real rho_s, const Real p_s, const Real e_s, const Real* xn,
             Real* rhostar_hist, Real* Ws_hist,
             eos_t& eos_state, bool& converged) {

    // Newton iterations -- we are zeroing the energy R-H jump condition
    // W^2 [e] = 1/2 [p^2]
    //
    // we write f(W) = W^2 (e(pstar, rhostar_s) - e_s) - (1/2)(pstar^2 - p_s)
    //
    // and then compute f'

    converged = false;

    int iter = 1;
    while (! converged && iter < max_iters) {

        Real rhostar_s, f, fprime;

        W_s_shock(W_s, pstar, rho_s, p_s, e_s, xn,
                  rhostar_s, eos_state, f, fprime);

        Real dW = -f / fprime;

        if (std::abs(dW) < tol * W_s) {
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
shock(const Real pstar, const Real rho_s, const Real u_s, const Real p_s, const Real* xn,
      const Real gammaE_bar, const Real gammaC_bar,
      Real& Z_s, Real& W_s) {


    const Real tol_p = 1.e-6_rt;

    Real rhostar_hist[max_iters], Ws_hist[max_iters];

    // compute the Z_s function for a shock following C&G Eq. 20 and
    // 23.  Here the "_s" variables are the state (L or R) that we are
    // connecting to the star region through a shock.

    // first we need to compute W_s -- this is iterative because of
    // the nonlinear EOS.  We use the R-H jump conditions + the EOS

    // get the s-state energy, e_s

    eos_t eos_state;
    eos_state.rho = rho_s;
    eos_state.p = p_s;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn[n];
    }
    eos_state.T = problem::initial_temp_guess;

    eos(eos_input_rp, eos_state);

    Real e_s = eos_state.e;

    // to kick things off, we need a guess for W_s.  We'll use the
    // approximation from Colella & Glaz (Eq. 34), which in turn
    // makes an approximation for gammaE_star, using equation 31.

    Real gammaE_s = p_s / (rho_s * e_s) + 1.0_rt;

    Real gammaE_star = gammaE_s +
        2.0_rt * (1.0_rt - gammaE_bar / gammaC_bar) * (gammaE_bar - 1.0_rt) *
        (pstar - p_s) / (pstar + p_s);

    std::cout << "pstar, ps = " << pstar << " " << p_s << " " << gammaE_s << " " << gammaE_star << std::endl;
    std::cout << (pstar/rho_s - (gammaE_star - 1.0_rt)/(gammaE_s - 1.0_rt) * p_s/rho_s);
    std::cout << (pstar + 0.5_rt * (gammaE_star - 1.0_rt) * (pstar + p_s));

    // there is a pathological case that if p_s - pstar ~ 0, the root finding
    // just doesn't work.  In this case, we use the ideas from CG, Eq. 35, and
    // take W = sqrt(gamma p rho)

    if (pstar - p_s < tol_p * p_s) {
        W_s = std::sqrt(eos_state.gam1 * p_s * rho_s);
    } else {
       W_s = std::sqrt((pstar - p_s) *
                       (pstar + 0.5_rt * (gammaE_star - 1.0_rt) * (pstar + p_s)) /
                       (pstar / rho_s - (gammaE_star - 1.0_rt) / (gammaE_s - 1.0_rt) * p_s / rho_s));
    }

    // we need rhostar -- get it from the R-H conditions

    Real taustar_s = (1.0_rt / rho_s) - (pstar - p_s) / (W_s * W_s);
    Real rhostar_s;

    if (taustar_s < 0.0_rt) {
        rhostar_s = smallrho;
        W_s = std::sqrt((pstar - p_s) / (1.0_rt / rho_s - 1.0_rt / rhostar_s));
    }

    Real W_s_guess = W_s;

    // newton

    bool converged;
    newton_shock(W_s, pstar, rho_s, p_s, e_s, xn,
                 rhostar_hist, Ws_hist,
                 eos_state, converged);


    // now did we converge?

    if (! converged) {
        for (int i = 0; i < max_iters; ++i) {
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

    Real C = std::sqrt(eos_state.gam1 * pstar * rhostar_s);

    Real p_e = eos_state.dpdT / eos_state.dedT;
    Real p_rho = eos_state.dpdr - eos_state.dpdT * eos_state.dedr / eos_state.dedT;

    Real p_tau = -std::pow(rhostar_s, 2) * p_rho;

    Real dW2dpstar = (C*C - W_s*W_s) * W_s * W_s /
        ((0.5_rt * (pstar + p_s) * p_e - p_tau) * (pstar - p_s));

    Real dWdpstar = 0.5_rt * dW2dpstar / W_s;

    // finally compute Z_s

    Z_s = W_s * W_s / (W_s - dWdpstar * (pstar - p_s));

}


AMREX_INLINE
void
riemann_invariant_rhs(const Real p, const Real tau, const Real u, const Real* xn, const int iwave,
                      Real& dtaudp, Real& dudp) {

    // here, p is out independent variable, and tau, u are the
    // dependent variables.  We return the derivatives of these
    // wrt p for integration.

    // get the thermodynamics

    eos_t eos_state;
    eos_state.rho = 1.0_rt / tau;
    eos_state.p = p;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn[n];
    }
    eos_state.T = problem::initial_temp_guess;

    eos(eos_input_rp, eos_state);

    Real C = std::sqrt(eos_state.gam1 * p / tau);

    dtaudp = -1.0_rt / (C * C);

    if (iwave == 1) {
        dudp = -1.0_rt / C;
    } else if (iwave == 3) {
        dudp = 1.0_rt / C;
    }

}


AMREX_INLINE
void
riemann_invariant_rhs2(const Real u, const Real tau, const Real p, const Real* xn, const int iwave,
                       Real& dtaudu, Real& dpdu) {

    // here, u is out independent variable, and tau, p are the
    // dependent variables.  We return the derivatives of these
    // wrt u for integration.

    // get the thermodynamics

    eos_t eos_state;
    eos_state.rho = 1.0_rt / tau;
    eos_state.p = p;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn[n];
    }
    eos_state.T = problem::initial_temp_guess;

    eos(eos_input_rp, eos_state);

    Real C = std::sqrt(eos_state.gam1 * p / tau);

    if (iwave == 3) {
       dpdu = C;
       dtaudu = -1.0_rt / C;

    } else if (iwave == 1) {
        dpdu = -C;
        dtaudu = 1.0_rt / C;
    }
}


AMREX_INLINE
void
rarefaction(const Real pstar, const Real rho_s, const Real u_s, const Real p_s,
            const Real* xn, const int iwave, Real& Z_s, Real& W_s, Real& rhostar) {

    const int npts = 1000;

    // Compute Z_s = C for a rarefaction connecting the state to the star
    // region by integrating the Riemann invariant from p_s to pstar.
    // This means solving a system of ODEs.  We use 4th-order R-K.

    Real tau = 1.0_rt / rho_s;
    Real u = u_s;
    Real p = p_s;

    Real dp = (pstar - p_s) / static_cast<Real>(npts);
    Real dp2 = 0.5_rt * dp;

    for (int i = 1; i <= npts; ++i) {

        // do 4th-order RT

        Real dtaudp1, dudp1;
        riemann_invariant_rhs(p, tau, u, xn, iwave, dtaudp1, dudp1);

        Real dtaudp2, dudp2;
        riemann_invariant_rhs(p+dp2, tau+dp2*dtaudp1, u+dp2*dudp1, xn, iwave, dtaudp2, dudp2);

        Real dtaudp3, dudp3;
        riemann_invariant_rhs(p+dp2, tau+dp2*dtaudp2, u+dp2*dudp2, xn, iwave, dtaudp3, dudp3);

        Real dtaudp4, dudp4;
        riemann_invariant_rhs(p+dp, tau+dp*dtaudp3, u+dp*dudp3, xn, iwave, dtaudp4, dudp4);

        p += dp;
        u += (1.0_rt/6.0_rt) * dp * (dudp1 + 2.0_rt * dudp2 + 2.0_rt * dudp3 + dudp4);
        tau += (1.0_rt/6.0_rt) * dp * (dtaudp1 + 2.0_rt * dtaudp2 + 2.0_rt * dtaudp3 + dtaudp4);

    }

    // Z_s is just the Lagrangian sound speed

    eos_t eos_state;
    eos_state.rho = 1.0_rt / tau;
    eos_state.p = p;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn[n];
    }
    eos_state.T = problem::initial_temp_guess;

    eos(eos_input_rp, eos_state);

    Z_s = std::sqrt(eos_state.gam1 * p / tau);

    // also need W_s -- this is C&G Eq. 16.  u above is ustar_s.

    if (u == u_s) {
        W_s = Z_s;
    } else {
        W_s = std::abs(pstar - p_s) / std::abs(u - u_s);
    }

    rhostar = 1.0_rt / tau;
}


AMREX_INLINE
void
rarefaction_to_u(const Real rho_s, const Real u_s, const Real p_s, const Real* xn, const int iwave, const Real xi,
                 Real& rho, Real& p, Real& u) {

    const int npts = 1000;

    // here we integrate the Riemann invariants for a rarefaction up to
    // some intermediate u (between u_s and ustar).  This accounts for
    // the fact that we are inside the rarefaction.

    // We reformulate the system of ODEs from C&G Eq. 13 to make u the
    // dependent variable.  Now we solve:

    // we actually don't know the stopping point.  For the 1-wave, we
    // stop at u = xi + c, for the 3-wave, we stop at u = xi - c, where
    // c is computed as we step.

    // dp/du =  C; dtau/du = -1/C   for the 1-wave
    // dp/du = -C; dtau/du =  1/C   for the 3-wave

    Real tau = 1.0_rt / rho_s;
    u = u_s;
    p = p_s;

    // estimate
    // compute c

    eos_t eos_state;
    eos_state.rho = 1.0_rt / tau;
    eos_state.p = p;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn[n];
    }
    eos_state.T = problem::initial_temp_guess;

    eos(eos_input_rp, eos_state);

    Real c = std::sqrt(eos_state.gam1 * p * tau);

    Real ustop;
    if (iwave == 1) {
        ustop = xi + c;
    } else if (iwave == 3) {
        ustop = xi - c;
    }

    Real du = (ustop - u_s) / static_cast<Real>(npts);

    bool finished = false;

    std::cout << "integrating from u: " << u << " " << ustop << " " << xi << " " << c << std::endl;

    Real du2 = 0.5_rt * du;

    while (! finished) {

        // do 4th-order RT

        Real dtaudu1, dpdu1;
        riemann_invariant_rhs2(u, tau, p, xn, iwave, dtaudu1, dpdu1);

        Real dtaudu2, dpdu2;
        riemann_invariant_rhs2(u+du2, tau+du2*dtaudu1, p+du2*dpdu1, xn, iwave, dtaudu2, dpdu2);

        Real dtaudu3, dpdu3;
        riemann_invariant_rhs2(u+du2, tau+du2*dtaudu2, p+du2*dpdu2, xn, iwave, dtaudu3, dpdu3);

        Real dtaudu4, dpdu4;
        riemann_invariant_rhs2(u+du, tau+du*dtaudu3, p+du*dpdu3, xn, iwave, dtaudu4, dpdu4);

        u += du;
        p += (1.0_rt/6.0_rt) * du * (dpdu1 + 2.0_rt * dpdu2 + 2.0_rt * dpdu3 + dpdu4);
        tau += (1.0_rt/6.0_rt) * du * (dtaudu1 + 2.0_rt * dtaudu2 + 2.0_rt * dtaudu3 + dtaudu4);

        // compute c

        eos_t eos_state;
        eos_state.rho = 1.0_rt/tau;
        eos_state.p = p;
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = xn[n];
        }
        eos_state.T = problem::initial_temp_guess;

        eos(eos_input_rp, eos_state);

        c = std::sqrt(eos_state.gam1 * p * tau);

        // check the step size

        if (iwave == 1) {
            ustop = xi + c;
        } else if (iwave == 3) {
            ustop = xi - c;
        }

        if (du * u > 0.0_rt) {
            while (std::abs(u + du) > std::abs(ustop) && du != 0.0_rt) {
                du = 0.5_rt * du;
            }
        } else {
            if (u > 0.0_rt) {
                while (u + du < ustop && du != 0.0_rt) {
                    du = 0.5_rt * du;
                }

            } else {
                while (u + du > ustop && du != 0.0_rt) {
                    du = 0.5_rt * du;
                }
            }
        }

        if (std::abs(du) < tol * std::abs(u)) {
            finished = true;
        }

    }

    rho = 1.0_rt / tau;
}

#endif
