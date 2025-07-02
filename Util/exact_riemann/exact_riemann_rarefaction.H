#ifndef RIEMANN_RAREFACTION_H
#define RIEMANN_RAREFACTION_H

#include <utility>

#include <AMReX_REAL.H>

#include <eos_type.H>
#include <eos.H>


AMREX_INLINE
void
riemann_invariant_rhs(const amrex::Real p, const amrex::Real tau, const amrex::Real u,
                      const amrex::Real* xn, const int iwave,
                      amrex::Real& T,
                      amrex::Real& dtaudp, amrex::Real& dudp) {

    amrex::ignore_unused(u);

    // here, p is out independent variable, and tau, u are the
    // dependent variables.  We return the derivatives of these
    // wrt p for integration.

    // T should be an initial guess coming in

    // get the thermodynamics

    eos_rep_t eos_state;
    eos_state.rho = 1.0_rt / tau;
    eos_state.p = p;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn[n];
    }
    eos_state.T = T;

    eos(eos_input_rp, eos_state);

    T = eos_state.T;

    amrex::Real C = std::sqrt(eos_state.gam1 * p / tau);

    dtaudp = -1.0_rt / (C * C);

    if (iwave == 1) {
        dudp = -1.0_rt / C;
    } else if (iwave == 3) {
        dudp = 1.0_rt / C;
    }

}


AMREX_INLINE
void
riemann_invariant_rhs2(const amrex::Real u, const amrex::Real tau, const amrex::Real p,
                       const amrex::Real* xn, const int iwave,
                       amrex::Real& T,
                       amrex::Real& dtaudu, amrex::Real& dpdu) {

    // here, u is out independent variable, and tau, p are the
    // dependent variables.  We return the derivatives of these
    // wrt u for integration.

    // here T is an initial guess

    // get the thermodynamics

    amrex::ignore_unused(u);

    eos_rep_t eos_state;
    eos_state.rho = 1.0_rt / tau;
    eos_state.p = p;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn[n];
    }
    eos_state.T = T;

    eos(eos_input_rp, eos_state);

    T = eos_state.T;

    amrex::Real C = std::sqrt(eos_state.gam1 * p / tau);

    if (iwave == 3) {
       dpdu = C;
       dtaudu = -1.0_rt / C;

    } else if (iwave == 1) {
        dpdu = -C;
        dtaudu = 1.0_rt / C;
    }
}


AMREX_INLINE
std::pair<amrex::Real, amrex::Real>
single_step_p(const amrex::Real pstar0, const amrex::Real dp,
              const amrex::Real u0, const amrex::Real tau0,
              const amrex::Real* xn, const int iwave, amrex::Real& T) {

    // this takes a single step of the Riemann invariant system where
    // p is the independent variable

    amrex::Real dtaudp1, dudp1;
    riemann_invariant_rhs(pstar0, tau0, u0, xn, iwave, T, dtaudp1, dudp1);

    amrex::Real dtaudp2, dudp2;
    riemann_invariant_rhs(pstar0+0.5*dp, tau0+0.5*dp*dtaudp1, u0+0.5*dp*dudp1,
                          xn, iwave, T, dtaudp2, dudp2);

    amrex::Real dtaudp3, dudp3;
    riemann_invariant_rhs(pstar0+0.5*dp, tau0+0.5*dp*dtaudp2, u0+0.5*dp*dudp2,
                          xn, iwave, T, dtaudp3, dudp3);

    amrex::Real dtaudp4, dudp4;
    riemann_invariant_rhs(pstar0+dp, tau0+dp*dtaudp3, u0+dp*dudp3,
                          xn, iwave, T, dtaudp4, dudp4);

    amrex::Real u = u0 + (1.0_rt/6.0_rt) * dp * (dudp1 + 2.0_rt * dudp2 + 2.0_rt * dudp3 + dudp4);
    amrex::Real tau = tau0 + (1.0_rt/6.0_rt) * dp * (dtaudp1 + 2.0_rt * dtaudp2 + 2.0_rt * dtaudp3 + dtaudp4);

    return {u, tau};
}


AMREX_INLINE
std::pair<amrex::Real, amrex::Real>
single_step_u(const amrex::Real u0, const amrex::Real du,
              const amrex::Real tau0, const amrex::Real p0,
              const amrex::Real* xn, const int iwave, amrex::Real& T) {

    // this takes a single step of the Riemann invariant system where
    // u is the independent variable

    amrex::Real dtaudu1, dpdu1;
    riemann_invariant_rhs2(u0, tau0, p0, xn, iwave, T, dtaudu1, dpdu1);

    amrex::Real dtaudu2, dpdu2;
    riemann_invariant_rhs2(u0+0.5*du, tau0+0.5*du*dtaudu1, p0+0.5*du*dpdu1,
                           xn, iwave, T, dtaudu2, dpdu2);

    amrex::Real dtaudu3, dpdu3;
    riemann_invariant_rhs2(u0+0.5*du, tau0+0.5*du*dtaudu2, p0+0.5*du*dpdu2,
                           xn, iwave, T, dtaudu3, dpdu3);

    amrex::Real dtaudu4, dpdu4;
    riemann_invariant_rhs2(u0+du, tau0+du*dtaudu3, p0+du*dpdu3,
                           xn, iwave, T, dtaudu4, dpdu4);

    amrex::Real p = p0 + (1.0_rt/6.0_rt) * du * (dpdu1 + 2.0_rt * dpdu2 + 2.0_rt * dpdu3 + dpdu4);
    amrex::Real tau = tau0 + (1.0_rt/6.0_rt) * du * (dtaudu1 + 2.0_rt * dtaudu2 + 2.0_rt * dtaudu3 + dtaudu4);

    return {tau, p};

}

AMREX_INLINE
void
rarefaction(const amrex::Real pstar,
            const amrex::Real rho_s, const amrex::Real u_s, const amrex::Real p_s,
            const amrex::Real* xn, const int iwave,
            amrex::Real& Z_s, amrex::Real& W_s, amrex::Real& rhostar) {

    // Compute Z_s = C(p*, rho*) for a rarefaction connecting the
    // state to the star region by integrating the Riemann invariant
    // from p_s to pstar.
    //
    // The system is:
    //
    //    dtau/dp = -1/C**2
    //    du/dp = +/- 1/C    (+ for 1-wave, - for 3-wave)

    const double S1{0.9};
    const double S2{4.0};

    // initial conditions

    amrex::Real tau = 1.0_rt / rho_s;
    amrex::Real u = u_s;
    amrex::Real p = p_s;

    // initial guess at integration step
    amrex::Real dp = (pstar - p_s) / static_cast<amrex::Real>(100);

    amrex::Real T = castro::T_guess;

    // adaptive RK4 loop

    // note: we should only be in here if p > pstar, since otherwise we'd have a shock

    amrex::Real dp_new{dp};
    int nstep{0};

    while (p > pstar) {

        amrex::Real u0{u};
        amrex::Real tau0{tau};

        amrex::Real u_new;
        amrex::Real tau_new;

        amrex::Real rel_error = std::numeric_limits<double>::max();

        // take a step

        while (rel_error > riemann_constants::riemann_integral_tol) {
            dp = dp_new;
            if (p + dp < pstar) {
                dp = pstar - p;
            }

            // take 2 half steps
            amrex::Real u_tmp;
            amrex::Real tau_tmp;
            std::tie(u_tmp, tau_tmp) = single_step_p(p, 0.5*dp, u0, tau0, xn, iwave, T);
            std::tie(u_new, tau_new) = single_step_p(p + 0.5*dp, 0.5*dp, u_tmp, tau_tmp, xn, iwave, T);

            // now take a single step to cover dp
            auto [u_single, tau_single] = single_step_p(p, dp, u0, tau0, xn, iwave, T);

            // estimate the relative error
            amrex::Real u_err = std::abs(u_new - u_single);
            if (u0 != 0.0) {
                u_err /= u0;
            }
            amrex::Real tau_err = std::abs((tau_new - tau_single) / tau0);

            rel_error = std::max(u_err, tau_err);

            // adaptive step estimate

            double dp_est = S1 * dp * std::pow(std::abs(riemann_constants::riemann_integral_tol/rel_error), 0.2);
            dp_new = dp_est; //std::clamp(S1*dp_est, dp/S2, S2*dp);

        }

        // success
        p += dp;
        u = u_new;
        tau = tau_new;

        nstep++;
    }

    // Z_s is just the Lagrangian sound speed

    eos_rep_t eos_state;
    eos_state.rho = 1.0_rt / tau;
    eos_state.p = p;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn[n];
    }
    eos_state.T = castro::T_guess;

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
rarefaction_to_u(const amrex::Real rho_s, const amrex::Real u_s, const amrex::Real p_s,
                 const amrex::Real* xn, const int iwave, const amrex::Real xi,
                 amrex::Real& rho, amrex::Real& p, amrex::Real& u) {

    const int npts = 1000;

    // here we integrate the Riemann invariants for a rarefaction up to
    // some intermediate u (between u_s and ustar).  This accounts for
    // the fact that we are inside the rarefaction.
    //
    // We reformulate the system of ODEs from C&G Eq. 13 to make u the
    // dependent variable.  Now we solve:
    //
    //   dp/du =  C; dtau/du = -1/C   for the 1-wave
    //   dp/du = -C; dtau/du =  1/C   for the 3-wave
    //
    // we actually don't know the stopping point.  For the 1-wave, we
    // stop at u = xi + c, for the 3-wave, we stop at u = xi - c, where
    // c is computed as we step.

    const double S1{0.9};

    // initial conditions

    amrex::Real tau = 1.0_rt / rho_s;
    u = u_s;
    p = p_s;

    // compute c and estimate the velocity for the endpoint of integration

    eos_rep_t eos_state;
    eos_state.rho = 1.0_rt / tau;
    eos_state.p = p;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn[n];
    }
    eos_state.T = castro::T_guess;

    eos(eos_input_rp, eos_state);

    amrex::Real c = std::sqrt(eos_state.gam1 * p * tau);

    amrex::Real ustop = (iwave == 1) ? xi + c : xi - c;

    // initial guess at the step size -- we'll adapt below

    amrex::Real du = (ustop - u_s) / static_cast<amrex::Real>(10);

    bool finished = false;

    // this will be used as an initial guess to accelerate the EOS inversions

    amrex::Real T = eos_state.T;

    amrex::Real du_new{du};

    while (! finished) {

        amrex::Real p0{p};
        amrex::Real tau0{tau};

        amrex::Real p_new;
        amrex::Real tau_new;

        amrex::Real rel_error = std::numeric_limits<double>::max();

        // take a step

        while (rel_error > riemann_constants::riemann_integral_tol) {
            du = du_new;

            // take 2 half steps
            amrex::Real p_tmp;
            amrex::Real tau_tmp;
            std::tie(tau_tmp, p_tmp) = single_step_u(u, 0.5*du, tau0, p0, xn, iwave, T);
            std::tie(tau_new, p_new) = single_step_u(u + 0.5*du, 0.5*du, tau_tmp, p_tmp, xn, iwave, T);

            // now take a single step to cover du
            auto [tau_single, p_single] = single_step_u(u, du, tau0, p0, xn, iwave, T);

            // estimate the relative error
            amrex::Real p_err = std::abs((p_new - p_single) / p0);
            amrex::Real tau_err = std::abs((tau_new - tau_single) / tau0);

            rel_error = std::max(p_err, tau_err);

            // adaptive step estimate

            // if our original step was too small, the error might be zero, in which case,
            // just bump up the step for next time

            amrex::Real du_est{};
            if (rel_error > 0) {
                du_est = S1 * du * std::pow(std::abs(riemann_constants::riemann_integral_tol / rel_error), 0.2);
            } else {
                du_est = 10.0 * du;
            }
            du_new = du_est;

        }

        // success
        u += du;
        p = p_new;
        tau = tau_new;

        // compute c and re-estimate the endpoint of integration

        eos_rep_t eos_state;
        eos_state.rho = 1.0_rt/tau;
        eos_state.p = p;
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = xn[n];
        }
        eos_state.T = castro::T_guess;

        eos(eos_input_rp, eos_state);

        c = std::sqrt(eos_state.gam1 * p * tau);

        ustop = (iwave == 1) ? xi + c : xi - c;

        // check the step size for the next step to ensure we don't
        // go past the stopping velocity, ustop

        if (du_new * u > 0.0_rt) {
            while (std::abs(u + du_new) > std::abs(ustop) && du_new != 0.0_rt) {
                du_new *= 0.5_rt;
            }
        } else {
            if (u > 0.0_rt) {
                while (u + du_new < ustop && du_new != 0.0_rt) {
                    du_new *= 0.5_rt;
                }

            } else {
                while (u + du_new > ustop && du_new != 0.0_rt) {
                    du_new *= 0.5_rt;
                }
            }
        }

        if (std::abs(du_new) < riemann_constants::riemann_u_tol * std::abs(u) + riemann_constants::riemann_u_tol) {
            finished = true;
        }

    }

    rho = 1.0_rt / tau;
}

#endif
