#ifndef RIEMANN_STAR_STATE_H
#define RIEMANN_STAR_STATE_H

#include <AMReX_REAL.H>

#include <riemann_support.H>

using namespace amrex::literals;

AMREX_INLINE
void
riemann_star_state(const amrex::Real rho_l, const amrex::Real u_l, const amrex::Real p_l, const amrex::Real* xn_l,
                   const amrex::Real rho_r, const amrex::Real u_r, const amrex::Real p_r, const amrex::Real* xn_r,
                   amrex::Real& ustar, amrex::Real& pstar, amrex::Real& W_l, amrex::Real& W_r) {

    const amrex::Real tol = 1.e-10_rt;
    const amrex::Real smallp = 1.e-8_rt;
    const amrex::Real SMALL = 1.e-13_rt;


    // get the initial sound speeds

    eos_t eos_state;
    eos_state.rho = rho_l;
    eos_state.p = p_l;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn_l[n];
    }
    eos_state.T = problem::initial_temp_guess;

    eos(eos_input_rp, eos_state);

    amrex::Real cs_l = std::sqrt(eos_state.gam1 * p_l / rho_l);

    amrex::Real gammaE_l = p_l / (rho_l * eos_state.e) + 1.0_rt;
    amrex::Real gammaC_l = eos_state.gam1;

    eos_state.rho = rho_r;
    eos_state.p = p_r;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn_r[n];
    }
    eos_state.T = problem::initial_temp_guess;

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
        pstar = ((W_r*p_l + W_l*p_r) + W_l*W_r*(u_l - u_r))/(W_l + W_r);
    }

    pstar = amrex::max(pstar, smallp);


    // find the exact pstar and ustar

    std::cout << "solving for star state: " << rho_l << " " << u_l << " " << p_l << " " << rho_r << " " << u_r << " " << p_r << std::endl;


    // this procedure follows directly from Colella & Glaz 1985, section 1

    amrex::Real ustar_l, ustar_r;
    bool converged = false;

    int iter = 1;
    while (! converged && iter < problem::riemann_max_iter) {

        // compute Z_l and Z_r -- the form of these depend on whether the
        // wave is a shock or a rarefaction

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

        std::cout << "done with iteration " << iter << std::endl;
        std::cout << "ustar_l/r, pstar: " << ustar_l << " " << ustar_r << " " << pstar_new << std::endl;

        // estimate the error in the current star solution
        amrex::Real err1 = std::abs(ustar_r - ustar_l);
        amrex::Real err2 = pstar_new - pstar;

        if (err1 < tol * amrex::max(std::abs(ustar_l), std::abs(ustar_r)) &&
            err2 < tol * pstar) {
            converged = true;
        }

        // get ready for the next iteration
        pstar = pstar_new;

        iter++;
    }

    ustar = 0.5_rt * (ustar_l + ustar_r);

    std::cout << "found pstar, ustar: " <<  pstar << " " << ustar << std::endl;
}
#endif
