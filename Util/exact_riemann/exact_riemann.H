#ifndef EXACT_RIEMANN_H
#define EXACT_RIEMANN_H

AMREX_INLINE
void
exact_riemann() {

    // we need a composition to interface with our EOS, but we are not
    // exploring composition jumps here.  We'll take a constant
    // composition.

    Real xn[NumSpec] = {0.0};
    xn[0] = 1.0_rt;


    // if we are using T as the independent variable (rather than p), then
    // get p now
    if (use_Tinit) {

        eos_t eos_state;
        eos_state.rho = rho_l;
        eos_state.T = T_l;
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = xn[n];
        }

        eos(eos_input_rt, eos_state);

        p_l = eos_state.p;

        eos_state.rho = rho_r;
        eos_state.T = T_r;

        eos(eos_input_rt, eos_state);

        p_r = eos_state.p;

    }

    if (co_moving_frame) {
        Real W_avg = 0.5_rt * (u_l + u_r);
        u_l -= W_avg;
        u_r -= W_avg;
    }

    riemann_star_state(rho_l, u_l, p_l, xn,
                       rho_r, u_r, p_r, xn,
                       ustar, pstar, W_l, W_r);


    // find the solution as a function of xi = x/t
    ofstream ofile;
    ofile.open("riemann.out");

    // This follows from the discussion around C&G Eq. 15

    Real dx = (xmax - xmin) / static_cast<Real>(npts);

    // loop over xi space

    for (int i = 1; i <= npts; i++) {

        // compute xi = x/t -- this is the similarity variable for the
        // solution
        Real x  = xmin + (static_cast<Real>(i) - 0.5_rt) * dx;

        call riemann_sample(rho_l, u_l, p_l, rho_r, u_r, p_r, xn, xn,
                            ustar, pstar, W_l, W_r,
                            x, xjump, t,
                            rho, u, p, xn_s);

        if (co_moving_frame) {
            u += W_avg;
            x += t * W_avg;
        }

        // get the thermodynamics for this state for output

        eos_t eos_state;
        eos_state.rho = rho;
        eos_state.p = p;
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn(n) = xn_s(n);
        }
        eos_state.T = initial_temp_guess;

        eos(eos_input_rp, eos_state);

        if (i == 1) {
            ofile << "# i ";
            ofile << std::setw(26) << "x" << "rho" << "u" << "p" << "T" << "e" << "gamma_1" << std::endl;
        }

        ofile << std::setw(4) << i;
        ofile << std::setw(26) << std::setprecision(12) << x << rho << u << p << eos_state.T << eos_state.e << eos_state.gam1 << std::endl;

    }

    ofile.close();
}
#endif
