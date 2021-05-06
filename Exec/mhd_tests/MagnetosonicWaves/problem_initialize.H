#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_INLINE
void problem_initialize ()
{

    // get unitary vectors for k and b
    // such that b is 45 degrees from k

    problem::u_kx = problem::k_x / std::sqrt(problem::k_x * problem::k_x +
                                             problem::k_y * problem::k_y);
    problem::u_ky = problem::k_y / std::sqrt(problem::k_x * problem::k_x +
                                             problem::k_y * problem::k_y);

    problem::u_bx = (problem::u_kx - problem::u_ky) * 0.5_rt * std::sqrt(2.0_rt);
    problem::u_by = (problem::u_kx + problem::u_ky) * 0.5_rt * std::sqrt(2.0_rt);

    problem::B_x = problem::B_0 * problem::u_bx;
    problem::B_y = problem::B_0 * problem::u_by;
    problem::B_z = 0.0_rt;

    Real xn[NumSpec] = {0.0};
    xn[0] = 1.0_rt;

    eos_t eos_state;

    // compute sound speed

    eos_state.rho = problem::rho_0;
    eos_state.p = problem::p_0;
    eos_state.T = 100.0_rt;  // initial guess
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = xn[n];
    }

    eos(eos_input_rp, eos_state);

    problem::c_s = eos_state.cs;
    problem::c_s = problem::c_s * problem::c_s;

    if (problem::type_wave == 1 ) {
        // Fast magnetosonic wave speed
        problem::cw = std::sqrt(0.5_rt * (problem::c_s + 1.0_rt + std::sqrt(1.0_rt + problem::c_s * problem::c_s)));
    } else {
        // Slow magnetosonic wave speed
        problem::cw = std::sqrt(0.5_rt * (problem::c_s + 1.0_rt - std::sqrt(1.0_rt + problem::c_s * problem::c_s)));
    }

}
#endif
