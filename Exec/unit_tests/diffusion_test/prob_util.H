#ifndef PROB_UTIL_H
#define PROB_UTIL_H

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real analytic(const Real* r, const Real time, const int coord_type) {

    Real exponent;
    if (coord_type == 2) {
        // Handle spherical coordinates
        exponent = 3.0_rt / 2.0_rt;
    } else if (AMREX_SPACEDIM == 2 && coord_type == 1) {
        // Handle cylindrical coordinates
        exponent = 3.0_rt / 2.0_rt;
    } else {
        exponent = AMREX_SPACEDIM / 2.0_rt;
    }

    Real dist2 = 0.0;
    for (int d = 0; d < AMREX_SPACEDIM; d++) {
        dist2 += r[d] * r[d];
    }

    // Real temp = problem::T1 + (problem::T2 - problem::T1) *
    //     std::pow(problem::t_0 / (time + problem::t_0), exponent) *
    //     std::exp(-0.25_rt * dist2 / (problem::diff_coeff * (time + problem::t_0)));

    Real kr = std::sqrt(1.0_rt / (problem::diff_coeff * problem::t_0)) * r[0];
    Real temp = problem::T1 + (problem::T2 - problem::T1) *
        (std::sin(kr) / (kr * kr) - std::cos(kr) / kr) * std::cos(r[1]) *
        std::exp(-1.0_rt * time / problem::t_0);

    return temp;
}
#endif
