#ifndef PROB_UTIL_H
#define PROB_UTIL_H

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real analytic(const Real* r, const Real time, const int coord_type) {

    Real temp = 0.0_rt;

#if AMREX_SPACEDIM == 2
    if ((!problem::gaussian_init) && (coord_type == 2)) {
        // For spherical geometry with azimuthal symmetry,
        // diffusion equation has analytic form of:
        // T(r, theta, t) = A_l j_l(kr) P_l(cos(theta)) exp(-t/t_0)
        // where k^2 = 1/(D*t_0)
        // Here we initialize problem to a single mode, l = 1.

        Real kr = std::sqrt(1.0_rt / (problem::diff_coeff * problem::t_0)) * r[0];
        temp = problem::T1 + (problem::T2 - problem::T1) *
            (std::sin(kr) / (kr * kr) - std::cos(kr) / kr) * std::cos(r[1]) *
            std::exp(-1.0_rt * time / problem::t_0);
    } else
#endif
    {
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
        if (coord_type == 2) {
            dist2 = r[0] * r[0];
#if AMREX_SPACEDIM == 2
            // problem::center[0] = r_0 and problem::center[1] = theta_0
            // dist2 = r^2 + r_0^2 - 2 r r_0 cos(theta - theta_0)
            dist2 += problem::center[0]*problem::center[0] -
                2.0_rt * r[0] * problem::center[0] * std::cos(r[1] - problem::center[1]);
#endif
        } else {
            for (int d = 0; d < AMREX_SPACEDIM; d++) {
                dist2 += r[d] * r[d];
            }
        }

        temp = problem::T1 + (problem::T2 - problem::T1) *
            std::pow(problem::t_0 / (time + problem::t_0), exponent) *
            std::exp(-0.25_rt * dist2 / (problem::diff_coeff * (time + problem::t_0)));
    }

    return temp;
}
#endif
