#ifndef CASTRO_FLUXLIMITER_H
#define CASTRO_FLUXLIMITER_H

#include <AMReX_REAL.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE inline
amrex::Real Edd_factor(Real lambda)
{

    Real f;

    if (radiation::closure == 0) {
        f = lambda;
    }
    else if (radiation::closure == 1) {
        f = 1.0_rt / 3.0_rt;
    }
    else if (radiation::closure == 2) {
        f = 1.e0_rt - 2.e0_rt * lambda;
    }
    else if (radiation::closure == 3) { // lambda + (lambda*R)**2
        if (radiation::limiter == 0) { // no limiter
            f = 1.0_rt / 3.0_rt;
        }
        else if (radiation::limiter < 10) { // approximate LP, [123]
            f = 0.5e0_rt * amrex::max(0.0_rt, (1.e0_rt - 3.e0_rt*  lambda)) +
                std::sqrt(amrex::max(0.0_rt, (1.e0_rt - 3.e0_rt * lambda)) *
                          (1.e0_rt + 5.e0_rt * lambda));
            f = lambda + f * f;
        }
        else if (radiation::limiter < 20) { // Bruenn, 1[123]
            f = 1.0e0_rt - 5.e0_rt * lambda + 9.e0_rt * lambda * lambda;
        }
        else if (radiation::limiter < 30) { // Larsen's square root, 2[123]
            f = 1.0e0_rt + lambda - 9.e0_rt * lambda * lambda;
        }
        else if (radiation::limiter < 40) { // Minerbo
            if (lambda > 2.0_rt / 9.0_rt) {
                f = 1.0_rt / 3.0_rt;
            }
            else {
                f = 1.e0_rt + 3.e0_rt * lambda - 2.e0_rt * std::sqrt(2.e0_rt * lambda);
            }
        }
#ifndef AMREX_USE_GPU
        else {
            amrex::Error("Unknown limiter");
        }
#endif
    }
    else if (radiation::closure == 4) { // 1/3 + 2/3*(lambda*R)**2
        if (radiation::limiter == 0) { // no limiter
            f = 1.0_rt / 3.0_rt;
        }
        else if (radiation::limiter < 10) { // approximate LP, [123]
            f = amrex::max(0.0_rt, (1.e0_rt - 3.e0_rt * lambda)) +
                 std::sqrt(amrex::max(0.0_rt, (1.e0_rt - 3.e0_rt * lambda)) *
                           (1.e0_rt + 5.e0_rt * lambda));
            f = (1.0_rt / 3.0_rt) + (1.0_rt / 6.0_rt) * f * f;
        }
        else if (radiation::limiter < 20) { // Bruenn, 1[123]
            f = (1.0_rt / 3.0_rt) + (2.0_rt / 3.0_rt) * (1.0e0_rt - 6.e0_rt * lambda + 9.e0_rt * lambda * lambda);
        }
        else if (radiation::limiter < 30) { // Larsen's square root, 2[123]
            f = (1.0_rt / 3.0_rt) + (2.0_rt / 3.0_rt) * (1.0e0_rt - 9.e0_rt * lambda * lambda);
        }
        else if (radiation::limiter < 40) { // Minerbo
            if (lambda > 2.0_rt / 9.0_rt) {
                f = (5.0_rt / 9.0_rt) - (2.0_rt / 3.0_rt) * lambda;
            }
            else {
                f = (1.0_rt / 3.0_rt) + (2.0_rt / 3.0_rt) * (1.e0_rt + 2.e0_rt * lambda - 2.e0_rt * std::sqrt(2.e0_rt * lambda));
            }
        }
#ifndef AMREX_USE_GPU
        else {
            amrex::Error("Unknown limiter");
        }
#endif
    }

    return f;

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real FLDalpha(Real lam)
{
    Real omtl = amrex::max(0.0_rt, 1.0_rt - 3.0_rt * lam);

    Real R;

    if (radiation::limiter == 0) { // no limiter
        R = 0.0_rt;
    }
    else if (radiation::limiter < 10) { // approximate LP, [123]
        R = (omtl + std::sqrt(omtl * (1.0_rt + 5.0_rt * lam))) / (2.0_rt * lam + 1.e-50_rt);
    }
    else if (radiation::limiter < 20) { // Bruenn, 1[123]
        R = omtl / (lam + 1.e-50_rt);
    }
    else if (radiation::limiter < 30) { // Larsen's square root, 2[123]
        R = std::sqrt(omtl * (1.0_rt + 3.0_rt * lam)) / (lam + 1.e-50_rt);
    }
    else if (radiation::limiter < 40) { // Minerbo
        if (lam > 2.0_rt / 9.0_rt) {
            R = std::sqrt(omtl / 3.0_rt) / (lam + 1.e-50_rt);
        }
        else {
            R = 1.0_rt / (lam + 1.e-50_rt) - std::sqrt(2.0_rt / (lam + 1.e-50_rt));
        }
    }
    else {
#ifndef AMREX_USE_GPU
        amrex::Error("Unknown limiter");
#endif
    }

    Real alpha;

    if (R < 1.e-6_rt) {
        alpha = 0.25_rt;
    }
    else if (R > 300.0_rt) {
        alpha = 0.5_rt;
    }
    else {
        Real cr = std::cosh(R);
        alpha = cr * std::log(cr) / (2.0_rt * R * std::sinh(R));
    }

    return alpha;
}

#endif
