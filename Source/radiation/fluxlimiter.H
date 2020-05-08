#ifndef CASTRO_FLUXLIMITER_H
#define CASTRO_FLUXLIMITER_H

#include <AMReX_REAL.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE inline
amrex::Real Edd_factor(Real lambda, int limiter, int closure)
{

    Real f;

    if (closure == 0) {
        f = lambda;
    }
    else if (closure == 1) {
        f = 1.0_rt / 3.0_rt;
    }
    else if (closure == 2) {
        f = 1.e0_rt - 2.e0_rt * lambda;
    }
    else if (closure == 3) { // lambda + (lambda*R)**2
        if (limiter == 0) { // no limiter
            f = 1.0_rt / 3.0_rt;
        }
        else if (limiter < 10) { // approximate LP, [123]
            f = 0.5e0_rt * amrex::max(0.0_rt, (1.e0_rt - 3.e0_rt*  lambda)) + 
                std::sqrt(amrex::max(0.0_rt, (1.e0_rt - 3.e0_rt * lambda)) *
                          (1.e0_rt + 5.e0_rt * lambda));
            f = lambda + f * f;
        }
        else if (limiter < 20) { // Bruenn, 1[123]
            f = 1.0e0_rt - 5.e0_rt * lambda + 9.e0_rt * lambda * lambda;
        }
        else if (limiter < 30) { // Larsen's square root, 2[123]
            f = 1.0e0_rt + lambda - 9.e0_rt * lambda * lambda;
        }
        else if (limiter < 40) { // Minerbo
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
    else if (closure == 4) { // 1/3 + 2/3*(lambda*R)**2
        if (limiter == 0) { // no limiter
            f = 1.0_rt / 3.0_rt;
        }
        else if (limiter < 10) { // approximate LP, [123]
            f = amrex::max(0.0_rt, (1.e0_rt - 3.e0_rt * lambda)) +
                 std::sqrt(amrex::max(0.0_rt, (1.e0_rt - 3.e0_rt * lambda)) *
                           (1.e0_rt + 5.e0_rt * lambda));
            f = (1.0_rt / 3.0_rt) + (1.0_rt / 6.0_rt) * f * f;
        }
        else if (limiter < 20) { // Bruenn, 1[123]
            f = (1.0_rt / 3.0_rt) + (2.0_rt / 3.0_rt) * (1.0e0_rt - 6.e0_rt * lambda + 9.e0_rt * lambda * lambda);
        }
        else if (limiter < 30) { // Larsen's square root, 2[123]
            f = (1.0_rt / 3.0_rt) + (2.0_rt / 3.0_rt) * (1.0e0_rt - 9.e0_rt * lambda * lambda);
        }
        else if (limiter < 40) { // Minerbo
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

#endif
