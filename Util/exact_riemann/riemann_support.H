#ifndef RIEMANN_SUPPORT_H
#define RIEMANN_SUPPORT_H

#include <AMReX_REAL.H>

namespace riemann_solver {
    const int max_iters = 100;

    const amrex::Real smallrho = 1.e-5_rt;
    const amrex::Real tol = 1.e-6_rt;
    const amrex::Real tol_p = 1.e-6_rt;
}

#endif
