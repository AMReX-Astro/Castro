#ifndef RIEMANN_CONSTANTS
#define RIEMANN_CONSTANTS

#include <AMReX_REAL.H>

using namespace amrex::literals;

namespace riemann_constants {
    constexpr amrex::Real smlp1 = 1.e-10_rt;
    constexpr amrex::Real small = 1.e-8_rt;
    constexpr amrex::Real smallu = 1.e-12_rt;
    constexpr amrex::Real riemann_integral_tol = 1.e-8_rt;
    constexpr amrex::Real riemann_u_tol = 1.e-6_rt;
    constexpr amrex::Real riemann_p_tol = 1.e-8_rt;
    constexpr int HISTORY_SIZE=40;
    constexpr int PSTAR_BISECT_FACTOR = 5;
}

#endif
