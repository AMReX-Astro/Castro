#ifndef RIEMANN_CONSTANTS_H
#define RIEMANN_CONSTANTS_H

#include <AMReX_REAL.H>

using namespace amrex::literals;

namespace riemann_constants {
    constexpr amrex::Real smlp1 = 1.e-10_rt;
    constexpr amrex::Real small = 1.e-8_rt;
    constexpr amrex::Real smallu = 1.e-12_rt;
    constexpr int HISTORY_SIZE=40;
    constexpr int PSTAR_BISECT_FACTOR = 5;
}

#endif
