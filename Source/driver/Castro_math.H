#ifndef CASTRO_MATH_H
#define CASTRO_MATH_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_BLassert.H>
#include <cmath>

using namespace amrex::literals;


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
cross_product(amrex::GpuArray<amrex::Real, 3> const& a,
              amrex::GpuArray<amrex::Real, 3> const& b,
              amrex::GpuArray<amrex::Real, 3>& c) {

  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];

}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real cot(amrex::Real x) {

    [[maybe_unused]] constexpr amrex::Real eps = 1.e-14_rt;
    AMREX_ASSERT(std::abs(std::sin(x)) > eps);

    // we want to reflect outside of x in [0, pi]
    // i.e. if x = -0.1 -> x = 0.1
    // and  if x = 5π/4 -> x = 3π/4

    if (x < 0.0_rt || x > M_PI) {
        x = M_PI - x;
    }
    return std::cos(x) / std::sin(x);
}
#endif
