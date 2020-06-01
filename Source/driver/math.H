#ifndef CASTRO_MATH_H
#define CASTRO_MATH_H

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
cross_product(GpuArray<Real, 3> const& a, GpuArray<Real, 3> const& b, GpuArray<Real, 3>& c) {

  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];

}

#endif
