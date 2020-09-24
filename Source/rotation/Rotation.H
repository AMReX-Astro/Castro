#ifndef ROTATION_H
#define ROTATION_H

#include <math.H>
#include <Castro.H>
#include <Castro_util.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real
rotational_potential(GpuArray<Real, 3>& r, GpuArray<Real, 3> const& omega) {

  // Construct rotational potential, phi_R = -1/2 | omega x r |**2
  //

  Real phi = 0.0_rt;

  if (state_in_rotating_frame == 1) {

    if (rotation_include_centrifugal == 1) {

      GpuArray<Real, 3> omega_cross_r;
      cross_product(omega, r, omega_cross_r);

      for (int idir = 0; idir < 3; idir++) {
        phi -= 0.5_rt * omega_cross_r[idir] * omega_cross_r[idir];
      }

    }
  }

  return phi;

}

AMREX_GPU_HOST_DEVICE 
void
inertial_to_rotational_velocity_c(const int i, const int j, const int k,
                                    const GeometryData& geomdata,
                                    const Real* omega,
                                    const Real time, Real* v);

AMREX_GPU_HOST_DEVICE 
void
inertial_to_rotational_velocity(const int i, const int j, const int k,
                                const GeometryData& geomdata,
                                const GpuArray<Real, 3>& omega,
                                const Real time, GpuArray<Real, 3>& v);

#endif
