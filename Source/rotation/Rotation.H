#include "math.H"

#ifndef CASTRO_ROTATION_H
#define CASTRO_ROTATION_H

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
get_omega(const int coord, GpuArray<Real, 3>& omega) {

  // If rot_period is less than zero, that means rotation is disabled,
  // and so we should effectively shut off the source term by setting
  // omega = 0. Note that by default rot_axis == 2 for Cartesian
  // coordinates and rot_axis == 1 for cylindrical coordinates.

  omega[0] = 0.0_rt;
  omega[1] = 0.0_rt;
  omega[2] = 0.0_rt;

  if (coord == 0 || coord == 1) {

    if (rotational_period > 0.0_rt) {

      // If we have a time rate of change of the rotational period,
      // adjust it accordingly in the calculation of omega. We assume
      // that the change has been linear and started at t == 0.

      Real curr_period = rotational_period;

      omega[rot_axis] = 2.0_rt * M_PI / curr_period;

    }

  } else { 

#ifndef AMREX_USE_GPU
    amrex::Error("Error:: rotation_nd.f90 :: invalid coord_type");
#endif
  }
}




AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
set_rot_period(const Real period) {
  rotational_period = period;
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
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



#endif
