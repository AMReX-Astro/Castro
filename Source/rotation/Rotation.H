#ifndef CASTRO_ROTATION_H
#define CASTRO_ROTATION_H

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
Castro::get_omega(const Real time, GpuArray<Real, 3>& omega) {

  // If rot_period is less than zero, that means rotation is disabled,
  // and so we should effectively shut off the source term by setting
  // omega = 0. Note that by default rot_axis == 3 for Cartesian
  // coordinates and rot_axis == 2 for cylindrical coordinates.

  omega[0] = 0.0_rt;
  omega[1] = 0.0_rt;
  omega[2] = 0.0_rt;

  coord = geom.Coord();

  if (coord_type == 0 || coord_type == 1) {

    if (rotational_period > 0.0_rt) {

      // If we have a time rate of change of the rotational period,
      // adjust it accordingly in the calculation of omega. We assume
      // that the change has been linear and started at t == 0.

      Real curr_period = rotation_period + rotational_dPdt * time;

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
Castro::get_domegadt(const Real time, GpuArray<Real, 3>& domegadt) {

  domegadt[0] = 0.0_rt;
  domegadt[1] = 0.0_rt;
  domegadt[2] = 0.0_rt;

  if (coord_type == 0 || coord_type == 1) {

    if (rotational_period > 0.0_rt) {

      // Rate of change of the rotational frequency is given by
      // d( ln(period) ) / dt = - d( ln(omega) ) / dt

      Real curr_period = rotational_period + rotational_dPdt * time;

      GpuArray<Real, 3> curr_omega;
      get_omega(time, curr_omega);

      for (int idir = 0; idir < 3; idir++) {
        domegadt[idir] = -curr_omega[idir] * (rot_period_dot / curr_period);
      }
    }

  } else {
#ifndef AMREX_USE_GPU
    amrex::Error("Error:: rotation_nd.f90 :: unknown coord_type");
#endif
  }
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
Castro::set_rot_period(const Real period) {
  rotational_period = period;
}

#endif
