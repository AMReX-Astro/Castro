#ifndef ROTATION_H
#define ROTATION_H

#include <math.H>
#include <AMReX_Array.H>
#include <Castro.H>
#include <Castro_util.H>

///
/// Return the omega vector corresponding to the current rotational period
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
GpuArray<Real, 3> get_omega()
{
    GpuArray<Real, 3> omega = {0.0_rt};

    // If rotational_period is less than zero, that means rotation is disabled, and so we should effectively
    // shut off the source term by setting omega = 0.

    if (castro::rotational_period > 0.0_rt) {
        omega[castro::rot_axis - 1] = 2.0_rt * M_PI / castro::rotational_period;
    }

    return omega;
}

///
/// Compute the rotational acceleration for a single zone
/// (Coriolis and centrifugal)
///
/// @param r         distance from origin of rotation vector
/// @param v         velocity
/// @param coriolis  do we include the Coriolis force
/// @param Sr        rotational acceleration
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
rotational_acceleration(GpuArray<Real, 3>& r, GpuArray<Real, 3>& v,
                        const bool coriolis, Real* Sr) {

  // Given a position and velocity, calculate
  // the rotational acceleration. This is the sum of:
  // the Coriolis force (-2 omega x v),
  // the centrifugal force (- omega x ( omega x r)),
  // and a changing rotation rate (-d(omega)/dt x r).

  Sr[0] = 0.0;
  Sr[1] = 0.0;
  Sr[2] = 0.0;

  auto omega = get_omega();

  if (castro::state_in_rotating_frame == 1) {

    // Allow the various terms to be turned off.  This is often used
    // for diagnostic purposes, but there are genuine science cases
    // for disabling certain terms in some cases (in particular, when
    // obtaining a system in rotational equilibrium through a
    // relaxation process involving damping or sponging, one may want
    // to turn off the Coriolis force during the relaxation process,
    // on the basis that the true equilibrium state will have zero
    // velocity anyway).

    bool c1 = (castro::rotation_include_centrifugal == 1) ? true : false;

    bool c2 = (castro::rotation_include_coriolis == 1 && coriolis) ? true : false;

    GpuArray<Real, 3> omega_cross_v;
    cross_product(omega, v, omega_cross_v);

    if (c1) {
      GpuArray<Real, 3> omega_cross_r;
      cross_product(omega, r, omega_cross_r);

      GpuArray<Real, 3> omega_cross_omega_cross_r;
      cross_product(omega, omega_cross_r, omega_cross_omega_cross_r);

      for (int idir = 0; idir < 3; idir++) {
        Sr[idir] -= omega_cross_omega_cross_r[idir];
      }
    }

    if (c2) {
      for (int idir = 0; idir < 3; idir++) {
        Sr[idir] -= 2.0_rt * omega_cross_v[idir];
      }
    }

  } else {

    // The source term for the momenta when we're not measuring state
    // variables in the rotating frame is not strictly the traditional
    // Coriolis force, but we'll still allow it to be disabled with
    // the same parameter.

    bool c2 = (castro::rotation_include_coriolis == 1 && coriolis) ? true : false;

    if (c2) {
      GpuArray<Real, 3> omega_cross_v;
      cross_product(omega, v, omega_cross_v);

      for (int idir = 0; idir < 3; idir++) {
        Sr[idir] -= omega_cross_v[idir];
      }
    }

  }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real
rotational_potential(GpuArray<Real, 3>& r) {

  // Construct rotational potential, phi_R = -1/2 | omega x r |**2
  //

  Real phi = 0.0_rt;

  auto omega = get_omega();

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
                                    const Real time, Real* v);

AMREX_GPU_HOST_DEVICE 
void
inertial_to_rotational_velocity(const int i, const int j, const int k,
                                const GeometryData& geomdata,
                                const Real time, GpuArray<Real, 3>& v);

#endif
