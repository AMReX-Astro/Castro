#ifndef ROTATION_H
#define ROTATION_H

#include <math.H>
#include <AMReX_Array.H>
#include <Castro.H>
#include <Castro_util.H>


///
/// Return the magnitude of omega along the rotating axis.
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real get_omega()
{
    Real omega = 0.0_rt;

    // If rotational_period is less than zero, that means rotation is disabled, and so we should effectively
    // shut off the source term by setting omega = 0.

    if (castro::rotational_period > 0.0_rt) {
        omega = 2.0_rt * M_PI / castro::rotational_period;
    }

    return omega;
}

///
/// Return the omega vector corresponding to the current rotational period
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
GpuArray<Real, 3> get_omega_vec(const GeometryData& geomdata,
                                const int j)
{
    GpuArray<Real, 3> omega_vec = {0.0_rt};

    const auto coord = geomdata.Coord();
    const auto problo = geomdata.ProbLo(1);
    const auto dtheta = geomdata.CellSize(1);

    auto omega = get_omega();

    // If coord == 2, i.e. Spherical2D, then it is assumed that rotational axis is in z-axis
    // then we convert back to r and theta direction.
    if (coord == 2) {
        Real theta = problo + (static_cast<Real>(j) + 0.5_rt) * dtheta;
        omega_vec[0] = std::cos(theta) * omega;
        omega_vec[1] = -std::sin(theta) * omega;
    } else {
        omega_vec[castro::rot_axis - 1] = omega;
    }

    return omega_vec;
}

///
/// Compute the rotational acceleration for a single zone
/// (Coriolis and centrifugal)
///
/// @param r         distance from origin of rotation vector
/// @param v         velocity
/// @param omega     angular velocity vector array
/// @param coord     coordinate system
/// @param coriolis  do we include the Coriolis force
/// @param Sr        rotational acceleration
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
rotational_acceleration(const GpuArray<Real, 3>& r, const GpuArray<Real, 3>& v,
                        const GpuArray<Real, 3>& omega, const int coord,
                        const bool coriolis, Real* Sr) {

  // Given a position and velocity, calculate
  // the rotational acceleration. This is the sum of:
  // the Coriolis force (-2 omega x v),
  // the centrifugal force (- omega x ( omega x r)),
  // and a changing rotation rate (-d(omega)/dt x r).

  Sr[0] = 0.0;
  Sr[1] = 0.0;
  Sr[2] = 0.0;

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

  auto r_vec = r;

  GpuArray<Real, 3> omega_cross_v;
  cross_product(omega, v, omega_cross_v);

  if (c1) {
      // For Spherical coordinate, the physical position vector is just r rhat,
      // but the input r will be in the format of (r,theta,0). So change it to (r,0,0)

      if (coord == 2) {
          r_vec[1] = 0.0_rt;
          r_vec[2] = 0.0_rt;
      }
      GpuArray<Real, 3> omega_cross_r;
      cross_product(omega, r_vec, omega_cross_r);

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

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real
rotational_potential(const GpuArray<Real, 3>& r, const GpuArray<Real, 3>& omega,
                     const int coord) {

  // Construct rotational potential, phi_R = -1/2 | omega x r |**2
  //

  Real phi = 0.0_rt;
  auto r_vec = r;

  if (coord == 2) {
      r_vec[1] = 0.0_rt;
      r_vec[2] = 0.0_rt;
  }

  if (rotation_include_centrifugal == 1) {

      GpuArray<Real, 3> omega_cross_r;
      cross_product(omega, r, omega_cross_r);

      for (int idir = 0; idir < 3; idir++) {
          phi -= 0.5_rt * omega_cross_r[idir] * omega_cross_r[idir];
      }

  }

  return phi;

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
inertial_to_rotational_velocity (const int i, const int j, const int k,
                                 const GeometryData& geomdata,
                                 const Real time, GpuArray<Real, 3>& v) {

  amrex::ignore_unused(time);

  // Given a velocity vector in the inertial frame, transform it to a
  // velocity vector in the rotating frame.

  // Note: this version assumes all cell-centers

  GpuArray<Real, 3> loc;

  position(i, j, k, geomdata, loc);

  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    loc[dir] -= problem::center[dir];
  }

  if (geomdata.Coord() == 2) {
      loc[1] = 0.0_rt;
      loc[2] = 0.0_rt;
  }

  auto omega = get_omega_vec(geomdata, j);

  // do the cross product Omega x loc
  v[0] += -(omega[1]*loc[2] - omega[2]*loc[1]);
  v[1] += -(omega[2]*loc[0] - omega[0]*loc[2]);
  v[2] += -(omega[0]*loc[1] - omega[1]*loc[0]);

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
rotational_to_inertial_velocity (const int i, const int j, const int k,
                                 const GeometryData& geomdata,
                                 const Real time, GpuArray<Real, 3>& v) {

  amrex::ignore_unused(time);

  // Given a velocity vector in the rotating frame, transform it to a
  // velocity vector in the inertial frame.

  // Note: this version assumes all cell-centers

  GpuArray<Real, 3> loc;

  position(i, j, k, geomdata, loc);

  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    loc[dir] -= problem::center[dir];
  }

  if (geomdata.Coord() == 2) {
      loc[1] = 0.0_rt;
      loc[2] = 0.0_rt;
  }

  auto omega = get_omega_vec(geomdata, j);

  // do the cross product Omega x loc
  v[0] += (omega[1]*loc[2] - omega[2]*loc[1]);
  v[1] += (omega[2]*loc[0] - omega[0]*loc[2]);
  v[2] += (omega[0]*loc[1] - omega[1]*loc[0]);

}

// If we are in a rotating reference frame, then rotate a vector
// by an amount corresponding to the time that has passed
// since the beginning of the simulation.

AMREX_GPU_HOST_DEVICE AMREX_INLINE
GpuArray<Real, 3> inertial_rotation(const GpuArray<Real, 3>& vec,
                                    const GpuArray<Real, 3>& omega, Real time)
{
    GpuArray<Real, 3> vec_i{};

    Array1D<Real, 0, 2> theta{};
    Array2D<Real, 0, 2, 0, 2> rot_matrix{};

    if (castro::do_rotation == 1) {
        for (int n = 0; n < 3; ++n) {
            theta(n) = omega[n] * time;
        }
    }

    // This is the 3D rotation matrix for converting between reference frames.
    // It is the composition of rotations along the x, y, and z axes. Therefore
    // it allows for the case where we are rotating about multiple axes. Normally
    // we use the right-hand convention for constructing the usual rotation matrix,
    // but this is the transpose of that rotation matrix to account for the fact
    // that we are rotating *back* to the inertial frame, rather than from the
    // inertial frame to the rotating frame.

    rot_matrix(0,0) =  cos(theta(1)) * cos(theta(2));
    rot_matrix(0,1) = -cos(theta(1)) * sin(theta(2));
    rot_matrix(0,2) =  sin(theta(1));
    rot_matrix(1,0) =  cos(theta(0)) * sin(theta(2)) + sin(theta(0)) * sin(theta(1)) * cos(theta(2));
    rot_matrix(1,1) =  cos(theta(0)) * cos(theta(2)) - sin(theta(0)) * sin(theta(1)) * sin(theta(2));
    rot_matrix(1,2) = -sin(theta(0)) * cos(theta(1));
    rot_matrix(2,0) =  sin(theta(0)) * sin(theta(2)) - cos(theta(0)) * sin(theta(1)) * cos(theta(2));
    rot_matrix(2,1) =  sin(theta(0)) * cos(theta(2)) + cos(theta(0)) * sin(theta(1)) * sin(theta(2));
    rot_matrix(2,2) =  cos(theta(0)) * cos(theta(1));

    for (int l = 0; l < 3; ++l) {
        for (int m = 0; m < 3; ++m) {
            vec_i[l] += rot_matrix(l,m) * vec[m];
        }
    }

    return vec_i;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void fill_dt_omega_matrix(const Real dt, const GpuArray<Real, 3>& omega,
                          Array2D<Real, 0, 2, 0, 2>& dt_omega_matrix)
{
    // Fill in dt_omega_matrix. This is used in function corrrsrc

    Real dt_omega[3];

    // Don't do anything here if we've got the Coriolis force disabled.

    if (castro::rotation_include_coriolis == 1) {

        for (int idir = 0; idir < 3; idir++) {
          dt_omega[idir] = dt * omega[idir];
        }

    } else {

      for (auto& e : dt_omega) {
        e = 0.0_rt;
      }

    }

    dt_omega_matrix(0, 0) = 1.0_rt + dt_omega[0] * dt_omega[0];
    dt_omega_matrix(0, 1) = dt_omega[0] * dt_omega[1] + dt_omega[2];
    dt_omega_matrix(0, 2) = dt_omega[0] * dt_omega[2] - dt_omega[1];

    dt_omega_matrix(1, 0) = dt_omega[1] * dt_omega[0] - dt_omega[2];
    dt_omega_matrix(1, 1) = 1.0_rt + dt_omega[1] * dt_omega[1];
    dt_omega_matrix(1, 2) = dt_omega[1] * dt_omega[2] + dt_omega[0];

    dt_omega_matrix(2, 0) = dt_omega[2] * dt_omega[0] + dt_omega[1];
    dt_omega_matrix(2, 1) = dt_omega[2] * dt_omega[1] - dt_omega[0];
    dt_omega_matrix(2, 2) = 1.0_rt + dt_omega[2] * dt_omega[2];

    for (int l = 0; l < 3; l++) {
      for (int m = 0; m < 3; m++) {
        dt_omega_matrix(l, m) /= (1.0_rt + dt_omega[0] * dt_omega[0] +
                                           dt_omega[1] * dt_omega[1] +
                                           dt_omega[2] * dt_omega[2]);
      }
    }

}

#endif
