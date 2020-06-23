#include "Castro.H"
#include "Castro_F.H"
#include "math.H"

void
Castro::rotational_acceleration(GpuArray<Real, 3>& r, GpuArray<Real, 3>& v,
                                GpuArray<Real, 3> const& omega,
                                const bool coriolis, Real* Sr) {

  // Given a position and velocity, calculate
  // the rotational acceleration. This is the sum of:
  // the Coriolis force (-2 omega x v),
  // the centrifugal force (- omega x ( omega x r)),
  // and a changing rotation rate (-d(omega)/dt x r).

  Sr[0] = 0.0;
  Sr[1] = 0.0;
  Sr[2] = 0.0;

  if (state_in_rotating_frame == 1) {

    // Allow the various terms to be turned off.  This is often used
    // for diagnostic purposes, but there are genuine science cases
    // for disabling certain terms in some cases (in particular, when
    // obtaining a system in rotational equilibrium through a
    // relaxation process involving damping or sponging, one may want
    // to turn off the Coriolis force during the relaxation process,
    // on the basis that the true equilibrium state will have zero
    // velocity anyway).

    bool c1 = (rotation_include_centrifugal == 1) ? true : false;

    bool c2 = (rotation_include_coriolis == 1 && coriolis) ? true : false;

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

    bool c2 = (rotation_include_coriolis == 1 && coriolis) ? true : false;

    if (c2) {
      GpuArray<Real, 3> omega_cross_v;
      cross_product(omega, v, omega_cross_v);

      for (int idir = 0; idir < 3; idir++) {
        Sr[idir] -= omega_cross_v[idir];
      }
    }

  }
}


void
Castro::fill_rotational_acceleration(const Box& bx,
                                     Array4<Real> const& rot,
                                     Array4<Real const> const& state) {

  GpuArray<Real, 3> center;
  ca_get_center(center.begin());

  auto problo = geom.ProbLoArray();

  auto dx = geom.CellSizeArray();

  auto coord_type = geom.Coord();

  GpuArray<Real, 3> omega;
  get_omega(omega.begin());

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    GpuArray<Real, 3> r;

    r[0] = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt) - center[0];
#if AMREX_SPACEDIM >= 2
    r[1] = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt) - center[1];
#else
    r[1] = 0.0_rt;
#endif
#if AMREX_SPACEDIM == 3
    r[2] = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt) - center[2];
#else
    r[2] = 0.0_rt;
#endif

    GpuArray<Real, 3> v;

    v[0] = state(i,j,k,UMX) / state(i,j,k,URHO);
    v[1] = state(i,j,k,UMY) / state(i,j,k,URHO);
    v[2] = state(i,j,k,UMZ) / state(i,j,k,URHO);

    bool coriolis = true;
    Real Sr[3];
    rotational_acceleration(r, v, omega, coriolis, Sr);

    for (int idir = 0; idir < 3; idir++) {
      rot(i,j,k,idir) = Sr[idir];
    }
  });
}

