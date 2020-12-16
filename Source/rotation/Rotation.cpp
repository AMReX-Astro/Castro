#include <Castro.H>
#include <Castro_F.H>
#include <math.H>
#include <Rotation.H>

void
Castro::fill_rotational_potential(const Box& bx,
                                  Array4<Real> const& phi,
                                  const Real time) {

  auto coord_type = geom.Coord();

  auto problo = geom.ProbLoArray();

  auto dx = geom.CellSizeArray();

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    GpuArray<Real, 3> r;

    r[0] = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt) - problem::center[0];
#if AMREX_SPACEDIM >= 2
    r[1] = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt) - problem::center[1];
#else
    r[1] = 0.0_rt;
#endif
#if AMREX_SPACEDIM == 3
    r[2] = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt) - problem::center[2];
#else
    r[2] = 0.0_rt;
#endif

    phi(i,j,k) = rotational_potential(r);

  });

}

void
Castro::fill_rotational_psi(const Box& bx,
                            Array4<Real> const& psi,
                            const Real time) {

  // Construct psi, which is the distance-related part of the rotation
  // law. See e.g. Hachisu 1986a, Equation 15.  For rigid-body
  // rotation, psi = -R^2 / 2, where R is the distance orthogonal to
  // the rotation axis. There are also v-constant and j-constant
  // rotation laws that we do not implement here. We will construct
  // this as potential / omega**2, so that the rotational_potential
  // routine uniquely determines the rotation law. For the other
  // rotation laws, we would simply divide by v_0^2 or j_0^2 instead.

  auto coord_type = geom.Coord();

  auto omega = get_omega();
  Real denom = omega[0] * omega[0] + omega[1] * omega[1] + omega[2] * omega[2];

  auto problo = geom.ProbLoArray();

  auto dx = geom.CellSizeArray();

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    GpuArray<Real, 3> r;

    r[0] = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt) - problem::center[0];
#if AMREX_SPACEDIM >= 2
    r[1] = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt) - problem::center[1];
#else
    r[1] = 0.0_rt;
#endif
#if AMREX_SPACEDIM == 3
    r[2] = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt) - problem::center[2];
#else
    r[2] = 0.0_rt;
#endif


    psi(i,j,k) = rotational_potential(r) / denom;

  });
}

AMREX_GPU_HOST_DEVICE 
void
inertial_to_rotational_velocity_c(const int i, const int j, const int k,
                                    const GeometryData& geomdata,
                                    const Real time, Real* v) {

    // Given a velocity vector in the inertial frame, transform it to a
    // velocity vector in the rotating frame.

    // Note: this version assumes all cell-centers

    GpuArray<Real, 3> loc;

    position(i, j, k, geomdata, loc);

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        loc[dir] -= problem::center[dir];
    }

    auto omega = get_omega();

    // do the cross product Omega x loc
    v[0] += -(omega[1]*loc[2] - omega[2]*loc[1]);
    v[1] += -(omega[2]*loc[0] - omega[0]*loc[2]);
    v[2] += -(omega[0]*loc[1] - omega[1]*loc[0]);

}

AMREX_GPU_HOST_DEVICE 
void
inertial_to_rotational_velocity(const int i, const int j, const int k,
                                const GeometryData& geomdata,
                                const Real time, GpuArray<Real, 3>& v) {

  // Given a velocity vector in the inertial frame, transform it to a
  // velocity vector in the rotating frame.

  // Note: this version assumes all cell-centers

  GpuArray<Real, 3> loc;

  position(i, j, k, geomdata, loc);

  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    loc[dir] -= problem::center[dir];
  }

  auto omega = get_omega();

  // do the cross product Omega x loc
  v[0] += -(omega[1]*loc[2] - omega[2]*loc[1]);
  v[1] += -(omega[2]*loc[0] - omega[0]*loc[2]);
  v[2] += -(omega[0]*loc[1] - omega[1]*loc[0]);

}
