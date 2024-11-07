#include <Castro.H>
#include <math.H>
#include <Rotation.H>

void
Castro::fill_rotational_psi(const Box& bx,
                            Array4<Real> const& psi,
                            const Real time) {

    amrex::ignore_unused(time);

  // Construct psi, which is the distance-related part of the rotation
  // law. See e.g. Hachisu 1986a, Equation 15.  For rigid-body
  // rotation, psi = -R^2 / 2, where R is the distance orthogonal to
  // the rotation axis. There are also v-constant and j-constant
  // rotation laws that we do not implement here. We will construct
  // this as potential / omega**2, so that the rotational_potential
  // routine uniquely determines the rotation law. For the other
  // rotation laws, we would simply divide by v_0^2 or j_0^2 instead.

  auto omega = get_omega();
  Real denom = omega[0] * omega[0] + omega[1] * omega[1] + omega[2] * omega[2];

  auto problo = geom.ProbLoArray();

  auto dx = geom.CellSizeArray();

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    GpuArray<Real, 3> r;
    position(i, j, k, geomdata, r);

    if (denom != 0.0) {
        psi(i,j,k) = rotational_potential(r) / denom;
    }
    else {
        psi(i,j,k) = 0.0;
    }

  });
}
