#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::estdt_cfl(const Box& bx,
                  Array4<Real const> const u,
                  Real time, Real* dt)
{

  // Courant-condition limited timestep

  GpuArray<Real, 3> center;
  ca_get_center(center.begin());

#ifdef ROTATION
  GpuArray<Real, 3> omega;
  get_omega(time, omega.begin());
#endif

  const auto dx = geom.CellSizeArray();

  AMREX_PARALLEL_FOR_3D(bx, i, j, k,
  {

    Real rhoInv = 1.0_rt / u(i,j,k,URHO);

    eos_t eos_state;
    eos_state.rho = u(i,j,k,URHO);
    eos_state.T = u(i,j,k,UTEMP);
    eos_state.e = u(i,j,k,UEINT) * rhoInv;
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = u(i,j,k,UFS+n) * rhoInv;
    }
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = u(i,j,k,UFX+n) * rhoInv;
    }

    eos(eos_input_re, eos_state);

    // Compute velocity and then calculate CFL timestep.

    Real ux = u(i,j,k,UMX) * rhoInv;
    Real uy = u(i,j,k,UMY) * rhoInv;
    Real uz = u(i,j,k,UMZ) * rhoInv;

#ifdef ROTATION
    if (do_rotation == 1 && state_in_rotating_frame != 1) {
      Real vel[3];
      vel[0] = ux;
      vel[1] = uy;
      vel[2] = uz;

      GeometryData geomdata = geom.data();

      inertial_to_rotational_velocity_c(i, j, k, geomdata, center.begin(), omega.begin(), time, vel);

      ux = vel[0];
      uy = vel[1];
      uz = vel[2];
    }
#endif

    Real c = eos_state.cs;

    Real dt1 = dx[0]/(c + std::abs(ux));

    Real dt2;
#if AMREX_SPACEDIM >= 2
    dt2 = dx[1]/(c + std::abs(uy));
#else
    dt2 = dt1;
#endif

    Real dt3;
#if AMREX_SPACEDIM == 3
    dt3 = dx[2]/(c + std::abs(uz));
#else
    dt3 = dt1;
#endif

    // The CTU method has a less restrictive timestep than MOL-based
    // schemes (including the true SDC).  Since the simplified SDC
    // solver is based on CTU, we can use its timestep.
    if (time_integration_method == 0 || time_integration_method == 3) {
      Gpu::Atomic::Min(dt, amrex::min(dt1, dt2, dt3));

    } else {
      // method of lines-style constraint is tougher
      Real dt_tmp = 1.0_rt/dt1;
#if AMREX_SPACEIM >= 2
      dt_tmp += 1.0_rt/dt2;
#endif
#if AMREX_SPACEDIM == 3
      dt_tmp += 1.0_rt/dt3;
#endif

      Gpu::Atomic::Min(dt, 1.0_rt/dt_tmp);
    }

  });

}
