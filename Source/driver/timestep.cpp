#include "Castro.H"
#include "Castro_F.H"

#ifdef DIFFUSION
#include "conductivity.H"
#endif

using namespace amrex;

Real
Castro::estdt_cfl(const Real time)
{

  // Courant-condition limited timestep

  GpuArray<Real, 3> center;
  ca_get_center(center.begin());

#ifdef ROTATION
  GpuArray<Real, 3> omega;
  get_omega(time, omega.begin());
#endif

  const auto dx = geom.CellSizeArray();

  ReduceOps<ReduceOpMin> reduce_op;
  ReduceData<Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  const MultiFab& stateMF = get_new_data(State_Type);

  const int ltime_integration_method = time_integration_method;
#ifdef ROTATION
  const int ldo_rotation = do_rotation;
  const int lstate_in_rotating_frame = state_in_rotating_frame;
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(stateMF, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const Box& box = mfi.tilebox();

    auto u = stateMF.array(mfi);

    reduce_op.eval(box, reduce_data,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
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
      if (ldo_rotation == 1 && lstate_in_rotating_frame != 1) {
        Real vel[3];
        vel[0] = ux;
        vel[1] = uy;
        vel[2] = uz;

        GeometryData geomdata = geom.data();

        inertial_to_rotational_velocity_c(i, j, k, geomdata,
                                          center.begin(), omega.begin(), time, vel);

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
      if (ltime_integration_method == 0 || ltime_integration_method == 3) {
        return {amrex::min(dt1, dt2, dt3)};

      } else {
        // method of lines-style constraint is tougher
        Real dt_tmp = 1.0_rt/dt1;
#if AMREX_SPACEIM >= 2
        dt_tmp += 1.0_rt/dt2;
#endif
#if AMREX_SPACEDIM == 3
        dt_tmp += 1.0_rt/dt3;
#endif

        return 1.0_rt/dt_tmp;
      }

    });

  }

  ReduceTuple hv = reduce_data.value();
  Real estdt_hydro = amrex::get<0>(hv);

  return estdt_hydro;

}

#ifdef DIFFUSION

Real
Castro::estdt_temp_diffusion(void)
{

  // Diffusion-limited timestep
  //
  // dt < 0.5 dx**2 / D
  // where D = k/(rho c_v), and k is the conductivity

  const auto dx = geom.CellSizeArray();

  ReduceOps<ReduceOpMin> reduce_op;
  ReduceData<Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  const MultiFab& stateMF = get_new_data(State_Type);

  const Real ldiffuse_cutoff_density = diffuse_cutoff_density;
  const Real lmax_dt = max_dt;
  const Real lcfl = cfl;

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(stateMF, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const Box& box = mfi.tilebox();

    auto ustate = stateMF.array(mfi);

    reduce_op.eval(box, reduce_data,
                   [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
                   {

                     if (ustate(i,j,k,URHO) > ldiffuse_cutoff_density) {

                       Real rho_inv = 1.0_rt/ustate(i,j,k,URHO);

                       // we need cv
                       eos_t eos_state;
                       eos_state.rho = ustate(i,j,k,URHO);
                       eos_state.T = ustate(i,j,k,UTEMP);
                       eos_state.e = ustate(i,j,k,UEINT) * rho_inv;
                       for (int n = 0; n < NumSpec; n++) {
                         eos_state.xn[n]  = ustate(i,j,k,UFS+n) * rho_inv;
                       }
                       for (int n = 0; n < NumAux; n++) {
                         eos_state.aux[n] = ustate(i,j,k,UFX+n) * rho_inv;
                       }

                       eos(eos_input_re, eos_state);

                       // we also need the conductivity
                       conductivity(eos_state);

                       // maybe we should check (and take action) on negative cv here?
                       Real D = eos_state.conductivity * rho_inv / eos_state.cv;

                       Real dt1 = 0.5_rt * dx[0]*dx[0] / D;

                       Real dt2;
#if AMREX_SPACEDIM >= 2
                       dt2 = 0.5_rt * dx[1]*dx[1] / D;
#else
                       dt2 = dt1;
#endif

                       Real dt3;
#if AMREX_SPACEDIM >= 3
                       dt3 = 0.5_rt * dx[2]*dx[2] / D;
#else
                       dt3 = dt1;
#endif

                       return {amrex::min(dt1, dt2, dt3)};

                     } else {
                       return lmax_dt/lcfl;
      }
    });
  }

  ReduceTuple hv = reduce_data.value();
  Real estdt_diff = amrex::get<0>(hv);

  return estdt_diff;
}
#endif

bool
Castro::check_timestep(Real dt)
{
    BL_PROFILE("Castro::check_timestep()");

    MultiFab& S_new = get_new_data(State_Type);

#ifdef REACTIONS
    MultiFab& R_new = get_new_data(Reactions_Type);
#endif

    auto dx = geom.CellSizeArray();

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        auto U = S_new.array(mfi);
#ifdef REACTIONS
        auto R = R_new.array(mfi);
#endif

        reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
        {
            Real zone_failed = 0.0_rt;

            const Real derivative_floor = 1.e-50_rt;

            Real dt_tmp = dt;

            Real rhoInv = 1.0_rt / U(i,j,k,URHO);

            // CFL hydrodynamic stability criterion

            // If the timestep created a velocity v and sound speed at the new time
            // such that (v+c) * dt / dx < CFL / change_max, where CFL is the user's
            // chosen timestep constraint and change_max is the factor that determines
            // how much the timestep can change during an advance, consider the advance
            // to have failed. This prevents the timestep from shrinking too much, whereas
            // in estdt change_max prevents the timestep from growing too much.

            if (castro::do_hydro == 1) {

                eos_t eos_state;

                eos_state.rho = U(i,j,k,URHO );
                eos_state.T   = U(i,j,k,UTEMP);
                eos_state.e   = U(i,j,k,UEINT) * rhoInv;
                for (int n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] = U(i,j,k,UFS+n) * rhoInv;
                }
                for (int n = 0; n < NumAux; ++n) {
                    eos_state.aux[n] = U(i,j,k,UFX+n) * rhoInv;
                }

                eos(eos_input_re, eos_state);

                Real u = std::abs(U(i,j,k,UMX)) * rhoInv;
                Real v = std::abs(U(i,j,k,UMY)) * rhoInv;
                Real w = std::abs(U(i,j,k,UMZ)) * rhoInv;

                Real c = eos_state.cs;

                Real tau_CFL = amrex::min(AMREX_D_DECL(dx[0] / (c + u), dx[1] / (c + v), dx[2] / (c + w)));

                dt_tmp = castro::cfl * tau_CFL;

            }

#ifdef REACTIONS
            // Burning stability criterion
            // See ca_estdt_burning for an explanation of these limiters.

            if (castro::do_react == 1) {

                Real tau_X = std::numeric_limits<Real>::max();

                for (int n = 0; n < NumSpec; ++n) {
                    Real X = amrex::max(small_x, U(i,j,k,UFS+n) * rhoInv);

                    Real X_dot = R(i,j,k,n);

                    if (X >= castro::dtnuc_X_threshold) {
                        X_dot = amrex::max(std::abs(X_dot), derivative_floor);
                    }
                    else {
                        X_dot = derivative_floor;
                    }

                    tau_X = amrex::min(tau_X, X / X_dot);
                }

                Real e = U(i,j,k,UEINT) * rhoInv;
                Real e_dot = R(i,j,k,NumSpec+1);
                e_dot = amrex::max(std::abs(e_dot), derivative_floor);

                Real tau_e = e / e_dot;

                dt_tmp = amrex::min(dt_tmp, dtnuc_e * tau_e);
                dt_tmp = amrex::min(dt_tmp, dtnuc_X * tau_X);

            }
#endif

            if (castro::change_max * dt_tmp < dt) {
                zone_failed = 1.0_rt;
            }

            return zone_failed;
        });

    }

    ReduceTuple hv = reduce_data.value();
    Real check_timestep_failure = amrex::get<0>(hv);

    ParallelDescriptor::ReduceRealSum(check_timestep_failure);

    return (check_timestep_failure == 0.0_rt);
}
