#ifdef SPONGE
#include <Castro.H>

#ifdef HYBRID_MOMENTUM
#include <hybrid.H>
#endif

using namespace amrex;

void
Castro::construct_old_sponge_source(MultiFab& source, MultiFab& state_in,
                                    Real time, Real dt)
{
    // We do not apply any sponge at the old time.

    amrex::ignore_unused(source);
    amrex::ignore_unused(state_in);
    amrex::ignore_unused(time);
    amrex::ignore_unused(dt);

}

void
Castro::construct_new_sponge_source(MultiFab& source, const MultiFab& state_old, MultiFab& state_new,
                                    Real time, Real dt)
{

    amrex::ignore_unused(state_old);
    amrex::ignore_unused(time);

    const Real strt_time = ParallelDescriptor::second();

    if (!do_sponge) {
        return;
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        apply_sponge(bx, state_new.array(mfi), source.array(mfi), dt);
    }

    if (verbose > 1)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "Castro::construct_new_sponge_source() time = " << run_time << "\n" << "\n";
        }
#ifdef BL_LAZY
        });
#endif
    }
}

void
Castro::apply_sponge(const Box& bx,
                     Array4<Real const> const& state_in,
                     Array4<Real> const& source,
                     Real dt) {

  // alpha is a dimensionless measure of the timestep size; if
  // sponge_timescale < dt, then the sponge will have a larger effect,
  // and if sponge_timescale > dt, then the sponge will have a diminished effect.
  Real alpha;
  if (sponge_timescale > 0.0_rt) {
    alpha = dt / sponge_timescale;
  } else {
    alpha = 0.0_rt;
  }

  auto dx = geom.CellSizeArray();
  auto problo = geom.ProbLoArray();
  auto geomdata = geom.data();

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    Real src[NSRC] = {0.0};

    GpuArray<Real, 3> r;

    r[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0];

#if AMREX_SPACEDIM >= 2
    r[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1];
#else
    r[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
    r[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - problem::center[2];
#else
    r[2] = 0.0_rt;
#endif

    Real rho = state_in(i,j,k,URHO);
    Real rhoInv = 1.0_rt / rho;

    // compute the update factor

    // Radial distance between upper and lower boundaries.
    Real delta_r = sponge_upper_radius - sponge_lower_radius;

    // Density difference between upper and lower cutoffs.
    Real delta_rho = sponge_lower_density - sponge_upper_density;

    // Pressure difference between upper and lower cutoffs.
    Real delta_p = sponge_lower_pressure - sponge_upper_pressure;


    // Apply radial sponge. By default sponge_lower_radius will be zero
    // so this sponge is applied only if set by the user.
    Real sponge_factor = 0.0_rt;

    if (sponge_lower_radius >= 0.0_rt && sponge_upper_radius > sponge_lower_radius) {
      Real rad = distance(geomdata, r);

      if (rad < sponge_lower_radius) {
        sponge_factor = sponge_lower_factor;

      } else if (rad >= sponge_lower_radius && rad <= sponge_upper_radius) {
        sponge_factor = sponge_lower_factor +
          0.5_rt * (sponge_upper_factor - sponge_lower_factor) *
          (1.0_rt - std::cos(M_PI * (rad - sponge_lower_radius) / delta_r));

      } else {
        sponge_factor = sponge_upper_factor;
      }
    }

    // Apply density sponge. This sponge is applied only if set by the user.

    // Note that because we do this second, the density sponge gets priority
    // over the radial sponge in cases where the two would overlap.

    if (sponge_upper_density > 0.0_rt && sponge_lower_density > 0.0_rt) {
      if (rho > sponge_upper_density) {
        sponge_factor = sponge_lower_factor;

      } else if (rho <= sponge_upper_density && rho >= sponge_lower_density) {
        sponge_factor = sponge_lower_factor +
          0.5_rt * (sponge_upper_factor - sponge_lower_factor) *
          (1.0_rt - std::cos(M_PI * (rho - sponge_upper_density) / delta_rho));

      } else {
        sponge_factor = sponge_upper_factor;
      }
    }

    // Apply pressure sponge. This sponge is applied only if set by the user.

    // Note that because we do this third, the pressure sponge gets priority
    // over the radial and density sponges in cases where the two would overlap.

    if (sponge_upper_pressure > 0.0_rt && sponge_lower_pressure >= 0.0_rt) {

      eos_rep_t eos_state;

      eos_state.rho = state_in(i,j,k,URHO);
      eos_state.T = state_in(i,j,k,UTEMP);
      for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state_in(i,j,k,UFS+n) * rhoInv;
      }
#if NAUX_NET > 0
      for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = state_in(i,j,k,UFX+n) * rhoInv;
      }
#endif

      eos(eos_input_rt, eos_state);

      Real p = eos_state.p;

      if (p > sponge_upper_pressure) {
        sponge_factor = sponge_lower_factor;

      } else if (p <= sponge_upper_pressure && p >= sponge_lower_pressure) {
        sponge_factor = sponge_lower_factor +
          0.5_rt * (sponge_upper_factor - sponge_lower_factor) *
          (1.0_rt - std::cos(M_PI * (p - sponge_upper_pressure) / delta_p));

      } else {
        sponge_factor = sponge_upper_factor;

      }
    }

    // For an explicit update (sponge_implicit /= 1), the source term is given by
    // -(rho v) * alpha * sponge_factor. We simply add this directly by using the
    // current value of the momentum.

    // For an implicit update (sponge_implicit == 1), we choose the (rho v) to be
    // the momentum after the update. This then leads to an update of the form
    // (rho v) --> (rho v) * ONE / (ONE + alpha * sponge_factor). To get an equivalent
    // explicit form of this source term, we can then solve
    //    (rho v) + Sr == (rho v) / (ONE + alpha * sponge_factor),
    // which yields Sr = - (rho v) * (ONE - ONE / (ONE + alpha * sponge_factor)).

    Real fac;
    if (sponge_implicit == 1) {
       fac = -(1.0_rt - 1.0_rt / (1.0_rt + alpha * sponge_factor));

    } else {
       fac = -alpha * sponge_factor;

    }


    // now compute the source
    GpuArray<Real, 3> Sr;
    GpuArray<Real, 3> sponge_target_velocity = {sponge_target_x_velocity,
                                                sponge_target_y_velocity,
                                                sponge_target_z_velocity};
    for (int n = 0; n < 3; n++) {
      Sr[n] = (state_in(i,j,k,UMX+n) - rho * sponge_target_velocity[n]) * fac / dt;
      src[UMX+n] = Sr[n];
    }

    // Kinetic energy is 1/2 rho u**2, or (rho u)**2 / (2 rho). This means
    // that d(KE)/dt = u d(rho u)/dt - 1/2 u**2 d(rho)/dt. In this case
    // the sponge has no contribution to rho, so the kinetic energy source
    // term, and thus the total energy source term, is u * momentum source.

    Real SrE = 0.0;
    for (int n = 0; n < 3; n++) {
      SrE += state_in(i,j,k,UMX+n) * rhoInv * Sr[n];
    }

    src[UEDEN] = SrE;

#ifdef HYBRID_MOMENTUM
    GpuArray<Real, 3> Sr_hybrid;
    set_hybrid_momentum_source(r, Sr, Sr_hybrid);
    for (int n = 0; n < 3; n++) {
      src[UMR+n] = Sr_hybrid[n];
    }
#endif

    // Add terms to the source array.
    for (int n = 0; n < NSRC; n++) {
      source(i,j,k,n) += src[n];
    }

  });
}


#endif
