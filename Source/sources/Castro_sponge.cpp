#ifdef SPONGE
#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::construct_old_sponge_source(MultiFab& source, MultiFab& state_in, Real time, Real dt)
{
    const Real strt_time = ParallelDescriptor::second();

    if (!do_sponge) return;

    update_sponge_params(&time);

    const Real mult_factor = 1.0;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state_in, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        do_sponge(bx, state_in.array(mfi), source.array(mfi), dt, mult_factor);

    }

    if (verbose > 1)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Castro::construct_old_thermo_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }
}

void
Castro::construct_new_sponge_source(MultiFab& source, MultiFab& state_old, MultiFab& state_new, Real time, Real dt)
{
    const Real strt_time = ParallelDescriptor::second();

    if (!do_sponge) return;

    const Real mult_factor_old = -0.5;
    const Real mult_factor_new =  0.5;

    // First, subtract half of the old-time source.
    // Note that the sponge parameters are still current
    // at this point from their evaluation at the old time.

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state_old, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        do_sponge(bx, state_old.array(mfi), source.array(mfi), dt, mult_factor_old);
    }

    // Now update to the new-time sponge parameter values
    // and then evaluate the new-time part of the corrector.

    update_sponge_params(&time);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        do_sponge(bx, state_new.array(mfi), source.array(mfi), dt, mult_factor_new);

    }

    if (verbose > 1)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Castro::construct_new_thermo_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }
}

void
Castro::sponge_init()
{
    ca_allocate_sponge_params();
}

void
Castro::sponge_finalize()
{
    ca_deallocate_sponge_params();
}

void
Castro::do_sponge(const Box& bx,
                  Array4<Real const> const state,
                  Array4<Real> const source,
                  Real dt, Real mult_factor) {

  // alpha is a dimensionless measure of the timestep size; if
  // sponge_timescale < dt, then the sponge will have a larger effect,
  // and if sponge_timescale > dt, then the sponge will have a diminished effect.
  Real alpha;
  if (sponge_timescale > 0.0_rt) {
    alpha = dt / sponge_timescale;
  } else {
    alpha = 0.0_rt;
  }

  auto dx = geomdata.CellSizeArray();

  AMREX_PARALLEL_FOR_3D(bx, i, j, k,
  {

    Real src(NSRC);

    for (int n = 0; n < NSRC; n++) {
      src[n] = 0.0;
    }

    GpuArray<Real, 3> r;

    r[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - center[0];
    r[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - center[1];
    r[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - center[2];

    Real rho = state(i,j,k,URHO);
    Real rhoInv = ONE / rho;

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
      Real radius = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

      if (radius < sponge_lower_radius) {
        sponge_factor = sponge_lower_factor;

      } else if (radius >= sponge_lower_radius && radius <= sponge_upper_radius) {
        sponge_factor = sponge_lower_factor +
          0.5_rt * (sponge_upper_factor - sponge_lower_factor) *
          (1.0_rt - std::cos(M_PI * (radius - sponge_lower_radius) / delta_r));

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

      eos_t eos_state;

      eos_state.rho = state(i,j,k,URHO);
      eos_state.T = state(i,j,k,UTEMP);
      for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state(i,j,k,UFS+n) * rhoInv;
      }
      for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = state(i,j,k,UFX+n) * rhoInv;
      }

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

    if (sponge_implicit == 1) {
       fac = -(1.0_rt - 1.0_rt / (1.0_rt + alpha * sponge_factor));

    } else {
       fac = -alpha * sponge_factor;

    }


    // now compute the source
    GpuArray<Real, 3> Sr;
    for (int i = 0; i < 3; i++) {
      Sr[i] = (state(i,j,k,UMX+i) - rho * sponge_target_velocity) * fac * mult_factor / dt;
      src[UMX+i] = Sr[i];
    }

    // Kinetic energy is 1/2 rho u**2, or (rho u)**2 / (2 rho). This means
    // that d(KE)/dt = u d(rho u)/dt - 1/2 u**2 d(rho)/dt. In this case
    // the sponge has no contribution to rho, so the kinetic energy source
    // term, and thus the total energy source term, is u * momentum source.

    Real SrE = 0.0;
    for (int i = 0; i < 3; i++) {
      SrE += state(i,j,k,UMX+i) * rhoInv * Sr[i];
    }

    src[UEDEN] = SrE;

#ifdef HYBRID_MOMENTUM
    GpuArray<Real, 3> Sr_hybrid;
    set_hybrid_momentum_source(r, Sr, Sr_hybrid);
    for (int i = 0; i < 3; i++) {
      src[UMR+i] = Sr_hybrid[i];
    }
#endif

    // Add terms to the source array.
    for (int n = 0; n < NSRC; n++) {
      source(i,j,k,n) += src[n];
    }
  });
}


#endif
