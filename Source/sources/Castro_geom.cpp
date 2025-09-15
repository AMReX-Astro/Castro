#include <Castro.H>
#include <Castro_math.H>

using namespace amrex;

///
/// this adds the geometric source terms for non-Cartesian Coordinates
/// This includes 2D Cylindrical (R-Z) coordinate as described in Bernand-Champmartin
/// and 2D Spherical (R-THETA) coordinate.
///

void
Castro::construct_old_geom_source(MultiFab& source, MultiFab& state_in, Real time, Real dt)
{

  amrex::ignore_unused(source);
  amrex::ignore_unused(state_in);
  amrex::ignore_unused(time);
  amrex::ignore_unused(dt);

  if (geom.Coord() == 0) {
      return;
  }

  if (use_geom_source == 0) {
      return;
  }

#if AMREX_SPACEDIM <= 2
  const Real strt_time = ParallelDescriptor::second();

  MultiFab geom_src(grids, dmap, source.nComp(), 0);

  geom_src.setVal(0.0);

  if (geom.Coord() == 1) {
      fill_RZ_geom_source(time, dt, state_in, geom_src);
  } else {
      fill_RTheta_geom_source(time, dt, state_in, geom_src);
  }

  Real mult_factor = 1.0;

  MultiFab::Saxpy(source, mult_factor, geom_src, 0, 0, source.nComp(), 0);  // NOLINT(readability-suspicious-call-argument)

  if (verbose > 1)
  {
      const int IOProc = ParallelDescriptor::IOProcessorNumber();
      amrex::Real run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
      ParallelDescriptor::ReduceRealMax(run_time,IOProc);

      amrex::Print() << "Castro::construct_old_geom_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
      });
#endif
    }

#endif
}



void
Castro::construct_new_geom_source(MultiFab& source, MultiFab& state_old,
                                  MultiFab& state_new, Real time, Real dt)
{

  amrex::ignore_unused(source);
  amrex::ignore_unused(state_old);
  amrex::ignore_unused(state_new);
  amrex::ignore_unused(time);
  amrex::ignore_unused(dt);

  if (geom.Coord() == 0) {
      return;
  }

  if (use_geom_source == 0) {
      return;
  }

#if AMREX_SPACEDIM <= 2
  const Real strt_time = ParallelDescriptor::second();

  MultiFab geom_src(grids, dmap, source.nComp(), 0);

  geom_src.setVal(0.0);

  // Subtract off the old-time value first
  Real old_time = time - dt;

  if (geom.Coord() == 1) {
      fill_RZ_geom_source(old_time, dt, state_old, geom_src);
  } else {
      fill_RTheta_geom_source(old_time, dt, state_old, geom_src);
  }

  Real mult_factor = -0.5;

  MultiFab::Saxpy(source, mult_factor, geom_src, 0, 0, source.nComp(), 0);   // NOLINT(readability-suspicious-call-argument)

  // Time center with the new data

  geom_src.setVal(0.0);

  mult_factor = 0.5;

  if (geom.Coord() == 1) {
      fill_RZ_geom_source(time, dt, state_new, geom_src);
  } else {
      fill_RTheta_geom_source(time, dt, state_new, geom_src);
  }

  MultiFab::Saxpy(source, mult_factor, geom_src, 0, 0, source.nComp(), 0);   // NOLINT(readability-suspicious-call-argument)

  if (verbose > 1)
  {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
      ParallelDescriptor::ReduceRealMax(run_time,IOProc);

      amrex::Print() << "Castro::construct_new_geom_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
      });
#endif
    }

#endif
}



void
Castro::fill_RZ_geom_source (Real time, Real dt, MultiFab& cons_state, MultiFab& geom_src)
{

  // Compute the geometric source for axisymmetric coordinates (R-Z)
  // resulting from taking the divergence of (rho U U) in cylindrical
  // coordinates.  See the paper by Bernard-Champmartin

  amrex::ignore_unused(time);
  amrex::ignore_unused(dt);

  auto dx = geom.CellSizeArray();
  auto prob_lo = geom.ProbLoArray();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(geom_src, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.tilebox();

    Array4<Real const> const U_arr = cons_state.array(mfi);

    Array4<Real> const src = geom_src.array(mfi);

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {

      // radius for non-Cartesian
      Real r = prob_lo[0] + (static_cast<Real>(i) + 0.5_rt)*dx[0];

      // radial momentum: F = rho v_phi**2 / r
      src(i,j,k,UMX) = U_arr(i,j,k,UMZ) * U_arr(i,j,k,UMZ) / (U_arr(i,j,k,URHO) * r);

      // azimuthal momentum: F = - rho v_r v_phi / r
      src(i,j,k,UMZ) = - U_arr(i,j,k,UMX) * U_arr(i,j,k,UMZ) / (U_arr(i,j,k,URHO) * r);

    });
  }
}



void
Castro::fill_RTheta_geom_source (Real time, Real dt, MultiFab& cons_state, MultiFab& geom_src)
{

  // Compute the geometric source resulting from taking the divergence of (rho U U)
  // in Spherical 2D (R-Theta) coordinate.

  amrex::ignore_unused(time);
  amrex::ignore_unused(dt);

  auto dx = geom.CellSizeArray();
  auto prob_lo = geom.ProbLoArray();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(geom_src, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.tilebox();

    Array4<Real const> const U_arr = cons_state.array(mfi);

    Array4<Real> const src = geom_src.array(mfi);

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {

      // Cell-centered Spherical Radius and Theta
      Real r = prob_lo[0] + (static_cast<Real>(i) + 0.5_rt)*dx[0];
      Real theta = prob_lo[1] + (static_cast<Real>(j) + 0.5_rt)*dx[1];

      // radial momentum: F = rho (v_theta**2 + v_phi**2) / r
      src(i,j,k,UMX) = (U_arr(i,j,k,UMY) * U_arr(i,j,k,UMY) +
                        U_arr(i,j,k,UMZ) * U_arr(i,j,k,UMZ)) / (U_arr(i,j,k,URHO) * r);

      // Theta momentum F = rho v_phi**2 cot(theta) / r - rho v_r v_theta / r
      src(i,j,k,UMY) = (U_arr(i,j,k,UMZ) * U_arr(i,j,k,UMZ) * cot(theta) -
                        U_arr(i,j,k,UMX) * U_arr(i,j,k,UMY)) / (U_arr(i,j,k,URHO) * r);

      // Phi momentum: F = - rho v_r v_phi / r - rho v_theta v_phi cot(theta) / r
      src(i,j,k,UMZ) = (- U_arr(i,j,k,UMY) * U_arr(i,j,k,UMZ) * cot(theta) -
                        U_arr(i,j,k,UMX) * U_arr(i,j,k,UMZ)) / (U_arr(i,j,k,URHO) * r);

    });
  }
}
