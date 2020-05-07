#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::construct_old_geom_source(MultiFab& source, MultiFab& state_in, Real time, Real dt)
{

  if (geom.Coord() != 1) {
    amrex::Abort("this is only defined for axisymmetric geometry");
  }

  const Real strt_time = ParallelDescriptor::second();

  MultiFab geom_src(grids, dmap, source.nComp(), 0);

  geom_src.setVal(0.0);

  fill_geom_source(time, dt, state_in, state_in, geom_src);

  Real mult_factor = 1.0;

  MultiFab::Saxpy(source, mult_factor, geom_src, 0, 0, source.nComp(), 0);

  if (verbose > 1)
  {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
      ParallelDescriptor::ReduceRealMax(run_time,IOProc);

      if (ParallelDescriptor::IOProcessor())
          std::cout << "Castro::construct_old_geom_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
      });
#endif
    }

}



void
Castro::construct_new_geom_source(MultiFab& source, MultiFab& state_old, MultiFab& state_new, Real time, Real dt)
{

  if (geom.Coord() != 1) {
    amrex::Abort("this is only defined for axisymmetric geometry");
  }

  const Real strt_time = ParallelDescriptor::second();
 
  MultiFab geom_src(grids, dmap, source.nComp(), 0);

  geom_src.setVal(0.0);

  //Substract off the old-time value first
  Real old_time = time - dt;

  fill_geom_source(old_time, dt, state_old, state_old, geom_src);

  Real mult_factor = -0.5;

  MultiFab::Saxpy(source, mult_factor, geom_src, 0, 0, source.nComp(), 0);

  //Time center with the new data
  
  geom_src.setVal(0.0);

  mult_factor = 0.5;
  
  fill_geom_source(time, dt, state_old, state_new, geom_src);
  
  MultiFab::Saxpy(source, mult_factor, geom_src, 0, 0, source.nComp(), 0);

  if (verbose > 1)
  {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
      ParallelDescriptor::ReduceRealMax(run_time,IOProc);

      if (ParallelDescriptor::IOProcessor())
          std::cout << "Castro::construct_new_geom_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
      });
#endif
    }


}



void
Castro::fill_geom_source (Real time, Real dt,
                          MultiFab& state_old, MultiFab& state_new,
                          MultiFab& geom_src)
{

  // Compute the geometric source for axisymmetric coordinates
  // resulting from taking the divergence of (rho U U) in cylindrical
  // coordinates.  See the paper by Bernard-Champmartin

  auto dx = geom.CellSizeArray();
  auto prob_lo = geom.ProbLoArray();

  auto coord = geom.Coord(); 


#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(geom_src, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.tilebox();

    Array4<Real const> const old_state = state_old.array(mfi);
    Array4<Real const> const new_state = state_new.array(mfi);
    Array4<Real> const src = geom_src.array(mfi);

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
    {

      // radius for non-Cartesian
      Real r = prob_lo[0] + (static_cast<Real>(i) + 0.5_rt)*dx[0];

      Real rho_half = 0.5_rt*(old_state(i,j,k,URHO) + new_state(i,j,k,URHO));

      // radial momentum: F = rho v_phi**2 / r
      src(i,j,k,UMX) = 0.25_rt * (old_state(i,j,k,UMZ) + new_state(i,j,k,UMZ)) *
                                 (old_state(i,j,k,UMZ) + new_state(i,j,k,UMZ)) / (rho_half * r);

      // azimuthal momentum: F = - rho v_r v_phi / r
      src(i,j,k,UMZ) = -0.25_rt * (old_state(i,j,k,UMX) + new_state(i,j,k,UMX)) *
                                  (old_state(i,j,k,UMZ) + new_state(i,j,k,UMZ)) / (rho_half * r);

    });
  }
}
