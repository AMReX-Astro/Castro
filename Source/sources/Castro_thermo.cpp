#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::construct_old_thermo_source(MultiFab& source, MultiFab& state_in, Real time, Real dt)
{
  if (!(time_integration_method == SpectralDeferredCorrections)) return;

  const Real strt_time = ParallelDescriptor::second();

  MultiFab thermo_src(grids, dmap, source.nComp(), 0);

  thermo_src.setVal(0.0);

  fill_thermo_source(time, dt, state_in, state_in, thermo_src);

  Real mult_factor = 1.0;

  MultiFab::Saxpy(source, mult_factor, thermo_src, 0, 0, source.nComp(), 0);

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
Castro::construct_new_thermo_source(MultiFab& source, MultiFab& state_old, MultiFab& state_new, Real time, Real dt)
{
  if (!(time_integration_method == SpectralDeferredCorrections)) return;

  amrex::Abort("you should not get here!");
}



void
Castro::fill_thermo_source (Real time, Real dt,
                            MultiFab& state_old, MultiFab& state_new,
                            MultiFab& thermo_src)
{
  const Real* dx = geom.CellSize();
  const Real* prob_lo = geom.ProbLo();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(thermo_src, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

      const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
      ca_thermo_src(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                    BL_TO_FORTRAN_ANYD(state_old[mfi]),
                    BL_TO_FORTRAN_ANYD(state_new[mfi]),
                    BL_TO_FORTRAN_ANYD(thermo_src[mfi]),
                    AMREX_REAL_ANYD(prob_lo),AMREX_REAL_ANYD(dx),time,dt);
    }
}
