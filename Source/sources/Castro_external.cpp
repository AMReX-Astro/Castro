#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::construct_old_ext_source(MultiFab& source, MultiFab& state, Real time, Real dt)
{
    const Real strt_time = ParallelDescriptor::second();

    if (!add_ext_src) return;

    MultiFab ext_src(grids, dmap, NUM_STATE, 0);

    ext_src.setVal(0.0);

    fill_ext_source(time, dt, state, state, ext_src);

    Real mult_factor = 1.0;

    MultiFab::Saxpy(source, mult_factor, ext_src, 0, 0, NUM_STATE, 0);

    if (verbose > 1)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Castro::construct_old_ext_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }
}



void
Castro::construct_new_ext_source(MultiFab& source, MultiFab& state_old, MultiFab& state_new, Real time, Real dt)
{
    const Real strt_time = ParallelDescriptor::second();

    if (!add_ext_src) return;

    MultiFab ext_src(grids, dmap, NUM_STATE, 0);

    ext_src.setVal(0.0);

    // Subtract off the old-time value first.

    Real old_time = time - dt;
    Real mult_factor = -0.5;

    fill_ext_source(old_time, dt, state_old, state_old, ext_src);

    MultiFab::Saxpy(source, mult_factor, ext_src, 0, 0, NUM_STATE, 0);

    // Time center with the new data.

    ext_src.setVal(0.0);

    mult_factor = 0.5;

    fill_ext_source(time, dt, state_old, state_new, ext_src);

    MultiFab::Saxpy(source, mult_factor, ext_src, 0, 0, NUM_STATE, 0);

    if (verbose > 1)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Castro::construct_new_ext_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }
}



void
Castro::fill_ext_source (Real time, Real dt, MultiFab& state_old, MultiFab& state_new, MultiFab& ext_src)
{
    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(ext_src,true); mfi.isValid(); ++mfi)
    {

        const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
        ca_ext_src
	  (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
	   BL_TO_FORTRAN_ANYD(state_old[mfi]),
	   BL_TO_FORTRAN_ANYD(state_new[mfi]),
	   BL_TO_FORTRAN_ANYD(ext_src[mfi]),
	   AMREX_REAL_ANYD(prob_lo), AMREX_REAL_ANYD(dx), time, dt);
    }
}
