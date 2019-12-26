#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::construct_old_hybrid_source(MultiFab& source, MultiFab& state, Real time, Real dt)
{
    BL_PROFILE("Castro::construct_old_hybrid_source()");

    const Real strt_time = ParallelDescriptor::second();

    Real mult_factor = 1.0;

    fill_hybrid_hydro_source(source, state, mult_factor);

    if (verbose > 1)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Castro::construct_old_hybrid_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }
}



void
Castro::construct_new_hybrid_source(MultiFab& source, MultiFab& state_old, MultiFab& state_new, Real time, Real dt)
{
    BL_PROFILE("Castro::construct_new_hybrid_source()");

    const Real strt_time = ParallelDescriptor::second();

    // Start by subtracting off the old-time data.

    Real mult_factor = -0.5;

    fill_hybrid_hydro_source(source, state_old, mult_factor);

    // Time center with the new data.

    mult_factor = 0.5;

    fill_hybrid_hydro_source(source, state_new, mult_factor);

    if (verbose > 1)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Castro::construct_new_hybrid_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }
}



void
Castro::fill_hybrid_hydro_source(MultiFab& sources, MultiFab& state, Real mult_factor)
{
    BL_PROFILE("Castro::fill_hybrid_hydro_source()");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
        ca_hybrid_hydro_source(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                               BL_TO_FORTRAN_ANYD(state[mfi]),
                               BL_TO_FORTRAN_ANYD(sources[mfi]),
                               mult_factor);

    }
}



void
Castro::linear_to_hybrid_momentum(MultiFab& state, int ng)
{
    BL_PROFILE("Castro::linear_to_hybrid_momentum()");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(ng);

#pragma gpu box(bx)
        ca_linear_to_hybrid_momentum(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()), BL_TO_FORTRAN_ANYD(state[mfi]));
    }
}



void
Castro::hybrid_to_linear_momentum(MultiFab& state, int ng)
{
    BL_PROFILE("Castro::hybrid_to_linear_momentum()");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(ng);

#pragma gpu box(bx)
        ca_hybrid_to_linear_momentum(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()), BL_TO_FORTRAN_ANYD(state[mfi]));
    }
}
