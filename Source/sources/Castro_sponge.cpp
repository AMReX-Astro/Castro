#ifdef SPONGE
#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::construct_old_sponge_source(MultiFab& source, MultiFab& state, Real time, Real dt)
{
    const Real strt_time = ParallelDescriptor::second();

    if (!do_sponge) return;

    update_sponge_params(&time);

    const Real *dx = geom.CellSize();

    const Real mult_factor = 1.0;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
        ca_sponge(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                  BL_TO_FORTRAN_ANYD(state[mfi]),
                  BL_TO_FORTRAN_ANYD(source[mfi]),
                  BL_TO_FORTRAN_ANYD(volume[mfi]),
                  AMREX_REAL_ANYD(dx), dt, time, mult_factor);
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

    const Real *dx = geom.CellSize();

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

#pragma gpu box(bx)
        ca_sponge(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                  BL_TO_FORTRAN_ANYD(state_old[mfi]),
                  BL_TO_FORTRAN_ANYD(source[mfi]),
                  BL_TO_FORTRAN_ANYD(volume[mfi]),
                  AMREX_REAL_ANYD(dx), dt, time, mult_factor_old);
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

#pragma gpu box(bx)
        ca_sponge(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                  BL_TO_FORTRAN_ANYD(state_new[mfi]),
                  BL_TO_FORTRAN_ANYD(source[mfi]),
                  BL_TO_FORTRAN_ANYD(volume[mfi]),
                  AMREX_REAL_ANYD(dx), dt, time, mult_factor_new);
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

#endif
