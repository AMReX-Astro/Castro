#include <Castro.H>
#include <Castro_util.H>

using namespace amrex;

void
Castro::construct_old_rotation_source(MultiFab& source, MultiFab& state_in, Real time, Real dt)
{

    amrex::ignore_unused(time);

    BL_PROFILE("Castro::construct_old_rotation_source()");

    const Real strt_time = ParallelDescriptor::second();

    // Fill the rotation data.

    if (!do_rotation) {
        return;

    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state_in, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        rsrc(bx, state_in.array(mfi), source.array(mfi), dt);

    }

    if (verbose > 1)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        amrex::Real run_time = ParallelDescriptor::second() - strt_time;
        amrex::Real llevel = level;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        amrex::Print() << "Castro::construct_old_rotation_source() time = " << run_time
                       << " on level " << llevel << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }
}



void
Castro::construct_new_rotation_source(MultiFab& source, MultiFab& state_old, MultiFab& state_new, Real time, Real dt)
{

    amrex::ignore_unused(time);

    BL_PROFILE("Castro::construct_new_rotation_source()");

    const Real strt_time = ParallelDescriptor::second();

    // Fill the rotation data.

    if (!do_rotation) {
        return;

    }

    // Now do corrector part of rotation source term update

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        for (MFIter mfi(state_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

            corrrsrc(bx,
                     state_old.array(mfi), state_new.array(mfi),
                     source.array(mfi),
                     (*mass_fluxes[0]).array(mfi), (*mass_fluxes[1]).array(mfi), (*mass_fluxes[2]).array(mfi),
                     dt, volume.array(mfi));
        }
    }

    if (verbose > 1)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;
        Real llevel = level;
#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        amrex::Print() << "Castro::construct_new_rotation_source() time = " << run_time << " on level " << llevel << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }
}
