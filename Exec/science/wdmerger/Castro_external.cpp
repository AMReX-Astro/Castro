#include <Castro.H>
#include <Castro_F.H>
#include <Castro_ext_src.H>

using namespace amrex;

void
Castro::construct_old_ext_source(MultiFab& source, MultiFab& state_in, Real time, Real dt)
{
    const Real strt_time = ParallelDescriptor::second();

    if (!add_ext_src) return;

    MultiFab ext_src(grids, dmap, source.nComp(), 0);

    ext_src.setVal(0.0);

    fill_ext_source(time, dt, state_in, state_in, ext_src);

    Real mult_factor = 1.0;

    MultiFab::Saxpy(source, mult_factor, ext_src, 0, 0, source.nComp(), 0);

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

    // In this routine, we have two options: we can either do an
    // explicit predictor-corrector solve, or an implicit solve.
    // If we are doing predictor-corrector, the goal is to have
    // 0.5 * source_old + 0.5 * source_new. If we are doing an
    // implicit solve, the goal is just to have 1.0 * source_new.
    // This choice informs how we select mult_factor below.
    // It is up to the implementer of the external source term
    // to get it right if they choose to do an implicit solve.

    Real mult_factor;

    MultiFab ext_src(grids, dmap, source.nComp(), 0);

    ext_src.setVal(0.0);

    // Subtract off the old-time value first.

    Real old_time = time - dt;

    if (ext_src_implicit) {
        mult_factor = -1.0;
    } else {
        mult_factor = -0.5;
    }

    fill_ext_source(old_time, dt, state_old, state_old, ext_src);

    MultiFab::Saxpy(source, mult_factor, ext_src, 0, 0, source.nComp(), 0);

    // Time center with the new data.

    ext_src.setVal(0.0);

    if (ext_src_implicit) {
        mult_factor = 1.0;
    } else {
        mult_factor = 0.5;
    }

    fill_ext_source(time, dt, state_old, state_new, ext_src);

    MultiFab::Saxpy(source, mult_factor, ext_src, 0, 0, source.nComp(), 0);

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
Castro::fill_ext_source (const Real time, const Real dt, const MultiFab& state_old, const MultiFab& state_new, MultiFab& ext_src)
{
    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    GeometryData geomdata = geom.data();

    GpuArray<Real, 3> center;
    ca_get_center(center.begin());

    GpuArray<Real, 3> omega;
    get_omega(omega.begin());

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(ext_src, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        Array4<Real const> const sold = state_old.array(mfi);
        Array4<Real const> const snew = state_new.array(mfi);
        Array4<Real> const src = ext_src.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
        {
            do_ext_src(i, j, k, geomdata, snew, src, center, omega, dt, time);
        });

#pragma gpu box(bx)
        ca_ext_src
          (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
           BL_TO_FORTRAN_ANYD(state_old[mfi]),
           BL_TO_FORTRAN_ANYD(state_new[mfi]),
           BL_TO_FORTRAN_ANYD(ext_src[mfi]),
           AMREX_REAL_ANYD(prob_lo), AMREX_REAL_ANYD(dx), time, dt);
    }
}

