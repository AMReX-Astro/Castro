#include <Castro.H>
#include <Castro_util.H>

#include <hybrid.H>

using namespace amrex;

void
Castro::construct_old_hybrid_source(MultiFab& source, MultiFab& state_old, Real time, Real dt)
{

    amrex::ignore_unused(time);
    amrex::ignore_unused(dt);

    BL_PROFILE("Castro::construct_old_hybrid_source()");

    const Real strt_time = ParallelDescriptor::second();

    Real mult_factor = 1.0;

    fill_hybrid_hydro_source(source, state_old, mult_factor);

    if (verbose > 1)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        amrex::Real run_time = ParallelDescriptor::second() - strt_time;

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

    amrex::ignore_unused(time);
    amrex::ignore_unused(dt);

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
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        amrex::Real run_time = ParallelDescriptor::second() - strt_time;

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
Castro::fill_hybrid_hydro_source(MultiFab& sources, const MultiFab& state_in, Real mult_factor)
{
    BL_PROFILE("Castro::fill_hybrid_hydro_source()");

    GeometryData geomdata = geom.data();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state_in, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        const auto u = state_in.array(mfi);
        auto src = sources.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            GpuArray<Real, 3> loc;

            position(i, j, k, geomdata, loc);

            Real R = amrex::max(std::sqrt(loc[0] * loc[0] + loc[1] * loc[1]),
                                std::numeric_limits<Real>::min());

            Real rhoInv = 1.0_rt / u(i,j,k,URHO);
            Real RInv = 1.0_rt / R;

            src(i,j,k,UMR) = src(i,j,k,UMR) + mult_factor * (rhoInv * RInv * RInv * RInv) *
                                              u(i,j,k,UML) * u(i,j,k,UML);

        });
    }
}



void
Castro::linear_to_hybrid_momentum(MultiFab& state_in, int ng)
{
    BL_PROFILE("Castro::linear_to_hybrid_momentum()");

    GeometryData geomdata = geom.data();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state_in, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(ng);

        auto u = state_in.array(mfi);

        // Convert linear momentum to hybrid momentum.

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            GpuArray<Real, 3> loc;

            position(i, j, k, geomdata, loc);

            GpuArray<Real, 3> linear_mom;

            for (int dir = 0; dir < 3; ++dir) {
                linear_mom[dir] = u(i,j,k,UMX+dir);
            }

            GpuArray<Real, 3> hybrid_mom;

            linear_to_hybrid(loc, linear_mom, hybrid_mom);

            for (int dir = 0; dir < 3; ++dir) {
                u(i,j,k,UMR+dir) = hybrid_mom[dir];
            }
        });
    }
}



void
Castro::hybrid_to_linear_momentum(MultiFab& state_in, int ng)
{
    BL_PROFILE("Castro::hybrid_to_linear_momentum()");

    GeometryData geomdata = geom.data();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state_in, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(ng);

        auto u = state_in.array(mfi);

        // Convert hybrid momentum to linear momentum.

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            GpuArray<Real, 3> loc;

            position(i, j, k, geomdata, loc);

            GpuArray<Real, 3> hybrid_mom;

            for (int dir = 0; dir < 3; ++dir) {
                hybrid_mom[dir] = u(i,j,k,UMR+dir);
            }

            GpuArray<Real, 3> linear_mom;

            hybrid_to_linear(loc, hybrid_mom, linear_mom);

            for (int dir = 0; dir < 3; ++dir) {
                u(i,j,k,UMX+dir) = linear_mom[dir];
            }
        });
    }
}
