#include <Castro.H>
#include <Castro_util.H>
#include <Castro_F.H>

#include <hybrid.H>

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

    GeometryData geomdata = geom.data();

    GpuArray<Real, 3> center;
    ca_get_center(center.begin());

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto u = state.array(mfi);
        auto src = sources.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            GpuArray<Real, 3> loc;

            position(i, j, k, geomdata, loc);

            loc[0] -= center[0];
            loc[1] -= center[1];

            Real R = amrex::max(std::sqrt(loc[0] * loc[0] + loc[1] * loc[1]), R_min);

            Real rhoInv = 1.0_rt / u(i,j,k,URHO);
            Real RInv = 1.0_rt / R;

            src(i,j,k,UMR) = src(i,j,k,UMR) + mult_factor * (rhoInv * RInv * RInv * RInv) *
                                              u(i,j,k,UML) * u(i,j,k,UML);

        });
    }
}



void
Castro::linear_to_hybrid_momentum(MultiFab& state, int ng)
{
    BL_PROFILE("Castro::linear_to_hybrid_momentum()");

    GeometryData geomdata = geom.data();

    GpuArray<Real, 3> center;
    ca_get_center(center.begin());

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(ng);

        auto u = state.array(mfi);

        // Convert linear momentum to hybrid momentum.

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            GpuArray<Real, 3> loc;

            position(i, j, k, geomdata, loc);

            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
                loc[dir] -= center[dir];

            GpuArray<Real, 3> linear_mom;

            for (int dir = 0; dir < 3; ++dir)
                linear_mom[dir] = u(i,j,k,UMX+dir);

            GpuArray<Real, 3> hybrid_mom;

            linear_to_hybrid(loc, linear_mom, hybrid_mom);

            for (int dir = 0; dir < 3; ++dir)
                u(i,j,k,UMR+dir) = hybrid_mom[dir];

        });
    }
}



void
Castro::hybrid_to_linear_momentum(MultiFab& state, int ng)
{
    BL_PROFILE("Castro::hybrid_to_linear_momentum()");

    GeometryData geomdata = geom.data();

    GpuArray<Real, 3> center;
    ca_get_center(center.begin());

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(ng);

        auto u = state.array(mfi);

        // Convert hybrid momentum to linear momentum.

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            GpuArray<Real, 3> loc;

            position(i, j, k, geomdata, loc);

            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
                loc[dir] -= center[dir];

            GpuArray<Real, 3> hybrid_mom;

            for (int dir = 0; dir < 3; ++dir)
                hybrid_mom[dir] = u(i,j,k,UMR+dir);

            GpuArray<Real, 3> linear_mom;

            hybrid_to_linear(loc, hybrid_mom, linear_mom);

            for (int dir = 0; dir < 3; ++dir)
                u(i,j,k,UMX+dir) = linear_mom[dir];

        });
    }
}
