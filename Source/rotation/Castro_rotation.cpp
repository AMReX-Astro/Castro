
#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::construct_old_rotation_source(MultiFab& source, MultiFab& state, Real time, Real dt)
{

    BL_PROFILE("Castro::construct_old_rotation_source()");

    const Real strt_time = ParallelDescriptor::second();

    MultiFab& phirot_old = get_old_data(PhiRot_Type);
    MultiFab& rot_old = get_old_data(Rotation_Type);

    // Fill the rotation data.

    if (!do_rotation) {

        phirot_old.setVal(0.0);
        rot_old.setVal(0.0);

        return;

    }

    fill_rotation_field(phirot_old, rot_old, state, time);

    const Real *dx = geom.CellSize();
    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
#pragma gpu box(bx)
        ca_rsrc(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
                BL_TO_FORTRAN_ANYD(phirot_old[mfi]),
                BL_TO_FORTRAN_ANYD(rot_old[mfi]),
                BL_TO_FORTRAN_ANYD(state[mfi]),
                BL_TO_FORTRAN_ANYD(source[mfi]),
                BL_TO_FORTRAN_ANYD(volume[mfi]),
                AMREX_REAL_ANYD(dx),dt,time);

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
            std::cout << "Castro::construct_old_rotation_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }
}



void
Castro::construct_new_rotation_source(MultiFab& source, MultiFab& state_old, MultiFab& state_new, Real time, Real dt)
{
    BL_PROFILE("Castro::construct_new_rotation_source()");

    const Real strt_time = ParallelDescriptor::second();

    MultiFab& phirot_old = get_old_data(PhiRot_Type);
    MultiFab& rot_old = get_old_data(Rotation_Type);

    MultiFab& phirot_new = get_new_data(PhiRot_Type);
    MultiFab& rot_new = get_new_data(Rotation_Type);

    // Fill the rotation data.

    if (!do_rotation) {

        phirot_new.setVal(0.);
        rot_new.setVal(0.);

        return;

    }

    fill_rotation_field(phirot_new, rot_new, state_new, time);

    // Now do corrector part of rotation source term update

    const Real *dx = geom.CellSize();
    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        for (MFIter mfi(state_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
            ca_corrrsrc(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                        AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
                        BL_TO_FORTRAN_ANYD(phirot_old[mfi]),
                        BL_TO_FORTRAN_ANYD(phirot_new[mfi]),
                        BL_TO_FORTRAN_ANYD(rot_old[mfi]),
                        BL_TO_FORTRAN_ANYD(rot_new[mfi]),
                        BL_TO_FORTRAN_ANYD(state_old[mfi]),
                        BL_TO_FORTRAN_ANYD(state_new[mfi]),
                        BL_TO_FORTRAN_ANYD(source[mfi]),
                        BL_TO_FORTRAN_ANYD((*mass_fluxes[0])[mfi]),
                        BL_TO_FORTRAN_ANYD((*mass_fluxes[1])[mfi]),
                        BL_TO_FORTRAN_ANYD((*mass_fluxes[2])[mfi]),
                        AMREX_REAL_ANYD(dx),dt,time,
                        BL_TO_FORTRAN_ANYD(volume[mfi]));
        }
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
            std::cout << "Castro::construct_new_rotation_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }
}



void Castro::fill_rotation_field(MultiFab& phi, MultiFab& rot, MultiFab& state, Real time)
{

    BL_PROFILE("Castro::fill_rotation_field()");

    const Real* dx = geom.CellSize();

    phi.setVal(0.0);

    int ng = phi.nGrow();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        const Box& bx = mfi.growntilebox(ng);
#pragma gpu box(bx)
        ca_fill_rotational_potential(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                                     BL_TO_FORTRAN_ANYD(phi[mfi]),
                                     AMREX_REAL_ANYD(dx),time);

    }

    rot.setVal(0.0);

    ng = state.nGrow();

    if (ng > rot.nGrow())
        amrex::Error("State MF has more ghost cells than rotation MF.");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        const Box& bx = mfi.growntilebox(ng);
#pragma gpu box(bx)
        ca_fill_rotational_acceleration(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                                        BL_TO_FORTRAN_ANYD(rot[mfi]),
                                        BL_TO_FORTRAN_ANYD(state[mfi]),
                                        AMREX_REAL_ANYD(dx),time);

    }

}
