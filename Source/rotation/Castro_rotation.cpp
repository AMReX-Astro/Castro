
#include "Castro.H"
#include "Castro_F.H"
#include "Castro_util.H"

using namespace amrex;

void
Castro::construct_old_rotation_source(MultiFab& source, MultiFab& state_in, Real time, Real dt)
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

    fill_rotation_field(phirot_old, rot_old, state_in, time);

    const Real *dx = geom.CellSize();
    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state_in, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        rsrc(bx,
             phirot_old.array(mfi),
             rot_old.array(mfi),
             state_in.array(mfi),
             source.array(mfi),
             volume.array(mfi),
             dt, time);

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

            corrrsrc(bx,
                     phirot_old.array(mfi),
                     phirot_new.array(mfi),
                     rot_old.array(mfi),
                     rot_new.array(mfi),
                     state_old.array(mfi),
                     state_new.array(mfi),
                     source.array(mfi),
                     (*mass_fluxes[0]).array(mfi),
                     (*mass_fluxes[1]).array(mfi),
                     (*mass_fluxes[2]).array(mfi),
                     dt, time,
                     volume.array(mfi));
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



void Castro::fill_rotation_field(MultiFab& phi, MultiFab& rot, MultiFab& state_in, Real time)
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

        fill_rotational_potential(bx, phi.array(mfi), time);

    }

    rot.setVal(0.0);

    ng = state_in.nGrow();

    if (ng > rot.nGrow())
        amrex::Error("State MF has more ghost cells than rotation MF.");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state_in, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        const Box& bx = mfi.growntilebox(ng);

        fill_rotational_acceleration(bx, rot.array(mfi), state_in.array(mfi), time);

    }

}


AMREX_GPU_HOST_DEVICE 
void
Castro::inertial_to_rotational_velocity_c(const int i, const int j, const int k,
                                          const GeometryData& geomdata,
                                          const Real* center,
                                          const Real* omega,
                                          const Real time, Real* v) {

  // Given a velocity vector in the inertial frame, transform it to a
  // velocity vector in the rotating frame.

  // Note: this version assumes all cell-centers

  GpuArray<Real, 3> loc;

  position(i, j, k, geomdata, loc);

  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    loc[dir] -= center[dir];
  }

  // do the cross product Omega x loc
  v[0] += -(omega[1]*loc[2] - omega[2]*loc[1]);
  v[1] += -(omega[2]*loc[0] - omega[0]*loc[2]);
  v[2] += -(omega[0]*loc[1] - omega[1]*loc[0]);

}
