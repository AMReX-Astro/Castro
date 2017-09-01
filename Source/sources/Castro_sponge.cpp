#ifdef SPONGE
#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::construct_old_sponge_source(Real time, Real dt)
{
    int ng = Sborder.nGrow();

    old_sources[sponge_src]->setVal(0.0);

    if (!time_center_sponge || !do_sponge) return;

    Real mult_factor = 1.0;

    update_sponge_params(&time);

    const Real *dx = geom.CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(Sborder,true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();

	ca_sponge(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		  BL_TO_FORTRAN_3D(Sborder[mfi]),
		  BL_TO_FORTRAN_3D((*old_sources[sponge_src])[mfi]),
		  BL_TO_FORTRAN_3D(volume[mfi]),
		  ZFILL(dx), dt, time, mult_factor);

    }

}

void
Castro::construct_new_sponge_source(Real time, Real dt)
{
    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    int ng = 0;

    new_sources[sponge_src]->setVal(0.0);

    if (!do_sponge) return;

    const Real *dx = geom.CellSize();

    // For the time centered version, do the old-time update first.
    // This way we can rely on the sponge parameter values which have
    // already been computed earlier in the step for the old-time source.

    Real mult_factor;

    if (time_center_sponge) {

        Real old_time = time - dt;
        mult_factor = -0.5;

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(S_old,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

            ca_sponge(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                      BL_TO_FORTRAN_3D(S_old[mfi]),
                      BL_TO_FORTRAN_3D((*new_sources[sponge_src])[mfi]),
                      BL_TO_FORTRAN_3D(volume[mfi]),
                      ZFILL(dx), dt, old_time, mult_factor);

        }

    }

    // Now update the sponge parameters in anticipation of the new-time source.

    update_sponge_params(&time);

    if (time_center_sponge) {
        mult_factor = 0.5;
    } else {
        mult_factor = 1.0;
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        ca_sponge(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                  BL_TO_FORTRAN_3D(S_new[mfi]),
                  BL_TO_FORTRAN_3D((*new_sources[sponge_src])[mfi]),
                  BL_TO_FORTRAN_3D(volume[mfi]),
                  ZFILL(dx), dt, time, mult_factor);

    }

}
#endif
