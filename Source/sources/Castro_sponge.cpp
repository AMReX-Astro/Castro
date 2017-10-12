#ifdef SPONGE
#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::construct_old_sponge_source(Real time, Real dt)
{
    int ng = Sborder.nGrow();

    old_sources[sponge_src]->setVal(0.0);

    if (!do_sponge) return;

    update_sponge_params(&time);

    const Real *dx = geom.CellSize();

    const Real mult_factor = 1.0;

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
		  ZFILL(dx), dt, &time, mult_factor);

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

    const Real mult_factor_old = -0.5;
    const Real mult_factor_new =  0.5;

    // First, subtract half of the old-time source.
    // Note that the sponge parameters are still current
    // at this point from their evaluation at the old time.

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();

	ca_sponge(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                  BL_TO_FORTRAN_3D(S_old[mfi]),
		  BL_TO_FORTRAN_3D((*new_sources[sponge_src])[mfi]),
		  BL_TO_FORTRAN_3D(volume[mfi]),
		  ZFILL(dx), dt, &time, mult_factor_old);

    }

    // Now update to the new-time sponge parameter values
    // and then evaluate the new-time part of the corrector.

    update_sponge_params(&time);

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
		  ZFILL(dx), dt, &time, mult_factor_new);

    }

}
#endif
