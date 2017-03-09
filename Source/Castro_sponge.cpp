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
		  ZFILL(dx), dt, &time);

    }

}

void
Castro::construct_new_sponge_source(Real time, Real dt)
{
    MultiFab& S_new = get_new_data(State_Type);

    int ng = 0;

    new_sources[sponge_src]->setVal(0.0);

    if (!do_sponge) return;

    update_sponge_params(&time);

    const Real *dx = geom.CellSize();

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
		  ZFILL(dx), dt, &time);

    }

    // Time center the source term.

    if (time_center_sponge) {
	new_sources[sponge_src]->mult(0.5);

	MultiFab::Saxpy(*new_sources[sponge_src],-0.5,*old_sources[sponge_src],0,0,NUM_STATE,ng);
    }

}
#endif
