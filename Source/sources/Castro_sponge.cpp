#ifdef SPONGE
#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::construct_old_sponge_source(Real time, Real dt)
{
    MultiFab& S_new = get_new_data(State_Type);

    int ng = Sborder.nGrow();

    old_sources[sponge_src]->setVal(0.0);

    update_sponge_params(&time);

    const Real *dx = geom.CellSize();

    int is_new = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(Sborder,true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();

	ca_sponge(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		  BL_TO_FORTRAN_3D(Sborder[mfi]),
                  BL_TO_FORTRAN_3D(S_new[mfi]),
		  BL_TO_FORTRAN_3D((*old_sources[sponge_src])[mfi]),
		  BL_TO_FORTRAN_3D(volume[mfi]),
		  ZFILL(dx), dt, &time, is_new);

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

    update_sponge_params(&time);

    const Real *dx = geom.CellSize();

    int is_new = 1;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();

	ca_sponge(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                  BL_TO_FORTRAN_3D(S_old[mfi]),
		  BL_TO_FORTRAN_3D(S_new[mfi]),
		  BL_TO_FORTRAN_3D((*new_sources[sponge_src])[mfi]),
		  BL_TO_FORTRAN_3D(volume[mfi]),
		  ZFILL(dx), dt, &time, is_new);

    }

}
#endif
