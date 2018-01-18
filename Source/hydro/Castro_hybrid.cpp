#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::construct_old_hybrid_source(MultiFab& source, MultiFab& state, Real time, Real dt)
{
    Real mult_factor = 1.0;

    fill_hybrid_hydro_source(source, state, mult_factor);
}



void
Castro::construct_new_hybrid_source(MultiFab& source, MultiFab& state_old, MultiFab& state_new, Real time, Real dt)
{
    // Start by subtracting off the old-time data.

    Real mult_factor = -0.5;

    fill_hybrid_hydro_source(source, state_old, mult_factor);

    // Time center with the new data.

    mult_factor = 0.5;

    fill_hybrid_hydro_source(source, state_new, mult_factor);
}



void
Castro::fill_hybrid_hydro_source(MultiFab& sources, MultiFab& state, Real mult_factor)
{

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(state, true); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.tilebox();

    ca_hybrid_hydro_source(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			   BL_TO_FORTRAN_3D(state[mfi]),
			   BL_TO_FORTRAN_3D(sources[mfi]),
                           mult_factor);

  }

}



void
Castro::hybrid_sync(MultiFab& state)
{

    if (hybrid_hydro) {

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(state, true); mfi.isValid(); ++mfi) {

	    const Box& bx = mfi.tilebox();

	    ca_hybrid_update(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()), BL_TO_FORTRAN_3D(state[mfi]));

	}

    }

}
