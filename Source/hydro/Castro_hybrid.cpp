#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::construct_old_hybrid_source(Real time, Real dt)
{
    int ng = Sborder.nGrow();

    Real mult_factor = 1.0;

    fill_hybrid_hydro_source(*old_sources[unified_src], Sborder, ng, mult_factor);
}



void
Castro::construct_new_hybrid_source(Real time, Real dt)
{
    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    int ng = 0;

    // Start by subtracting off the old-time data.

    Real mult_factor = -0.5;

    fill_hybrid_hydro_source(*new_sources[unified_src], S_old, ng, mult_factor);

    // Time center with the new data.

    mult_factor = 0.5;

    fill_hybrid_hydro_source(*new_sources[unified_src], S_new, ng, mult_factor);

}



void
Castro::fill_hybrid_hydro_source(MultiFab& sources, MultiFab& state, int ng, Real mult_factor)
{
  BL_ASSERT(sources.nGrow() >= ng);
  BL_ASSERT(state.nGrow() >= ng);

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(state, true); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.growntilebox(ng);

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
