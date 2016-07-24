#include "Castro.H"
#include "Castro_F.H"

void
Castro::construct_old_hybrid_source(Real time, Real dt)
{
    int ng = S_old.nGrow();

    fill_hybrid_hydro_source(old_sources[hybrid_src], (*Sborder));

    // Add to the hydro source terms.

    MultiFab::Add(*sources_for_hydro,old_sources[hybrid_src],0,0,NUM_STATE,ng);

}



void
Castro::construct_new_hybrid_source(Real time, Real dt)
{
    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    int ng = S_new.nGrow();

    fill_hybrid_hydro_source(new_sources[hybrid_src], S_new);

    // Time center the source term.

    old_sources[hybrid_src].mult(-0.5);
    new_sources[hybrid_src].mult( 0.5);

    MultiFab::Add(new_sources[hybrid_src],old_sources[hybrid_src],0,0,NUM_STATE,ng);

    old_sources[hybrid_src].mult(-2.0);

    // Add to the hydro source terms.

    if (source_term_predictor == 1)
        MultiFab::Add(*sources_for_hydro,old_sources[hybrid_src],0,0,NUM_STATE,ng);

}



void
Castro::fill_hybrid_hydro_source(MultiFab& sources, MultiFab& state)
{
  int ng = state.nGrow();

  BL_ASSERT(sources.nGrow() >= ng);

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(state, true); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.growntilebox(ng);

    ca_hybrid_hydro_source(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			   BL_TO_FORTRAN_3D(state[mfi]),
			   BL_TO_FORTRAN_3D(sources[mfi]));

  }

}
