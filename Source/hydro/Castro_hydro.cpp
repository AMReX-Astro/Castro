#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

using namespace amrex;

void
Castro::cons_to_prim(const Real time)
{

#ifdef RADIATION
    AmrLevel::FillPatch(*this, Erborder, NUM_GROW, time, Rad_Type, 0, Radiation::nGroups);

    MultiFab lamborder(grids, dmap, Radiation::nGroups, NUM_GROW);
    if (radiation->pure_hydro) {
      lamborder.setVal(0.0, NUM_GROW);
    }
    else {
      radiation->compute_limiter(level, grids, Sborder, Erborder, lamborder);
    }
#endif

    MultiFab& S_new = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

        const Box& qbx = mfi.growntilebox(NUM_GROW);

        // Convert the conservative state to the primitive variable state.
        // This fills both q and qaux.

#pragma gpu
        ca_ctoprim(AMREX_INT_ANYD(qbx.loVect()), AMREX_INT_ANYD(qbx.hiVect()),
                   BL_TO_FORTRAN_ANYD(Sborder[mfi]),
#ifdef RADIATION
                   BL_TO_FORTRAN_ANYD(Erborder[mfi]),
                   BL_TO_FORTRAN_ANYD(lamborder[mfi]),
#endif
                   BL_TO_FORTRAN_ANYD(q[mfi]),
                   BL_TO_FORTRAN_ANYD(qaux[mfi]));

        // Convert the source terms expressed as sources to the conserved state to those
        // expressed as sources for the primitive state.
        if (time_integration_method == CornerTransportUpwind || time_integration_method == SimplifiedSpectralDeferredCorrections) {
#pragma gpu
            ca_srctoprim(BL_TO_FORTRAN_BOX(qbx),
                         BL_TO_FORTRAN_ANYD(q[mfi]),
                         BL_TO_FORTRAN_ANYD(qaux[mfi]),
                         BL_TO_FORTRAN_ANYD(sources_for_hydro[mfi]),
                         BL_TO_FORTRAN_ANYD(src_q[mfi]));
        }

#ifndef RADIATION
#ifdef REACTIONS
        // Add in the reactions source term; only done in simplified SDC.

        if (time_integration_method == SimplifiedSpectralDeferredCorrections) {

            MultiFab& SDC_react_source = get_new_data(Simplified_SDC_React_Type);

            if (do_react)
                src_q[mfi].plus(SDC_react_source[mfi],qbx,qbx,0,0,NQSRC);

        }
#endif
#endif

    }

}

// Convert a MultiFab with conservative state data u to a primitive MultiFab q.
void
Castro::cons_to_prim(MultiFab& u, MultiFab& q, MultiFab& qaux, Real time)
{

    BL_PROFILE("Castro::cons_to_prim()");

    BL_ASSERT(u.nComp() == NUM_STATE);
    BL_ASSERT(q.nComp() == NQ);
    BL_ASSERT(u.nGrow() >= q.nGrow());

    int ng = q.nGrow();

#ifdef RADIATION
    AmrLevel::FillPatch(*this, Erborder, NUM_GROW, time, Rad_Type, 0, Radiation::nGroups);

    MultiFab lamborder(grids, dmap, Radiation::nGroups, NUM_GROW);
    if (radiation->pure_hydro) {
      lamborder.setVal(0.0, NUM_GROW);
    }
    else {
      radiation->compute_limiter(level, grids, Sborder, Erborder, lamborder);
    }
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(u, true); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.growntilebox(ng);

	ca_ctoprim(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		   BL_TO_FORTRAN_ANYD(u[mfi]),
#ifdef RADIATION
                   BL_TO_FORTRAN_ANYD(Erborder[mfi]),
                   BL_TO_FORTRAN_ANYD(lamborder[mfi]),
#endif
		   BL_TO_FORTRAN_ANYD(q[mfi]),
		   BL_TO_FORTRAN_ANYD(qaux[mfi]));

    }

}

void
Castro::cons_to_prim_fourth(const Real time)
{
  // convert the conservative state cell averages to primitive cell
  // averages with 4th order accuracy

    MultiFab& S_new = get_new_data(State_Type);

    // we don't support radiation here
#ifdef RADIATION
    amrex::Abort("radiation not supported to fourth order");
#else
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      const Box& qbx = mfi.growntilebox(NUM_GROW);
      const Box& qbxm1 = mfi.growntilebox(NUM_GROW-1);

      // note: these conversions are using a growntilebox, so it
      // will include ghost cells

      // convert U_avg to U_cc -- this will use a Laplacian
      // operation and will result in U_cc defined only on
      // NUM_GROW-1 ghost cells at the end.
      FArrayBox U_cc;
      U_cc.resize(qbx, NUM_STATE);

      ca_make_cell_center(BL_TO_FORTRAN_BOX(qbxm1),
                          BL_TO_FORTRAN_FAB(Sborder[mfi]),
                          BL_TO_FORTRAN_FAB(U_cc));

      // convert U_avg to q_bar -- this will be done on all NUM_GROW
      // ghost cells.
      ca_ctoprim(BL_TO_FORTRAN_BOX(qbx),
                 BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                 BL_TO_FORTRAN_ANYD(q_bar[mfi]),
                 BL_TO_FORTRAN_ANYD(qaux_bar[mfi]));

      // this is what we should construct the flattening coefficient
      // from

      // convert U_cc to q_cc (we'll store this temporarily in q,
      // qaux).  This will remain valid only on the NUM_GROW-1 ghost
      // cells.
      ca_ctoprim(BL_TO_FORTRAN_BOX(qbxm1),
                 BL_TO_FORTRAN_ANYD(U_cc),
                 BL_TO_FORTRAN_ANYD(q[mfi]),
                 BL_TO_FORTRAN_ANYD(qaux[mfi]));
    }


    // check for NaNs
#ifndef AMREX_USE_CUDA
    check_for_nan(q);
    check_for_nan(q_bar);
#endif


#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      const Box& qbxm1 = mfi.growntilebox(NUM_GROW-1);

      // now convert q, qaux into 4th order accurate averages
      // this will create q, qaux in NUM_GROW-1 ghost cells, but that's
      // we need here

      ca_make_fourth_average(BL_TO_FORTRAN_BOX(qbxm1),
                             BL_TO_FORTRAN_FAB(q[mfi]),
                             BL_TO_FORTRAN_FAB(q_bar[mfi]));

      // not sure if we need to convert qaux this way, or if we can
      // just evaluate it (we may not need qaux at all actually)
      ca_make_fourth_average(BL_TO_FORTRAN_BOX(qbxm1),
                             BL_TO_FORTRAN_FAB(qaux[mfi]),
                             BL_TO_FORTRAN_FAB(qaux_bar[mfi]));

    }

#endif // RADIATION
}

void
Castro::check_for_cfl_violation(const Real dt)
{

    Real courno = -1.0e+200;

    const Real *dx = geom.CellSize();

    MultiFab& S_new = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel reduction(max:courno)
#endif
    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

#pragma gpu
        ca_compute_cfl(BL_TO_FORTRAN_BOX(bx),
                       BL_TO_FORTRAN_ANYD(q[mfi]),
                       BL_TO_FORTRAN_ANYD(qaux[mfi]),
                       dt, AMREX_REAL_ANYD(dx), AMREX_MFITER_REDUCE_MAX(&courno), print_fortran_warnings);

    }

    ParallelDescriptor::ReduceRealMax(courno);

    if (courno > 1.0) {
        amrex::Print() << "WARNING -- EFFECTIVE CFL AT LEVEL " << level << " IS " << courno << std::endl << std::endl;

        cfl_violation = 1;
    }

}
