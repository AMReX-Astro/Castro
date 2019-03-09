#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

using namespace amrex;

void
Castro::do_sdc_update(int m_start, int m_end, Real dt_m) {

  // this routine needs to do the update from time node m to m+1
  //
  // We come in with:
  //   A_new[m_start] : this is the advective update at node m_start
  //   A_old[:] : this is the advective source for all nodes at the old iterate
  //
  //   R_old[:] : this is the reaction source for all nodes at the old iterate

  // If we do advection only, then the update is explicit.  If we do
  // reactions, then the update is implicit within a zone.

  // for 4th order reactive SDC, we need to first compute the source, C
  // and do a ghost cell fill on it

#ifdef REACTIONS
  // SDC_Source_Type is only defined for 4th order
  MultiFab tmp;
  MultiFab& C_source = (fourth_order == 1) ? get_new_data(SDC_Source_Type) : tmp;

  if (fourth_order == 1) {

    // for 4th order reacting flow, we need to create the "source" C
    // as averages and then convert it to cell centers.  The cell-center
    // version needs to have 2 ghost cells
    for (MFIter mfi(*k_new[0]); mfi.isValid(); ++mfi) {

      const Box& bx = mfi.tilebox();

      ca_sdc_compute_C4(BL_TO_FORTRAN_BOX(bx),
                        BL_TO_FORTRAN_3D((*A_new[m_start])[mfi]),
                        BL_TO_FORTRAN_3D((*A_old[0])[mfi]),
                        BL_TO_FORTRAN_3D((*A_old[1])[mfi]),
                        BL_TO_FORTRAN_3D((*A_old[2])[mfi]),
                        BL_TO_FORTRAN_3D((*R_old[0])[mfi]),
                        BL_TO_FORTRAN_3D((*R_old[1])[mfi]),
                        BL_TO_FORTRAN_3D((*R_old[2])[mfi]),
                        BL_TO_FORTRAN_3D(C_source[mfi]),
                        &m_start);
    }

    // need to construct the time for this stage -- but it is not really
    // at a single instance in time.  For single level this does not matter,
    Real time = state[SDC_Source_Type].curTime();
    AmrLevel::FillPatch(*this, C_source, C_source.nGrow(), time,
                        SDC_Source_Type, 0, NUM_STATE);

  }
#endif

  // main update loop -- we are updating k_new[m_start] to
  // k_new[m_end]

  FArrayBox U_center;
  FArrayBox C_center;
  FArrayBox U_new_center;
  FArrayBox R_new;

  for (MFIter mfi(*k_new[0]); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.tilebox();
    const Box& bx1 = mfi.growntilebox(1);

#ifdef REACTIONS
    // advection + reactions
    if (sdc_order == 2) {
      ca_sdc_update_o2(BL_TO_FORTRAN_BOX(bx), &dt_m,
                       BL_TO_FORTRAN_3D((*k_new[m_start])[mfi]),
                       BL_TO_FORTRAN_3D((*k_new[m_end])[mfi]),
                       BL_TO_FORTRAN_3D((*A_new[m_start])[mfi]),
                       BL_TO_FORTRAN_3D((*A_old[0])[mfi]),
                       BL_TO_FORTRAN_3D((*A_old[1])[mfi]),
                       BL_TO_FORTRAN_3D((*R_old[0])[mfi]),
                       BL_TO_FORTRAN_3D((*R_old[1])[mfi]),
                       &sdc_iteration,
                       &m_start);
    } else {

      // convert the starting U to cell-centered on a fab-by-fab basis
      // -- including one ghost cell
      U_center.resize(bx1, NUM_STATE);
      ca_make_cell_center(BL_TO_FORTRAN_BOX(bx1),
                          BL_TO_FORTRAN_FAB(Sborder[mfi]),
                          BL_TO_FORTRAN_FAB(U_center));

      // convert the C source to cell-centers
      C_center.resize(bx1, NUM_STATE);
      ca_make_cell_center(BL_TO_FORTRAN_BOX(bx1),
                          BL_TO_FORTRAN_FAB(C_source[mfi]),
                          BL_TO_FORTRAN_FAB(C_center));

      // solve for the updated cell-center U using our cell-centered C -- we
      // need to do this with one ghost cell
      U_new_center.resize(bx1, NUM_STATE);
      ca_sdc_update_centers_o4(BL_TO_FORTRAN_BOX(bx1), &dt_m,
                               BL_TO_FORTRAN_3D(U_center),
                               BL_TO_FORTRAN_3D(U_new_center),
                               BL_TO_FORTRAN_3D(C_center),
                               &sdc_iteration);

      // compute R_i and in 1 ghost cell and then convert to <R> in
      // place (only for the interior)
      R_new.resize(bx1, NUM_STATE);
      ca_instantaneous_react(BL_TO_FORTRAN_BOX(bx1),
                             BL_TO_FORTRAN_3D(U_new_center),
                             BL_TO_FORTRAN_3D(R_new));

      ca_make_cell_center_in_place(BL_TO_FORTRAN_BOX(bx),
                                   BL_TO_FORTRAN_FAB(R_new));

      // now do the conservative update using this <R> to get <U>
      // We'll also need to pass in <C>
      ca_sdc_conservative_update(BL_TO_FORTRAN_BOX(bx), &dt_m,
                                 BL_TO_FORTRAN_3D((*k_new[m_start])[mfi]),
                                 BL_TO_FORTRAN_3D((*k_new[m_end])[mfi]),
                                 BL_TO_FORTRAN_3D(C_source[mfi]),
                                 BL_TO_FORTRAN_3D(R_new));

    }
#else
    // pure advection
    if (sdc_order == 2) {
      ca_sdc_update_advection_o2(BL_TO_FORTRAN_BOX(bx), &dt_m,
                                 BL_TO_FORTRAN_3D((*k_new[m_start])[mfi]),
                                 BL_TO_FORTRAN_3D((*k_new[m_end])[mfi]),
                                 BL_TO_FORTRAN_3D((*A_new[m_start])[mfi]),
                                 BL_TO_FORTRAN_3D((*A_old[0])[mfi]),
                                 BL_TO_FORTRAN_3D((*A_old[1])[mfi]),
                                 &m_start);
    } else {
      ca_sdc_update_advection_o4(BL_TO_FORTRAN_BOX(bx), &dt_m,
                                 BL_TO_FORTRAN_3D((*k_new[m_start])[mfi]),
                                 BL_TO_FORTRAN_3D((*k_new[m_end])[mfi]),
                                 BL_TO_FORTRAN_3D((*A_new[m_start])[mfi]),
                                 BL_TO_FORTRAN_3D((*A_old[0])[mfi]),
                                 BL_TO_FORTRAN_3D((*A_old[1])[mfi]),
                                 BL_TO_FORTRAN_3D((*A_old[2])[mfi]),
                                 &m_start);
    }
#endif

  }
}


#ifdef REACTIONS
void
Castro::construct_old_react_source(amrex::MultiFab& U_state,
                                   amrex::MultiFab& R_source) {

  // this routine simply fills R_source with the reactive source from
  // state U_state.  Note: it is required that U_state have atleast 2
  // valid ghost cells for 4th order.

  // at this point, k_new has not yet been updated, so it represents
  // the state at the SDC nodes from the previous iteration
  for (MFIter mfi(U_state); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.tilebox();

    // construct the reactive source term
    ca_instantaneous_react(BL_TO_FORTRAN_BOX(bx),
                           BL_TO_FORTRAN_3D(U_state[mfi]),
                           BL_TO_FORTRAN_3D(R_source[mfi]));

  }
}
#endif
