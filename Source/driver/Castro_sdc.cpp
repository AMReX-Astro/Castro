#include "Castro.H"
#include "Castro_F.H"

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
  if (fourth_order == 1) {
    MultiFab& C_source = get_new_data(SDC_Source_Type);

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

  for (MFIter mfi(*k_new[0]); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.tilebox();

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
      amrex::Abort("sdc_order != 2 not implemented");

      // convert the starting U to cell-centered on a fab-by-fab basis
      // -- including one ghost cell

      // solve for the cell-center U using our cell-centered C -- we
      // need to do this with one ghost cell

      // do the conservative correction -- we need to compute <R>
      // first, then get <U>.  We'll also need to pass in <C> too


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
Castro::construct_old_react_source() {

  // this routine simply fills R_old with the old-iteration reactive
  // source

  // at this point, k_new has not yet been updated, so it represents
  // the state at the SDC nodes from the previous iteration
  for (MFIter mfi(*k_new[0]); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.tilebox();

    for (int m=0; m < SDC_NODES; m++) {

      // if we are not the first iteration, then we don't need to
      // recompute the reactive source at the first node -- this is
      // the old time, and hence doesn't change from one iteration to
      // the next.
      if (sdc_iteration > 0 && m == 0)
        continue;

      // construct the reactive source term
      ca_instantaneous_react(BL_TO_FORTRAN_BOX(bx),
                             BL_TO_FORTRAN_3D((*k_new[m])[mfi]),
                             BL_TO_FORTRAN_3D((*R_old[m])[mfi]));

      // if we are the very first iteration, then there is no old state
      // at all the time nodes, so we just copy R_old[0] into the other
      // nodes
      if (sdc_iteration == 0) {

        for (int n=1; n < SDC_NODES; n++) {
          MultiFab::Copy(*(R_old[n]), *(R_old[0]), 0, 0, R_old[0]->nComp(), 0);
        }

        break;
      }

    }

  }
}
#endif
