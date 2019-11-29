
#include "Castro.H"
#include "Castro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

#ifdef SELF_GRAVITY
#include "Gravity.H"
#endif

#ifdef DIFFUSION
#include "Diffusion.H"
#endif

#include <cmath>
#include <climits>

using std::string;
using namespace amrex;

#ifndef AMREX_USE_CUDA
Real
Castro::do_advance_sdc (Real time,
                        Real dt,
                        int  amr_iteration,
                        int  amr_ncycle)
{

  // this is the new "formal" SDC integration routine.

  // this does the entire update in time for 1 SDC iteration.

  BL_PROFILE("Castro::do_advance_sdc()");

  const Real prev_time = state[State_Type].prevTime();
  const Real  cur_time = state[State_Type].curTime();

  MultiFab& S_old = get_old_data(State_Type);
  MultiFab& S_new = get_new_data(State_Type);

  const int* domain_lo = geom.Domain().loVect();
  const int* domain_hi = geom.Domain().hiVect();

  // Perform initialization steps.

  initialize_do_advance(time, dt, amr_iteration, amr_ncycle);

  // Check for NaN's.

  check_for_nan(S_old);

  MultiFab& old_source = get_old_data(Source_Type);
  MultiFab& new_source = get_new_data(Source_Type);

  // we loop over all nodes, even the last, since we need to compute
  // the advective update source at each node

  for (int m=0; m < SDC_NODES; m++) {

    current_sdc_node = m;

    // k_new represents carries the solution.  Coming into here, it
    // will be entirely the old state, but we update it on each time
    // node in place.

    Real node_time = time + dt_sdc[m]*dt;

    // fill Sborder with the starting node's info -- we use S_new as
    // our staging area.  Note we need to pass new_time here to the
    // FillPatch so it only pulls from the new MF -- this will not
    // work for multilevel.
    MultiFab::Copy(S_new, *(k_new[m]), 0, 0, S_new.nComp(), 0);
    clean_state(S_new, cur_time, 0);
    expand_state(Sborder, cur_time, NUM_GROW);


    // the next chunk of code constructs the advective term for the
    // current node, m.  First we get the sources, then full hydro
    // source term from the MOL driver.  Note, for m = 0, we only have
    // to do all of this the first iteration, since that state never
    // changes
    if (!(sdc_iteration > 0 && m == 0) &&
        !(sdc_iteration == sdc_order+sdc_extra-1 && m == SDC_NODES-1)) {

      // Construct the "old-time" sources from Sborder.  Since we are
      // working from Sborder, this will actually evaluate the sources
      // using the current stage's starting point.

      // TODO: this is not using the density at the current stage
#ifdef SELF_GRAVITY
      construct_old_gravity(amr_iteration, amr_ncycle, prev_time);
#endif

      if (apply_sources()) {
#ifndef AMREX_USE_CUDA
        if (sdc_order == 4) {
          // if we are 4th order, convert to cell-center Sborder -> Sborder_cc
          // we'll use Sburn for this memory buffer at the moment

          for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
            const Box& gbx = mfi.growntilebox(1);
            ca_make_cell_center(BL_TO_FORTRAN_BOX(gbx),
                                BL_TO_FORTRAN_FAB(Sborder[mfi]),
                                BL_TO_FORTRAN_FAB(Sburn[mfi]),
                                AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

          }

          // we pass in the stage time here
          do_old_sources(old_source, Sburn, node_time, dt, amr_iteration, amr_ncycle);

          // fill the ghost cells for the sources -- note since we have
          // not defined the new_source yet, we either need to copy this
          // into new_source for the time-interpolation in the ghost
          // fill to make sense, or so long as we are not multilevel,
          // just use the old time (prev_time) in the fill instead of
          // the node time (time)
          AmrLevel::FillPatch(*this, old_source, old_source.nGrow(), prev_time, Source_Type, 0, NSRC);

          // Now convert to cell averages.  This loop cannot be tiled.
          for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            ca_make_fourth_in_place(BL_TO_FORTRAN_BOX(bx),
                                    BL_TO_FORTRAN_FAB(old_source[mfi]),
                                    AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));
          }

        } else {
          // there is a ghost cell fill hidden in diffusion, so we need
          // to pass in the time associate with Sborder
          do_old_sources(old_source, Sborder, cur_time, dt, amr_iteration, amr_ncycle);
        }

        // note: we don't need a FillPatch on the sources, since they
        // are only used in the valid box in the conservative flux
        // update construction

#endif

        // store the result in sources_for_hydro -- this is what will
        // be used in the final conservative update
        MultiFab::Copy(sources_for_hydro, old_source, 0, 0, NSRC, 0);

      } else {
        sources_for_hydro.setVal(0.0, 0);
      }

      // Now compute the advective term for the current node -- this
      // will be used to advance us to the next node the new time


      // Construct the primitive variables.
      if (sdc_order == 4) {
        cons_to_prim_fourth(time);
      } else {
        cons_to_prim(time);
      }

      if (do_hydro) {
        // Check for CFL violations.
        check_for_cfl_violation(dt);

        // If we detect one, return immediately.
        if (cfl_violation)
          return dt;
      }

      // construct the update for the current stage -- this fills
      // A_new[m] with the righthand side for this stage.
      A_new[m]->setVal(0.0);
      construct_mol_hydro_source(time, dt, *A_new[m]);

    } // end of the m = 0 sdc_iter > 0 check

    // also, if we are the first SDC iteration, we haven't yet stored
    // any old advective terms, so we cannot yet do the quadrature
    // over nodes.  Initialize those now.  Recall, A_new[0] and A_old[0]
    // are aliased.
    if (sdc_iteration == 0 && m == 0) {
      for (int n=1; n < SDC_NODES; n++) {
        MultiFab::Copy(*(A_old[n]), *(A_new[0]), 0, 0, NUM_STATE, 0);
      }

#ifdef REACTIONS
      // if this is the first node of a new iteration, then we need
      // to compute and store the old reactive source

      // we already have the node state with ghost cells in Sborder,
      // so we can just use that as the starting point
      bool input_is_average = true;
      construct_old_react_source(Sborder, *(R_old[0]), input_is_average);

      // copy to the other nodes -- since the state is the same on all
      // nodes for sdc_iteration == 0
      for (int n = 1; n < SDC_NODES; n++) {
        MultiFab::Copy(*(R_old[n]), *(R_old[0]), 0, 0, R_old[0]->nComp(), 0);
      }
#endif
    }

    // update to the next stage -- this involves computing the
    // integral over the k-1 iteration data.  Note we don't do
    // this if we are on the final node (since there is nothing to
    // update to
    if (m < SDC_NODES-1) {

      amrex::Print() << "... doing the SDC update, iteration = " << sdc_iteration << " from node " << m << " to " << m+1 << std::endl;

      do_sdc_update(m, m+1, dt); //(dt_sdc[m+1] - dt_sdc[m])*dt);

      // we now have a new value of k_new[m+1], do a clean_state on it
      clean_state(S_new, cur_time, 0);

    }


  } // node iteration

  // we are done with the integration over all nodes for this iteration

  // note: as of right now, S_new is not updated with the state from
  // the final time node.  This means we can still use S_new as
  // "scratch" until we finally set it.

  if (sdc_iteration != sdc_order+sdc_extra-1) {
    // store A_old for the next SDC iteration -- don't need to do n=0,
    // since that is unchanged
    for (int n=1; n < SDC_NODES; n++) {
      MultiFab::Copy(*(A_old[n]), *(A_new[n]), 0, 0, NUM_STATE, 0);
    }
  }

#ifdef REACTIONS
  // we just did the update, so now recompute the "old" reactive
  // source for the next SDC iteration.  We don't need to do this for
  // m = 0, since that state never changes.

  for (int m = 1; m < SDC_NODES; ++m) {
    // TODO: do we need a clean state here?
    MultiFab::Copy(S_new, *(k_new[m]), 0, 0, S_new.nComp(), 0);
    expand_state(Sburn, cur_time, 2);
    bool input_is_average = true;
    construct_old_react_source(Sburn, *(R_old[m]), input_is_average);
  }
#endif

  if (sdc_iteration == sdc_order+sdc_extra-1) {

    // store the new solution
    MultiFab::Copy(S_new, *(k_new[SDC_NODES-1]), 0, 0, S_new.nComp(), 0);

    // We need to make source_old and source_new be the source terms at
    // the old and new time.  we never actually evaluated the sources
    // using the new time state (since we just constructed it).  Note:
    // we always use do_old_sources here, since we want the actual
    // source and not a correction.

    // note: we need to have ghost cells here cause some sources (in
    // particular pdivU) need them.  Perhaps it would be easier to just
    // always require State_Type to have 1 ghost cell?

    // TODO: we also need to make these 4th order!
    clean_state(S_old, prev_time, 0);
    expand_state(Sborder, prev_time, Sborder.nGrow());
    do_old_sources(old_source, Sborder, prev_time, dt, amr_iteration, amr_ncycle);
    AmrLevel::FillPatch(*this, old_source, old_source.nGrow(), prev_time, Source_Type, 0, NSRC);

    clean_state(S_new, cur_time, 0);
    expand_state(Sborder, cur_time, Sborder.nGrow());
    do_old_sources(new_source, Sborder, cur_time, dt, amr_iteration, amr_ncycle);
    AmrLevel::FillPatch(*this, new_source, new_source.nGrow(), cur_time, Source_Type, 0, NSRC);
  }

  finalize_do_advance(time, dt, amr_iteration, amr_ncycle);

  return dt;
}
#endif
