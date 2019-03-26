
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

Real
Castro::do_advance_sdc (Real time,
                        Real dt,
                        int  amr_iteration,
                        int  amr_ncycle)
{

  // this is the new "formal" SDC integration routine.

  // unlike the MOL version which just operates on a single stage,
  // this does the entire update in time for 1 SDC iteration.

  BL_PROFILE("Castro::do_advance_sdc()");

  const Real prev_time = state[State_Type].prevTime();
  const Real  cur_time = state[State_Type].curTime();

  MultiFab& S_old = get_old_data(State_Type);
  MultiFab& S_new = get_new_data(State_Type);

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
    expand_state(Sborder, cur_time, 1, NUM_GROW);

    // Construct the "old-time" sources from Sborder.  Since we are
    // working from Sborder, this will actually evaluate the sources
    // using the current stage's starting point.

    // TODO: this is not using the density at the current stage
#ifdef SELF_GRAVITY
    construct_old_gravity(amr_iteration, amr_ncycle, prev_time);
#endif

    if (apply_sources()) {

      // we pass in the stage time here
      do_old_sources(old_source, Sborder, node_time, dt, amr_iteration, amr_ncycle);

      // hack: copy the source to the new data too, so fillpatch doesn't have to
      // worry about time
      MultiFab::Copy(new_source, old_source, 0, 0, NUM_STATE, 0);

      // Apply the old sources to the sources for the hydro.  Note that
      // we are doing an fill here, not an add (like we do for CTU --
      // this is because the source term predictor doesn't make sense
      // here).

      // we only need a fill here if sources_for_hydro has more ghost
      // cells than Source_Type, because otherwise, do_old_sources
      // already did the fill for us
      AmrLevel::FillPatch(*this, sources_for_hydro, NUM_GROW, time, Source_Type, 0, NUM_STATE);

    } else {
      old_source.setVal(0.0, old_source.nGrow());
    }

    // Now compute the advective term for the current node -- this
    // will be used to advance us to the next node the new time

    // Construct the primitive variables.
    if (fourth_order) {
      cons_to_prim_fourth(time);
    } else {
      cons_to_prim(time);
    }

    // Check for CFL violations.
    check_for_cfl_violation(dt);

    // If we detect one, return immediately.
    if (cfl_violation)
      return dt;

    // construct the update for the current stage -- this fills
    // A_new[m] with the righthand side for this stage.  Note, for m =
    // 0, the starting state is S_old and never changes with SDC
    // iteration, so we only do this once.
    if (!(sdc_iteration > 0 && m == 0)) {
      A_new[m]->setVal(0.0);
      construct_mol_hydro_source(time, dt, *A_new[m]);
    }

    // also, if we are the first SDC iteration, we haven't yet stored
    // any old advective terms, so we cannot yet do the quadrature
    // over nodes.  Initialize those now.  Recall, A_new[0] and A_old[0]
    // are aliased.
    if (sdc_iteration == 0 && m == 0) {
      for (int n=1; n < SDC_NODES; n++) {
        MultiFab::Copy(*(A_old[n]), *(A_new[0]), 0, 0, NUM_STATE, 0);
      }
    }

#ifdef REACTIONS
    // if this is the first node of a new iteration, then we need
    // to compute and store the old reactive source
    if (m == 0 && sdc_iteration == 0) {
      construct_old_react_source(Sborder, *(R_old[0]));

      // copy to the other nodes -- since the state is the same on all
      // nodes for sdc_iteration == 0
      for (int n = 1; n < SDC_NODES; n++) {
        MultiFab::Copy(*(R_old[n]), *(R_old[0]), 0, 0, R_old[0]->nComp(), 0);
      }
    }
#endif

    // update to the next stage -- this involves computing the
    // integral over the k-1 iteration data.  Note we don't do
    // this if we are on the final node (since there is nothing to
    // update to
    if (m < SDC_NODES-1) {
      do_sdc_update(m, m+1, dt); //(dt_sdc[m+1] - dt_sdc[m])*dt);
    }

  } // node iteration

  // we are done with the integration over all nodes for this iteration

  // note: as of right now, S_new is not updated with the state from
  // the final time node.  This means we can still use S_new as
  // "scratch" until we finally set it.

  // store A_old for the next SDC iteration -- don't need to do n=0,
  // since that is unchanged
  for (int n=1; n < SDC_NODES; n++) {
    MultiFab::Copy(*(A_old[n]), *(A_new[n]), 0, 0, NUM_STATE, 0);
  }

#ifdef REACTIONS
  // we just did the update, so now recompute the "old" reactive
  // source for the next SDC iteration.  We don't need to do this for
  // m = 0, since that state never changes.

  for (int m = 1; m < SDC_NODES; ++m) {
    // use a temporary storage
    // TODO: do we need a clean state here?
    MultiFab::Copy(S_new, *(k_new[m]), 0, 0, S_new.nComp(), 0);
    expand_state(Sborder, cur_time, -1, Sborder.nGrow());
    construct_old_react_source(Sborder, *(R_old[m]));
  }
#endif

  // store the new solution
  MultiFab::Copy(S_new, *(k_new[SDC_NODES-1]), 0, 0, S_new.nComp(), 0);


  // I think this bit only needs to be done for the last iteration...

  // We need to make source_old and source_new be the source terms at
  // the old and new time.  we never actually evaluate the sources
  // using the new time state (since we just constructed it).  Note:
  // we always use do_old_sources here, since we want the actual
  // source and not a correction.

  // note: we need to have ghost cells here cause some sources (in
  // particular pdivU) need them.  Perhaps it would be easier to just
  // always require State_Type to have 1 ghost cell?
  expand_state(Sborder, prev_time, 0, Sborder.nGrow());
  do_old_sources(old_source, Sborder, prev_time, dt, amr_iteration, amr_ncycle);

  expand_state(Sborder, cur_time, 1, Sborder.nGrow());
  do_old_sources(new_source, Sborder, cur_time, dt, amr_iteration, amr_ncycle);


  finalize_do_advance(time, dt, amr_iteration, amr_ncycle);

  return dt;
}

