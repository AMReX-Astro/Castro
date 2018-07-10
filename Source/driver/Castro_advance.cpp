
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
Castro::advance (Real time,
                 Real dt,
                 int  amr_iteration,
                 int  amr_ncycle)

  // the main driver for a single level.  This will do either the SDC
  // algorithm or the Strang-split reactions algorithm.
  //
  // arguments:
  //    time          : the current simulation time
  //    dt            : the timestep to advance (e.g., go from time to
  //                    time + dt)
  //    amr_iteration : where we are in the current AMR subcycle.  Each
  //                    level will take a number of steps to reach the
  //                    final time of the coarser level below it.  This
  //                    counter starts at 1
  //    amr_ncycle    : the number of subcycles at this level

{
    BL_PROFILE("Castro::advance()");

    // Save the wall time when we started the step.

    wall_time_start = ParallelDescriptor::second();

    Real dt_new = dt;

    initialize_advance(time, dt, amr_iteration, amr_ncycle);

    // Do the advance.

#ifdef SDC
    // this is the old SDC methodology

    for (int n = 0; n < sdc_iters; ++n) {

        sdc_iteration = n;

        if (ParallelDescriptor::IOProcessor())
	    std::cout << "\nBeginning SDC iteration " << n + 1 << " of " << sdc_iters << ".\n\n";

	// First do the non-reacting advance and construct the relevant source terms.

	dt_new = do_advance(time, dt, amr_iteration, amr_ncycle);

#ifdef REACTIONS
	if (do_react) {

            // Do the ODE integration to capture the reaction source terms.

	    react_state(time, dt);

	    MultiFab& S_new = get_new_data(State_Type);

            int is_new=1;
	    clean_state(is_new, S_new.nGrow());

	    // Compute the reactive source term for use in the next iteration.

	    MultiFab& SDC_react_new = get_new_data(SDC_React_Type);
	    get_react_source_prim(SDC_react_new, dt);

	    // Check for NaN's.

	    check_for_nan(S_new);

        }
#endif

        if (ParallelDescriptor::IOProcessor())
	    std::cout << "\nEnding SDC iteration " << n + 1 << " of " << sdc_iters << ".\n\n";

    }

#else
    // we are either CTU, MOL, or the new SDC

#ifndef AMREX_USE_CUDA
    if (time_integration_method == CTU) {

        dt_new = std::min(dt_new, subcycle_advance(time, dt, amr_iteration, amr_ncycle));

    } else if (time_integration_method == MOL) {
#endif

      for (int iter = 0; iter < MOL_STAGES; ++iter) {
	mol_iteration = iter;
	dt_new = do_advance_mol(time + c_mol[iter]*dt, dt, amr_iteration, amr_ncycle);
      }

#ifndef AMREX_USE_CUDA
    } else if (time_integration_method == SDC) {

      for (int iter = 0; iter < sdc_order; ++iter) {
	sdc_iteration = iter;
	dt_new = do_advance_sdc(time, dt, amr_iteration, amr_ncycle);
      }

      // store the new solution
      MultiFab& S_new = get_new_data(State_Type);
      MultiFab::Copy(S_new, *(k_new[SDC_NODES-1]), 0, 0, S_new.nComp(), 0);

#ifdef REACTIONS
      // store the reaction information as well -- note: this will be
      // the instantaneous reactive source.  In the future, we might
      // want to do a quadrature over R_new[]
      MultiFab& R_new = get_new_data(Reactions_Type);

      for (MFIter mfi(R_new, hydro_tile_size); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const int idx = mfi.tileIndex();

        ca_store_reaction_state(BL_TO_FORTRAN_BOX(bx),
                                BL_TO_FORTRAN_3D((*R_old[SDC_NODES-1])[mfi]),
                                BL_TO_FORTRAN_3D(S_new[mfi]),
                                BL_TO_FORTRAN_3D(R_new[mfi]));

      }
#endif
    }
#endif

    // Optionally kill the job at this point, if we've detected a violation.

    if (cfl_violation && hard_cfl_limit && !use_retry)
        amrex::Abort("CFL is too high at this level -- go back to a checkpoint and restart with lower cfl number");

    // If we didn't kill the job, reset the violation counter.

    cfl_violation = 0;

#endif

    if (use_post_step_regrid)
	check_for_post_regrid(time + dt);

#ifdef AUX_UPDATE
    advance_aux(time, dt);
#endif

#ifdef SELF_GRAVITY
#if (BL_SPACEDIM > 1)
    // We do this again here because the solution will have changed
    if ( (level == 0) && (spherical_star == 1) ) {
       int is_new = 1;
       make_radial_data(is_new);
    }
#endif
#endif

#ifdef POINTMASS
    // Update the point mass.
    if (use_point_mass)
        pointmass_update(time, dt);
#endif

#ifdef RADIATION
    MultiFab& S_new = get_new_data(State_Type);
    final_radiation_call(S_new, amr_iteration, amr_ncycle);
#endif

#ifdef AMREX_PARTICLES
    advance_particles(amr_iteration, time, dt);
#endif

    finalize_advance(time, dt, amr_iteration, amr_ncycle);

    return dt_new;
}



#ifndef AMREX_USE_CUDA
Real
Castro::do_advance (Real time,
                    Real dt,
                    int  amr_iteration,
                    int  amr_ncycle)
{

  // this routine will advance the old state data (called S_old here)
  // to the new time, for a single level.  The new data is called
  // S_new here.  The update includes reactions (if we are not doing
  // SDC), hydro, and the source terms.

    BL_PROFILE("Castro::do_advance()");

    const Real prev_time = state[State_Type].prevTime();
    const Real  cur_time = state[State_Type].curTime();

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    // Perform initialization steps.

    initialize_do_advance(time, dt, amr_iteration, amr_ncycle);

    // Check for NaN's.

    check_for_nan(S_old);

    // Since we are Strang splitting the reactions, do them now


#ifdef REACTIONS
#ifndef SDC
    // this operates on Sborder (which is initially S_old).  The result
    // of the reactions is added directly back to Sborder.
    strang_react_first_half(prev_time, 0.5 * dt);
#endif
#endif

    // Initialize the new-time data. This copy needs to come after the
    // reactions.

    MultiFab::Copy(S_new, Sborder, 0, 0, NUM_STATE, S_new.nGrow());

#ifdef REACTIONS
#ifndef SDC
    // Do this for the reactions as well, in case we cut the timestep
    // short due to it being rejected.

    MultiFab& R_old = get_old_data(Reactions_Type);
    MultiFab& R_new = get_new_data(Reactions_Type);
    MultiFab::Copy(R_new, R_old, 0, 0, R_new.nComp(), R_new.nGrow());

    // Skip the rest of the advance if the burn was unsuccessful.

    if (burn_success != 1)
        return dt;
#endif
#endif

    // Construct the old-time sources from Sborder.  This will already
    // be applied to S_new (with full dt weighting), to be correctly
    // later.  Note -- this does not affect the prediction of the
    // interface state, an explict source will be traced there as
    // needed.

#ifdef SELF_GRAVITY
    construct_old_gravity(amr_iteration, amr_ncycle, prev_time);
#endif

    MultiFab& old_source = get_old_data(Source_Type);

    if (apply_sources()) {

      do_old_sources(old_source, Sborder, prev_time, dt, amr_iteration, amr_ncycle);

      int is_new=1;
      apply_source_to_state(is_new, S_new, old_source, dt, S_new.nGrow());

      // Apply the old sources to the sources for the hydro.
      // Note that we are doing an add here, not a copy,
      // in case we have already started with some source
      // terms (e.g. the source term predictor, or the SDC source).

      AmrLevel::FillPatchAdd(*this, sources_for_hydro, NUM_GROW, time, Source_Type, 0, NUM_STATE);

    } else {
      old_source.setVal(0.0, NUM_GROW);

    }


    // Do the hydro update.  We build directly off of Sborder, which
    // is the state that has already seen the burn

    if (do_hydro)
    {
      // Construct the primitive variables.
      cons_to_prim(time);

      // Check for CFL violations.
      check_for_cfl_violation(dt);

      // If we detect one, return immediately.
      if (cfl_violation && hard_cfl_limit)
          return dt;

      construct_hydro_source(time, dt);
      int is_new=1;
      apply_source_to_state(is_new, S_new, hydro_source, dt);

    }


    // Sync up state after old sources and hydro source.
    int is_new=1;
    frac_change = clean_state(is_new, Sborder, S_new.nGrow());

    // If the state has ghost zones, sync them up now
    // since the hydro source only works on the valid zones.

    if (S_new.nGrow() > 0) {
      expand_state(S_new, cur_time, 1, S_new.nGrow());
    }

    // Check for NaN's.

    check_for_nan(S_new);

    // if we are done with the update do the source correction and
    // then the second half of the reactions

#ifdef SELF_GRAVITY
    // Must define new value of "center" before we call new gravity
    // solve or external source routine
    if (moving_center == 1)
      define_new_center(S_new, time);
#endif

#ifdef SELF_GRAVITY
    // We need to make the new radial data now so that we can use it when we
    // FillPatch in creating the new source.

#if (BL_SPACEDIM > 1)
    if ( (level == 0) && (spherical_star == 1) ) {
      int is_new = 1;
      make_radial_data(is_new);
    }
#endif
#endif

    // Construct and apply new-time source terms.

#ifdef SELF_GRAVITY
    construct_new_gravity(amr_iteration, amr_ncycle, cur_time);
#endif

    MultiFab& new_source = get_new_data(Source_Type);

    if (apply_sources()) {

      do_new_sources(new_source, Sborder, S_new, cur_time, dt, amr_iteration, amr_ncycle);

      int is_new=1;
      apply_source_to_state(is_new, S_new, new_source, dt, S_new.nGrow());

    } else {

      new_source.setVal(0.0, NUM_GROW);

    }

    // Do the second half of the reactions.

#ifdef REACTIONS
#ifndef SDC
    strang_react_second_half(cur_time - 0.5 * dt, 0.5 * dt);

    // Skip the rest of the advance if the burn was unsuccessful.

    if (burn_success != 1)
        return dt;
#endif
#endif

    finalize_do_advance(time, dt, amr_iteration, amr_ncycle);

    return dt;
}
#endif


Real
Castro::do_advance_mol (Real time,
                        Real dt,
                        int  amr_iteration,
                        int  amr_ncycle)
{

  // this routine will advance the old state data (called S_old here)
  // to the new time, for a single level.  The new data is called
  // S_new here.  The update includes reactions, hydro, and the source
  // terms.

  // NOTE: the time that passes through here is the time for the
  // current stage

  BL_PROFILE("Castro::do_advance_mol()");

  const Real prev_time = state[State_Type].prevTime();
  const Real  cur_time = state[State_Type].curTime();

  MultiFab& S_old = get_old_data(State_Type);
  MultiFab& S_new = get_new_data(State_Type);

  // Perform initialization steps.

  initialize_do_advance(time, dt, amr_iteration, amr_ncycle);

#ifndef AMREX_USE_CUDA
  // Check for NaN's.

  check_for_nan(S_old);
#endif

  // Since we are Strang splitting the reactions, do them now (only
  // for first stage of MOL)

  if (mol_iteration == 0) {

#ifndef SDC
#ifdef REACTIONS
    // this operates on Sborder (which is initially S_old).  The result
    // of the reactions is added directly back to Sborder.
    strang_react_first_half(prev_time, 0.5 * dt);
#endif
#endif
    // store the result of the burn in Sburn for later stages
    MultiFab::Copy(Sburn, Sborder, 0, 0, NUM_STATE, 0);
  }


  // Construct the "old-time" sources from Sborder.  Since we are
  // working from Sborder, this will actually evaluate the sources
  // using the current stage's starting point.

  // We do not apply the sources here -- they will be integrated
  // in the RK integration to come

  // TODO: this is not using the density at the current stage
#ifdef SELF_GRAVITY
  construct_old_gravity(amr_iteration, amr_ncycle, prev_time);
#endif

  MultiFab& old_source = get_old_data(Source_Type);
  MultiFab& new_source = get_new_data(Source_Type);

  if (apply_sources()) {

#ifndef AMREX_USE_CUDA
    if (fourth_order) {
      // if we are 4th order, convert to cell-center Sborder -> Sborder_cc
      // we'll reuse sources_for_hydro for this memory buffer at the moment

      for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const int idx = mfi.tileIndex();
        ca_make_cell_center(BL_TO_FORTRAN_BOX(bx),
                            BL_TO_FORTRAN_FAB(Sborder[mfi]),
                            BL_TO_FORTRAN_FAB(sources_for_hydro[mfi]));

      }
    }

    // we pass in the stage time here
    if (fourth_order) {
      do_old_sources(old_source, sources_for_hydro, time, dt, amr_iteration, amr_ncycle);

      // Note: this filled the ghost cells for us, so we can now convert to
      // cell averages.  This loop cannot be tiled.
      for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const int idx = mfi.tileIndex();
        ca_make_fourth_in_place(BL_TO_FORTRAN_BOX(bx),
                                BL_TO_FORTRAN_FAB(old_source[mfi]));

      }

      // now that we redid these, redo the ghost fill
      AmrLevel::FillPatch(*this, old_source, old_source.nGrow(), time, Source_Type, 0, NUM_STATE);

    } else {
      do_old_sources(old_source, Sborder, time, dt, amr_iteration, amr_ncycle);
    }
#endif

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


  // Do the hydro update.  We build directly off of Sborder, which
  // is the state that has already seen the burn

  if (do_hydro)
    {
      // Construct the primitive variables.
      if (fourth_order) {
#ifndef AMREX_USE_CUDA
        cons_to_prim_fourth(time);
#endif
      } else {
        cons_to_prim(time);
      }

      // Check for CFL violations.
      check_for_cfl_violation(dt);

      // If we detect one, return immediately.
      if (cfl_violation)
        return dt;

      // construct the update for the current stage -- this fills k_mol
      // with the righthand side for this stage
      construct_mol_hydro_source(time, dt, *k_mol[mol_iteration]);
    }

  // For MOL integration, we are done with this stage, unless it is
  // the last stage
  if (mol_iteration < MOL_STAGES-1) {
    finalize_do_advance(time, dt, amr_iteration, amr_ncycle);
    return dt;
  }

  // we just finished the last stage of the MOL integration.
  // Construct S_new now using the weighted sum of the k_mol
  // updates -- this will include both the advective and
  // source terms

  // Apply the update -- we need to build on Sburn, so
  // start with that state
  MultiFab::Copy(S_new, Sburn, 0, 0, S_new.nComp(), 0);
  for (int i = 0; i < MOL_STAGES; ++i)
    MultiFab::Saxpy(S_new, dt*b_mol[i], *k_mol[i], 0, 0, S_new.nComp(), 0);

  // define the temperature now
  int is_new=1;
  clean_state(is_new, S_new.nGrow());

  // If the state has ghost zones, sync them up now
  // since the hydro source only works on the valid zones.

  if (S_new.nGrow() > 0) {
    expand_state(S_new, cur_time, 1, S_new.nGrow());
  }

#ifndef AMREX_USE_CUDA
  // Check for NaN's.
  check_for_nan(S_new);
#endif

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


  // Do the second half of the reactions.

#ifndef SDC
#ifdef REACTIONS
  strang_react_second_half(cur_time - 0.5 * dt, 0.5 * dt);
#endif
#endif

  finalize_do_advance(time, dt, amr_iteration, amr_ncycle);

  return dt;
}


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
    if (m == 0) {
      construct_old_react_source();
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

  // store A_old for the next SDC iteration -- don't need to do n=0,
  // since that is unchanged
  for (int n=1; n < SDC_NODES; n++) {
    MultiFab::Copy(*(A_old[n]), *(A_new[n]), 0, 0, NUM_STATE, 0);
  }

#ifdef REACTIONS
  // we just did the update, so now recompute the "old" reactive source
  // for the next SDC iteration
  construct_old_react_source();
#endif


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



void
Castro::initialize_do_advance(Real time, Real dt, int amr_iteration, int amr_ncycle)
{

    // Reset the change from density resets

    frac_change = 1.e0;

    // Reset the CFL violation flag.

    cfl_violation = 0;

    // Reset the burn success flag.

    burn_success = 1;

    int finest_level = parent->finestLevel();

#ifdef RADIATION
    // make sure these are filled to avoid check/plot file errors:
    if (do_radiation) {
      get_old_data(Rad_Type).setBndry(0.0);
      get_new_data(Rad_Type).setBndry(0.0);
    }
    else {
      get_old_data(Rad_Type).setVal(0.0);
      get_new_data(Rad_Type).setVal(0.0);
    }
    get_old_data(State_Type).setBndry(0.0);
    get_new_data(State_Type).setBndry(0.0);
#endif

    // Reset the grid loss tracking.

    if (track_grid_losses)
      for (int i = 0; i < n_lost; i++)
	material_lost_through_boundary_temp[i] = 0.0;

#ifdef SELF_GRAVITY
    if (moving_center == 1)
        define_new_center(get_old_data(State_Type), time);

#if (BL_SPACEDIM > 1)
    if ( (level == 0) && (spherical_star == 1) ) {
       swap_outflow_data();
       int is_new = 0;
       make_radial_data(is_new);
    }
#endif
#endif

    // Scale the source term predictor by the current timestep.

#ifndef SDC
    if (time_integration_method == CTU && source_term_predictor == 1) {
        sources_for_hydro.mult(0.5 * dt, NUM_GROW);
    }
#endif

    // For the hydrodynamics update we need to have NUM_GROW ghost
    // zones available, but the state data does not carry ghost
    // zones. So we use a FillPatch using the state data to give us
    // Sborder, which does have ghost zones.

    if (time_integration_method == CTU) {
      // for the CTU unsplit method, we always start with the old state
      Sborder.define(grids, dmap, NUM_STATE, NUM_GROW);
      const Real prev_time = state[State_Type].prevTime();
      expand_state(Sborder, prev_time, 0, NUM_GROW);

    } else if (time_integration_method == MOL)  {
      // for Method of lines, our initialization of Sborder depends on
      // which stage in the RK update we are working on

      if (mol_iteration == 0) {

	// first MOL stage
	Sborder.define(grids, dmap, NUM_STATE, NUM_GROW);
	const Real prev_time = state[State_Type].prevTime();
	expand_state(Sborder, prev_time, 0, NUM_GROW);

      } else {

	// the initial state for the kth stage follows the Butcher
	// tableau.  We need to create the proper state starting with
	// the result after the first dt/2 burn (which we copied into
	// Sburn) and we need to fill ghost cells.

	// We'll overwrite S_new with this information, since we don't
	// need it anymorebuild this state temporarily in S_new (which
	// is State_Data) to allow for ghost filling.
	MultiFab& S_new = get_new_data(State_Type);

	MultiFab::Copy(S_new, Sburn, 0, 0, S_new.nComp(), 0);
	for (int i = 0; i < mol_iteration; ++i)
	  MultiFab::Saxpy(S_new, dt*a_mol[mol_iteration][i], *k_mol[i], 0, 0, S_new.nComp(), 0);

        // not sure if this is needed
        int is_new=1;
        clean_state(is_new, S_new.nGrow());

	Sborder.define(grids, dmap, NUM_STATE, NUM_GROW);
	const Real new_time = state[State_Type].curTime();
	expand_state(Sborder, new_time, 1, NUM_GROW);

      }

    } else if (time_integration_method == SDC) {

      // we'll handle the filling inside of do_advance_sdc 
      Sborder.define(grids, dmap, NUM_STATE, NUM_GROW);

    } else {
      amrex::Abort("invalid time_integration_method");
    }

}



void
Castro::finalize_do_advance(Real time, Real dt, int amr_iteration, int amr_ncycle)
{

#ifdef RADIATION
    if (!do_hydro && Radiation::rad_hydro_combined) {
	MultiFab& Er_old = get_old_data(Rad_Type);
	MultiFab& Er_new = get_new_data(Rad_Type);
	Er_new.copy(Er_old);
    }
#endif

    Sborder.clear();

}



void
Castro::initialize_advance(Real time, Real dt, int amr_iteration, int amr_ncycle)
{
    // Save the current iteration.

    iteration = amr_iteration;

    do_subcycle = false;
    sub_iteration = 0;
    sub_ncycle = 0;
    dt_subcycle = 1.e200;
    dt_advance = dt;

    keep_prev_state = false;

    // Reset the retry timestep information.

    lastDtRetryLimited = 0;
    lastDtFromRetry = 1.e200;

    if (use_post_step_regrid && level > 0) {

	if (getLevel(level-1).post_step_regrid && amr_iteration == 1) {

            // If the level below this just triggered a special regrid,
            // the coarse contribution to this level's FluxRegister
            // is no longer valid because the grids have, in general, changed.
            // Zero it out, and add them back using the saved copy of the fluxes.

	    getLevel(level-1).FluxRegCrseInit();

            // If we're coming off a new regrid at the end of the last coarse
            // timestep, then we want to subcycle this timestep at the timestep
            // suggested by this level, since the data on this level will not
            // have been taken into account when calculating the timestep
            // constraint using the coarser data. This is true even if the level
            // previously existed, because in general there can be new data at this
            // level as a result of the regrid.

            // This step MUST be done before the time level swap because estTimeStep
            // looks at the "new" time data for calculating the timestep constraint.
            // It should also be done before the call to ca_set_amr_info since estTimeStep
            // temporarily resets the level data.

            dt_subcycle = estTimeStep(dt);

            if (dt_subcycle < dt) {

                sub_ncycle = ceil(dt / dt_subcycle);

                if (ParallelDescriptor::IOProcessor()) {
                    std::cout << std::endl;
                    std::cout << "  Subcycling with maximum dt = " << dt_subcycle << " at level " << level
                              << " to avoid timestep constraint violations after a post-timestep regrid."
                              << std::endl << std::endl;
                }

                do_subcycle = true;

            }

        }

    }

    // Pass some information about the state of the simulation to a Fortran module.

    ca_set_amr_info(level, amr_iteration, amr_ncycle, time, dt);

    // The option of whether to do a multilevel initialization is
    // controlled within the radiation class.  This step belongs
    // before the swap.

#ifdef RADIATION
    if (do_radiation)
        radiation->pre_timestep(level);

    Erborder.define(grids, dmap, Radiation::nGroups, NUM_GROW);
    lamborder.define(grids, dmap, Radiation::nGroups, NUM_GROW);
#endif

#ifdef SELF_GRAVITY
    // If we're on level 0, update the maximum density used in the gravity solver
    // for setting the tolerances. This will be used in all level solves to follow.
    // This must be done before the swap because it relies on the new data.

    if (level == 0 && gravity->get_gravity_type() == "PoissonGrav") {
	gravity->update_max_rhs();
    }
#endif

    // If we're going to do a retry, or more generally if we're about to
    // subcycle the advance, save the simulation times of the
    // previous state data. This must happen before the swap.

    if (use_retry || do_subcycle) {

        prev_state_old_time = get_state_data(State_Type).prevTime();
        prev_state_new_time = get_state_data(State_Type).curTime();

        prev_state_had_old_data = get_state_data(State_Type).hasOldData();

    }

    // This array holds the sum of all source terms that affect the
    // hydrodynamics.  If we are doing the source term predictor,
    // we'll also use this after the hydro update to store the sum of
    // the new-time sources, so that we can compute the time
    // derivative of the source terms.

    sources_for_hydro.define(grids, dmap, NUM_STATE, NUM_GROW);
    sources_for_hydro.setVal(0.0, NUM_GROW);

#ifndef SDC

    // Add the source term predictor.
    // This must happen before the swap.

    if (time_integration_method == CTU && source_term_predictor == 1) {
        apply_source_term_predictor();
    }

#else

    // If we're doing SDC, time-center the source term (using the
    // current iteration's old sources and the last iteration's new
    // sources). Since the "new-time" sources are just the corrector step
    // of the predictor-corrector formalism, we want to add the full
    // value of the "new-time" sources to the old-time sources to get a
    // time-centered value.

    AmrLevel::FillPatch(*this, sources_for_hydro, NUM_GROW, time, Source_Type, 0, NUM_STATE);

#endif


    // Swap the new data from the last timestep into the old state data.

    swap_state_time_levels(dt);

#ifdef SELF_GRAVITY
    if (do_grav)
	gravity->swapTimeLevels(level);
#endif

    // Ensure data is valid before beginning advance. This addresses
    // the fact that we may have new data on this level that was interpolated
    // from a coarser level, and the interpolation in general cannot be
    // trusted to respect the consistency between certain state variables
    // (e.g. UEINT and UEDEN) that we demand in every zone.

    int is_new=0;
    MultiFab& S_old = get_old_data(State_Type);
    clean_state(is_new, S_old.nGrow());

    // Initialize the previous state data container now, so that we can
    // always ask if it has valid data.

    for (int k = 0; k < num_state_type; ++k)
        prev_state[k].reset(new StateData());

    // Make a copy of the MultiFabs in the old and new state data in case we may do a retry.

    if (use_retry || do_subcycle) {

      // Store the old and new time levels.

      for (int k = 0; k < num_state_type; k++) {

	StateData::Initialize(*prev_state[k], state[k]);

      }

    }

    // This array holds the hydrodynamics update.

    hydro_source.define(grids,dmap,NUM_STATE,0);



    // Allocate space for the primitive variables.

    q.define(grids, dmap, NQ, NUM_GROW);
    q.setVal(0.0);
    qaux.define(grids, dmap, NQAUX, NUM_GROW);
    if (time_integration_method == CTU)
      src_q.define(grids, dmap, QVAR, NUM_GROW);
    if (fourth_order)
      q_bar.define(grids, dmap, NQ, NUM_GROW);
      qaux_bar.define(grids, dmap, NQAUX, NUM_GROW);

    if (time_integration_method == MOL) {
      // if we are not doing CTU advection, then we are doing a method
      // of lines, and need storage for hte intermediate stages
      k_mol.resize(MOL_STAGES);
      for (int n = 0; n < MOL_STAGES; ++n) {
	k_mol[n].reset(new MultiFab(grids, dmap, NUM_STATE, 0));
	k_mol[n]->setVal(0.0);
      }

      // for the post-burn state
      Sburn.define(grids, dmap, NUM_STATE, 0);
    }

    if (time_integration_method == SDC) {

      MultiFab& S_old = get_old_data(State_Type);
      k_new.resize(SDC_NODES);
      k_new[0].reset(new MultiFab(S_old, amrex::make_alias, 0, NUM_STATE));
      for (int n = 1; n < SDC_NODES; ++n) {
	k_new[n].reset(new MultiFab(grids, dmap, NUM_STATE, 0));
	k_new[n]->setVal(0.0);
      }

      A_old.resize(SDC_NODES);
      for (int n = 0; n < SDC_NODES; ++n) {
	A_old[n].reset(new MultiFab(grids, dmap, NUM_STATE, 0));
	A_old[n]->setVal(0.0);
      }

      A_new.resize(SDC_NODES);
      A_new[0].reset(new MultiFab(*A_old[0], amrex::make_alias, 0, NUM_STATE));
      for (int n = 1; n < SDC_NODES; ++n) {
	A_new[n].reset(new MultiFab(grids, dmap, NUM_STATE, 0));
        A_new[n]->setVal(0.0);
      }

#ifdef REACTIONS
      R_old.resize(SDC_NODES);
      for (int n = 0; n < SDC_NODES; ++n) {
	R_old[n].reset(new MultiFab(grids, dmap, NUM_STATE, 0));
        R_old[n]->setVal(0.0);
      }
#endif
    }

    // Zero out the current fluxes.

    for (int dir = 0; dir < 3; ++dir)
	fluxes[dir]->setVal(0.0);

    for (int dir = 0; dir < 3; ++dir)
        mass_fluxes[dir]->setVal(0.0);

#if (BL_SPACEDIM <= 2)
    if (!Geometry::IsCartesian())
	P_radial.setVal(0.0);
#endif

#ifdef RADIATION
    if (Radiation::rad_hydro_combined)
	for (int dir = 0; dir < BL_SPACEDIM; ++dir)
	    rad_fluxes[dir]->setVal(0.0);
#endif

}



void
Castro::finalize_advance(Real time, Real dt, int amr_iteration, int amr_ncycle)
{

    // Add the material lost in this timestep to the cumulative losses.

    if (track_grid_losses) {

      ParallelDescriptor::ReduceRealSum(material_lost_through_boundary_temp, n_lost);

      for (int i = 0; i < n_lost; i++)
	material_lost_through_boundary_cumulative[i] += material_lost_through_boundary_temp[i];

    }

    if (do_reflux) {
	FluxRegCrseInit();
	FluxRegFineAdd();
    }

    Real cur_time = state[State_Type].curTime();

    hydro_source.clear();

    q.clear();
    qaux.clear();
    if (time_integration_method == CTU)
      src_q.clear();
    if (fourth_order) {
      q_bar.clear();
      qaux_bar.clear();
    }

#ifdef RADIATION
    Erborder.clear();
    lamborder.clear();
#endif

    sources_for_hydro.clear();

    if (!keep_prev_state)
        amrex::FillNull(prev_state);

    if (time_integration_method == MOL) {
      k_mol.clear();
      Sburn.clear();
    }

    if (time_integration_method == SDC) {
      k_new.clear();
      A_new.clear();
      A_old.clear();
#ifdef REACTIONS
      R_old.clear();
#endif
    }

    // Record how many zones we have advanced.

    num_zones_advanced += grids.numPts() / getLevel(0).grids.numPts();

}



#ifndef AMREX_USE_CUDA
bool
Castro::retry_advance(Real& time, Real dt, int amr_iteration, int amr_ncycle)
{

    Real dt_new = 1.e200;
    Real dt_sub = 1.e200;

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

#ifdef REACTIONS
    MultiFab& R_old = get_old_data(Reactions_Type);
    MultiFab& R_new = get_new_data(Reactions_Type);
#endif

    const Real* dx = geom.CellSize();

    bool do_retry = false;

    // By default, we don't do a retry unless the criteria are violated.

#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_sub)
#endif
    for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const int* lo = bx.loVect();
        const int* hi = bx.hiVect();

        ca_check_timestep(ARLIM_3D(lo), ARLIM_3D(hi),
                          BL_TO_FORTRAN_3D(S_old[mfi]),
                          BL_TO_FORTRAN_3D(S_new[mfi]),
#ifdef REACTIONS
                          BL_TO_FORTRAN_3D(R_old[mfi]),
                          BL_TO_FORTRAN_3D(R_new[mfi]),
#endif
                          ZFILL(dx),
                          &dt, &dt_sub);

    }

    if (retry_neg_dens_factor > 0.0) {

        // Negative density criterion
        // Reset so that the desired maximum fractional change in density
        // is not larger than retry_neg_dens_factor.

        ParallelDescriptor::ReduceRealMin(frac_change);

        if (frac_change < 0.0)
            dt_sub = std::min(dt_sub, dt * -(retry_neg_dens_factor / frac_change));

    }

    ParallelDescriptor::ReduceRealMin(dt_sub);

    // Do the retry if the suggested timestep is smaller than the actual one.
    // A user-specified tolerance parameter can be used here to prevent
    // retries that are caused by small differences.

    if (dt_sub * (1.0 + retry_tolerance) < std::min(dt, dt_subcycle) || burn_success != 1) {

        do_retry = true;

        dt_subcycle = std::min(dt, dt_subcycle) * retry_subcycle_factor;

        if (verbose && ParallelDescriptor::IOProcessor()) {
            std::cout << std::endl;
            std::cout << "  Timestep " << dt << " rejected at level " << level << "." << std::endl;
            std::cout << "  Performing a retry, with subcycled timesteps of maximum length dt = " << dt_subcycle << std::endl;
            std::cout << std::endl;
        }

        // Restore the original values of the state data.        

        for (int k = 0; k < num_state_type; k++) {

            if (prev_state[k]->hasOldData())
                state[k].copyOld(*prev_state[k]);

            if (prev_state[k]->hasNewData())
                state[k].copyNew(*prev_state[k]);

        }

        // Reset the source term predictor.

        sources_for_hydro.setVal(0.0, NUM_GROW);

        // Clear the contribution to the fluxes from this step.

        for (int dir = 0; dir < 3; ++dir)
            fluxes[dir]->setVal(0.0);

        for (int dir = 0; dir < 3; ++dir)
            mass_fluxes[dir]->setVal(0.0);

#if (BL_SPACEDIM <= 2)
        if (!Geometry::IsCartesian())
            P_radial.setVal(0.0);
#endif

#ifdef RADIATION
        if (Radiation::rad_hydro_combined)
            for (int dir = 0; dir < BL_SPACEDIM; ++dir)
                rad_fluxes[dir]->setVal(0.0);
#endif

#ifndef SDC
        if (source_term_predictor == 1) {

            // Normally the source term predictor is done before the swap,
            // but the prev_state data is saved after the initial swap had
            // been done. So we will temporarily swap the state data back,
            // and reset the time levels.

            // Note that unlike the initial application of the source term
            // predictor before the swap, the old data will have already
            // been allocated when we get to this point. So we want to skip
            // this step if we didn't have old data initially.

            if (prev_state_had_old_data) {

                swap_state_time_levels(0.0);

                const Real dt_old = prev_state_new_time - prev_state_old_time;

                for (int k = 0; k < num_state_type; k++)
                    state[k].setTimeLevel(prev_state_new_time, dt_old, 0.0);

                apply_source_term_predictor();

                swap_state_time_levels(0.0);

                for (int k = 0; k < num_state_type; k++)
                    state[k].setTimeLevel(time + dt_subcycle, dt_subcycle, 0.0);

            }

        }
#endif

        if (track_grid_losses)
            for (int i = 0; i < n_lost; i++)
                material_lost_through_boundary_temp[i] = 0.0;

    }

    return do_retry;

}



Real
Castro::subcycle_advance(const Real time, const Real dt, int amr_iteration, int amr_ncycle)
{

    // Start the subcycle time off with the main dt,
    // unless we already came in here with an estimate
    // that is different from the initial value we assigned,
    // for example from the post-step regrid algorithm.

    if (dt_subcycle == 1.e200)
        dt_subcycle = dt;

    Real subcycle_time = time;

    Real dt_new = 1.e200;

    sub_iteration = 0;

    // Subcycle until we've reached the target time.
    // Compare against a slightly smaller number to
    // avoid roundoff concerns.

    Real eps = 1.0e-14;

    bool do_swap = false;

    while (subcycle_time < (1.0 - eps) * (time + dt)) {

        sub_iteration += 1;

        if (dt_subcycle < dt_cutoff) {
            if (ParallelDescriptor::IOProcessor()) {
                std::cout << std::endl;
                std::cout << "  The subcycle mechanism requested subcycled timesteps of maximum length dt = " << dt_subcycle << "," << std::endl
                          << "  but this timestep is shorter than the user-defined minimum, " << std::endl
                          << "  castro.dt_cutoff = " << dt_cutoff << ". Aborting." << std::endl;
            }
            amrex::Abort("Error: subcycled timesteps too short.");
        }

        // Shorten the last timestep so that we don't overshoot
        // the ending time.

        if (subcycle_time + dt_subcycle > (time + dt))
            dt_subcycle = (time + dt) - subcycle_time;

        // Check on whether we are going to take too many subcycles.

        int num_subcycles_remaining = int(round(((time + dt) - subcycle_time) / dt_subcycle));

        if (num_subcycles_remaining > max_subcycles) {
            amrex::Print() << std::endl
                           << "  The subcycle mechanism requested " << num_subcycles_remaining << " subcycled timesteps, which is larger than the maximum of " << max_subcycles << "." << std::endl
                           << "  If you would like to override this, increase the parameter castro.max_subcycles." << std::endl;
            amrex::Abort("Error: too many subcycles.");
        }

        // If we get to this point, we survived the sanity checks. Print out the current subcycle iteration.

        if (verbose && ParallelDescriptor::IOProcessor()) {
            std::cout << std::endl;
            std::cout << "  Beginning subcycle " << sub_iteration << " starting at time " << subcycle_time
                      << " with dt = " << dt_subcycle << std::endl;
            std::cout << "  Estimated number of subcycles remaining: " << num_subcycles_remaining << std::endl << std::endl;
        }

        // Swap the time levels. Only do this after the first iteration,
        // and when we are not doing a retry (which handles the swap).

        if (do_swap) {

            // Reset the source term predictor.
            // This must come before the swap.

            sources_for_hydro.setVal(0.0, NUM_GROW);

#ifndef SDC
            if (source_term_predictor == 1)
                apply_source_term_predictor();
#endif

            swap_state_time_levels(0.0);

#ifdef SELF_GRAVITY
            if (do_grav) {
                gravity->swapTimeLevels(level);
            }
#endif

        }

        // Assume we want to do a swap in the next iteration,
        // unless the retry tells us otherwise.

        do_swap = true;

        // Set the relevant time levels.

        for (int k = 0; k < num_state_type; k++)
            state[k].setTimeLevel(subcycle_time + dt_subcycle, dt_subcycle, 0.0);

        do_advance(subcycle_time, dt_subcycle, amr_iteration, amr_ncycle);

        if (verbose && ParallelDescriptor::IOProcessor()) {
            std::cout << std::endl;
            std::cout << "  Subcycle completed" << std::endl;
            std::cout << std::endl;
        }

        subcycle_time += dt_subcycle;

        // If we have hit a CFL violation during this subcycle, we must abort.

        if (cfl_violation && hard_cfl_limit && !use_retry)
            amrex::Abort("CFL is too high at this level, and we are already inside a retry -- go back to a checkpoint and restart with lower cfl number");

        if (burn_success != 1 && !use_retry)
            amrex::Abort("Burn was unsuccessful");

        // If we're allowing for retries, check for that here.

        if (use_retry) {

            // If we hit a retry, signal that we want to try again.
            // The retry function will handle resetting the state,
            // and updating dt_subcycle.

            if (retry_advance(subcycle_time, dt_subcycle, amr_iteration, amr_ncycle)) {
                do_swap = false;
                sub_iteration = 0;
                subcycle_time = time;
                lastDtRetryLimited = true;
                lastDtFromRetry = dt_subcycle;
            }

        }

    }

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "  Subcycling complete" << std::endl << std::endl;

    if (sub_iteration > 1) {

        // Finally, copy the original data back to the old state
        // data so that externally it appears like we took only
        // a single timestep. We'll do this as a swap so that
        // we still have the last iteration's old data if we need
        // it later.

        for (int k = 0; k < num_state_type; k++) {

            if (prev_state[k]->hasOldData())
                state[k].replaceOldData(*prev_state[k]);

            state[k].setTimeLevel(time + dt, dt, 0.0);
            prev_state[k]->setTimeLevel(time + dt, dt_subcycle, 0.0);

        }

        // If we took more than one step and are going to do a reflux,
        // keep the data past the end of the step.

        if (do_reflux && update_sources_after_reflux) {

            // Note that since we only want to do this if there's actually a
            // reflux immediately following this, skip this if we're on the
            // finest level and this is not the last iteration.

            if (!(amr_iteration < amr_ncycle && level == parent->finestLevel()))
                keep_prev_state = true;

        }

    }

    // We want to return the subcycled timestep as a suggestion.

    dt_new = std::min(dt_new, dt_subcycle);

    return dt_new;

}
#endif
