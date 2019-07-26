
#include "Castro.H"
#include "Castro_F.H"

#ifdef SELF_GRAVITY
#include "Gravity.H"
#endif

using namespace amrex;

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

  const int* domain_lo = geom.Domain().loVect();
  const int* domain_hi = geom.Domain().hiVect();

  // Perform initialization steps.

  initialize_do_advance(time, dt, amr_iteration, amr_ncycle);

#ifndef AMREX_USE_CUDA
  // Check for NaN's.

  check_for_nan(S_old);
#endif

  // Since we are Strang splitting the reactions, do them now (only
  // for first stage of MOL)

  if (mol_iteration == 0) {

#ifdef REACTIONS
    // this operates on Sborder (which is initially S_old).  The result
    // of the reactions is added directly back to Sborder.
    strang_react_first_half(prev_time, 0.5 * dt);
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
    if (mol_order == 4) {
      // if we are 4th order, convert to cell-center Sborder -> Sborder_cc
      // we'll reuse sources_for_hydro for this memory buffer at the moment

      for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {
        const Box& gbx = mfi.growntilebox(1);
        ca_make_cell_center(BL_TO_FORTRAN_BOX(gbx),
                            BL_TO_FORTRAN_FAB(Sborder[mfi]),
                            BL_TO_FORTRAN_FAB(sources_for_hydro[mfi]),
                            AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      }
    }

    // we pass in the stage time here
    if (mol_order == 4) {
      // time here is the stage time
      do_old_sources(old_source, sources_for_hydro, time, dt, amr_iteration, amr_ncycle);

      // fill the ghost cells for the sources -- we are storing these
      // in the "old" time slot of Source_Type, so we should only use
      // that date -- note this is not multilevel save
      AmrLevel::FillPatch(*this, old_source, old_source.nGrow(), prev_time, Source_Type, 0, NUM_STATE);

      // Convert to cell averages.  This loop cannot be tiled.
      for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        ca_make_fourth_in_place(BL_TO_FORTRAN_BOX(bx),
                                BL_TO_FORTRAN_FAB(old_source[mfi]),
                                AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      }

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
      if (mol_order == 4) {
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
  clean_state(S_new, cur_time, S_new.nGrow());

  // If the state has ghost zones, sync them up now
  // since the hydro source only works on the valid zones.

  if (S_new.nGrow() > 0) {
      clean_state(S_new, cur_time, 0);
      expand_state(S_new, cur_time, S_new.nGrow());
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
  clean_state(S_old, prev_time, 0);
  expand_state(Sborder, prev_time, Sborder.nGrow());
  do_old_sources(old_source, Sborder, prev_time, dt, amr_iteration, amr_ncycle);
  AmrLevel::FillPatch(*this, old_source, old_source.nGrow(), prev_time, Source_Type, 0, NUM_STATE);

  clean_state(S_new, cur_time, 0);
  expand_state(Sborder, cur_time, Sborder.nGrow());
  do_old_sources(new_source, Sborder, cur_time, dt, amr_iteration, amr_ncycle);
  AmrLevel::FillPatch(*this, old_source, old_source.nGrow(), cur_time, Source_Type, 0, NUM_STATE);

  // Do the second half of the reactions.

#ifdef REACTIONS
  strang_react_second_half(cur_time - 0.5 * dt, 0.5 * dt);
#endif

  finalize_do_advance(time, dt, amr_iteration, amr_ncycle);

  return dt;
}
