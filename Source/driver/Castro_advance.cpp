
#include "Castro.H"
#include "Castro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

#ifdef SELF_GRAVITY
#include "Gravity.H"
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

    MultiFab::RegionTag amrlevel_tag("AmrLevel_Level_" + std::to_string(level));

    Real dt_new = dt;

    initialize_advance(time, dt, amr_iteration, amr_ncycle);

    // Do the advance.

    if (time_integration_method == CornerTransportUpwind) {

        dt_new = std::min(dt_new, subcycle_advance_ctu(time, dt, amr_iteration, amr_ncycle));

#ifndef AMREX_USE_CUDA
    } else if (time_integration_method == SpectralDeferredCorrections) {

      for (int iter = 0; iter < sdc_order+sdc_extra; ++iter) {
	sdc_iteration = iter;
	dt_new = do_advance_sdc(time, dt, amr_iteration, amr_ncycle);
      }

#ifdef REACTIONS
      // store the reaction information as well -- note: this will be
      // the instantaneous reactive source.  In the future, we might
      // want to do a quadrature over R_new[]

      // this is done only for the plotfile
      MultiFab& R_new = get_new_data(Reactions_Type);
      MultiFab& S_new = get_new_data(State_Type);

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
    else if (time_integration_method == SimplifiedSpectralDeferredCorrections) {

        for (int n = 0; n < sdc_iters; ++n) {

            sdc_iteration = n;

	    amrex::Print() << "Beginning SDC iteration " << n + 1 << " of " << sdc_iters << "." << std::endl << std::endl;

            // First do the non-reacting advance and construct the relevant source terms.
            // We use the CTU advance here, with the Strang-split reactions skipped,
            // but we call do_advance_ctu directly rather than subcycle_advance_ctu,
            // as the simplified SDC logic is not compatible with the subcycling.

            dt_new = do_advance_ctu(time, dt, amr_iteration, amr_ncycle);

#ifdef REACTIONS
            if (do_react) {

                // Do the ODE integration to capture the reaction source terms.

                react_state(time, dt);

                MultiFab& S_new = get_new_data(State_Type);

                clean_state(S_new, state[State_Type].curTime(), S_new.nGrow());

                // Compute the reactive source term for use in the next iteration.

                MultiFab& SDC_react_new = get_new_data(Simplified_SDC_React_Type);
                get_react_source_prim(SDC_react_new, time, dt);

                // Check for NaN's.

                check_for_nan(S_new);

            }
#endif

            amrex::Print() << "Ending SDC iteration " << n + 1 << " of " << sdc_iters << "." << std::endl << std::endl;

        }

#endif // AMREX_USE_CUDA
    }

    // Optionally kill the job at this point, if we've detected a violation.

    if (cfl_violation && hard_cfl_limit && !use_retry)
        amrex::Abort("CFL is too high at this level -- go back to a checkpoint and restart with lower cfl number");

    // If we didn't kill the job, reset the violation counter.

    cfl_violation = 0;

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

#ifdef GRAVITY
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



void
Castro::initialize_do_advance(Real time, Real dt, int amr_iteration, int amr_ncycle)
{
    BL_PROFILE("Castro::initialize_do_advance()");

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

    if (time_integration_method == CornerTransportUpwind && source_term_predictor == 1) {
        sources_for_hydro.mult(0.5 * dt, NUM_GROW);
    }

    // For the hydrodynamics update we need to have NUM_GROW ghost
    // zones available, but the state data does not carry ghost
    // zones. So we use a FillPatch using the state data to give us
    // Sborder, which does have ghost zones.

    MultiFab& S_old = get_old_data(State_Type);

    if (time_integration_method == CornerTransportUpwind || time_integration_method == SimplifiedSpectralDeferredCorrections) {
      // for the CTU unsplit method, we always start with the old state
      Sborder.define(grids, dmap, NUM_STATE, NUM_GROW, MFInfo().SetTag("Sborder"));
      const Real prev_time = state[State_Type].prevTime();
      clean_state(S_old, prev_time, 0);
      expand_state(Sborder, prev_time, NUM_GROW);

    } else if (time_integration_method == SpectralDeferredCorrections) {

      // we'll handle the filling inside of do_advance_sdc 
      Sborder.define(grids, dmap, NUM_STATE, NUM_GROW, MFInfo().SetTag("Sborder"));

    } else {
      amrex::Abort("invalid time_integration_method");
    }

#ifdef SHOCK_VAR
    // Zero out the shock data, and fill it during the advance.
    // For subcycling cases this will always give the shock
    // variable for the latest subcycle, rather than averaging.

    Sborder.setVal(0.0, Shock, 1, Sborder.nGrow());
#endif

}



void
Castro::finalize_do_advance(Real time, Real dt, int amr_iteration, int amr_ncycle)
{
    BL_PROFILE("Castro::finalize_do_advance()");

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
    BL_PROFILE("Castro::initialize_advance()");

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

    // Add the source term predictor.
    // This must happen before the swap.

    if (time_integration_method == CornerTransportUpwind && source_term_predictor == 1) {
        apply_source_term_predictor();
    }

    // If we're doing simplified SDC, time-center the source term (using the
    // current iteration's old sources and the last iteration's new
    // sources). Since the "new-time" sources are just the corrector step
    // of the predictor-corrector formalism, we want to add the full
    // value of the "new-time" sources to the old-time sources to get a
    // time-centered value.

    if (time_integration_method == SimplifiedSpectralDeferredCorrections) {
        AmrLevel::FillPatch(*this, sources_for_hydro, NUM_GROW, time, Source_Type, 0, NUM_STATE);
    }

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

    MultiFab& S_old = get_old_data(State_Type);
    clean_state(S_old, time, S_old.nGrow());

    // Initialize the previous state data container now, so that we can
    // always ask if it has valid data.

    for (int k = 0; k < num_state_type; ++k)
        prev_state[k].reset(new StateData());

    // Make a copy of the MultiFabs in the old and new state data in case we may do a retry.

    if (use_retry || do_subcycle) {

      // Store the old and new time levels.

      for (int k = 0; k < num_state_type; k++) {
        *prev_state[k] = state[k];
      }

    }

    // This array holds the hydrodynamics update.
    if (time_integration_method == CornerTransportUpwind || time_integration_method == SimplifiedSpectralDeferredCorrections) {
      hydro_source.define(grids,dmap,NUM_STATE,0);
    }


    // Allocate space for the primitive variables.

    q.define(grids, dmap, NQ, NUM_GROW);
    q.setVal(0.0);
    qaux.define(grids, dmap, NQAUX, NUM_GROW);

    if (time_integration_method == CornerTransportUpwind || time_integration_method == SimplifiedSpectralDeferredCorrections) {
      src_q.define(grids, dmap, NQSRC, NUM_GROW);
    }

    if (sdc_order == 4) {
      q_bar.define(grids, dmap, NQ, NUM_GROW);
      qaux_bar.define(grids, dmap, NQAUX, NUM_GROW);
#ifdef DIFFUSION
      T_cc.define(grids, dmap, 1, NUM_GROW);
#endif
    }


    if (time_integration_method == SpectralDeferredCorrections) {

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
      // for the temporary storage of the reaction terms
      Sburn.define(grids, dmap, NUM_STATE, 2);

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
    if (!Geom().IsCartesian())
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
    BL_PROFILE("Castro::finalize_advance()");

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

    if (time_integration_method == CornerTransportUpwind || time_integration_method == SimplifiedSpectralDeferredCorrections) {
      hydro_source.clear();
    }

    q.clear();
    qaux.clear();

    if (time_integration_method == CornerTransportUpwind || time_integration_method == SimplifiedSpectralDeferredCorrections) {
      src_q.clear();
    }

    if (sdc_order == 4) {
      q_bar.clear();
      qaux_bar.clear();
#ifdef DIFFUSION
      T_cc.clear();
#endif
    }

#ifdef RADIATION
    Erborder.clear();
    lamborder.clear();
#endif

    sources_for_hydro.clear();

    if (!keep_prev_state)
        amrex::FillNull(prev_state);

    if (time_integration_method == SpectralDeferredCorrections) {
      k_new.clear();
      A_new.clear();
      A_old.clear();
#ifdef REACTIONS
      R_old.clear();
      Sburn.clear();
#endif
    }

    // Record how many zones we have advanced.

    num_zones_advanced += grids.numPts() / getLevel(0).grids.numPts();

}
