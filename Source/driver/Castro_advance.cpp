
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

    Real dt_new = dt;

    initialize_advance(time, dt, amr_iteration, amr_ncycle);

    // Do the advance.

#ifdef SDC

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

	    clean_state(S_new);

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
    // no SDC

    if (do_ctu) {

        if (do_subcycle) {

            dt_new = std::min(dt_new, subcycle_advance(time, dt, amr_iteration, amr_ncycle));
            
        }
        else {

            dt_new = do_advance(time, dt, amr_iteration, amr_ncycle);

        }

    } else {
      for (int iter = 0; iter < MOL_STAGES; ++iter) {
	mol_iteration = iter;
	dt_new = do_advance(time, dt, amr_iteration, amr_ncycle);
      }
    }

    // Check to see if this advance violated certain stability criteria.
    // If so, get a new timestep and do subcycled advances until we reach
    // t = time + dt.

    if (use_retry)
        dt_new = std::min(dt_new, retry_advance(time, dt, amr_iteration, amr_ncycle));
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
    pointmass_update(time, dt);
#endif

#ifdef RADIATION
    MultiFab& S_new = get_new_data(State_Type);
    final_radiation_call(S_new, amr_iteration, amr_ncycle);
#endif

#ifdef PARTICLES
    advance_particles(amr_iteration, time, dt);
#endif

    finalize_advance(time, dt, amr_iteration, amr_ncycle);

    return dt_new;
}

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

    // Since we are Strang splitting the reactions, do them now (only
    // for first stage of MOL)

    if (do_ctu || (!do_ctu && mol_iteration == 0)) {

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


      // Construct the old-time sources from Sborder.  For both CTU
      // and MOL integration, this will already be applied to S_new
      // (with full dt weighting), to be correctly later.  Note --
      // this does not affect the prediction of the interface state,
      // an explict source will be traced there as needed.

#ifdef SELF_GRAVITY
      construct_old_gravity(amr_iteration, amr_ncycle, prev_time);
#endif

      do_old_sources(prev_time, dt, amr_iteration, amr_ncycle);

      if (!do_ctu) {
	// store the result of the burn and old-time sources in Sburn for later stages
	MultiFab::Copy(Sburn, S_new, 0, 0, NUM_STATE, S_new.nGrow());
      }

    }

    // Do the hydro update.  We build directly off of Sborder, which
    // is the state that has already seen the burn

    if (do_hydro)
    {
      if (do_ctu) {
        construct_hydro_source(time, dt);
	apply_source_to_state(S_new, hydro_source, dt);
      } else {
        construct_mol_hydro_source(time, dt);
      }
    }

    // For MOL integration, we are done with this stage, unless it is
    // the last stage

    if (!do_ctu && mol_iteration == MOL_STAGES-1) {
      // we just finished the last stage of the MOL integration.
      // Construct S_new now using the weighted sum of the k_mol
      // updates

      // Apply the update -- we need to build on Sburn, so
      // start with that state
      MultiFab::Copy(S_new, Sburn, 0, 0, S_new.nComp(), 0);
      MultiFab::Saxpy(S_new, dt, hydro_source, 0, 0, S_new.nComp(), 0);

      // define the temperature now
      clean_state(S_new);

    } else if (do_ctu) {

      // Sync up state after old sources and hydro source.

      frac_change = clean_state(S_new, Sborder);

      // Check for NaN's.

      check_for_nan(S_new);
    }


    // if we are done with the update do the source correction and
    // then the second half of the reactions

    if (do_ctu || mol_iteration == MOL_STAGES-1) {

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

      do_new_sources(cur_time, dt, amr_iteration, amr_ncycle);


      // Do the second half of the reactions.

      // last part of reactions for CTU and if we are done with the
      // MOL stages

#ifdef REACTIONS
#ifndef SDC
      strang_react_second_half(cur_time - 0.5 * dt, 0.5 * dt);
#endif
#endif
      }

    finalize_do_advance(time, dt, amr_iteration, amr_ncycle);

    return dt;

}



void
Castro::initialize_do_advance(Real time, Real dt, int amr_iteration, int amr_ncycle)
{

    // Reset the change from density resets

    frac_change = 1.e0;

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

    // For the hydrodynamics update we need to have NUM_GROW ghost zones available,
    // but the state data does not carry ghost zones. So we use a FillPatch
    // using the state data to give us Sborder, which does have ghost zones.

    if (do_ctu) {
      // for the CTU unsplit method, we always start with the old state
      Sborder.define(grids, dmap, NUM_STATE, NUM_GROW);
      const Real prev_time = state[State_Type].prevTime();
      expand_state(Sborder, prev_time, NUM_GROW);

    } else {
      // for Method of lines, our initialization of Sborder depends on
      // which stage in the RK update we are working on

      if (mol_iteration == 0) {

	// first MOL stage
	Sborder.define(grids, dmap, NUM_STATE, NUM_GROW);
	const Real prev_time = state[State_Type].prevTime();
	expand_state(Sborder, prev_time, NUM_GROW);

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

	Sborder.define(grids, dmap, NUM_STATE, NUM_GROW);
	const Real new_time = state[State_Type].curTime();
	expand_state(Sborder, new_time, NUM_GROW);

      }
    }
}



void
Castro::finalize_do_advance(Real time, Real dt, int amr_iteration, int amr_ncycle)
{

#ifndef SDC
    // Update the dSdt MultiFab. Since we want (S^{n+1} - S^{n}) / dt,
    // we only need to take twice the new-time source term, since in
    // the predictor-corrector approach, the new-time source term is
    // 1/2 * S^{n+1} - 1/2 * S^{n}. This is untrue in general for the
    // non-momentum sources, but those don't appear in the hydro
    // anyway, and for safety we'll only do this on the momentum
    // terms.

    if (source_term_predictor == 1) {

        MultiFab& dSdt_new = get_new_data(Source_Type);

	dSdt_new.setVal(0.0, NUM_GROW);

        MultiFab::Add(dSdt_new, new_sources, Xmom, Xmom, 3, 0);

	dSdt_new.mult(2.0 / dt);

    }
#else
    // The new sources are broken into types (ext, diff, hybrid, grav,
    // ...) via an enum.  For SDC, store the sum of the new_sources
    // over these different physics types in the state data -- that's
    // what hydro really cares about.

    MultiFab& SDC_source_new = get_new_data(SDC_Source_Type);
    SDC_source_new.setVal(0.0, SDC_source_new.nGrow());
    MultiFab::Add(SDC_source_new, new_sources, 0, 0, NUM_STATE, new_sources.nGrow());
#endif

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
#endif

#ifdef SELF_GRAVITY
    // If we're on level 0, update the maximum density used in the gravity solver
    // for setting the tolerances. This will be used in all level solves to follow.
    // This must be done before the swap because it relies on the new data.

    if (level == 0 && gravity->get_gravity_type() == "PoissonGrav") {
	gravity->update_max_rhs();
    }
#endif

    // Swap the new data from the last timestep into the old state data.

    for (int k = 0; k < num_state_type; k++) {

	// The following is a hack to make sure that we only
	// ever have new data for a few state types that only
	// ever need new time data; by doing a swap now, we'll
	// guarantee that allocOldData() does nothing. We do
	// this because we never need the old data, so we
	// don't want to allocate memory for it.

	if (k == Source_Type)
	    state[k].swapTimeLevels(0.0);
#ifdef SDC
	else if (k == SDC_Source_Type)
	    state[k].swapTimeLevels(0.0);
#ifdef REACTIONS
	else if (k == SDC_React_Type)
	    state[k].swapTimeLevels(0.0);
#endif
#endif

	state[k].allocOldData();
	state[k].swapTimeLevels(dt);

    }

#ifdef SELF_GRAVITY
    if (do_grav)
	gravity->swapTimeLevels(level);
#endif

    // Ensure data is valid before beginning advance. This addresses
    // the fact that we may have new data on this level that was interpolated
    // from a coarser level, and the interpolation in general cannot be
    // trusted to respect the consistency between certain state variables
    // (e.g. UEINT and UEDEN) that we demand in every zone.

    clean_state(get_old_data(State_Type));

    // Make a copy of the MultiFabs in the old and new state data in case we may do a retry.

    if (use_retry) {

      // Store the old and new time levels.

      for (int k = 0; k < num_state_type; k++) {

	prev_state[k].reset(new StateData());

	StateData::Initialize(*prev_state[k], state[k]);

      }

    }

    if (!(do_reflux && update_sources_after_reflux)) {

      // These arrays hold all source terms that update the state.
      // Normally these are allocated at initialization because
      // of the source term update after a reflux, but if the user
      // chooses not to do that, we can save memory by only allocating
      // the data temporarily for the duration of the advance.

      old_sources.define(grids, dmap, NUM_STATE, get_new_data(State_Type).nGrow());
      new_sources.define(grids, dmap, NUM_STATE, get_new_data(State_Type).nGrow());

      // This array holds the hydrodynamics update.

      hydro_source.define(grids,dmap,NUM_STATE,0);

    }

    // This array holds the sum of all source terms that affect the
    // hydrodynamics.  If we are doing the source term predictor,
    // we'll also use this after the hydro update to store the sum of
    // the new-time sources, so that we can compute the time
    // derivative of the source terms.

    sources_for_hydro.define(grids,dmap,NUM_STATE,NUM_GROW);


    if (!do_ctu) {
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

    // Zero out the current fluxes.

    for (int dir = 0; dir < 3; ++dir)
	fluxes[dir]->setVal(0.0);

#if (BL_SPACEDIM <= 2)
    if (!Geometry::IsCartesian())
	P_radial.setVal(0.0);
#endif

#ifdef RADIATION
    if (Radiation::rad_hydro_combined)
	for (int dir = 0; dir < BL_SPACEDIM; ++dir)
	    rad_fluxes[dir]->setVal(0.0);
#endif

    mass_fluxes.resize(3);

    for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
	mass_fluxes[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, 1, 0));
        mass_fluxes[dir]->setVal(0.0);
    }

    for (int dir = BL_SPACEDIM; dir < 3; ++dir) {
	mass_fluxes[dir].reset(new MultiFab(get_new_data(State_Type).boxArray(), dmap, 1, 0));
        mass_fluxes[dir]->setVal(0.0);
    }

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

    if (!(do_reflux && update_sources_after_reflux)) {

	old_sources.clear();
	new_sources.clear();
	hydro_source.clear();

    }

    sources_for_hydro.clear();

    amrex::FillNull(prev_state);

    if (!do_ctu) {
      k_mol.clear();
      Sburn.clear();
    }

    for (int dir = 0; dir < 3; ++dir)
	mass_fluxes[dir].reset();

}



Real
Castro::retry_advance(Real time, Real dt, int amr_iteration, int amr_ncycle)
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

#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_sub)
#endif
    for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

	const int* lo = bx.loVect();
	const int* hi = bx.hiVect();

	ca_check_timestep(BL_TO_FORTRAN_3D(S_old[mfi]),
			  BL_TO_FORTRAN_3D(S_new[mfi]),
#ifdef REACTIONS
			  BL_TO_FORTRAN_3D(R_old[mfi]),
			  BL_TO_FORTRAN_3D(R_new[mfi]),
#endif
			  ARLIM_3D(lo), ARLIM_3D(hi), ZFILL(dx),
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

    if (dt_sub * (1.0 + retry_tolerance) < std::min(dt, dt_subcycle)) {

        dt_subcycle = dt_sub;

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

	if (track_grid_losses)
	  for (int i = 0; i < n_lost; i++)
	    material_lost_through_boundary_temp[i] = 0.0;

        dt_new = std::min(dt_new, subcycle_advance(time, dt, amr_iteration, amr_ncycle));

    }

    return dt_new;

}



Real
Castro::subcycle_advance(const Real time, const Real dt, int amr_iteration, int amr_ncycle)
{

    // Check the input number of subcycles against various safety conditions.

    // Do a basic sanity check to make sure we're not about to overflow.

    if (dt_subcycle * INT_MAX < dt) {
        if (ParallelDescriptor::IOProcessor()) {
	    std::cout << std::endl;
	    std::cout << "  The subcycle mechanism requested subcycled timesteps of maximum length dt = " << dt_subcycle << "," << std::endl
                      << "  but this would imply a number of timesteps that overflows an integer." << std::endl;
	    std::cout << "  The code will abort. Consider using a stronger stability criterion" << std::endl
                      << "  to avoid unstable timesteps." << std::endl;
        }
        amrex::Abort("Error: integer overflow in subcycling.");
    }

    sub_ncycle = ceil(dt / dt_subcycle);

    // Abort if we would take more subcycled timesteps than the user has permitted.

    if (max_subcycles > 0 && sub_ncycle > max_subcycles) {
        if (clamp_subcycles) {
            if (ParallelDescriptor::IOProcessor()) {
                std::cout << std::endl;
                std::cout << "  The subcycle mechanism requested " << sub_ncycle << " subcycled timesteps of maximum length dt = " << dt_subcycle << "," << std::endl
                          << "  but this is more than the maximum number of permitted subcycles, " << max_subcycles << "." << std::endl;
                std::cout << "  The code will decrease the number of subcycles to this maximum, and the subcycled timestep" << std::endl
                          << "  will be updated to dt = " << dt / max_subcycles << "." << std::endl;
                std::cout << "  If you want more subcycles you can update the max_subcycles parameter." << std::endl;
                std::cout << "  If you would prefer the simulation to stop when this maximum number is" << std::endl
                          << "  exceeded, set castro.clamp_subcycles = 0." << std::endl;
            }

            sub_ncycle = max_subcycles;
        } else {
            if (ParallelDescriptor::IOProcessor()) {
                std::cout << std::endl;
                std::cout << "  The subcycle mechanism requested " << sub_ncycle << " subcycled timesteps of maximum length dt = " << dt_subcycle << "," << std::endl
                          << "  but this is more than the maximum number of permitted subcycles, " << max_subcycles << "." << std::endl;
                std::cout << "  The code will abort. Consider using a stronger stability criterion" << std::endl
                          << "  to avoid unstable timesteps, or consider increasing the parameter " << std::endl
                          << "  castro.max_subcycles to permit more subcycled timesteps." << std::endl;
            }
            amrex::Abort("Error: too many subcycled timesteps.");
        }
    }

    // Abort if our subcycled timestep would be shorter than the minimum permitted timestep.

    if (dt_subcycle < dt_cutoff) {
        if (ParallelDescriptor::IOProcessor()) {
	    std::cout << std::endl;
	    std::cout << "  The subcycle mechanism requested " << sub_ncycle << " subcycled timesteps of maximum length dt = " << dt_subcycle << "," << std::endl
                      << "  but this timestep is shorter than the user-defined minimum, " << std::endl
                      << "  castro.dt_cutoff = " << dt_cutoff << ". Aborting." << std::endl;
	  }
	  amrex::Abort("Error: subcycled timesteps too short.");
    }

    // Subcycle until we've reached the target time.
    // Compare against a slightly smaller number to
    // avoid roundoff concerns.

    Real eps = 1.0e-14;

    Real subcycle_time = time;
    dt_subcycle = dt / sub_ncycle;
    Real dt_advance = dt_subcycle;

    Real dt_new = 1.e200;

    sub_iteration = 1;

    while (subcycle_time < (1.0 - eps) * (time + dt)) {

        // Shorten the last timestep so that we don't overshoot
        // the ending time.

        if (subcycle_time + dt_advance > (time + dt))
            dt_advance = (time + dt) - subcycle_time;

        if (verbose && ParallelDescriptor::IOProcessor()) {
            std::cout << std::endl;
            std::cout << "  Beginning subcycle " << sub_iteration << " of " << sub_ncycle
                      << ", starting at time " << subcycle_time
                      << " with dt = " << dt_advance << std::endl << std::endl;
        }

        // Swap the time levels. Only do this after the first iteration;
        // we assume that the values coming into this routine are
        // consistent in the sense that the "old" data is really the time
        // level n data. For a normal advance this is guaranteed by the
        // swap in initialize_advance and for a retry this is done prior
        // to entering this routine.

        if (sub_iteration > 1) {

            for (int k = 0; k < num_state_type; k++) {

                if (k == Source_Type)
                    state[k].swapTimeLevels(0.0);
#ifdef SDC
                else if (k == SDC_Source_Type)
                    state[k].swapTimeLevels(0.0);
#ifdef REACTIONS
                else if (k == SDC_React_Type)
                    state[k].swapTimeLevels(0.0);
#endif
#endif

                state[k].swapTimeLevels(0.0);

            }

#ifdef SELF_GRAVITY
            if (do_grav) {
                gravity->swapTimeLevels(level);
            }
#endif

        }

        // Set the relevant time levels.

        for (int k = 0; k < num_state_type; k++)
            state[k].setTimeLevel(subcycle_time + dt_advance, dt_advance, 0.0);

        do_advance(subcycle_time, dt_advance, amr_iteration, amr_ncycle);

        if (verbose && ParallelDescriptor::IOProcessor()) {
            std::cout << std::endl;
            std::cout << "  Subcycle " << sub_iteration << " of " << sub_ncycle << " completed" << std::endl;
            std::cout << std::endl;
        }

        subcycle_time += dt_advance;
        sub_iteration += 1;

    }

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "  Subcycling complete" << std::endl << std::endl;

    // Finally, copy the original data back to the old state
    // data so that externally it appears like we took only
    // a single timestep.

    for (int k = 0; k < num_state_type; k++) {

        if (prev_state[k]->hasOldData())
            state[k].copyOld(*prev_state[k]);

        state[k].setTimeLevel(time + dt, dt, 0.0);

    }

    // We want to return the subcycled timestep as a suggestion.

    dt_new = std::min(dt_new, dt_subcycle);

    return dt_new;

}
