#include <winstd.H>

#include "Castro.H"
#include "Castro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

#ifdef GRAVITY
#include "Gravity.H"
#endif

#ifdef DIFFUSION
#include "Diffusion.H"
#endif

#ifdef LEVELSET
#include "LevelSet_F.H"
#endif

#ifdef ROTATION
#include "Rotation.H"
#endif

#include <cmath>

using std::string;

Real
Castro::advance (Real time,
                 Real dt,
                 int  amr_iteration,
                 int  amr_ncycle)
{
    BL_PROFILE("Castro::advance()");

    Real dt_new = dt;

    initialize_advance(time, dt, amr_iteration, amr_ncycle);

    // Do the advance.

#ifdef SDC

    for (int n = 0; n < sdc_iters; ++n) {

        if (ParallelDescriptor::IOProcessor())
	    std::cout << "\nBeginning SDC iteration " << n + 1 << " of " << sdc_iters << ".\n\n";

	// First do the non-reacting advance and construct the relevant source terms.

	dt_new = do_advance(time, dt, amr_iteration, amr_ncycle, n, sdc_iters);

#ifdef REACTIONS
	if (do_react) {

            // Do the ODE integration to capture the reaction source terms.

	    react_state(time, dt);

	    MultiFab& S_new = get_new_data(State_Type);

	    clean_state(S_new);

	    // Compute the reactive source term for use in the next iteration.

	    get_react_source_prim(*react_src, dt);

	    // Check for NaN's.

	    check_for_nan(S_new);

        }
#endif

        if (ParallelDescriptor::IOProcessor())
	    std::cout << "\nEnding SDC iteration " << n + 1 << " of " << sdc_iters << ".\n\n";

    }

#else

    int sub_iteration = 0;
    int sub_ncycle = 0;

    dt_new = do_advance(time, dt, amr_iteration, amr_ncycle, sub_iteration, sub_ncycle);

#endif

    // Check to see if this advance violated certain stability criteria.
    // If so, get a new timestep and do subcycled advances until we reach
    // t = time + dt.

    if (use_retry)
        dt_new = std::min(dt_new, retry_advance(time, dt, amr_iteration, amr_ncycle));

#ifdef AUX_UPDATE
    advance_aux(time, dt);
#endif

#ifdef LEVELSET
    advance_levelset(time, dt);
#endif

#if (BL_SPACEDIM > 1)
    // We do this again here because the solution will have changed
    if ( (level == 0) && (spherical_star == 1) ) {
       int is_new = 1;
       make_radial_data(is_new);
    }
#endif

    // Update the point mass.

#ifdef POINTMASS
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
                    int  amr_ncycle,
                    int  sub_iteration,
                    int  sub_ncycle)
{
    BL_PROFILE("Castro::do_advance()");

    const Real prev_time = state[State_Type].prevTime();
    const Real  cur_time = state[State_Type].curTime();

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    // Perform initialization steps.

    initialize_do_advance(time, dt, amr_iteration, amr_ncycle, sub_iteration, sub_ncycle);

    // Check for NaN's.

    check_for_nan(S_old);

    // Since we are Strang splitting the reactions, do them now.

#ifdef REACTIONS
#ifndef SDC
    strang_react_first_half(prev_time, 0.5 * dt);
#endif
#endif

    // Initialize the new-time data. This copy needs to come after the reactions.

    MultiFab::Copy(S_new, *Sborder, 0, 0, NUM_STATE, S_new.nGrow());

    // Construct the old-time sources.

    for (int n = 0; n < num_src; ++n)
	construct_old_source(n, amr_iteration, amr_ncycle, sub_iteration, sub_ncycle, prev_time, dt);

    // Apply the old-time sources.

    for (int n = 0; n < num_src; ++n)
        apply_source_to_state(S_new, old_sources[n], dt);

    // Do the hydro update.

    if (do_hydro)
    {
        construct_hydro_source(time, dt);
	apply_source_to_state(S_new, *hydro_source, dt);
	frac_change = clean_state(S_new);
    }

    // Check for NaN's.

    check_for_nan(S_new);

    // Must define new value of "center" before we call new gravity solve or external source routine

#ifdef GRAVITY
    if (moving_center == 1)
        define_new_center(S_new, time);
#endif

    // We need to make the new radial data now so that we can use it when we
    // FillPatch in creating the new source.

#if (BL_SPACEDIM > 1)
    if ( (level == 0) && (spherical_star == 1) ) {
        int is_new = 1;
	make_radial_data(is_new);
    }
#endif

    // Construct the new-time source terms.

    for (int n = 0; n < num_src; ++n)
	construct_new_source(n, amr_iteration, amr_ncycle, sub_iteration, sub_ncycle, cur_time, dt);

    // Apply the new-time sources to the state.

    for (int n = 0; n < num_src; ++n)
        apply_source_to_state(S_new, new_sources[n], dt);

    // Sync up the temperature now that all sources have been applied.

    computeTemp(S_new);

    // Do the second half of the reactions.

#ifdef REACTIONS
#ifndef SDC
    strang_react_second_half(cur_time - 0.5 * dt, 0.5 * dt);
#endif
#endif

    finalize_do_advance(time, dt, amr_iteration, amr_ncycle, sub_iteration, sub_ncycle);

    return dt;

}



void
Castro::initialize_do_advance(Real time, Real dt, int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle)
{

    // Reset the change from density resets

    frac_change = 1.e0;

    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level && sub_iteration < 2) {
        //
        // Set reflux registers to zero.
        //
        getFluxReg(level+1).setVal(0.0);
#ifdef RADIATION
	if (Radiation::rad_hydro_combined) {
	  getRADFluxReg(level+1).setVal(0.0);
	}
#endif
    }

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

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
	u_gdnv[dir]->setVal(1.e40,1);
    }

    // Reset the grid loss tracking.

    if (track_grid_losses)
      for (int i = 0; i < n_lost; i++)
	material_lost_through_boundary_temp[i] = 0.0;

#ifdef GRAVITY
    if (moving_center == 1)
        define_new_center(get_old_data(State_Type), time);
#endif

#if (BL_SPACEDIM > 1)
    if ( (level == 0) && (spherical_star == 1) ) {
       swap_outflow_data();
       int is_new = 0;
       make_radial_data(is_new);
    }
#endif

#ifdef DIFFUSION
#ifdef TAU
    tau_diff->setVal(0.);
    define_tau(grav_old,time);
#endif
#endif

    for (int j = 0; j < 3; j++)
        fluxes[j]->setVal(0.0);

    hydro_source->setVal(0.0);

    // For the hydrodynamics update we need to have NUM_GROW ghost zones available,
    // but the state data does not carry ghost zones. So we use a FillPatch
    // using the state data to give us Sborder, which does have ghost zones.

    Sborder = new MultiFab(grids, NUM_STATE, NUM_GROW, Fab_allocate);
    const Real prev_time = state[State_Type].prevTime();
    expand_state(*Sborder, prev_time, NUM_GROW);

}



void
Castro::finalize_do_advance(Real time, Real dt, int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle)
{

#ifndef SDC
    // Update the dSdt MultiFab. Take the difference of the old and new
    // sources and then divide by dt.

    if (source_term_predictor == 1) {

        MultiFab& dSdt_new = get_new_data(Source_Type);

	dSdt_new.setVal(0.0, NUM_GROW);

	for (int n = 0; n < num_src; ++n) {
	    MultiFab::Add(dSdt_new, new_sources[n], 0, 0, NUM_STATE, 0);
	    MultiFab::Subtract(dSdt_new, old_sources[n], 0, 0, NUM_STATE, 0);
	}

	dSdt_new.mult(1.0 / dt);

    }
#endif

    // Sync up the hybrid and linear momenta.

#ifdef HYBRID_MOMENTUM
    if (hybrid_hydro) {

        MultiFab& S_new = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {

	    const Box& bx = mfi.tilebox();

	    hybrid_update(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()), BL_TO_FORTRAN_3D(S_new[mfi]));

	}

    }
#endif

}



void
Castro::initialize_advance(Real time, Real dt, int amr_iteration, int amr_ncycle)
{
    // Pass some information about the state of the simulation to a Fortran module.

    set_amr_info(level, amr_iteration, amr_ncycle, time, dt);

    // The option of whether to do a multilevel initialization is
    // controlled within the radiation class.  This step belongs
    // before the swap.

#ifdef RADIATION
    if (do_radiation)
        radiation->pre_timestep(level);
#endif

    // Swap the new data from the last timestep into the old state data.
    // If we're on level 0, do it for all levels below this one as well.
    // Or, if we're on a later iteration at a finer timestep, swap for all
    // lower time levels as well.

    if (level == 0 || amr_iteration > 1) {

        for (int lev = level; lev <= parent->finestLevel(); lev++) {

	    Real dt_lev = parent->dtLevel(lev);
            for (int k = 0; k < NUM_STATE_TYPE; k++) {

	        // The following is a hack to make sure that
	        // we only ever have new data for the Source_Type;
	        // by doing a swap now, we'll guarantee that
	        // allocOldData() does nothing. We do this because
	        // we never need the old data, so we don't want to
	        // allocate memory for it.

	        if (k == Source_Type) {
		  getLevel(lev).state[k].swapTimeLevels(0.0);
		}

	        getLevel(lev).state[k].allocOldData();
                getLevel(lev).state[k].swapTimeLevels(dt_lev);
            }

#ifdef GRAVITY
	    if (do_grav)
               gravity->swapTimeLevels(lev);
#endif

        }
    }

    // Make a copy of the MultiFabs in the old and new state data in case we may do a retry.

    if (use_retry) {

      // Store the old and new time levels.

      for (int k = 0; k < NUM_STATE_TYPE; k++) {

	prev_state.set(k, new StateData());

	StateData::Initialize(prev_state[k], state[k]);

      }

    }

    // These arrays hold all source terms that update the state.

#ifndef SDC
    for (int n = 0; n < num_src; ++n) {
        old_sources.set(n, new MultiFab(grids, NUM_STATE, NUM_GROW));
        new_sources.set(n, new MultiFab(grids, NUM_STATE, 0));
    }
#endif

    // This array holds the hydrodynamics update.

    hydro_source = new MultiFab(grids,NUM_STATE,0,Fab_allocate);

    // This array holds the sum of all source terms that affect the hydrodynamics.
    // If we are doing the source term predictor, we'll also use this after the
    // hydro update to store the sum of the new-time sources, so that we can
    // compute the time derivative of the source terms.

    sources_for_hydro = new MultiFab(grids,NUM_STATE,NUM_GROW,Fab_allocate);

    for (int j = 0; j < BL_SPACEDIM; j++)
    {
        fluxes[j] = new MultiFab(getEdgeBoxArray(j), NUM_STATE, 0, Fab_allocate);
    }

    for (int j = BL_SPACEDIM; j < 3; j++)
    {
        BoxArray ba = get_new_data(State_Type).boxArray();
	fluxes[j] = new MultiFab(ba, NUM_STATE, 0, Fab_allocate);
    }

#ifdef RADIATION
    MultiFab& Er_new = get_new_data(Rad_Type);
    if (Radiation::rad_hydro_combined) {
        for (int dir = 0; dir < BL_SPACEDIM; dir++) {
	    rad_fluxes[dir] = new MultiFab(getEdgeBoxArray(dir), Radiation::nGroups, 0, Fab_allocate);
	}
    }
#endif

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
	u_gdnv[dir] = new MultiFab(getEdgeBoxArray(dir),1,1,Fab_allocate);
    }

#ifdef DIFFUSION
#ifdef TAU
    tau_diff = new MultiFab(grids,1,NUM_GROW);
#endif
#endif

#ifdef DIFFUSION
    OldTempDiffTerm = new MultiFab(grids, 1, 1);
    OldSpecDiffTerm = new MultiFab(grids,NumSpec,1);
    OldViscousTermforMomentum = new MultiFab(grids,BL_SPACEDIM,1);
    OldViscousTermforEnergy = new MultiFab(grids,1,1);
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

    Real cur_time = state[State_Type].curTime();
    set_special_tagging_flag(cur_time);

    delete hydro_source;
    delete sources_for_hydro;

#ifndef SDC
    old_sources.clear();
    new_sources.clear();
#endif

    for (int n = 0; n < 3; ++n)
        delete fluxes[n];

#ifdef RADIATION
    for (int n = 0; n < BL_SPACEDIM; ++n)
        delete rad_fluxes[n];
#endif

#ifndef LEVELSET
    for (int n = 0; n < BL_SPACEDIM; ++n)
        delete u_gdnv[n];
#endif

    prev_state.clear();

#ifdef DIFFUSION
    delete OldTempDiffTerm;
    delete OldSpecDiffTerm;
    delete OldViscousTermforMomentum;
    delete OldViscousTermforEnergy;
#endif

#ifdef TAU
    delete tau_diff;
#endif

}



Real
Castro::retry_advance(Real time, Real dt, int amr_iteration, int amr_ncycle)
{

    Real dt_new = 1.e200;
    Real dt_subcycle = 1.e200;

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

#ifdef REACTIONS
    MultiFab& R_old = get_old_data(Reactions_Type);
    MultiFab& R_new = get_new_data(Reactions_Type);
#endif

    const Real* dx = geom.CellSize();

#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_subcycle)
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
			  &dt, &dt_subcycle);

    }

    if (retry_neg_dens_factor > 0.0) {

        // Negative density criterion
	// Reset so that the desired maximum fractional change in density
	// is not larger than retry_neg_dens_factor.

        ParallelDescriptor::ReduceRealMin(frac_change);

	if (frac_change < 0.0)
	  dt_subcycle = std::min(dt_subcycle, dt * -(retry_neg_dens_factor / frac_change));

    }

    ParallelDescriptor::ReduceRealMin(dt_subcycle);

    if (dt_subcycle < dt) {

        int sub_ncycle = ceil(dt / dt_subcycle);

	if (verbose && ParallelDescriptor::IOProcessor()) {
	  std::cout << std::endl;
	  std::cout << "  Timestep " << dt << " rejected at level " << level << "." << std::endl;
	  std::cout << "  Performing a retry, with " << sub_ncycle
		    << " subcycled timesteps of maximum length dt = " << dt_subcycle << std::endl;
	  std::cout << std::endl;
	}

	Real subcycle_time = time;
	int sub_iteration = 1;
	Real dt_advance = dt / sub_ncycle;

	// Restore the original values of the state data.

	for (int k = 0; k < NUM_STATE_TYPE; k++) {

	  if (prev_state[k].hasOldData())
	    state[k].copyOld(prev_state[k]);

	  if (prev_state[k].hasNewData())
	    state[k].copyNew(prev_state[k]);

	  // Anticipate the swapTimeLevels to come.

	  if (k == Source_Type)
	    state[k].swapTimeLevels(0.0);

	  state[k].swapTimeLevels(0.0);

	  state[k].setTimeLevel(time, 0.0, 0.0);

	}

	if (track_grid_losses)
	  for (int i = 0; i < n_lost; i++)
	    material_lost_through_boundary_temp[i] = 0.0;

	// Subcycle until we've reached the target time.

	while (subcycle_time < time + dt) {

	    // Shorten the last timestep so that we don't overshoot
	    // the ending time. We want to protect against taking
	    // a very small last timestep due to precision issues,
	    // so subtract a small number from that time.

	    Real eps = 1.0e-10 * dt;

	    if (subcycle_time + dt_advance > time + dt - eps)
	        dt_advance = (time + dt) - subcycle_time;

	    if (verbose && ParallelDescriptor::IOProcessor()) {
	        std::cout << "  Beginning retry subcycle " << sub_iteration << " of " << sub_ncycle
		          << ", starting at time " << subcycle_time
		         << " with dt = " << dt_advance << std::endl << std::endl;
	    }

	    for (int k = 0; k < NUM_STATE_TYPE; k++) {

	        if (k == Source_Type)
		    state[k].swapTimeLevels(0.0);

		state[k].swapTimeLevels(dt_advance);

	    }

#ifdef GRAVITY
	    if (do_grav)
	        gravity->swapTimeLevels(level);
#endif

	    do_advance(subcycle_time,dt_advance,amr_iteration,amr_ncycle,sub_iteration,sub_ncycle);

	    if (verbose && ParallelDescriptor::IOProcessor()) {
	        std::cout << std::endl;
	        std::cout << "  Retry subcycle " << sub_iteration << " of " << sub_ncycle << " completed" << std::endl;
	        std::cout << std::endl;
	    }

	  subcycle_time += dt_advance;
	  sub_iteration += 1;

	}

	// We want to return this subcycled timestep as a suggestion,
	// if it is smaller than what the hydro estimates.

	dt_new = std::min(dt_new, dt_subcycle);

	if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "  Retry subcycling complete" << std::endl << std::endl;

	// Finally, copy the original data back to the old state
	// data so that externally it appears like we took only
	// a single timestep.

	for (int k = 0; k < NUM_STATE_TYPE; k++) {

           if (prev_state[k].hasOldData())
	      state[k].copyOld(prev_state[k]);

	   state[k].setTimeLevel(time + dt, dt, 0.0);

	}

    }

    return dt_new;

}
