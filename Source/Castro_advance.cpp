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

    // Pass some information about the state of the simulation to a Fortran module.

    set_amr_info(level, amr_iteration, amr_ncycle, time, dt);

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

    PArray<StateData> prev_state(NUM_STATE_TYPE,PArrayManage);

    int sub_iteration = 0;
    int sub_ncycle = 0;

    if (use_retry) {

      // Store the old and new time levels.

      for (int k = 0; k < NUM_STATE_TYPE; k++) {

	prev_state.set(k, new StateData());

	StateData::Initialize(prev_state[k], state[k]);

      }

    }

    // Reset the grid loss tracking.

    if (track_grid_losses)
      for (int i = 0; i < n_lost; i++)
	material_lost_through_boundary_temp[i] = 0.0;

    // Do the advance.

    dt_new = advance_hydro(time,dt,amr_iteration,amr_ncycle,sub_iteration,sub_ncycle);

    // Check to see if this advance violated certain stability criteria.
    // If so, get a new timestep and do subcycled advances until we reach
    // t = time + dt.

    if (use_retry)
    {

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
	sub_iteration = 1;
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

	  advance_hydro(subcycle_time,dt_advance,amr_iteration,amr_ncycle,sub_iteration,sub_ncycle);

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

	if (verbose && ParallelDescriptor::IOProcessor()) {
	  std::cout << "  Retry subcycling complete" << std::endl << std::endl;
	}

	// Finally, copy the original data back to the old state
	// data so that externally it appears like we took only
	// a single timestep.

	for (int k = 0; k < NUM_STATE_TYPE; k++) {

	  if (prev_state[k].hasOldData())
	    state[k].copyOld(prev_state[k]);

	  state[k].setTimeLevel(time + dt, dt, 0.0);

	}

      }

    }

    // Add the material lost in this timestep to the cumulative losses.

    if (track_grid_losses) {

      ParallelDescriptor::ReduceRealSum(material_lost_through_boundary_temp, n_lost);

      for (int i = 0; i < n_lost; i++)
	material_lost_through_boundary_cumulative[i] += material_lost_through_boundary_temp[i];

    }

    Real cur_time = state[State_Type].curTime();
    set_special_tagging_flag(cur_time);

#ifdef AUX_UPDATE
    advance_aux(time,dt);
#endif

#ifdef LEVELSET
    advance_levelset(time,dt);
#endif

#if (BL_SPACEDIM > 1)
    // We do this again here because the solution will have changed
    if ( (level == 0) && (spherical_star == 1) ) {
       int is_new = 1;
       make_radial_data(is_new);
    }
#endif

#ifdef RADIATION
    MultiFab& S_new = get_new_data(State_Type);
    final_radiation_call(S_new,amr_iteration,amr_ncycle);
#endif

#ifdef PARTICLES
    UpdateParticles(amr_iteration, time, dt);
#endif

    return dt_new;
}

Real
Castro::advance_hydro (Real time,
                       Real dt,
                       int  amr_iteration,
                       int  amr_ncycle,
		       int  sub_iteration,
		       int  sub_ncycle)
{
    BL_PROFILE("Castro::advance_hydro()");

#ifdef RADIATION
    if (do_radiation) {
        // The option of whether to do a multilevel initialization is
        // controlled within the radiation class.  This step belongs
        // before the swap.

        radiation->pre_timestep(level);
    }
#endif

    // These arrays hold all source terms that update the state.

    for (int n = 0; n < num_src; ++n) {
        old_sources.set(n, new MultiFab(grids, NUM_STATE, NUM_GROW));
        new_sources.set(n, new MultiFab(grids, NUM_STATE, 0));

	old_sources[n].setVal(0.0);
	new_sources[n].setVal(0.0);
    }

    // This array holds the hydrodynamics update.

    hydro_source = new MultiFab(grids,NUM_STATE,0,Fab_allocate);

    hydro_source->setVal(0.0);

    // This array holds the sum of all source terms that affect the hydrodynamics.
    // If we are doing the source term predictor, we'll also use this after the
    // hydro update to store the sum of the new-time sources, so that we can
    // compute the time derivative of the source terms.

    sources_for_hydro = new MultiFab(grids,NUM_STATE,NUM_GROW,Fab_allocate);
    sources_for_hydro->setVal(0.0,NUM_GROW);

    // Reset the change from density resets

    frac_change = 1.e0;

    u_gdnv = new MultiFab[BL_SPACEDIM];
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
	u_gdnv[dir].define(getEdgeBoxArray(dir),1,1,Fab_allocate);
	u_gdnv[dir].setVal(1.e40,1);
    }

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

    const Real prev_time = state[State_Type].prevTime();
    const Real  cur_time = state[State_Type].curTime();

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    // Check for NaN's.

    check_for_nan(S_old);

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
    S_old.setBndry(0.0);
    S_new.setBndry(0.0);
#endif

    // We want to define this on every grid, because we may need the fluxes
    // when computing the source terms later.

    for (int j = 0; j < BL_SPACEDIM; j++)
       {
         fluxes[j] = new MultiFab(getEdgeBoxArray(j), NUM_STATE, 0, Fab_allocate);
         fluxes[j]->setVal(0.0);
       }

    for (int j = BL_SPACEDIM; j < 3; j++)
      {
	BoxArray ba = S_new.boxArray();
	fluxes[j] = new MultiFab(ba, NUM_STATE, 0, Fab_allocate);
	fluxes[j]->setVal(0.0);
      }

#ifdef RADIATION
    MultiFab& Er_new = get_new_data(Rad_Type);
    if (Radiation::rad_hydro_combined) {
      for (int dir = 0; dir < BL_SPACEDIM; dir++) {
	rad_fluxes[dir] = new MultiFab(getEdgeBoxArray(dir), Radiation::nGroups, 0, Fab_allocate);
      }
    }
#endif

    // Do preliminary steps for constructing old-time sources.

    set_up_for_old_sources(time);

    // For the hydrodynamics update we need to have NUM_GROW ghost zones available,
    // but the state data does not carry ghost zones. So we use a FillPatch
    // using the state data to give us Sborder, which does have ghost zones.

    Sborder = new MultiFab(grids, NUM_STATE, NUM_GROW, Fab_allocate);
    expand_state(*Sborder, prev_time, NUM_GROW);

    // Since we are Strang splitting the reactions, do them now.

    // Reactions are expensive and we would usually rather do a communication
    // step than burn on the ghost zones. So what we will do here is create a mask
    // that indicates that we want to turn on the valid interior zones but NOT
    // on the ghost zones that are interior to the level. However, we DO want to
    // burn on the ghost zones that are on the coarse-fine interfaces, since that
    // is going to be more accurate than interpolating from coarse zones. So we will
    // not mask out those zones, and the subsequent FillBoundary call will not
    // interfere with it.

#ifdef REACTIONS

    MultiFab& reactions_old = get_old_data(Reactions_Type);

    const int react_ngrow_first_half = NUM_GROW;
    const iMultiFab& interior_mask_first_half = build_interior_boundary_mask(react_ngrow_first_half);

    react_state(*Sborder,reactions_old,interior_mask_first_half,time,0.5*dt,react_ngrow_first_half);

    BoxLib::fill_boundary(*Sborder, geom);

#endif

    // Initialize the new-time data.

    MultiFab::Copy(S_new, *Sborder, 0, 0, NUM_STATE, S_new.nGrow());

    // Construct the old-time sources.

    construct_old_sources(amr_iteration, amr_ncycle, sub_iteration, sub_ncycle, prev_time, dt);

    // Apply the old-time sources.

    for (int n = 0; n < num_src; ++n)
        apply_source_to_state(S_new, old_sources[n], dt);

    // Optionally we can predict the source terms to t + dt/2,
    // which is the time-level n+1/2 value, To do this we use a
    // lagged predictor estimate: dS/dt_n = (S_n - S_{n-1}) / dt, so 
    // S_{n+1/2} = S_n + (dt / 2) * dS/dt_n.

    if (source_term_predictor == 1) {

      MultiFab& dSdt_new = get_new_data(Source_Type);

      AmrLevel::FillPatch(*this,dSdt_new,NUM_GROW,cur_time,Source_Type,0,NUM_STATE);       

      dSdt_new.mult(dt / 2.0, NUM_GROW);

      MultiFab::Add(*sources_for_hydro,dSdt_new,0,0,NUM_STATE,NUM_GROW);

    }

    // Do the hydro update.

    if (do_hydro)
        hydro_update(time, dt);

    // Update the point mass.

#ifdef POINTMASS
    pointmass_update(time, dt);
#endif

    // Check for NaN's.

    check_for_nan(S_new);

    // Preliminary steps before constructing the new sources.

    set_up_for_new_sources(cur_time);

    // Construct the new-time source terms.

    construct_new_sources(amr_iteration, amr_ncycle,
			  sub_iteration, sub_ncycle,
			  cur_time, dt);

    // Apply the new-time sources to the state.

    for (int n = 0; n < num_src; ++n)
        apply_source_to_state(S_new, new_sources[n], dt);

    // Sync up the temperature now that all sources have been applied.

    computeTemp(S_new);

    // Do the final update for dSdt.

    if (source_term_predictor == 1) {

      MultiFab& dSdt_new = get_new_data(Source_Type);

      // Calculate the time derivative of the source terms.

      MultiFab::Add(dSdt_new,*sources_for_hydro,0,0,NUM_STATE,0);

      dSdt_new.mult(1.0/dt);

    }

    // Do the second half of the reactions.

#ifdef REACTIONS

    MultiFab& reactions_new = get_new_data(Reactions_Type);

    const int react_ngrow_second_half = 0;
    const iMultiFab& interior_mask_second_half = build_interior_boundary_mask(react_ngrow_second_half);

    react_state(S_new,reactions_new,interior_mask_second_half,cur_time-0.5*dt,0.5*dt,react_ngrow_second_half);

    BoxLib::fill_boundary(S_new, geom);

#endif

    // Sync up the hybrid and linear momenta.

#ifdef HYBRID_MOMENTUM
    if (hybrid_hydro) {

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {

	    const Box& bx = mfi.tilebox();

	    hybrid_update(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()), BL_TO_FORTRAN_3D(S_new[mfi]));

	}

    }
#endif

    delete hydro_source;
    delete sources_for_hydro;

    old_sources.clear();
    new_sources.clear();

    for (int n = 0; n < 3; ++n)
        delete fluxes[n];

#ifdef RADIATION
    for (int n = 0; n < BL_SPACEDIM; ++n)
        delete rad_fluxes[n];
#endif

#ifndef LEVELSET
    delete [] u_gdnv;
#endif

#ifdef DIFFUSION
    OldTempDiffTerm = new MultiFab(grids, 1, 1);
    OldSpecDiffTerm = new MultiFab(grids,NumSpec,1);
    OldViscousTermforMomentum = new MultiFab(grids,BL_SPACEDIM,1);
    OldViscousTermforEnergy = new MultiFab(grids,1,1);
#endif

#ifdef TAU
    delete tau_diff;
#endif

    return dt;
}



void
Castro::hydro_update(Real time, Real dt)
{

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... Entering hydro advance" << std::endl << std::endl;

    int finest_level = parent->finestLevel();

    //
    // Get pointers to Flux registers, or set pointer to zero if not there.
    //
    FluxRegister *fine    = 0;
    FluxRegister *current = 0;

    if (do_reflux && level < finest_level)
      fine = &getFluxReg(level+1);
    if (do_reflux && level > 0)
      current = &getFluxReg(level);

#ifdef RADIATION
    FluxRegister *rad_fine    = 0;
    FluxRegister *rad_current = 0;
    if (Radiation::rad_hydro_combined && do_reflux && level < finest_level)
      rad_fine = &getRADFluxReg(level+1);
    if (Radiation::rad_hydro_combined && do_reflux && level > 0)
      rad_current = &getRADFluxReg(level);
#endif

    const Real *dx = geom.CellSize();
    Real courno    = -1.0e+200;

    MultiFab& S_new = get_new_data(State_Type);

#ifdef RADIATION
    if (Radiation::rad_hydro_combined) {

	FillPatchIterator fpi_rad(*this, Er_new, NUM_GROW, time, Rad_Type, 0, Radiation::nGroups);
	MultiFab& Erborder = fpi_rad.get_mf();

	MultiFab lamborder(grids, Radiation::nGroups, NUM_GROW);
	MultiFab kappa_s;
	if (radiation->do_inelastic_scattering) {
	    kappa_s.define(grids, 1, NUM_GROW, Fab_allocate);
	    kappa_s.setVal(0.0, NUM_GROW);
	}
	if (radiation->pure_hydro) {
	    lamborder.setVal(0.0, NUM_GROW);
	}
	else {
	    radiation->compute_limiter(level, grids, *Sborder, Erborder, lamborder, kappa_s);
	}

	int nstep_fsp = -1;

	BL_PROFILE_VAR("Castro::advance_hydro_ca_umdrv_rad()", CA_UMDRV_RAD);

#ifdef _OPENMP
#pragma omp parallel
#endif
	{
	    FArrayBox flux[BL_SPACEDIM], ugdn[BL_SPACEDIM], rad_flux[BL_SPACEDIM];

	    int priv_nstep_fsp = -1;
	    Real cflLoc = -1.0e+200;
	    int is_finest_level = (level == finest_level) ? 1 : 0;
	    const int*  domain_lo = geom.Domain().loVect();
	    const int*  domain_hi = geom.Domain().hiVect();

	    for (MFIter mfi(S_new,hydro_tile_size); mfi.isValid(); ++mfi)
	    {
		const Box &bx    = mfi.tilebox();

		const int* lo = bx.loVect();
		const int* hi = bx.hiVect();

		FArrayBox &statein  = (*Sborder)[mfi];
		FArrayBox &stateout = S_new[mfi];

		FArrayBox &Er = Erborder[mfi];
		FArrayBox &lam = lamborder[mfi];
		FArrayBox &Erout = Er_new[mfi];

		FArrayBox& vol      = volume[mfi];

		// Allocate fabs for fluxes and Godunov velocities.
		for (int i = 0; i < BL_SPACEDIM ; i++)  {
		    const Box& bxtmp = BoxLib::surroundingNodes(bx,i);
		    flux[i].resize(BoxLib::surroundingNodes(bx,i),NUM_STATE);
		    rad_flux[i].resize(BoxLib::surroundingNodes(bx,i),Radiation::nGroups);
		    ugdn[i].resize(BoxLib::grow(bxtmp,1),1);
		}

		ca_umdrv_rad
		    (&is_finest_level,&time,
		     bx.loVect(), bx.hiVect(),
		     domain_lo, domain_hi,
		     BL_TO_FORTRAN(statein), BL_TO_FORTRAN(stateout),
		     BL_TO_FORTRAN(Er), BL_TO_FORTRAN(lam),
		     BL_TO_FORTRAN(Erout),
		     D_DECL(BL_TO_FORTRAN(ugdn[0]),
			    BL_TO_FORTRAN(ugdn[1]),
			    BL_TO_FORTRAN(ugdn[2])),
		     BL_TO_FORTRAN((*sources_for_hydro)[mfi]),
		     dx, &dt,
		     D_DECL(BL_TO_FORTRAN(flux[0]),
			    BL_TO_FORTRAN(flux[1]),
			    BL_TO_FORTRAN(flux[2])),
		     D_DECL(BL_TO_FORTRAN(rad_flux[0]),
			    BL_TO_FORTRAN(rad_flux[1]),
			    BL_TO_FORTRAN(rad_flux[2])),
		     D_DECL(BL_TO_FORTRAN(area[0][mfi]),
			    BL_TO_FORTRAN(area[1][mfi]),
			    BL_TO_FORTRAN(area[2][mfi])),
#if (BL_SPACEDIM < 3)
		     BL_TO_FORTRAN(dLogArea[0][mfi]),
#endif
		     BL_TO_FORTRAN(volume[mfi]),
		     &cflLoc, verbose, &priv_nstep_fsp);

		for (int i = 0; i < BL_SPACEDIM ; i++) {
		    u_gdnv[i][mfi].copy(ugdn[i],mfi.nodaltilebox(i));
		}

		if (do_reflux) {
		    for (int i = 0; i < BL_SPACEDIM ; i++) {
		        (*fluxes    [i])[mfi].copy(    flux[i],mfi.nodaltilebox(i));
			(*rad_fluxes[i])[mfi].copy(rad_flux[i],mfi.nodaltilebox(i));
		    }
		}

		if (radiation->do_inelastic_scattering) {
		    ca_inelastic_sct(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
				     BL_TO_FORTRAN_3D(stateout),
				     BL_TO_FORTRAN_3D(Erout),
				     BL_TO_FORTRAN_3D(kappa_s[mfi]),
				     dt);
		}

	    }

#ifdef _OPENMP
#pragma omp critical (radhydro_courno)
#endif
	    {
		courno = std::max(courno,cflLoc);
		nstep_fsp = std::max(nstep_fsp, priv_nstep_fsp);
	    }
	}  // end of omp parallel region

	BL_PROFILE_VAR_STOP(CA_UMDRV_RAD);

	if (radiation->verbose>=1) {
#ifdef BL_LAZY
	    Lazy::QueueReduction( [=] () mutable {
#endif
	    ParallelDescriptor::ReduceIntMax(nstep_fsp, ParallelDescriptor::IOProcessorNumber());
	    if (ParallelDescriptor::IOProcessor() && nstep_fsp > 0) {
		std::cout << "Radiation f-space advection on level " << level
			  << " takes as many as " << nstep_fsp;
		if (nstep_fsp == 1) {
		    std::cout<< " substep.\n";
		}
		else {
		    std::cout<< " substeps.\n";
		}
	    }
#ifdef BL_LAZY
	    });
#endif
	}

    }
    else {
      BoxLib::Abort("Castro::advance -- we don't implement a mode where we have radiation, but it is not coupled to hydro");
    }
#else

	// pure hydro (no radiation)

	Real E_added_flux    = 0.;
	Real mass_added      = 0.;
	Real eint_added      = 0.;
	Real eden_added      = 0.;
	Real dens_change     = 1.e200;
	Real mass_added_flux = 0.;
	Real xmom_added_flux = 0.;
	Real ymom_added_flux = 0.;
	Real zmom_added_flux = 0.;
	Real mass_lost       = 0.;
	Real xmom_lost       = 0.;
	Real ymom_lost       = 0.;
	Real zmom_lost       = 0.;
	Real eden_lost       = 0.;
	Real xang_lost       = 0.;
	Real yang_lost       = 0.;
	Real zang_lost       = 0.;

	BL_PROFILE_VAR("Castro::advance_hydro_ca_umdrv()", CA_UMDRV);

#ifdef _OPENMP
#pragma omp parallel reduction(+:E_added_flux) \
                     reduction(+:mass_added,eint_added,eden_added,mass_added_flux) \
		     reduction(+:xmom_added_flux,ymom_added_flux,zmom_added_flux) \
		     reduction(+:mass_lost,xmom_lost,ymom_lost,zmom_lost) \
		     reduction(+:eden_lost,xang_lost,yang_lost,zang_lost) \
		     reduction(min:dens_change)
#endif
	{
	    FArrayBox flux[BL_SPACEDIM], ugdn[BL_SPACEDIM];

	    Real cflLoc = -1.0e+200;
	    int is_finest_level = (level == finest_level) ? 1 : 0;
	    const int* domain_lo = geom.Domain().loVect();
	    const int* domain_hi = geom.Domain().hiVect();

	    for (MFIter mfi(S_new,hydro_tile_size); mfi.isValid(); ++mfi)
	    {
		const Box& bx = mfi.tilebox();

		const int* lo = bx.loVect();
		const int* hi = bx.hiVect();

		FArrayBox &statein  = (*Sborder)[mfi];
		FArrayBox &stateout = S_new[mfi];

		FArrayBox &source = (*hydro_source)[mfi];

		FArrayBox &vol = volume[mfi];

		// Allocate fabs for fluxes and Godunov velocities.
		for (int i = 0; i < BL_SPACEDIM; i++) {
		    const Box& bxtmp = BoxLib::surroundingNodes(bx,i);
		    flux[i].resize(bxtmp,NUM_STATE);
		    ugdn[i].resize(BoxLib::grow(bxtmp,1),1);
		}

		ca_umdrv
		    (&is_finest_level,&time,
		     lo, hi, domain_lo, domain_hi,
		     BL_TO_FORTRAN(statein),
		     BL_TO_FORTRAN(stateout),
		     BL_TO_FORTRAN(source),
		     D_DECL(BL_TO_FORTRAN(ugdn[0]),
			    BL_TO_FORTRAN(ugdn[1]),
			    BL_TO_FORTRAN(ugdn[2])),
		     BL_TO_FORTRAN((*sources_for_hydro)[mfi]),
		     dx, &dt,
		     D_DECL(BL_TO_FORTRAN(flux[0]),
			    BL_TO_FORTRAN(flux[1]),
			    BL_TO_FORTRAN(flux[2])),
		     D_DECL(BL_TO_FORTRAN(area[0][mfi]),
			    BL_TO_FORTRAN(area[1][mfi]),
			    BL_TO_FORTRAN(area[2][mfi])),
#if (BL_SPACEDIM < 3)
		     BL_TO_FORTRAN(dLogArea[0][mfi]),
#endif
		     BL_TO_FORTRAN(vol),
		     &cflLoc, verbose,
		     mass_added_flux,
		     xmom_added_flux,
		     ymom_added_flux,
		     zmom_added_flux,
		     E_added_flux,
		     mass_lost, xmom_lost, ymom_lost, zmom_lost,
		     eden_lost, xang_lost, yang_lost, zang_lost);

		// Apply the hydro source term.

		stateout.saxpy(dt,source,bx,bx,0,0,NUM_STATE);

		// Enforce the density >= small_dens.

		enforce_minimum_density(statein.dataPtr(),ARLIM_3D(statein.loVect()),ARLIM_3D(statein.hiVect()),
					stateout.dataPtr(),ARLIM_3D(stateout.loVect()),ARLIM_3D(stateout.hiVect()),
					vol.dataPtr(),ARLIM_3D(vol.loVect()),ARLIM_3D(vol.hiVect()),
					ARLIM_3D(bx.loVect()),ARLIM_3D(bx.hiVect()),
					&mass_added,&eint_added,&eden_added,&dens_change,&verbose);

		// Renormalize species mass fractions

		ca_normalize_species(stateout.dataPtr(),ARLIM_3D(stateout.loVect()), ARLIM_3D(stateout.hiVect()),
				     ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()));

		// Copy the normal velocities from the Riemann solver

		for (int i = 0; i < BL_SPACEDIM ; i++) {
		    u_gdnv[i][mfi].copy(ugdn[i],mfi.nodaltilebox(i));
		}

		// Since we may need the fluxes later on, we'll copy them
		// to the fluxes MultiFAB even if we aren't on a fine grid.

		for (int i = 0; i < BL_SPACEDIM ; i++)
		    (*fluxes[i])[mfi].copy(flux[i],mfi.nodaltilebox(i));

	    }

#ifdef _OPENMP
#pragma omp critical (hydro_courno)
#endif
	    {
		courno = std::max(courno,cflLoc);
	    }
	}  // end of omp parallel region

	BL_PROFILE_VAR_STOP(CA_UMDRV);

	frac_change = dens_change;

	// Flush Fortran output

	if (verbose)
	  flush_output();

	if (track_grid_losses)
	{

	  material_lost_through_boundary_temp[0] += mass_lost;
	  material_lost_through_boundary_temp[1] += xmom_lost;
	  material_lost_through_boundary_temp[2] += ymom_lost;
	  material_lost_through_boundary_temp[3] += zmom_lost;
	  material_lost_through_boundary_temp[4] += eden_lost;
	  material_lost_through_boundary_temp[5] += xang_lost;
	  material_lost_through_boundary_temp[6] += yang_lost;
	  material_lost_through_boundary_temp[7] += zang_lost;

	}

	if (print_energy_diagnostics)
	{
	   Real foo[8] = {mass_added, eint_added, eden_added,
			  E_added_flux,
			  xmom_added_flux, ymom_added_flux, zmom_added_flux,
			  mass_added_flux};
#ifdef BL_LAZY
	   Lazy::QueueReduction( [=] () mutable {
#endif
	   ParallelDescriptor::ReduceRealSum(foo, 20, ParallelDescriptor::IOProcessorNumber());
	   if (ParallelDescriptor::IOProcessor())
	   {
	       mass_added = foo[0];
	       eint_added = foo[1];
	       eden_added = foo[2];
	       E_added_flux = foo[3];
	       xmom_added_flux = foo[4];
	       ymom_added_flux = foo[5];
	       zmom_added_flux = foo[6];
	       mass_added_flux    = foo[7];
	       if (std::abs(mass_added) != 0.0)
	       {
		  std::cout << "   Mass added from negative density correction : " <<
				mass_added << std::endl;
		  std::cout << "(rho e) added from negative density correction : " <<
				eint_added << std::endl;
		  std::cout << "(rho E) added from negative density correction : " <<
				eden_added << std::endl;
	       }

	       std::cout << "mass added from fluxes                      : " <<
			     mass_added_flux << std::endl;
	       std::cout << "xmom added from fluxes                      : " <<
			     xmom_added_flux << std::endl;
	       std::cout << "ymom added from fluxes                      : " <<
			     ymom_added_flux << std::endl;
	       std::cout << "zmom added from fluxes                      : " <<
			     zmom_added_flux << std::endl;
	       std::cout << "(rho E) added from fluxes                   : " <<
			     E_added_flux << std::endl;
	   }
#ifdef BL_LAZY
	   });
#endif
	}

#endif    // RADIATION

    if (use_retry)
      ParallelDescriptor::ReduceRealMin(frac_change);

    // Save the fluxes.

    if (do_reflux) {
	if (current) {
	    for (int i = 0; i < BL_SPACEDIM ; i++)
	        current->FineAdd(*fluxes[i],i,0,0,NUM_STATE,1.);
	}
	if (fine) {
	    for (int i = 0; i < BL_SPACEDIM ; i++)
                fine->CrseInit(*fluxes[i],i,0,0,NUM_STATE,-1.,FluxRegister::ADD);
	}
#ifdef RADIATION
	if (rad_current) {
	    for (int i = 0; i < BL_SPACEDIM ; i++)
                rad_current->FineAdd(*rad_fluxes[i],i,0,0,Radiation::nGroups,1.);
	}
	if (rad_fine) {
	    for (int i = 0; i < BL_SPACEDIM ; i++)
                rad_fine->CrseInit(*rad_fluxes[i],i,0,0,Radiation::nGroups,-1.,FluxRegister::ADD);
        }
#endif
    }

    if (courno > 1.0) {
	std::cout << "WARNING -- EFFECTIVE CFL AT THIS LEVEL " << level << " IS " << courno << '\n';
	if (hard_cfl_limit == 1)
	  BoxLib::Abort("CFL is too high at this level -- go back to a checkpoint and restart with lower cfl number");
    }

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << std::endl << "... Leaving hydro advance" << std::endl << std::endl;

}



void
Castro::set_up_for_old_sources(Real time)
{

    MultiFab& S_old = get_old_data(State_Type);

#ifdef GRAVITY
    if (moving_center == 1)
       define_new_center(S_old,time);
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
    tau_diff = new MultiFab(grids,1,NUM_GROW);
    tau_diff->setVal(0.);
    define_tau(grav_old,time);
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
Castro::set_up_for_new_sources(Real time)
{

    // Now we'll start updating the dSdt MultiFab. First,
    // get rid of the dt/2 * dS/dt that we added from the last
    // timestep, then subtract the old sources data to get the
    // first half of the update for the next calculation of dS/dt.

    if (source_term_predictor == 1) {
      MultiFab& dSdt_new = get_new_data(Source_Type);
      MultiFab::Subtract(*sources_for_hydro,dSdt_new,0,0,NUM_STATE,NUM_GROW);
      dSdt_new.setVal(0.0, NUM_GROW);
      MultiFab::Subtract(dSdt_new,*sources_for_hydro,0,0,NUM_STATE,0);
    }

    sources_for_hydro->setVal(0.0,NUM_GROW);

    MultiFab& S_new = get_new_data(State_Type);

#ifdef GRAVITY
    // Must define new value of "center" before we call new gravity solve or external source routine
    if (moving_center == 1)
       define_new_center(S_new,time);
#endif

#if (BL_SPACEDIM > 1)
      // We need to make the new radial data now so that we can use it when we
      //   FillPatch in creating the new source
      if ( (level == 0) && (spherical_star == 1) ) {
	  int is_new = 1;
	  make_radial_data(is_new);
      }
#endif

#ifdef DIFFUSION
    NewTempDiffTerm = OldTempDiffTerm;
    NewSpecDiffTerm = OldSpecDiffTerm;
    NewViscousTermforMomentum = OldViscousTermforMomentum;
    NewViscousTermforEnergy   = OldViscousTermforEnergy;
#endif

    // Compute the current temperature for use in the source term evaluation.

    computeTemp(S_new);

}
