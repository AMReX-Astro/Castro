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
                 int  iteration,
                 int  ncycle)
{
    BL_PROFILE("Castro::advance()");

    Real dt_new = dt;

    // Pass some information about the state of the simulation to a Fortran module.

    set_amr_info(level, iteration, ncycle, time, dt);

    // Swap the new data from the last timestep into the old state data.
    // If we're on level 0, do it for all levels below this one as well.
    // Or, if we're on a later iteration at a finer timestep, swap for all
    // lower time levels as well.

    if (level == 0 || iteration > 1) {

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

    PArray<MultiFab> old_data(NUM_STATE_TYPE,PArrayManage);
    PArray<MultiFab> new_data(NUM_STATE_TYPE,PArrayManage);

    Array<int> got_old_data(NUM_STATE_TYPE);
    Array<int> got_new_data(NUM_STATE_TYPE);

    Array<Real> t_old(NUM_STATE_TYPE);
    Array<Real> t_new(NUM_STATE_TYPE);

    int subcycle_iter = 0;

    if (use_retry) {

      for (int k = 0; k < NUM_STATE_TYPE; k++) {

	// Store the old and new time levels.

	t_new[k] = state[k].curTime();
	t_old[k] = state[k].prevTime();

	got_old_data[k] = state[k].hasOldData();

	if (got_old_data[k]) {

	  MultiFab& state_old = get_old_data(k);

	  int nc = state_old.nComp();
	  int ng = state_old.nGrow();

	  old_data.set(k, new MultiFab(grids, nc, ng, Fab_allocate));

	  MultiFab::Copy(old_data[k], state_old, 0, 0, nc, ng);

	}

	got_new_data[k] = state[k].hasNewData();

	if (got_new_data[k]) {

	  MultiFab& state_new = get_new_data(k);

	  int nc = state_new.nComp();
	  int ng = state_new.nGrow();

	  new_data.set(k, new MultiFab(grids, nc, ng, Fab_allocate));

	  MultiFab::Copy(new_data[k], state_new, 0, 0, nc, ng);

	}

      }

    }

    // Do the advance.

    if (do_hydro) 
    {
        dt_new = advance_hydro(time,dt,iteration,ncycle,subcycle_iter);
    } 
    else 
    {
#ifdef SGS
        BoxLib::Abort("Castro::advance -- doesn't make sense to have SGS defined but not do_hydro");
        return 0.;
#else
        dt_new = advance_no_hydro(time,dt,iteration,ncycle,subcycle_iter);
#endif
    } 

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

      ParallelDescriptor::ReduceRealMin(dt_subcycle);

      Real negative_density_reset_factor = 1.0e-2;

      // Negative density criterion
      // Reset so that the desired maximum fractional change in density
      // is not larger than negative_density_reset_factor.

      if (frac_change < 0.0)
	dt_subcycle = std::min(dt_subcycle, dt * -(negative_density_reset_factor / frac_change));



      if (dt_subcycle < dt) {

	int n_subcycle_iters = ceil(dt / dt_subcycle);

	if (verbose && ParallelDescriptor::IOProcessor()) {
	  std::cout << std::endl;
	  std::cout << "  Timestep " << dt << " rejected at level " << level << "." << std::endl;
	  std::cout << "  Performing a retry, with " << n_subcycle_iters
		    << " subcycled timesteps of maximum length dt = " << dt_subcycle << std::endl;
	  std::cout << std::endl;
	}

	Real subcycle_time = time;
	subcycle_iter = 1;
	Real dt_advance = dt_subcycle;

	// Restore the original values of the state data.

	for (int k = 0; k < NUM_STATE_TYPE; k++) {

	  if (got_old_data[k]) {

	    MultiFab& state_old = get_old_data(k);

	    int nc = state_old.nComp();
	    int ng = state_old.nGrow();

	    MultiFab::Copy(state_old, old_data[k], 0, 0, nc, ng);

	  }

	  if (got_new_data[k]) {

	    MultiFab& state_new = get_new_data(k);

	    int nc = state_new.nComp();
	    int ng = state_new.nGrow();

	    MultiFab::Copy(state_new, new_data[k], 0, 0, nc, ng);

	  }

	  // Anticipate the swapTimeLevels to come.

	  if (k == Source_Type)
	    state[k].swapTimeLevels(0.0);

	  state[k].swapTimeLevels(0.0);

	  state[k].setTimeLevel(time, 0.0, 0.0);

	}

	// Subcycle until we've reached the target time.

	while (subcycle_time < time + dt) {

	  if (subcycle_time + dt_advance > time + dt)
	    dt_advance = (time + dt) - subcycle_time;

	  if (verbose && ParallelDescriptor::IOProcessor()) {
	    std::cout << "  Beginning retry subcycle " << subcycle_iter << " starting at time " << subcycle_time
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

	  Real dt_temp;

	  if (do_hydro)
	  {
	    dt_temp = advance_hydro(subcycle_time,dt_advance,iteration,ncycle,subcycle_iter);
	  }
	  else
	  {
	    dt_temp = advance_no_hydro(subcycle_time,dt_advance,iteration,ncycle,subcycle_iter);
	  }



	  if (verbose && ParallelDescriptor::IOProcessor()) {
	    std::cout << std::endl;
	    std::cout << "  Retry subcycle " << subcycle_iter << " completed" << std::endl;
	    std::cout << std::endl;
	  }

	  subcycle_time += dt_advance;
	  subcycle_iter += 1;

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

	  state[k].setTimeLevel(time + dt, dt, 0.0);

	  if (got_old_data[k] && state[k].hasOldData()) {

	    MultiFab& state_old = get_old_data(k);

	    int nc = state_old.nComp();
	    int ng = state_old.nGrow();

	    MultiFab::Copy(state_old, old_data[k], 0, 0, nc, ng);

	  }

	}

      }

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
    final_radiation_call(S_new,iteration,ncycle);
#endif

#ifdef PARTICLES
    if (TracerPC)
    {
	int ng = iteration;
	Real t = time + 0.5*dt;

	MultiFab Ucc(grids,BL_SPACEDIM,ng); // cell centered velocity

	{
	    FillPatchIterator fpi(*this, Ucc, ng, t, State_Type, 0, BL_SPACEDIM+1);
	    MultiFab& S = fpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel
#endif
	    for (MFIter mfi(Ucc,true); mfi.isValid(); ++mfi)
	    {
		const Box& bx = mfi.growntilebox();
		S[mfi].invert(1.0, bx, 0, 1);
		for (int dir=0; dir < BL_SPACEDIM; ++dir) {
		    Ucc[mfi].copy(S[mfi], bx, dir+1, bx, dir, 1);
		    Ucc[mfi].mult(S[mfi], bx, 0, dir);
		}
	    }
	}

	TracerPC->AdvectWithUcc(Ucc, level, dt);
    }
#endif

    return dt_new;
}

Real
Castro::advance_hydro (Real time,
                       Real dt,
                       int  iteration,
                       int  ncycle,
		       int  subcycle_iter)
{
    BL_PROFILE("Castro::advance_hydro()");

    if (!do_hydro)  
       BoxLib::Abort("In advance_hydro but do_hydro not true");

#ifdef RADIATION
    if (do_radiation) {
        // The option of whether to do a multilevel initialization is
        // controlled within the radiation class.  This step belongs
        // before the swap.

        radiation->pre_timestep(level);
    } 
#endif

    // Reset the change from density resets

    frac_change = 1.e0;

    u_gdnv = new MultiFab[BL_SPACEDIM];
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
	u_gdnv[dir].define(getEdgeBoxArray(dir),1,1,Fab_allocate);
	u_gdnv[dir].setVal(1.e40,1);
    }
    
    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level && subcycle_iter < 2) {
        //
        // Set reflux registers to zero.
        //
        getFluxReg(level+1).setVal(0.0);
#ifdef SGS
        getSGSFluxReg(level+1).setVal(0.0);
#endif
#ifdef RADIATION
	if (Radiation::rad_hydro_combined) {
	  getRADFluxReg(level+1).setVal(0.0);
	}
#endif
    }

#ifdef SGS
    // Make sure this is zero in case we turn off source terms
    MultiFab& SGS_new = get_new_data(SGS_Type);
    SGS_new.setVal(0.);
#endif

    const Real prev_time = state[State_Type].prevTime();
    const Real  cur_time = state[State_Type].curTime();

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    if (S_old.contains_nan(Density,S_old.nComp(),0,true))
    {
        for (int i = 0; i < S_old.nComp(); i++)
        {
            if (S_old.contains_nan(Density+i,1,0,true))
            {
	        std::cout << "Testing " << desc_lst[State_Type].name(i) << " for NaNs" << std::endl;
                BoxLib::Abort("S_old has NaNs in this component::advance_hydro()");
            }
        }
    }

#ifdef GRAVITY
    if (moving_center == 1)
       define_new_center(S_old,time);
#endif

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

#if (BL_SPACEDIM > 1)
    if ( (level == 0) && (spherical_star == 1) ) {
       swap_outflow_data();
       int is_new = 0;
       make_radial_data(is_new);
    }
#endif

#ifdef GRAVITY
    // Old and new gravitational potential.
       
    MultiFab& phi_old = get_old_data(PhiGrav_Type);
    MultiFab& phi_new = get_new_data(PhiGrav_Type);

    MultiFab& grav_old = get_old_data(Gravity_Type);
    MultiFab& grav_new = get_new_data(Gravity_Type);
    
    if (do_grav) {

       // Swap the old and new data at this level and all finer levels,
       // and zero out the flux registers. 

       if (level == 0 || iteration > 1) {

	 if (subcycle_iter < 2)
	   for (int lev = level; lev < finest_level; lev++) {
               if (do_reflux && gravity->get_gravity_type() == "PoissonGrav")
                   gravity->zeroPhiFluxReg(lev+1);
           }

       }

       // Define the old gravity vector (aka grad_phi on cell centers)
       //   Note that this is based on the multilevel solve when doing "PoissonGrav".

       gravity->get_old_grav_vector(level,grav_old,time);
       
       if (gravity->get_gravity_type() == "PoissonGrav" &&
           gravity->test_results_of_solves() == 1)
          gravity->test_level_grad_phi_prev(level);
    }
    else
    {
       grav_new.setVal(0.0);
       phi_new.setVal(0.0);
    }
#endif

#ifdef REACTIONS
    MultiFab& reactions_old = get_old_data(Reactions_Type);
    MultiFab& reactions_new = get_new_data(Reactions_Type);
#endif
    
    // It's possible for interpolation to create very small negative values for
    //   species so we make sure here that all species are non-negative after this point
    enforce_nonnegative_species(S_old);

#ifdef DIFFUSION
#ifdef TAU
    MultiFab tau_diff(grids,1,NUM_GROW);
    tau_diff.setVal(0.);
    define_tau(tau_diff,grav_old,time);
#endif
#endif

#ifdef GRAVITY
    MultiFab comp_minus_level_phi;
    PArray<MultiFab> comp_minus_level_grad_phi(BL_SPACEDIM,PArrayManage);

    // Do level solve at beginning of time step in order to compute the
    //   difference between the multilevel and the single level solutions.
    if (do_grav && gravity->get_gravity_type() == "PoissonGrav")
    {
        if (gravity->NoComposite() != 1 && level < parent->finestLevel()) {
	    comp_minus_level_phi.define(grids,1,0,Fab_allocate);
	    for (int n=0; n<BL_SPACEDIM; ++n) {
		  comp_minus_level_grad_phi.set(n, new MultiFab(getEdgeBoxArray(n),1,0));
	    }
	    gravity->create_comp_minus_level_grad_phi(level,
                            comp_minus_level_phi,comp_minus_level_grad_phi);
        } else {
           if (verbose && ParallelDescriptor::IOProcessor()) {
              std::cout << " " << '\n';
              std::cout << "... old-time level solve at level " << level << '\n';
           }
	   int is_new = 0;
           gravity->solve_for_phi(level, 
				  phi_old,
				  gravity->get_grad_phi_prev(level),
				  is_new);
        }
    }
#endif

    //
    // Get pointers to Flux registers, or set pointer to zero if not there.
    //
    FluxRegister *fine    = 0;
    FluxRegister *current = 0;
    
    if (do_reflux && level < finest_level)
      fine = &getFluxReg(level+1);
    if (do_reflux && level > 0)
      current = &getFluxReg(level);
    
#ifdef SGS
    FluxRegister *sgs_fine    = 0;
    FluxRegister *sgs_current = 0;
    if (do_reflux && level < finest_level)
      sgs_fine = &getSGSFluxReg(level+1);
    if (do_reflux && level > 0)
      sgs_current = &getSGSFluxReg(level);
#endif
    
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
    
    MultiFab fluxes[3];

    // We want to define this on every grid, because we may need the fluxes
    // when computing the source terms later.

    for (int j = 0; j < BL_SPACEDIM; j++)
       {
         fluxes[j].define(getEdgeBoxArray(j), NUM_STATE, 0, Fab_allocate);
         fluxes[j].setVal(0.0);
       }

    for (int j = BL_SPACEDIM; j < 3; j++)
      {
	BoxArray ba = S_new.boxArray();
	fluxes[j].define(ba, NUM_STATE, 0, Fab_allocate);
	fluxes[j].setVal(0.0);
      }

#ifdef SGS
    // We need these even if we are single-level because they are used in the source construction.
    MultiFab sgs_fluxes[BL_SPACEDIM];
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
      {
	sgs_fluxes[dir].define(getEdgeBoxArray(dir), NUM_STATE, 0, Fab_allocate);
      }
#endif
    
#ifdef RADIATION
    MultiFab& Er_new = get_new_data(Rad_Type);
    MultiFab rad_fluxes[BL_SPACEDIM];
    if (Radiation::rad_hydro_combined) {
      for (int dir = 0; dir < BL_SPACEDIM; dir++) {
	rad_fluxes[dir].define(getEdgeBoxArray(dir), Radiation::nGroups, 0, Fab_allocate);
      }
    }
#endif


    
    // For the hydrodynamics update we need to have NUM_GROW ghost zones available,
    // but the state data does not carry ghost zones. So we use a FillPatch
    // using the state data to give us Sborder, which does have ghost zones. 

    MultiFab Sborder(grids,NUM_STATE,NUM_GROW,Fab_allocate);

    AmrLevel::FillPatch(*this,Sborder,NUM_GROW,prev_time,State_Type,0,NUM_STATE);

    // This array will hold the source terms that go into the hydro update through umdrv.
    
    MultiFab sources(grids,NUM_STATE,NUM_GROW,Fab_allocate);
    sources.setVal(0.0,NUM_GROW);

    // Set up external source terms
    
    MultiFab ext_src_old(grids,NUM_STATE,NUM_GROW,Fab_allocate);
    ext_src_old.setVal(0.0,NUM_GROW);

    MultiFab ext_src_new(grids,NUM_STATE,0,Fab_allocate);
    ext_src_new.setVal(0.0);

#ifdef SGS
    if (add_ext_src) {
      reset_old_sgs(dt);
      MultiFab& sgs_old = get_old_data(SGS_Type);
      getSource(prev_time,dt,Sborder,ext_src_old,sgs_old,sgs_fluxes);
    }
#else
    if (add_ext_src)
      getSource(prev_time,dt,Sborder,Sborder,ext_src_old,NUM_GROW);
#endif

#ifdef DIFFUSION
    MultiFab OldTempDiffTerm(grids,1,1);
    MultiFab OldSpecDiffTerm(grids,NumSpec,1);
    MultiFab OldViscousTermforMomentum(grids,BL_SPACEDIM,1);
    MultiFab OldViscousTermforEnergy(grids,1,1);
#ifdef TAU
    add_temp_diffusion_to_source(ext_src_old,OldTempDiffTerm,prev_time,tau_diff);
#else
    add_temp_diffusion_to_source(ext_src_old,OldTempDiffTerm,prev_time);
#endif
#if (BL_SPACEDIM == 1) 
    add_spec_diffusion_to_source(ext_src_old,OldSpecDiffTerm,prev_time);
    add_viscous_term_to_source(ext_src_old,OldViscousTermforMomentum,OldViscousTermforEnergy,prev_time);
#endif
#endif
    
    BoxLib::fill_boundary(ext_src_old, geom);    

    MultiFab::Add(sources,ext_src_old,0,0,NUM_STATE,NUM_GROW);    

#ifdef GRAVITY
    if (do_grav)
      add_force_to_sources(grav_old, sources, Sborder);
#endif

#ifdef ROTATION
    MultiFab& phirot_old = get_old_data(PhiRot_Type);
    MultiFab& rot_old = get_old_data(Rotation_Type);
    MultiFab& phirot_new = get_new_data(PhiRot_Type);
    MultiFab& rot_new = get_new_data(Rotation_Type);

    // Fill the old rotation data.    
    
    if (do_rotation) {
      
      fill_rotation_field(phirot_old, rot_old, Sborder, prev_time);
      
    } else {
      
      phirot_old.setVal(0.);
      rot_old.setVal(0.);
      
    } 

    if (do_rotation)
      add_force_to_sources(rot_old, sources, Sborder);
#endif
    
    
#ifdef POINTMASS
    Real mass_change_at_center = 0.;
#endif

    // Permit the user to update the sponge parameters as a function of time.

    if (do_sponge)
      update_sponge_params(&time);

    // Set up the time-rate of change of the source terms.
    
    MultiFab& dSdt_new = get_new_data(Source_Type);    
    
    // Optionally we can predict the source terms to t + dt/2,
    // which is the time-level n+1/2 value, To do this we use a
    // lagged predictor estimate: dS/dt_n = (S_n - S_{n-1}) / dt, so 
    // S_{n+1/2} = S_n + (dt / 2) * dS/dt_n.
      
    if (source_term_predictor == 1) {

      AmrLevel::FillPatch(*this,dSdt_new,NUM_GROW,cur_time,Source_Type,0,NUM_STATE);       
      
      dSdt_new.mult(dt / 2.0, NUM_GROW);

      MultiFab::Add(sources,dSdt_new,0,0,NUM_STATE,NUM_GROW);

    }

    {

      // Note that we do the react_half_dt on Sborder because of our Strang
      // splitting approach -- the "old" data sent to the hydro,
      // which Sborder represents, has already had a half-timestep of burning.      
      
#ifdef REACTIONS
#ifdef TAU
	react_half_dt(Sborder,reactions_old,tau_diff,time,dt,NUM_GROW);
#else
	react_half_dt(Sborder,reactions_old,time,dt,NUM_GROW);
#endif
#endif

	if (verbose && ParallelDescriptor::IOProcessor())
	  std::cout << "... Entering hydro advance" << std::endl << std::endl;

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
		radiation->compute_limiter(level, grids, Sborder, Erborder, lamborder, kappa_s);
	    }

	    int nstep_fsp = -1;

	    BL_PROFILE_VAR("Castro::advance_hydro_ca_umdrv_rad()", CA_UMDRV_RAD);

#ifdef _OPENMP
#ifdef POINTMASS
#pragma omp parallel reduction(+:mass_change_at_center)
#else
#pragma omp parallel
#endif
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
		    
		    FArrayBox &stateold = S_old[mfi];		    
		    FArrayBox &statein  = Sborder[mfi];
		    FArrayBox &stateout = S_new[mfi];
		
		    FArrayBox &Er = Erborder[mfi];
		    FArrayBox &lam = lamborder[mfi];
		    FArrayBox &Erout = Er_new[mfi];
		
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
			 BL_TO_FORTRAN(sources[mfi]),
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

		    // Add dt * old-time external source terms

		    stateout.saxpy(dt,ext_src_old[mfi],bx,bx,0,0,NUM_STATE);		    

		    // Gravitational source term for the time-level n data.

		    Real E_added_grav = 0.0;
		    Real mom_added[3] = { 0.0 };

#ifdef GRAVITY
		    if (do_grav)
		      ca_gsrc(ARLIM_3D(lo), ARLIM_3D(hi),
			      ARLIM_3D(domain_lo), ARLIM_3D(domain_hi),
			      BL_TO_FORTRAN_3D(phi_old[mfi]),
			      BL_TO_FORTRAN_3D(grav_old[mfi]),
			      BL_TO_FORTRAN_3D(stateold),
			      BL_TO_FORTRAN_3D(stateout),
			      ZFILL(dx),dt,&time,
			      E_added_grav,mom_added);
#endif		    

		    for (int dir = 0; dir < 3; dir++)
			 mom_added[dir] = 0.0;

		    // Rotational source term for the time-level n data.

		    Real E_added_rot = 0.0;
		    
#ifdef ROTATION		    
		    if (do_rotation)
		      ca_rsrc(ARLIM_3D(lo), ARLIM_3D(hi),
			      ARLIM_3D(domain_lo), ARLIM_3D(domain_hi),
			      BL_TO_FORTRAN_3D(phirot_old[mfi]),
			      BL_TO_FORTRAN_3D(rot_old[mfi]),
			      BL_TO_FORTRAN_3D(stateold),
			      BL_TO_FORTRAN_3D(stateout),
			      ZFILL(dx),dt,&time,
			      E_added_rot,mom_added);		    
#endif

		    for (int dir = 0; dir < 3; dir++)
			 mom_added[dir] = 0.0;

		    Real E_added_sponge = 0.0;
		    
		    if (do_sponge)
		      ca_sponge(ARLIM_3D(lo), ARLIM_3D(hi),
				BL_TO_FORTRAN_3D(stateout), ZFILL(dx), dt, &time,
				E_added_sponge,mom_added);

		    if (radiation->do_inelastic_scattering) {
			ca_inelastic_sct(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
					 BL_TO_FORTRAN_3D(stateout),
					 BL_TO_FORTRAN_3D(Erout),
					 BL_TO_FORTRAN_3D(kappa_s[mfi]),
					 dt);
		    }

		    for (int i = 0; i < BL_SPACEDIM ; i++) {
			u_gdnv[i][mfi].copy(ugdn[i],mfi.nodaltilebox(i));
		    }
		    
		    if (do_reflux) {
			for (int i = 0; i < BL_SPACEDIM ; i++) {
			    fluxes    [i][mfi].copy(    flux[i],mfi.nodaltilebox(i));
			    rad_fluxes[i][mfi].copy(rad_flux[i],mfi.nodaltilebox(i));
			}
		    }

#ifdef POINTMASS
		    if (level == finest_level)
			pm_compute_delta_mass
			    (&mass_change_at_center,
			     ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			     BL_TO_FORTRAN_3D(stateold),
			     BL_TO_FORTRAN_3D(stateout),
			     BL_TO_FORTRAN_3D(volume[mfi]),
			     ZFILL(geom.ProbLo()), ZFILL(dx), &time, &dt);
#endif	
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
	  
	    Real E_added_grav    = 0.;
	    Real E_added_flux    = 0.;
	    Real E_added_rot     = 0.;
	    Real E_added_sponge  = 0.;
	    Real mass_added      = 0.;
	    Real eint_added      = 0.;
	    Real eden_added      = 0.;
	    Real dens_change     = 1.e200;
	    Real xmom_added_flux = 0.;
	    Real ymom_added_flux = 0.;
	    Real zmom_added_flux = 0.;
	    Real xmom_added_grav = 0.;
	    Real ymom_added_grav = 0.;
	    Real zmom_added_grav = 0.;
	    Real xmom_added_rot  = 0.;
	    Real ymom_added_rot  = 0.;
	    Real zmom_added_rot  = 0.;
	    Real xmom_added_sponge = 0.;
	    Real ymom_added_sponge = 0.;
	    Real zmom_added_sponge = 0.;

	    BL_PROFILE_VAR("Castro::advance_hydro_ca_umdrv()", CA_UMDRV);
	    
#ifdef _OPENMP
#ifdef POINTMASS
#pragma omp parallel reduction(+:E_added_grav,E_added_flux,E_added_rot,E_added_sponge) \
                     reduction(+:mass_added,eint_added,eden_added) \
                     reduction(+:xmom_added_flux,ymom_added_flux,zmom_added_flux) \
                     reduction(+:xmom_added_grav,ymom_added_grav,zmom_added_grav) \
                     reduction(+:xmom_added_rot,ymom_added_rot,zmom_added_rot) \
                     reduction(+:xmom_added_sponge,ymom_added_sponge,zmom_added_sponge) \
                     reduction(+:mass_change_at_center) \
                     reduction(min:dens_change)
#else
#pragma omp parallel reduction(+:E_added_grav,E_added_flux,E_added_rot,E_added_sponge) \
                     reduction(+:mass_added,eint_added,eden_added) \
                     reduction(+:xmom_added_flux,ymom_added_flux,zmom_added_flux) \
                     reduction(+:xmom_added_grav,ymom_added_grav,zmom_added_grav) \
                     reduction(+:xmom_added_rot,ymom_added_rot,zmom_added_rot) \
                     reduction(+:xmom_added_sponge,ymom_added_sponge,zmom_added_sponge) \
                     reduction(min:dens_change)
#endif
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
		    
		    FArrayBox &stateold = S_old[mfi];
		    FArrayBox &statein  = Sborder[mfi];
		    FArrayBox &stateout = S_new[mfi];

		    // Allocate fabs for fluxes and Godunov velocities.
		    for (int i = 0; i < BL_SPACEDIM; i++) {
			const Box& bxtmp = BoxLib::surroundingNodes(bx,i);
			flux[i].resize(bxtmp,NUM_STATE);
			ugdn[i].resize(BoxLib::grow(bxtmp,1),1);
		    }
	  
		    ca_umdrv
			(&is_finest_level,&time,
			 lo, hi, domain_lo, domain_hi,
			 BL_TO_FORTRAN(statein), BL_TO_FORTRAN(stateout),
			 D_DECL(BL_TO_FORTRAN(ugdn[0]), 
				BL_TO_FORTRAN(ugdn[1]), 
				BL_TO_FORTRAN(ugdn[2])), 
			 BL_TO_FORTRAN(sources[mfi]),
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
			 BL_TO_FORTRAN(volume[mfi]), 
			 &cflLoc, verbose, 
			 mass_added, eint_added, eden_added,
			 dens_change,
			 xmom_added_flux, 
                  	 ymom_added_flux, 
	                 zmom_added_flux,
                         E_added_flux);

		    // Add dt * old-time external source terms

		    stateout.saxpy(dt,ext_src_old[mfi],bx,bx,0,0,NUM_STATE);

		    // Copy the normal velocities from the Riemann solver
		    
		    for (int i = 0; i < BL_SPACEDIM ; i++) {
			u_gdnv[i][mfi].copy(ugdn[i],mfi.nodaltilebox(i));
		    }

		    // Since we may need the fluxes later on, we'll copy them
		    // to the fluxes MultiFAB even if we aren't on a fine grid.

		    for (int i = 0; i < BL_SPACEDIM ; i++)
		      fluxes[i][mfi].copy(flux[i],mfi.nodaltilebox(i));	  

		    // Gravitational source term for the time-level n data.

		    Real mom_added[3] = { 0.0 };

#ifdef GRAVITY
		    if (do_grav)
		      ca_gsrc(ARLIM_3D(lo), ARLIM_3D(hi),
			      ARLIM_3D(domain_lo), ARLIM_3D(domain_hi),
			      BL_TO_FORTRAN_3D(phi_old[mfi]),
			      BL_TO_FORTRAN_3D(grav_old[mfi]),
			      BL_TO_FORTRAN_3D(stateold),
			      BL_TO_FORTRAN_3D(stateout),
			      ZFILL(dx),dt,&time,
			      E_added_grav,mom_added);
#endif

		    xmom_added_grav += mom_added[0];
		    ymom_added_grav += mom_added[1];
		    zmom_added_grav += mom_added[2];

		    for (int dir = 0; dir < 3; dir++)
			 mom_added[dir] = 0.0;

		    // Rotational source term for the time-level n data.

#ifdef ROTATION
		    if (do_rotation)
		      ca_rsrc(ARLIM_3D(lo), ARLIM_3D(hi),
			      ARLIM_3D(domain_lo), ARLIM_3D(domain_hi),
			      BL_TO_FORTRAN_3D(phirot_old[mfi]),
			      BL_TO_FORTRAN_3D(rot_old[mfi]),
			      BL_TO_FORTRAN_3D(stateold),
			      BL_TO_FORTRAN_3D(stateout),
			      ZFILL(dx),dt,&time,
			      E_added_rot,mom_added);
#endif

		    xmom_added_rot += mom_added[0];
		    ymom_added_rot += mom_added[1];
		    zmom_added_rot += mom_added[2];

		    for (int dir = 0; dir < 3; dir++)
			 mom_added[dir] = 0.0;

		    if (do_sponge)
		      ca_sponge(ARLIM_3D(lo), ARLIM_3D(hi),
				BL_TO_FORTRAN_3D(stateout), ZFILL(dx), dt, &time,
				E_added_sponge,mom_added);

		    xmom_added_sponge += mom_added[0];
		    ymom_added_sponge += mom_added[1];
		    zmom_added_sponge += mom_added[2];
			 
		    
#ifdef POINTMASS
		    if (level == finest_level)
			pm_compute_delta_mass
			    (&mass_change_at_center, 
			     ARLIM_3D(lo), ARLIM_3D(hi),
			     BL_TO_FORTRAN_3D(stateold),
			     BL_TO_FORTRAN_3D(stateout),
			     BL_TO_FORTRAN_3D(volume[mfi]),
			     ZFILL(geom.ProbLo()), ZFILL(dx), &time, &dt);
#endif
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

	    if (print_energy_diagnostics)
	    {
	       const Real cell_vol = D_TERM(dx[0], *dx[1], *dx[2]);
	       Real foo[19] = {mass_added, eint_added, eden_added, 
			       E_added_flux, E_added_grav, E_added_rot, E_added_sponge,
			       xmom_added_flux, ymom_added_flux, zmom_added_flux,
			       xmom_added_grav, ymom_added_grav, zmom_added_grav,
			       xmom_added_rot,  ymom_added_rot,  zmom_added_rot,
                               xmom_added_sponge, ymom_added_sponge, zmom_added_sponge};
#ifdef BL_LAZY
	       Lazy::QueueReduction( [=] () mutable {
#endif
	       ParallelDescriptor::ReduceRealSum(foo, 19, ParallelDescriptor::IOProcessorNumber());
	       if (ParallelDescriptor::IOProcessor()) 
	       {
		   mass_added = foo[0];
		   eint_added = foo[1];
		   eden_added = foo[2];
		   E_added_flux = foo[3];
		   E_added_grav = foo[4];
		   E_added_rot  = foo[5];
		   E_added_sponge = foo[6];
		   xmom_added_flux = foo[7];
		   ymom_added_flux = foo[8];
		   zmom_added_flux = foo[9];
		   xmom_added_grav = foo[10];
		   ymom_added_grav = foo[11];
		   zmom_added_grav = foo[12];
		   xmom_added_rot  = foo[13];
		   ymom_added_rot  = foo[14];
		   zmom_added_rot  = foo[15];
		   xmom_added_sponge  = foo[16];
		   ymom_added_sponge  = foo[17];
		   zmom_added_sponge  = foo[18];		   
		   if (std::abs(mass_added) != 0)
		   {
		      std::cout << "   Mass added from negative density correction : " << 
				    mass_added*cell_vol << std::endl;
		      std::cout << "(rho e) added from negative density correction : " << 
				    eint_added*cell_vol << std::endl;
		      std::cout << "(rho E) added from negative density correction : " << 
				    eden_added*cell_vol << std::endl;
		   }

		   std::cout << "(rho E) added from fluxes                      : " << 
				 E_added_flux << std::endl;
		   std::cout << "xmom added from fluxes                      : " << 
				 xmom_added_flux << std::endl;
		   std::cout << "ymom added from fluxes                      : " << 
				 ymom_added_flux << std::endl;
		   std::cout << "zmom added from fluxes                      : " << 
				 zmom_added_flux << std::endl;
#ifdef GRAVITY
		   if (do_grav) 
		   {	 
		      std::cout << "(rho E) added from grav. source terms          : " << 
				    E_added_grav*cell_vol << std::endl;
		      std::cout << "xmom added from grav. source terms             : " << 
				    xmom_added_grav*cell_vol << std::endl;
		      std::cout << "ymom added from grav. source terms             : " << 
				    ymom_added_grav*cell_vol << std::endl;
		      std::cout << "zmom added from grav. source terms             : " << 
				    zmom_added_grav*cell_vol << std::endl;
		   }
#endif
#ifdef ROTATION
		   if (do_rotation) 
		   {	 
		      std::cout << "(rho E) added from rot. source terms          : " << 
				    E_added_rot*cell_vol << std::endl;
		      std::cout << "xmom added from rot. source terms             : " << 
				    xmom_added_rot*cell_vol << std::endl;
		      std::cout << "ymom added from rot. source terms             : " << 
				    ymom_added_rot*cell_vol << std::endl;
		      std::cout << "zmom added from rot. source terms             : " << 
				    zmom_added_rot*cell_vol << std::endl;
		   }
#endif

		   if (do_sponge) 
		   {	 
		      std::cout << "(rho E) added from sponge                     : " << 
				    E_added_sponge*cell_vol << std::endl;
		      std::cout << "xmom added from sponge                        : " << 
				    xmom_added_sponge*cell_vol << std::endl;
		      std::cout << "ymom added from sponge                        : " << 
				    ymom_added_sponge*cell_vol << std::endl;
		      std::cout << "zmom added from sponge                        : " << 
				    zmom_added_sponge*cell_vol << std::endl;
		   }
	       }
#ifdef BL_LAZY
	       });
#endif
	    }

#endif    // RADIATION

	    if (verbose && ParallelDescriptor::IOProcessor())
	      std::cout << std::endl << "... Leaving hydro advance" << std::endl << std::endl;

    }

    if (use_retry)
      ParallelDescriptor::ReduceRealMin(frac_change);

#ifdef POINTMASS
    if (level == finest_level)
    {
          ParallelDescriptor::ReduceRealSum(mass_change_at_center);
  	  if (mass_change_at_center > 0.)
          {
	     point_mass += mass_change_at_center;
	     for (MFIter mfi(S_old); mfi.isValid(); ++mfi)
             {
		const Box& bx = mfi.validbox();

		pm_fix_solution
		  (ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		   BL_TO_FORTRAN_3D(S_old[mfi]), BL_TO_FORTRAN_3D(S_new[mfi]),
		   ZFILL(geom.ProbLo()), ZFILL(dx), &time, &dt);
             }
          }
    }
#endif

    if (do_reflux) {
	if (current) {
	    for (int i = 0; i < BL_SPACEDIM ; i++)
		current->FineAdd(fluxes[i],i,0,0,NUM_STATE,1.);
	}
	if (fine) {
	    for (int i = 0; i < BL_SPACEDIM ; i++)
	        fine->CrseInit(fluxes[i],i,0,0,NUM_STATE,-1.,FluxRegister::ADD);
	}
#ifdef RADIATION
	if (rad_current) {
	    for (int i = 0; i < BL_SPACEDIM ; i++)
		rad_current->FineAdd(rad_fluxes[i],i,0,0,Radiation::nGroups,1.);
	}
	if (rad_fine) {
	    for (int i = 0; i < BL_SPACEDIM ; i++)
	        rad_fine->CrseInit(rad_fluxes[i],i,0,0,Radiation::nGroups,-1.,FluxRegister::ADD);
        }
#endif
    }

    if (courno > 1.0) {
	std::cout << "WARNING -- EFFECTIVE CFL AT THIS LEVEL " << level << " IS " << courno << '\n';
	if (hard_cfl_limit == 1)
	  BoxLib::Abort("CFL is too high at this level -- go back to a checkpoint and restart with lower cfl number");
    }

    if (S_new.contains_nan(Density,S_new.nComp(),0,true))
    {
        for (int i = 0; i < S_new.nComp(); i++)
        {
	if (S_new.contains_nan(Density + i, 1, 0,true))
            {
  	        std::cout << "Testing " << desc_lst[State_Type].name(i) << " for NaNs " << std::endl;
                BoxLib::Abort("S_new has NaNs in this component::advance_hydro()");
            }
        }
    }

    // Now we'll start updating the dSdt MultiFab. First,
    // get rid of the dt/2 * dS/dt that we added from the last
    // timestep, then subtract the old sources data to get the
    // first half of the update for the next calculation of dS/dt.
    
    if (source_term_predictor == 1) {
      MultiFab::Subtract(sources,dSdt_new,0,0,NUM_STATE,NUM_GROW);
      dSdt_new.setVal(0.0, NUM_GROW);
      MultiFab::Subtract(dSdt_new,sources,0,0,NUM_STATE,0);
    }
    
    sources.setVal(0.0,NUM_GROW);       
    
#ifdef GRAVITY
    // Must define new value of "center" before we call new gravity solve or external source routine
    if (moving_center == 1)
       define_new_center(S_new,cur_time);
#endif

    if (add_ext_src)
      {

#ifdef SGS
           // Re-compute source at old time because we may have added something to ext_src_old
	reset_old_sgs(dt);
	getOldSource(prev_time,dt,ext_src_old,sgs_fluxes);
	
	// Add half of old fluxes to the flux register
	if (do_reflux)
	  {
	    if (finest_level > 0) 
	      {
		for (int dir = 0; dir < BL_SPACEDIM ; dir++)
		  sgs_fluxes[dir].mult(0.5);
		if (sgs_current)
                  {
		    for (int dir = 0; dir < BL_SPACEDIM ; dir++)
		      sgs_current->FineAdd(sgs_fluxes[dir],area[dir],dir,0,0,NUM_STATE,dt);
                  }
		
		if (sgs_fine)
                  {
		    for (int dir = 0; dir < BL_SPACEDIM ; dir++)
		      sgs_fine->CrseInit(sgs_fluxes[dir],area[dir],dir,0,0,NUM_STATE,-dt);
                  }
              }
	  }
#endif
	// Must compute new temperature in case it is needed in the source term evaluation
	computeTemp(S_new);

	// Compute source at new time (no ghost cells needed)
	
#if (BL_SPACEDIM > 1)
	// We need to make the new radial data now so that we can use it when we
	//   FillPatch in creating the new source
	if ( (level == 0) && (spherical_star == 1) ) {
	  int is_new = 1;
	  make_radial_data(is_new);
	}
#endif

#ifdef SGS
           // Need to put this line here so that the state going into the source calculation
           //  satisfies K > energy_sgs_min
	reset_new_sgs(dt);

	sgs_new = get_new_data(SGS_Type);
	
	getSource(cur_time,dt,S_new,ext_src_new,sgs_new,sgs_fluxes);

	// Add half of new fluxes to the flux register
	if (do_reflux)
	  {
	    if (finest_level > 0) 
	      {
		for (int dir = 0; dir < BL_SPACEDIM ; dir++)
		  sgs_fluxes[dir].mult(0.5);
		
		if (sgs_current)
                  {
		    for (int dir = 0; dir < BL_SPACEDIM ; dir++)
		      sgs_current->FineAdd(sgs_fluxes[dir],area[dir],dir,0,0,NUM_STATE,dt);
                  }
		
		if (sgs_fine)
                  {
		    for (int dir = 0; dir < BL_SPACEDIM ; dir++)
		      sgs_fine->CrseInit(sgs_fluxes[dir],area[dir],dir,0,0,NUM_STATE,-dt);
                  }
	      }
	  }
#else
	ext_src_new.setVal(0.0);	
	getSource(cur_time,dt,S_old,S_new,ext_src_new,0);
#endif
      }

#ifdef SGS
    
// old way: time-centering for ext_src, diffusion are separated.
    if (add_ext_src) {
	time_center_source_terms(S_new,ext_src_old,ext_src_new,dt);
	reset_new_sgs(dt);
	computeTemp(S_new);
    }
    
#else
    
// New way for non-SGS: time-centering for ext_src, diffusion are merged.
#ifdef DIFFUSION
    MultiFab& NewTempDiffTerm = OldTempDiffTerm;
    MultiFab& NewSpecDiffTerm = OldSpecDiffTerm;
    MultiFab& NewViscousTermforMomentum = OldViscousTermforMomentum;
    MultiFab& NewViscousTermforEnergy   = OldViscousTermforEnergy;
#ifdef TAU
    add_temp_diffusion_to_source(ext_src_new,NewTempDiffTerm,cur_time,tau_diff);
#else
    add_temp_diffusion_to_source(ext_src_new,NewTempDiffTerm,cur_time);
#endif
#if (BL_SPACEDIM == 1) 
    add_spec_diffusion_to_source(ext_src_new,NewSpecDiffTerm,cur_time);
    add_viscous_term_to_source(ext_src_new,NewViscousTermforMomentum,NewViscousTermforEnergy,cur_time);
#endif
#endif

    if (add_ext_src) {
      time_center_source_terms(S_new,ext_src_old,ext_src_new,dt);
      computeTemp(S_new);
    }
    
#endif

    MultiFab::Add(sources,ext_src_new,0,0,NUM_STATE,0);    
    
#ifdef GRAVITY
    if (do_grav)
      {
	if (gravity->get_gravity_type() == "PoissonGrav")
	  {
            if (verbose && ParallelDescriptor::IOProcessor()) {
	      std::cout << " " << '\n';
	      std::cout << "... new-time level solve at level " << level << '\n';
            }
	    
            // Here we use the "old" phi from the current time step as a guess for this solve
	    MultiFab::Copy(phi_new,phi_old,0,0,1,phi_new.nGrow());
            if ( level < parent->finestLevel() && (gravity->NoComposite() != 1) ) {
		phi_new.minus(comp_minus_level_phi, 0, 1, 0);
	    }
            int is_new = 1;
            gravity->solve_for_phi(level,
				   phi_new,
				   gravity->get_grad_phi_curr(level),
				   is_new);
	    
            if (gravity->test_results_of_solves() == 1)
	      {
		if (verbose && ParallelDescriptor::IOProcessor()) 
		  {
		    std::cout << " " << '\n';
		    std::cout << "... testing grad_phi_curr before adding comp_minus_level_grad_phi " << '\n';
		  }
		gravity->test_level_grad_phi_curr(level);
	      }
	    
            if ( level < parent->finestLevel() && (gravity->NoComposite() != 1) ) {
	      phi_new.plus(comp_minus_level_phi, 0, 1, 0);
	      gravity->plus_grad_phi_curr(level,comp_minus_level_grad_phi);
            }
	    
            if (gravity->test_results_of_solves() == 1) {
	      if (level < parent->finestLevel()) {
		if (verbose && ParallelDescriptor::IOProcessor()) {
		  std::cout << " " << '\n';
		  std::cout << "... testing grad_phi_curr after adding comp_minus_level_grad_phi " << '\n';
		}
		gravity->test_level_grad_phi_curr(level);
	      }
            }
	    
            if (do_reflux)  gravity->add_to_fluxes(level,iteration,ncycle);
	  }
	else {
	    phi_new.setVal(0.0);  // so that plotfiles do not contain nans
	}

	// Now do corrector part of source term update
	gravity->get_new_grav_vector(level,grav_new,cur_time);

        Real E_added    = 0.;
	Real xmom_added = 0.;
	Real ymom_added = 0.;
	Real zmom_added = 0.;

	const int* domlo = geom.Domain().loVect();
	const int* domhi = geom.Domain().hiVect();

#ifdef _OPENMP
#pragma omp parallel reduction(+:E_added,xmom_added,ymom_added,zmom_added)
#endif
	{
	    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
	    {
		const Box& bx = mfi.tilebox();
		
		Real mom_added[3] = { 0.0 };

		ca_corrgsrc(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			    ARLIM_3D(domlo), ARLIM_3D(domhi),
			    BL_TO_FORTRAN_3D(phi_old[mfi]),
			    BL_TO_FORTRAN_3D(phi_new[mfi]),
			    BL_TO_FORTRAN_3D(grav_old[mfi]),
			    BL_TO_FORTRAN_3D(grav_new[mfi]),
			    BL_TO_FORTRAN_3D(S_old[mfi]),
			    BL_TO_FORTRAN_3D(S_new[mfi]),
			    BL_TO_FORTRAN_3D(fluxes[0][mfi]),
			    BL_TO_FORTRAN_3D(fluxes[1][mfi]),
			    BL_TO_FORTRAN_3D(fluxes[2][mfi]),
			    ZFILL(dx),dt,&cur_time,
			    BL_TO_FORTRAN_3D(volume[mfi]),
			    E_added, mom_added);

		xmom_added += mom_added[0];
		ymom_added += mom_added[1];
		zmom_added += mom_added[2];
	    }
	}

        if (print_energy_diagnostics)
        {
	    const Real cell_vol = D_TERM(dx[0], *dx[1], *dx[2]);
	    Real foo[1+BL_SPACEDIM] = {E_added, D_DECL(xmom_added, ymom_added, zmom_added)};
#ifdef BL_LAZY
            Lazy::QueueReduction( [=] () mutable {
#endif
	    ParallelDescriptor::ReduceRealSum(foo, 1+BL_SPACEDIM, ParallelDescriptor::IOProcessorNumber());
	    if (ParallelDescriptor::IOProcessor()) {
		E_added = foo[0];
		D_EXPR(xmom_added = foo[1],
		       ymom_added = foo[2],
		       zmom_added = foo[3]);

		std::cout << "(rho E) added from grav. corr.  terms          : " << E_added*cell_vol << std::endl;
		std::cout << "xmom added from grav. corr. terms              : " << xmom_added*cell_vol << std::endl;
#if (BL_SPACEDIM >= 2)
		std::cout << "ymom added from grav. corr. terms              : " << ymom_added*cell_vol << std::endl;
#endif
#if (BL_SPACEDIM == 3)
		std::cout << "zmom added from grav. corr. terms              : " << zmom_added*cell_vol << std::endl;
#endif
	    }
#ifdef BL_LAZY
	    });
#endif
        }	

	// Add this to the source term array if we're using the source term predictor.
	// If not, don't bother because sources isn't actually used in the update after this point.

	if (source_term_predictor == 1)
	  add_force_to_sources(grav_new, sources, S_new);

	computeTemp(S_new);
      }
#endif

#ifdef SGS  // for non-SGS, diffusion has been time-centered.
#ifdef DIFFUSION
#ifdef TAU
    time_center_temp_diffusion(S_new, OldTempDiffTerm, cur_time, dt, tau_diff);
#else
    time_center_temp_diffusion(S_new, OldTempDiffTerm, cur_time, dt);
#endif
    time_center_spec_diffusion(S_new, OldSpecDiffTerm, cur_time, dt);
    time_center_viscous_terms(S_new, OldViscousTermforMomentum, OldViscousTermforEnergy, cur_time, dt);
#endif
#endif

#ifdef ROTATION
    if (do_rotation) {

        // Fill the new rotation data.
      
        fill_rotation_field(phirot_new, rot_new, S_new, cur_time);

	// Now do corrector part of rotation source term update

        Real E_added    = 0.;
	Real xmom_added = 0.;
	Real ymom_added = 0.;
	Real zmom_added = 0.;

	const int* domlo = geom.Domain().loVect();
	const int* domhi = geom.Domain().hiVect();

#ifdef _OPENMP
#pragma omp parallel reduction(+:E_added,xmom_added,ymom_added,zmom_added)
#endif
	{
	    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
	    {
		const Box& bx = mfi.tilebox();

		Real mom_added[3] = { 0.0 };

		ca_corrrsrc(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			    ARLIM_3D(domlo), ARLIM_3D(domhi),
			    BL_TO_FORTRAN_3D(phirot_old[mfi]),
			    BL_TO_FORTRAN_3D(phirot_new[mfi]),
			    BL_TO_FORTRAN_3D(rot_old[mfi]),
			    BL_TO_FORTRAN_3D(rot_new[mfi]),
			    BL_TO_FORTRAN_3D(S_old[mfi]),
			    BL_TO_FORTRAN_3D(S_new[mfi]),
			    BL_TO_FORTRAN_3D(fluxes[0][mfi]),
			    BL_TO_FORTRAN_3D(fluxes[1][mfi]),
			    BL_TO_FORTRAN_3D(fluxes[2][mfi]),
			    ZFILL(dx),dt,&cur_time,
			    BL_TO_FORTRAN_3D(volume[mfi]),
			    E_added,mom_added);

		xmom_added += mom_added[0];
		ymom_added += mom_added[1];
		zmom_added += mom_added[2];
	    }
	} 

        if (print_energy_diagnostics)
        {
	    const Real cell_vol = D_TERM(dx[0], *dx[1], *dx[2]);
	    Real foo[4] = {E_added, xmom_added, ymom_added, zmom_added};
#ifdef BL_LAZY
            Lazy::QueueReduction( [=] () mutable {
#endif
	    ParallelDescriptor::ReduceRealSum(foo, 4, ParallelDescriptor::IOProcessorNumber());
	    if (ParallelDescriptor::IOProcessor()) {
		E_added = foo[0];
		xmom_added = foo[1],
		ymom_added = foo[2],
		zmom_added = foo[3];

		std::cout << "(rho E) added from rot. corr.  terms          : " << E_added*cell_vol << std::endl;
		std::cout << "xmom added from rot. corr. terms              : " << xmom_added*cell_vol << std::endl;
		std::cout << "ymom added from rot. corr. terms              : " << ymom_added*cell_vol << std::endl;
		std::cout << "zmom added from rot. corr. terms              : " << zmom_added*cell_vol << std::endl;
	    }
#ifdef BL_LAZY
	    });
#endif
        }

	computeTemp(S_new);

    } else {

        phirot_new.setVal(0.0);
        rot_new.setVal(0.0);

    }

    // Add this to the source term array if we're using the source term predictor.
    // If not, don't bother because sources isn't actually used in the update after this point.

    if (source_term_predictor == 1)
      add_force_to_sources(rot_new, sources, S_new);

#endif

    reset_internal_energy(S_new);

    if (source_term_predictor == 1) {
    
      // Calculate the time derivative of the source terms.

      MultiFab::Add(dSdt_new,sources,0,0,NUM_STATE,0);
      
      dSdt_new.mult(1.0/dt);

    } 
    
#ifdef REACTIONS
#ifdef TAU
    react_half_dt(S_new,reactions_new,tau_diff,cur_time,dt);
#else
    react_half_dt(S_new,reactions_new,cur_time,dt);
#endif
#endif

#ifndef LEVELSET
    delete [] u_gdnv;
#endif

    return dt;
}

#ifndef SGS
Real
Castro::advance_no_hydro (Real time,
                          Real dt,
                          int  iteration,
                          int  ncycle,
			  int  subcycle_iter)
{
    BL_PROFILE("Castro::advance_no_hydro()");

    if (do_hydro)  
       BoxLib::Abort("In advance_no_hydro but do_hydro is true");
    
#ifdef SGS
    BoxLib::Abort("In advance_no_hydro but SGS is defined");
#endif

#ifdef RADIATION
    if (do_radiation) {
        // The option of whether to do a multilevel initialization is
        // controlled within the radiation class.  This step belongs
        // before the swap.
        radiation->pre_timestep(level);
    }
#endif
    
    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level) {
        //
        // Set reflux registers to zero.
        //
        getFluxReg(level+1).setVal(0.0);
    }

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

#ifdef RADIATION
    // Make sure these are filled to avoid check/plot file errors:
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

#if (BL_SPACEDIM > 1)
    if ( (level == 0) && (spherical_star == 1) ) {
       swap_outflow_data();
       int is_new = 0;
       make_radial_data(is_new);
    }
#endif

#ifdef GRAVITY
    // Old and new gravitational potential.
       
    MultiFab& phi_old = get_old_data(PhiGrav_Type);
    MultiFab& phi_new = get_new_data(PhiGrav_Type);

    MultiFab& grav_old = get_old_data(Gravity_Type);
    MultiFab& grav_new = get_new_data(Gravity_Type);
    
    if (do_grav) {
       if (do_reflux && level < finest_level && gravity->get_gravity_type() == "PoissonGrav" && subcycle_iter < 2)
           gravity->zeroPhiFluxReg(level+1);

       // Define the old gravity vector (aka grad_phi on cell centers)
       //   Note that this is based on the multilevel solve when doing "PoissonGrav".

       gravity->get_old_grav_vector(level,grav_old,time);

       if (gravity->get_gravity_type() == "PoissonGrav" &&
           gravity->test_results_of_solves() == 1)
          gravity->test_level_grad_phi_prev(level);
    }
    else
    {
       grav_new.setVal(0.0);
       phi_new.setVal(0.0);
    }
#endif 

#ifdef REACTIONS
    MultiFab& reactions_old = get_old_data(Reactions_Type);
    MultiFab& reactions_new = get_new_data(Reactions_Type);
#endif    
    
#ifdef DIFFUSION
#ifdef TAU
    MultiFab tau_diff(grids,1,1);
    tau_diff.setVal(0.);
    define_tau(tau_diff,grav_old,time);
#endif
#endif

    const Real prev_time = state[State_Type].prevTime();

    // It's possible for interpolation to create very small negative values for
    //   species so we make sure here that all species are non-negative after this point
    enforce_nonnegative_species(S_old);

#ifdef REACTIONS
#ifdef TAU
    react_half_dt(S_old,reactions_old,tau_diff,time,dt);
#else
    react_half_dt(S_old,reactions_old,time,dt);
#endif
#endif
 
#ifdef GRAVITY
    MultiFab comp_minus_level_phi;
    PArray<MultiFab> comp_minus_level_grad_phi(BL_SPACEDIM,PArrayManage);

    // Do level solve at beginning of time step in order to compute the
    //   difference between the multilevel and the single level solutions.
    if (do_grav && gravity->get_gravity_type() == "PoissonGrav")
    {
        if (gravity->NoComposite() != 1 && level < parent->finestLevel()) {
	    comp_minus_level_phi.define(grids,1,0,Fab_allocate);
	    for (int n=0; n<BL_SPACEDIM; ++n) {
		comp_minus_level_grad_phi.set(n, new MultiFab(getEdgeBoxArray(n),1,0));
	    }
	    gravity->create_comp_minus_level_grad_phi(level,
                            comp_minus_level_phi,comp_minus_level_grad_phi);
        } else {
           if (verbose && ParallelDescriptor::IOProcessor()) {
              std::cout << " " << '\n';
              std::cout << "... old-time level solve at level " << level << '\n';
           }
           int is_new = 0;
           gravity->solve_for_phi(level,
				  phi_old,
				  gravity->get_grad_phi_prev(level),
				  is_new);
        }
    }
#endif

    Real cur_time = state[State_Type].curTime();
        
    // Copy old data into new data.
    MultiFab::Copy(S_new, S_old, 0, 0, NUM_STATE, S_new.nGrow());
        
    if (add_ext_src) {
           MultiFab ext_src_old(grids,NUM_STATE,0,Fab_allocate);
           getSource(prev_time,dt,S_old,S_old,ext_src_old,0);
           ext_src_old.mult(dt);
           MultiFab::Add(S_new,ext_src_old,0,0,NUM_STATE,0);

           // Must compute new temperature in case it is needed in the source term evaluation
           computeTemp(S_new);

           // Compute source at new time
           MultiFab ext_src_new(grids,NUM_STATE,0,Fab_allocate);
           getSource(cur_time,dt,S_old,S_new,ext_src_new,0);

           ext_src_old.mult(-0.5);
           ext_src_new.mult( 0.5*dt);

           // Subtract off half of the old source term, and add half of the new.
           MultiFab::Add(S_new,ext_src_old,0,0,S_new.nComp(),0);
           MultiFab::Add(S_new,ext_src_new,0,0,S_new.nComp(),0);
    }

    computeTemp(S_new);

#ifdef GRAVITY
    MultiFab::Copy(phi_new, phi_old, 0, 0, 1, phi_new.nGrow());
#endif

    // ******************
    //  Note: If add_ext_src changed the density, 
    //        we need to update gravity and the associated source terms and
    //        state variables here
    // ******************

#ifdef DIFFUSION
#ifdef TAU
        full_temp_diffusion_update(S_new,prev_time,cur_time,dt,tau_diff);
#else
        full_temp_diffusion_update(S_new,prev_time,cur_time,dt);
#endif
        full_spec_diffusion_update(S_new,prev_time,cur_time,dt);
#endif
        
#ifdef REACTIONS
#ifdef TAU
       react_half_dt(S_new,reactions_new,tau_diff,cur_time,dt);
#else
       react_half_dt(S_new,reactions_new,cur_time,dt);
#endif
#endif

#ifdef RADIATION
    if (Radiation::rad_hydro_combined) {
      MultiFab& Er_old = get_old_data(Rad_Type);
      MultiFab& Er_new = get_new_data(Rad_Type);
      Er_new.copy(Er_old);
    }
#endif

    return dt;
}
#endif
