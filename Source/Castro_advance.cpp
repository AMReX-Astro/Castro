#include <winstd.H>

#include "Castro.H"
#include "Castro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

#ifdef PARTICLES
#include <Particles_F.H>
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

    if (do_hydro) 
    {
        dt_new = advance_hydro(time,dt,iteration,ncycle);
    } 
    else 
    {
#ifdef SGS
        BoxLib::Abort("Castro::advance -- doesn't make sense to have SGS defined but not do_hydro");
        return 0.;
#else
        dt_new = advance_no_hydro(time,dt,iteration,ncycle);
#endif
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

    return dt_new;
}

Real
Castro::advance_hydro (Real time,
                       Real dt,
                       int  iteration,
                       int  ncycle)
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

    u_gdnv = new MultiFab[BL_SPACEDIM];
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
	BoxArray edge_grids(grids);
	edge_grids.surroundingNodes(dir);
	u_gdnv[dir].define(edge_grids,1,1,Fab_allocate);
	u_gdnv[dir].setVal(1.e40,1);
    }
    
    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level) {
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

    // Swap the new data from the last timestep into the old state data.
    // If we're on level 0, do it for all levels below this one as well.
    // Or, if we're on a later iteration at a finer timestep, swap for all
    // lower time levels as well.

    if (level == 0 || iteration > 1) {
        for (int lev = level; lev <= finest_level; lev++) {
            Real dt_lev = parent->dtLevel(lev);
            for (int k = 0; k < NUM_STATE_TYPE; k++) {
	        getLevel(lev).state[k].allocOldData();
                getLevel(lev).state[k].swapTimeLevels(dt_lev);
            }
        } 
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
                std::cout << "Testing component i for NaNs: " << i << std::endl;
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
       BL_FORT_PROC_CALL(SWAP_OUTFLOW_DATA,swap_outflow_data)();
       int is_new = 0;
       make_radial_data(is_new);
    }
#endif

    MultiFab grav_vec_old;

#ifdef GRAVITY
    // Old and new gravitational potential.
       
    MultiFab& phi_old = get_old_data(PhiGrav_Type);
    MultiFab& phi_new = get_new_data(PhiGrav_Type);
    
    if (do_grav) {

       // Swap the old and new data at this level and all finer levels,
       // and zero out the flux registers. 

       if (level == 0 || iteration > 1) {
	   for (int lev = level; lev < finest_level; lev++) {
               if (do_reflux && gravity->get_gravity_type() == "PoissonGrav")
                   gravity->zeroPhiFluxReg(lev+1);
           }

           for (int lev = level; lev <= finest_level; lev++)
               gravity->swapTimeLevels(lev);
       }

       grav_vec_old.define(grids,BL_SPACEDIM,NUM_GROW,Fab_allocate); 
       
       // Define the old gravity vector (aka grad_phi on cell centers)
       //   Note that this is based on the multilevel solve when doing "PoissonGrav".

       gravity->get_old_grav_vector(level,grav_vec_old,time);
       
       if (gravity->get_gravity_type() == "PoissonGrav" &&
           gravity->test_results_of_solves() == 1)
          gravity->test_level_grad_phi_prev(level);
    }
    else
    {
       MultiFab& new_grav_mf = get_new_data(Gravity_Type);
       new_grav_mf.setVal(0.0);
    }
#endif

    // It's possible for interpolation to create very small negative values for
    //   species so we make sure here that all species are non-negative after this point
    enforce_nonnegative_species(S_old);

#ifdef DIFFUSION
#ifdef TAU
    MultiFab tau_diff(grids,1,NUM_GROW);
    tau_diff.setVal(0.);
    define_tau(tau_diff,grav_vec_old,time);
#endif
#endif

#ifdef GRAVITY
    MultiFab comp_minus_level_phi(grids,1,0,Fab_allocate);
    PArray<MultiFab> comp_minus_level_grad_phi(BL_SPACEDIM,PArrayManage);
    for (int n=0; n<BL_SPACEDIM; ++n)
    {
        comp_minus_level_grad_phi.clear(n);
        comp_minus_level_grad_phi.set(n,new MultiFab(BoxArray(grids).surroundingNodes(n),1,0));
    }

    // Do level solve at beginning of time step in order to compute the
    //   difference between the multilevel and the single level solutions.
    if (do_grav)
    {
        if (gravity->NoComposite() != 1 && level < parent->finestLevel()) {
           gravity->create_comp_minus_level_grad_phi(level,
                           comp_minus_level_phi,comp_minus_level_grad_phi);
        } else {
           if (verbose && ParallelDescriptor::IOProcessor()) {
              std::cout << " " << '\n';
              std::cout << "... old-time level solve at level " << level << '\n';
           }
           int fill_interior = 0;
           gravity->solve_for_old_phi(level,phi_old,gravity->get_grad_phi_prev(level),fill_interior);
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
    
    MultiFab fluxes[BL_SPACEDIM];

    // We want to define this on every grid, because we may need the fluxes
    // when computing the source terms later.

    for (int j = 0; j < BL_SPACEDIM; j++)
       {
         BoxArray ba = S_new.boxArray();
         ba.surroundingNodes(j);
         fluxes[j].define(ba, NUM_STATE, 0, Fab_allocate);
         fluxes[j].setVal(0.0);
       }

#ifdef SGS
    // We need these even if we are single-level because they are used in the source construction.
    MultiFab sgs_fluxes[BL_SPACEDIM];
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
      {
	BoxArray ba = S_new.boxArray();
	ba.surroundingNodes(dir);
	sgs_fluxes[dir].define(ba, NUM_STATE, 0, Fab_allocate);
      }
#endif
    
#ifdef RADIATION
    MultiFab& Er_new = get_new_data(Rad_Type);
    MultiFab rad_fluxes[BL_SPACEDIM];
    if (Radiation::rad_hydro_combined) {
      for (int dir = 0; dir < BL_SPACEDIM; dir++) {
	BoxArray ba = Er_new.boxArray();
	ba.surroundingNodes(dir);
	rad_fluxes[dir].define(ba, Radiation::nGroups, 0, Fab_allocate);
      }
    }
#endif
    
    // Define the gravity vector so we can pass this to ca_umdrv.
    MultiFab grav_vector(grids,BL_SPACEDIM,NUM_GROW);
    grav_vector.setVal(0.);
    
#ifdef GRAVITY
    if (do_grav) {
      // Copy the gravity vector (including NUM_GROW ghost cells) for passing to umdrv.
      MultiFab::Copy(grav_vector,grav_vec_old,0,0,BL_SPACEDIM,NUM_GROW);
    }
#endif
    
#ifdef POINTMASS
    Real mass_change_at_center = 0.;
#endif
    
    MultiFab ext_src_old(grids,NUM_STATE,1,Fab_allocate);
    ext_src_old.setVal(0.0);
    
    MultiFab ext_src_new;
    bool have_source_terms = add_ext_src;
#ifdef DIFFUSION
    have_source_terms = true;
#endif
    if (have_source_terms) {
	ext_src_new.define(grids,NUM_STATE,0,Fab_allocate);
	ext_src_new.setVal(0.0);
    }

#ifdef SGS
    if (add_ext_src) {
      reset_old_sgs(dt);
      getOldSource(prev_time,dt,ext_src_old,sgs_fluxes);
    }
#else
    if (add_ext_src)
      getOldSource(prev_time,dt,ext_src_old);
#endif
    
#ifdef DIFFUSION
    MultiFab OldTempDiffTerm(grids,1,1);
#ifdef TAU
    add_diffusion_to_source(ext_src_old,OldTempDiffTerm,prev_time,tau_diff);
#else
    add_diffusion_to_source(ext_src_old,OldTempDiffTerm,prev_time);
#endif
#endif

    ext_src_old.FillBoundary();
    geom.FillPeriodicBoundary(ext_src_old,0,NUM_STATE,true);

    { // limit the scope of fpi
	FillPatchIterator fpi(*this, S_new, NUM_GROW, time, State_Type, 0, NUM_STATE);
	MultiFab& Sborder = fpi.get_mf();
    
#ifdef REACTIONS
	// Make sure to zero these even if do_react == 0.
	MultiFab& ReactMF_old = get_old_data(Reactions_Type);
	MultiFab& ReactMF     = get_new_data(Reactions_Type);
	ReactMF_old.setVal(0.);
	ReactMF.setVal(0.);
	
	for (MFIter mfi(Sborder); mfi.isValid(); ++mfi) {
#ifdef TAU
	    react_first_half_dt(Sborder[mfi],ReactMF[mfi],tau_diff[mfi],time,dt);
#else
	    react_first_half_dt(Sborder[mfi],ReactMF[mfi],time,dt);
#endif
	}
#endif
      
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

		    FArrayBox &state = Sborder[mfi];
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
		
		    BL_FORT_PROC_CALL(CA_UMDRV_RAD,ca_umdrv_rad)
			(&is_finest_level,&time,
			 bx.loVect(), bx.hiVect(),
			 domain_lo, domain_hi,
			 BL_TO_FORTRAN(state), BL_TO_FORTRAN(stateout),
			 BL_TO_FORTRAN(Er), BL_TO_FORTRAN(lam),
			 BL_TO_FORTRAN(Erout),
			 D_DECL(BL_TO_FORTRAN(ugdn[0]), 
				BL_TO_FORTRAN(ugdn[1]), 
				BL_TO_FORTRAN(ugdn[2])), 
			 BL_TO_FORTRAN(ext_src_old[mfi]),
			 BL_TO_FORTRAN(grav_vector[mfi]), 
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

		    if (radiation->do_inelastic_scattering) {
			BL_FORT_PROC_CALL(CA_INELASTIC_SCT, ca_inelastic_sct)
			    (ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
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
			BL_FORT_PROC_CALL(PM_COMPUTE_DELTA_MASS,pm_compute_delta_mass)
			    (&mass_change_at_center,
			     bx.loVect(), bx.hiVect(),
			     BL_TO_FORTRAN(state), BL_TO_FORTRAN(stateout),
			     BL_TO_FORTRAN(volume[mfi]),
			     geom.ProbLo(), dx, &time, &dt);
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
#endif
	    Real E_added_grav    = 0.;
	    Real E_added_flux    = 0.;
	    Real E_added_rot     = 0.;
	    Real E_added_sponge  = 0.;
	    Real mass_added      = 0.;
	    Real eint_added      = 0.;
	    Real eden_added      = 0.;
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
                     reduction(+:mass_change_at_center)
#else
#pragma omp parallel reduction(+:E_added_grav,E_added_flux,E_added_rot,E_added_sponge) \
                     reduction(+:mass_added,eint_added,eden_added) \
                     reduction(+:xmom_added_flux,ymom_added_flux,zmom_added_flux) \
                     reduction(+:xmom_added_grav,ymom_added_grav,zmom_added_grav) \
                     reduction(+:xmom_added_rot,ymom_added_rot,zmom_added_rot) \
                     reduction(+:xmom_added_sponge,ymom_added_sponge,zmom_added_sponge)
#endif
#endif
	    {
		FArrayBox flux[BL_SPACEDIM], ugdn[BL_SPACEDIM];
		
		Real cflLoc = -1.0e+200;
		int is_finest_level = (level == finest_level) ? 1 : 0;
		const int*  domain_lo = geom.Domain().loVect();
		const int*  domain_hi = geom.Domain().hiVect();
		
		for (MFIter mfi(S_new,hydro_tile_size); mfi.isValid(); ++mfi)
		{
		    const Box& bx  = mfi.tilebox();
		    
		    FArrayBox &state = Sborder[mfi];
		    FArrayBox &stateout = S_new[mfi];
		    
		    // Allocate fabs for fluxes and Godunov velocities.
		    for (int i = 0; i < BL_SPACEDIM ; i++) {
			const Box& bxtmp = BoxLib::surroundingNodes(bx,i);
			flux[i].resize(bxtmp,NUM_STATE);
			ugdn[i].resize(BoxLib::grow(bxtmp,1),1);
		    }
	  
		    BL_FORT_PROC_CALL(CA_UMDRV,ca_umdrv)
			(&is_finest_level,&time,
			 bx.loVect(), bx.hiVect(),
			 domain_lo, domain_hi,
			 BL_TO_FORTRAN(state), BL_TO_FORTRAN(stateout),
			 D_DECL(BL_TO_FORTRAN(ugdn[0]), 
				BL_TO_FORTRAN(ugdn[1]), 
				BL_TO_FORTRAN(ugdn[2])), 
			 BL_TO_FORTRAN(ext_src_old[mfi]),
			 BL_TO_FORTRAN(grav_vector[mfi]), 
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
			 xmom_added_flux, 
#if (BL_SPACEDIM >= 2)
                  	 ymom_added_flux, 
#endif
#if (BL_SPACEDIM == 3)
	                 zmom_added_flux,
#endif
	                 xmom_added_grav, 
#if (BL_SPACEDIM >= 2)
	                 ymom_added_grav, 
#endif
#if (BL_SPACEDIM == 3)
	                 zmom_added_grav,
#endif
#if (BL_SPACEDIM > 1)
	                 xmom_added_rot,
#endif
#if (BL_SPACEDIM >= 2)
	                 ymom_added_rot,  
#endif
#if (BL_SPACEDIM == 3)
	                 zmom_added_rot,
#endif
	                 xmom_added_sponge, 
#if (BL_SPACEDIM >= 2)
	                 ymom_added_sponge, 
#endif
#if (BL_SPACEDIM == 3)
	                 zmom_added_sponge,
#endif
#if (BL_SPACEDIM > 1)
	                 E_added_rot, 
#endif
                         E_added_flux, E_added_grav, E_added_sponge);
	    
		    for (int i = 0; i < BL_SPACEDIM ; i++) {
			u_gdnv[i][mfi].copy(ugdn[i],mfi.nodaltilebox(i));
		    }

		    // Since we may need the fluxes later on, we'll copy them
		    // to the fluxes MultiFAB even if we aren't on a fine grid.

		    for (int i = 0; i < BL_SPACEDIM ; i++)
		      fluxes[i][mfi].copy(flux[i],mfi.nodaltilebox(i));	  
		    
#ifdef POINTMASS
		    if (level == finest_level)
			BL_FORT_PROC_CALL(PM_COMPUTE_DELTA_MASS,pm_compute_delta_mass)
			    (&mass_change_at_center, 
			     bx.loVect(), bx.hiVect(),
			     BL_TO_FORTRAN(state), BL_TO_FORTRAN(stateout),
			     BL_TO_FORTRAN(volume[mfi]),
			     geom.ProbLo(), dx, &time, &dt);
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
				 E_added_flux*cell_vol << std::endl;
		   std::cout << "xmom added from fluxes                      : " << 
				 xmom_added_flux*cell_vol << std::endl;
#if (BL_SPACEDIM >= 2)
		   std::cout << "ymom added from fluxes                      : " << 
				 ymom_added_flux*cell_vol << std::endl;
#endif
#if (BL_SPACEDIM == 3)
		   std::cout << "zmom added from fluxes                      : " << 
				 zmom_added_flux*cell_vol << std::endl;
#endif
#ifdef GRAVITY
		   if (do_grav) 
		   {	 
		      std::cout << "(rho E) added from grav. source terms          : " << 
				    E_added_grav*cell_vol << std::endl;
		      std::cout << "xmom added from grav. source terms             : " << 
				    xmom_added_grav*cell_vol << std::endl;
#if (BL_SPACEDIM >= 2)
		      std::cout << "ymom added from grav. source terms             : " << 
				    ymom_added_grav*cell_vol << std::endl;
#endif
#if (BL_SPACEDIM == 3)
		      std::cout << "zmom added from grav. source terms             : " << 
				    zmom_added_grav*cell_vol << std::endl;
#endif
		   }
#endif
#ifdef ROTATION
		   if (do_rotation) 
		   {	 
		      std::cout << "(rho E) added from rot. source terms          : " << 
				    E_added_rot*cell_vol << std::endl;
#if (BL_SPACEDIM >= 2)
		      std::cout << "xmom added from rot. source terms             : " << 
				    xmom_added_rot*cell_vol << std::endl;
		      std::cout << "ymom added from rot. source terms             : " << 
				    ymom_added_rot*cell_vol << std::endl;
#endif
#if (BL_SPACEDIM == 3)
		      std::cout << "zmom added from rot. source terms             : " << 
				    zmom_added_rot*cell_vol << std::endl;
#endif
		   }
#endif

		   if (do_sponge) 
		   {	 
		      std::cout << "(rho E) added from sponge                     : " << 
				    E_added_sponge*cell_vol << std::endl;
		      std::cout << "xmom added from sponge                        : " << 
				    xmom_added_sponge*cell_vol << std::endl;
#if (BL_SPACEDIM >= 2)
		      std::cout << "ymom added from sponge                        : " << 
				    ymom_added_sponge*cell_vol << std::endl;
#endif
#if (BL_SPACEDIM == 3)
		      std::cout << "zmom added from sponge                        : " << 
				    zmom_added_sponge*cell_vol << std::endl;
#endif
		   }
	       }
#ifdef BL_LAZY
	       });
#endif
	    }

#ifdef RADIATION
	}
#endif
    }

#ifdef PARTICLES
    if (do_dm_particles && particle_move_type == "Gravitational")
    {
	BL_ASSERT(level == 0);
	Castro::theDMPC()->movePredict(grav_vector, level, dt);
    }
#endif
    
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
		BL_FORT_PROC_CALL(PM_FIX_SOLUTION,pm_fix_solution)
		  (bx.loVect(), bx.hiVect(),
		   BL_TO_FORTRAN(S_old[mfi]), BL_TO_FORTRAN(S_new[mfi]),
		   geom.ProbLo(), dx, &time, &dt);
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
		fine->CrseInit(fluxes[i],i,0,0,NUM_STATE,-1.);
	}
#ifdef RADIATION
	if (rad_current) {
	    for (int i = 0; i < BL_SPACEDIM ; i++)
		rad_current->FineAdd(rad_fluxes[i],i,0,0,Radiation::nGroups,1.);
	}
	if (rad_fine) {
	    for (int i = 0; i < BL_SPACEDIM ; i++)
	        rad_fine->CrseInit(rad_fluxes[i],i,0,0,Radiation::nGroups,-1.);
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
                std::cout << "Testing component i for NaNs: " << i << std::endl;
                BoxLib::Abort("S_new has NaNs in this component::advance_hydro()");
            }
        }
    }

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
	
	getNewSource(cur_time,dt,ext_src_new,sgs_fluxes);

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
	getNewSource(prev_time,cur_time,dt,ext_src_new);
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
#ifdef TAU
    add_diffusion_to_source(ext_src_new,NewTempDiffTerm,cur_time,tau_diff);
#else
    add_diffusion_to_source(ext_src_new,NewTempDiffTerm,cur_time);
#endif
#endif

#endif

    if (have_source_terms) {
	time_center_source_terms(S_new,ext_src_old,ext_src_new,dt);
	computeTemp(S_new);
    }

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
	    MultiFab::Copy(phi_new,phi_old,0,0,1,0);
            int fill_interior = 0;
            gravity->solve_for_new_phi(level,phi_new,gravity->get_grad_phi_curr(level),fill_interior);
	    
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

	// Now do corrector part of source term update
	MultiFab grav_vec_new(grids,BL_SPACEDIM,1,Fab_allocate);
	gravity->get_new_grav_vector(level,grav_vec_new,cur_time);

        Real E_added = 0.;
	Real xmom_added      = 0.;
	Real ymom_added      = 0.;
	Real zmom_added      = 0.;

#ifdef _OPENMP
#pragma omp parallel reduction(+:E_added,xmom_added,ymom_added,zmom_added)
#endif
	{
	    FArrayBox single_cell_fab(Box(IntVect::TheZeroVector(),IntVect::TheZeroVector()));
	    
	    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
	    {
		const Box& bx = mfi.tilebox();
		
		BL_FORT_PROC_CALL(CA_CORRGSRC,ca_corrgsrc)
		    (bx.loVect(), bx.hiVect(),
		     BL_TO_FORTRAN(grav_vec_old[mfi]),
		     BL_TO_FORTRAN(grav_vec_new[mfi]),
		     BL_TO_FORTRAN(S_old[mfi]),
		     BL_TO_FORTRAN(S_new[mfi]),
#if (BL_SPACEDIM == 3)
		     BL_TO_FORTRAN(phi_old[mfi]),
		     BL_TO_FORTRAN(phi_new[mfi]),
		     BL_TO_FORTRAN(fluxes[0][mfi]),
		     BL_TO_FORTRAN(fluxes[1][mfi]),
		     BL_TO_FORTRAN(fluxes[2][mfi]),
#endif
		     dx,dt,
		     BL_TO_FORTRAN(volume[mfi]),
		     xmom_added,
#if (BL_SPACEDIM >= 2)
		     ymom_added,
#endif
#if (BL_SPACEDIM == 3)
		     zmom_added,
#endif
		     E_added);
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

	computeTemp(S_new);
      }
#endif
    
#ifdef SGS  // for non-SGS, diffusion has been time-centered.
#ifdef DIFFUSION
#ifdef TAU
    time_center_diffusion(S_new, OldTempDiffTerm, cur_time, dt, tau_diff);
#else
    time_center_diffusion(S_new, OldTempDiffTerm, cur_time, dt);
#endif
#endif
#endif

#ifdef ROTATION
#if (BL_SPACEDIM > 1)
    if (do_rotation)
      {
	// Now do corrector part of rotation source term update

        Real E_added    = 0.;
	Real xmom_added = 0.;
	Real ymom_added = 0.;
	Real zmom_added = 0.;

#ifdef _OPENMP
#pragma omp parallel reduction(+:E_added,xmom_added,ymom_added,zmom_added)
#endif
	{
	    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
	    {
		const Box& bx = mfi.tilebox();
		
		BL_FORT_PROC_CALL(CA_CORRRSRC,ca_corrrsrc)
		    (bx.loVect(), bx.hiVect(),
		     BL_TO_FORTRAN(S_old[mfi]),
		     BL_TO_FORTRAN(S_new[mfi]),
		     BL_TO_FORTRAN(fluxes[0][mfi]),
#if (BL_SPACEDIM >= 2)
		     BL_TO_FORTRAN(fluxes[1][mfi]),
#endif
#if (BL_SPACEDIM == 3)
		     BL_TO_FORTRAN(fluxes[2][mfi]),
#endif
		     dx,dt,
		     BL_TO_FORTRAN(volume[mfi]),
		     xmom_added,
#if (BL_SPACEDIM >= 2)
		     ymom_added,
#endif
#if (BL_SPACEDIM == 3)
		     zmom_added,
#endif
		     E_added);
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

		std::cout << "(rho E) added from rot. corr.  terms          : " << E_added*cell_vol << std::endl;
		std::cout << "xmom added from rot. corr. terms              : " << xmom_added*cell_vol << std::endl;
#if (BL_SPACEDIM >= 2)	       
		std::cout << "ymom added from rot. corr. terms              : " << ymom_added*cell_vol << std::endl;
#endif
#if (BL_SPACEDIM == 3)
		std::cout << "zmom added from rot. corr. terms              : " << zmom_added*cell_vol << std::endl;
#endif
	    }
#ifdef BL_LAZY
	    });
#endif
        }	

	computeTemp(S_new);
      }
#endif
#endif

    reset_internal_energy(S_new);
    
#ifdef REACTIONS
#ifdef TAU
    react_second_half_dt(S_new,tau_diff,cur_time,dt,1);
#else
    react_second_half_dt(S_new,cur_time,dt,1);
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
                          int  ncycle)
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

    for (int k = 0; k < NUM_STATE_TYPE; k++) {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
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
       BL_FORT_PROC_CALL(SWAP_OUTFLOW_DATA,swap_outflow_data)();
       int is_new = 0;
       make_radial_data(is_new);
    }
#endif

    MultiFab grav_vec_old;

#ifdef GRAVITY
    // Old and new gravitational potential.
       
    MultiFab& phi_old = get_old_data(PhiGrav_Type);
    MultiFab& phi_new = get_new_data(PhiGrav_Type);
    
    if (do_grav) {
       if (do_reflux && level < finest_level && gravity->get_gravity_type() == "PoissonGrav")
           gravity->zeroPhiFluxReg(level+1);

       gravity->swapTimeLevels(level);

       grav_vec_old.define(grids,BL_SPACEDIM,NUM_GROW,Fab_allocate); 

       // Define the old gravity vector (aka grad_phi on cell centers)
       //   Note that this is based on the multilevel solve when doing "PoissonGrav".

       gravity->get_old_grav_vector(level,grav_vec_old,time);

       if (gravity->get_gravity_type() == "PoissonGrav" &&
           gravity->test_results_of_solves() == 1)
          gravity->test_level_grad_phi_prev(level);
    }
    else
    {
       MultiFab& new_grav_mf = get_new_data(Gravity_Type);
       new_grav_mf.setVal(0.0);
    }
#endif 

#ifdef DIFFUSION
#ifdef TAU
    MultiFab tau_diff(grids,1,1);
    tau_diff.setVal(0.);
    define_tau(tau_diff,grav_vec_old,time);
#endif
#endif

    const Real prev_time = state[State_Type].prevTime();

    // It's possible for interpolation to create very small negative values for
    //   species so we make sure here that all species are non-negative after this point
    enforce_nonnegative_species(S_old);

#ifdef REACTIONS
    // Make sure to zero these even if do_react == 0.
    MultiFab& ReactMF_old = get_old_data(Reactions_Type);
    MultiFab& ReactMF     = get_new_data(Reactions_Type);
    ReactMF_old.setVal(0.);
    ReactMF.setVal(0.);

    for (FillPatchIterator fpi(*this, S_old, 1, time, State_Type, 0, NUM_STATE);
	   fpi.isValid(); ++fpi)
    {
	  FArrayBox &state = fpi();
#ifdef TAU
          react_first_half_dt(state,ReactMF[fpi],tau_diff[fpi],time,dt);
#else
          react_first_half_dt(state,ReactMF[fpi],time,dt);
#endif
    }
#endif
 
#ifdef GRAVITY
    MultiFab comp_minus_level_phi(grids,1,0,Fab_allocate);
    PArray<MultiFab> comp_minus_level_grad_phi(BL_SPACEDIM,PArrayManage);
    for (int n=0; n<BL_SPACEDIM; ++n)
    {
        comp_minus_level_grad_phi.clear(n);
        comp_minus_level_grad_phi.set(n,new MultiFab(BoxArray(grids).surroundingNodes(n),1,0));
    }

    // Do level solve at beginning of time step in order to compute the
    //   difference between the multilevel and the single level solutions.
    if (do_grav && gravity->get_gravity_type() == "PoissonGrav")
    {
        if (gravity->NoComposite() != 1 && level < parent->finestLevel()) {
           gravity->create_comp_minus_level_grad_phi(level,
                           comp_minus_level_phi,comp_minus_level_grad_phi);
        } else {
           if (verbose && ParallelDescriptor::IOProcessor()) {
              std::cout << " " << '\n';
              std::cout << "... old-time level solve at level " << level << '\n';
           }
           int fill_interior = 0;
           gravity->solve_for_old_phi(level,phi_old,gravity->get_grad_phi_prev(level),fill_interior);
        }
    }
#endif

    Real cur_time = state[State_Type].curTime();
        
    // Copy old data into new data.
    MultiFab::Copy(S_new, S_old, 0, 0, NUM_STATE, S_new.nGrow());
        
    if (add_ext_src) {
           MultiFab ext_src_old(grids,NUM_STATE,0,Fab_allocate);
           getOldSource(prev_time,dt,ext_src_old);
           ext_src_old.mult(dt);
           MultiFab::Add(S_new,ext_src_old,0,0,NUM_STATE,0);

           // Must compute new temperature in case it is needed in the source term evaluation
           computeTemp(S_new);

           // Compute source at new time
           MultiFab ext_src_new(grids,NUM_STATE,0,Fab_allocate);
           getNewSource(prev_time,cur_time,dt,ext_src_new);

           ext_src_old.mult(-0.5);
           ext_src_new.mult( 0.5*dt);

           // Subtract off half of the old source term, and add half of the new.
           MultiFab::Add(S_new,ext_src_old,0,0,S_new.nComp(),0);
           MultiFab::Add(S_new,ext_src_new,0,0,S_new.nComp(),0);
    }

    computeTemp(S_new);

    // ******************
    //  Note: If add_ext_src changed the density, 
    //        we need to update gravity and the associated source terms and
    //        state variables here
    // ******************

#ifdef DIFFUSION
#ifdef TAU
        full_diffusion_update(S_new,prev_time,cur_time,dt,tau_diff);
#else
        full_diffusion_update(S_new,prev_time,cur_time,dt);
#endif
#endif
        
#ifdef REACTIONS
#ifdef TAU
       react_second_half_dt(S_new,tau_diff,cur_time,dt,1);
#else
       react_second_half_dt(S_new,cur_time,dt,1);
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
