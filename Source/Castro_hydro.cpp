#include "Castro.H"
#include "Castro_F.H"

void
Castro::hydro_update(Real time, Real dt)
{

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... Entering hydro advance" << std::endl << std::endl;

    const Real cur_time = state[State_Type].curTime();

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
		    (*u_gdnv[i])[mfi].copy(ugdn[i],mfi.nodaltilebox(i));
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
#pragma omp parallel reduction(+:E_added_flux,mass_added_flux) \
		     reduction(+:xmom_added_flux,ymom_added_flux,zmom_added_flux) \
		     reduction(+:mass_lost,xmom_lost,ymom_lost,zmom_lost) \
		     reduction(+:eden_lost,xang_lost,yang_lost,zang_lost)
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

		// Copy the normal velocities from the Riemann solver

		for (int i = 0; i < BL_SPACEDIM ; i++) {
		    (*u_gdnv[i])[mfi].copy(ugdn[i],mfi.nodaltilebox(i));
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

	    Real foo[5] = {E_added_flux, xmom_added_flux, ymom_added_flux, zmom_added_flux, mass_added_flux};

#ifdef BL_LAZY
	    Lazy::QueueReduction( [=] () mutable {
#endif
	    ParallelDescriptor::ReduceRealSum(foo, 5, ParallelDescriptor::IOProcessorNumber());

	    if (ParallelDescriptor::IOProcessor())
	    {
	        E_added_flux    = foo[0];
	        xmom_added_flux = foo[1];
	        ymom_added_flux = foo[2];
	        zmom_added_flux = foo[3];
	        mass_added_flux = foo[4];

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
