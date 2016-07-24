#include "Castro.H"
#include "Castro_F.H"
#include "Gravity.H"

void
Castro::construct_old_gravity(int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle, Real time)
{

    // Old and new gravitational potential.

    MultiFab& phi_old = get_old_data(PhiGrav_Type);
    MultiFab& grav_old = get_old_data(Gravity_Type);

    int finest_level = parent->finestLevel();

    if (do_grav) {

       // Swap the old and new data at this level and all finer levels,
       // and zero out the flux registers. 

       if (level == 0 || amr_iteration > 1) {

	 if (sub_iteration < 2)
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
       grav_old.setVal(0.0);
       phi_old.setVal(0.0);
    }

    // Do level solve at beginning of time step in order to compute the
    //   difference between the multilevel and the single level solutions.
    if (do_grav && gravity->get_gravity_type() == "PoissonGrav")
    {
        if (gravity->NoComposite() != 1 && level < parent->finestLevel()) {

	    comp_minus_level_phi = new MultiFab;
	    comp_minus_level_phi->define(grids,1,0,Fab_allocate);

	    for (int n = 0; n < BL_SPACEDIM; ++n) {
	        comp_minus_level_grad_phi[n] = new MultiFab;
	        comp_minus_level_grad_phi[n]->define(getEdgeBoxArray(n),1,0,Fab_allocate);
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

}

void
Castro::construct_new_gravity(int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle, Real time)
{

    MultiFab& phi_new = get_new_data(PhiGrav_Type);
    MultiFab& grav_new = get_new_data(Gravity_Type);

    if (do_grav)
    {
	if (gravity->get_gravity_type() == "PoissonGrav")
	  {
            if (verbose && ParallelDescriptor::IOProcessor()) {
	      std::cout << " " << '\n';
	      std::cout << "... new-time level solve at level " << level << '\n';
            }
	    
            // Here we use the "old" phi from the current time step as a guess for this solve
	    MultiFab& phi_old = get_old_data(PhiGrav_Type);
	    MultiFab::Copy(phi_new,phi_old,0,0,1,phi_new.nGrow());
            if ( level < parent->finestLevel() && (gravity->NoComposite() != 1) ) {
		phi_new.minus(*comp_minus_level_phi, 0, 1, 0);
	    }
            int is_new = 1;
            gravity->solve_for_phi(level,
				   phi_new,
				   gravity->get_grad_phi_curr(level),
				   is_new);

            if (gravity->test_results_of_solves() == 1) {
                if (verbose && ParallelDescriptor::IOProcessor()) {
		    std::cout << " " << '\n';
		    std::cout << "... testing grad_phi_curr before adding comp_minus_level_grad_phi " << '\n';
		}
		gravity->test_level_grad_phi_curr(level);
	    }

            if ( level < parent->finestLevel() && (gravity->NoComposite() != 1) ) {
	      phi_new.plus(*comp_minus_level_phi, 0, 1, 0);
	      gravity->plus_grad_phi_curr(level,comp_minus_level_grad_phi);

	      delete comp_minus_level_phi;
	      for (int n = 0; n < BL_SPACEDIM; ++n)
		delete comp_minus_level_grad_phi[n];
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

            if (do_reflux)  gravity->add_to_fluxes(level,amr_iteration,amr_ncycle);
	  }
	else {
	    phi_new.setVal(0.0);  // so that plotfiles do not contain nans
	}

	// Now do corrector part of source term update
	gravity->get_new_grav_vector(level,grav_new,time);
    }

}


void Castro::construct_old_gravity_source(Real time, Real dt)
{
    MultiFab& phi_old = get_old_data(PhiGrav_Type);
    MultiFab& grav_old = get_old_data(Gravity_Type);

    // Gravitational source term for the time-level n data.

    if (do_grav)
    {

        Real E_added    = 0.;
	Real xmom_added = 0.;
	Real ymom_added = 0.;
	Real zmom_added = 0.;

	const Real* dx = geom.CellSize();
	const int* domlo = geom.Domain().loVect();
	const int* domhi = geom.Domain().hiVect();

#ifdef _OPENMP
#pragma omp parallel reduction(+:E_added,xmom_added,ymom_added,zmom_added)
#endif
	for (MFIter mfi((*Sborder),true); mfi.isValid(); ++mfi)
	{
	    const Box& bx = mfi.tilebox();

	    Real mom_added[3] = { 0.0 };

	    ca_gsrc(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		    ARLIM_3D(domlo), ARLIM_3D(domhi),
		    BL_TO_FORTRAN_3D(phi_old[mfi]),
		    BL_TO_FORTRAN_3D(grav_old[mfi]),
		    BL_TO_FORTRAN_3D((*Sborder)[mfi]),
		    BL_TO_FORTRAN_3D(old_sources[grav_src][mfi]),
		    BL_TO_FORTRAN_3D(volume[mfi]),
		    ZFILL(dx),dt,&time,
		    E_added,mom_added);

	    xmom_added += mom_added[0];
	    ymom_added += mom_added[1];
	    zmom_added += mom_added[2];
	}

        if (print_energy_diagnostics)
        {
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

		std::cout << "(rho E) added from grav. source  terms          : " << E_added << std::endl;
		std::cout << "xmom added from grav. source terms              : " << xmom_added << std::endl;
#if (BL_SPACEDIM >= 2)
		std::cout << "ymom added from grav. source terms              : " << ymom_added << std::endl;
#endif
#if (BL_SPACEDIM == 3)
		std::cout << "zmom added from grav. source terms              : " << zmom_added << std::endl;
#endif
	    }
#ifdef BL_LAZY
	    });
#endif
        }

	add_force_to_sources(grav_old, *sources_for_hydro, (*Sborder));

    }

}



void Castro::construct_new_gravity_source(Real time, Real dt)
{
    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    MultiFab& phi_old = get_old_data(PhiGrav_Type);
    MultiFab& phi_new = get_new_data(PhiGrav_Type);

    MultiFab& grav_old = get_old_data(Gravity_Type);
    MultiFab& grav_new = get_new_data(Gravity_Type);

    if (do_grav) {

        Real E_added    = 0.;
	Real xmom_added = 0.;
	Real ymom_added = 0.;
	Real zmom_added = 0.;

	const Real *dx = geom.CellSize();
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
			    BL_TO_FORTRAN_3D(new_sources[grav_src][mfi]),
			    BL_TO_FORTRAN_3D((*fluxes[0])[mfi]),
			    BL_TO_FORTRAN_3D((*fluxes[1])[mfi]),
			    BL_TO_FORTRAN_3D((*fluxes[2])[mfi]),
			    ZFILL(dx),dt,&time,
			    BL_TO_FORTRAN_3D(volume[mfi]),
			    E_added, mom_added);

		xmom_added += mom_added[0];
		ymom_added += mom_added[1];
		zmom_added += mom_added[2];
	    }
	}

        if (print_energy_diagnostics)
        {
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

		std::cout << "(rho E) added from grav. corr.  terms          : " << E_added << std::endl;
		std::cout << "xmom added from grav. corr. terms              : " << xmom_added << std::endl;
#if (BL_SPACEDIM >= 2)
		std::cout << "ymom added from grav. corr. terms              : " << ymom_added << std::endl;
#endif
#if (BL_SPACEDIM == 3)
		std::cout << "zmom added from grav. corr. terms              : " << zmom_added << std::endl;
#endif
	    }
#ifdef BL_LAZY
	    });
#endif
        }

	// Add this to the source term array if we're using the source term predictor.
	// If not, don't bother because sources isn't actually used in the update after this point.

	if (source_term_predictor == 1)
	  add_force_to_sources(grav_new, *sources_for_hydro, S_new);

    }
}
