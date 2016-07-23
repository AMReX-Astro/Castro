#include "Castro.H"
#include "Gravity.H"

void
Castro::construct_old_gravity(int amr_iteration, int amr_ncycle,
			      int sub_iteration, int sub_ncycle,
			      Real time)
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
Castro::construct_new_gravity(int amr_iteration, int amr_ncycle,
			      int sub_iteration, int sub_ncycle,
			      Real time)
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
