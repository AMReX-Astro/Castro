#include "Castro.H"
#include "Castro_F.H"

#ifdef SELF_GRAVITY
#include "Gravity.H"

using namespace amrex;

void
Castro::construct_old_gravity(int amr_iteration, int amr_ncycle, Real time)
{
    MultiFab& grav_old = get_old_data(Gravity_Type);
    MultiFab& phi_old = get_old_data(PhiGrav_Type);

    // Always set phi to zero initially since some gravity modes
    // don't use it and we want to have valid data.

    if (gravity->get_gravity_type() != "PoissonGrav")
	phi_old.setVal(0.0);

    if (!do_grav) {

	grav_old.setVal(0.0);

	return;

    }

    // Do level solve at beginning of time step in order to compute the
    // difference between the multilevel and the single level solutions.

    if (gravity->get_gravity_type() == "PoissonGrav")
    {

	// Create a copy of the current (composite) data on this level.

	MultiFab comp_phi;
	Vector<std::unique_ptr<MultiFab> > comp_gphi(BL_SPACEDIM);

        if (gravity->NoComposite() != 1 && gravity->DoCompositeCorrection() && level < parent->finestLevel()) {

	    comp_phi.define(phi_old.boxArray(), phi_old.DistributionMap(), phi_old.nComp(), phi_old.nGrow());
	    MultiFab::Copy(comp_phi, phi_old, 0, 0, phi_old.nComp(), phi_old.nGrow());

	    for (int n = 0; n < BL_SPACEDIM; ++n) {
		comp_gphi[n].reset(new MultiFab(getEdgeBoxArray(n), dmap, 1, 0));
		comp_gphi[n]->copy(*gravity->get_grad_phi_prev(level)[n], 0, 0, 1);
	    }

	}

	if (verbose && ParallelDescriptor::IOProcessor()) {
	    std::cout << " " << '\n';
	    std::cout << "... old-time level solve at level " << level << '\n';
	}

	int is_new = 0;

	// If we are doing composite solves, then this is a placeholder solve
	// to get the difference between the composite and level solutions. If
	// we are only doing level solves, then this is the main result.

	gravity->solve_for_phi(level,
			       phi_old,
			       amrex::GetVecOfPtrs(gravity->get_grad_phi_prev(level)),
			       is_new);

        if (gravity->NoComposite() != 1 && gravity->DoCompositeCorrection() && level < parent->finestLevel()) {

	    // Subtract the level solve from the composite solution.

	    gravity->create_comp_minus_level_grad_phi(level,
						      comp_phi,
						      amrex::GetVecOfPtrs(comp_gphi),
						      comp_minus_level_phi,
						      comp_minus_level_grad_phi);

	    // Copy the composite data back. This way the forcing
	    // uses the most accurate data we have.

	    MultiFab::Copy(phi_old, comp_phi, 0, 0, phi_old.nComp(), phi_old.nGrow());

	    for (int n = 0; n < BL_SPACEDIM; ++n)
		gravity->get_grad_phi_prev(level)[n]->copy(*comp_gphi[n], 0, 0, 1);

        }

	if (gravity->test_results_of_solves() == 1) {

	    if (verbose && ParallelDescriptor::IOProcessor()) {
		std::cout << " " << '\n';
		std::cout << "... testing grad_phi_curr after doing single level solve " << '\n';
	    }

	    gravity->test_level_grad_phi_prev(level);

	}
 
    }

    // Define the old gravity vector.

    gravity->get_old_grav_vector(level, grav_old, time);

}

void
Castro::construct_new_gravity(int amr_iteration, int amr_ncycle, Real time)
{
    MultiFab& grav_new = get_new_data(Gravity_Type);
    MultiFab& phi_new = get_new_data(PhiGrav_Type);

    // Always set phi to zero initially since some gravity modes
    // don't use it and we want to have valid data.

    if (gravity->get_gravity_type() != "PoissonGrav")
	phi_new.setVal(0.0);

    if (!do_grav) {

	grav_new.setVal(0.0);

	return;

    }

    // If we're doing Poisson gravity, do the new-time level solve here.

    if (gravity->get_gravity_type() == "PoissonGrav")
    {

	// Use the "old" phi from the current time step as a guess for this solve.

	MultiFab& phi_old = get_old_data(PhiGrav_Type);

	MultiFab::Copy(phi_new, phi_old, 0, 0, 1, phi_new.nGrow());

	// Subtract off the (composite - level) contribution for the purposes
	// of the level solve. We'll add it back later.

	if (gravity->NoComposite() != 1 && gravity->DoCompositeCorrection() && level < parent->finestLevel())
	    phi_new.minus(comp_minus_level_phi, 0, 1, 0);

	if (verbose && ParallelDescriptor::IOProcessor()) {
	    std::cout << " " << '\n';
	    std::cout << "... new-time level solve at level " << level << '\n';
	}

	int is_new = 1;

	gravity->solve_for_phi(level,
			       phi_new,
			       amrex::GetVecOfPtrs(gravity->get_grad_phi_curr(level)),
			       is_new);

	if (gravity->NoComposite() != 1 && gravity->DoCompositeCorrection() == 1 && level < parent->finestLevel()) {

	    if (gravity->test_results_of_solves() == 1) {

		if (verbose && ParallelDescriptor::IOProcessor()) {
		    std::cout << " " << '\n';
		    std::cout << "... testing grad_phi_curr before adding comp_minus_level_grad_phi " << '\n';
		}

		gravity->test_level_grad_phi_curr(level);

	    }

	    // Add back the (composite - level) contribution. This ensures that
	    // if we are not doing a sync solve, then we still get the difference
	    // between the composite and level solves added to the force we
	    // calculate, so it is slightly more accurate than it would have been.

	    phi_new.plus(comp_minus_level_phi, 0, 1, 0);
	    for (int n = 0; n < BL_SPACEDIM; ++n)
		gravity->get_grad_phi_curr(level)[n]->plus(*comp_minus_level_grad_phi[n], 0, 1, 0);

	    if (gravity->test_results_of_solves() == 1) {

		if (verbose && ParallelDescriptor::IOProcessor()) {
		    std::cout << " " << '\n';
		    std::cout << "... testing grad_phi_curr after adding comp_minus_level_grad_phi " << '\n';
		}

		gravity->test_level_grad_phi_curr(level);

	    }

	}

    }

    // Define new gravity vector.

    gravity->get_new_grav_vector(level, grav_new, time);

    if (gravity->get_gravity_type() == "PoissonGrav") {

	if (gravity->NoComposite() != 1 && gravity->DoCompositeCorrection() == 1 && level < parent->finestLevel()) {

	    // Now that we have calculated the force, if we are going to do a sync
	    // solve then subtract off the (composite - level) contribution, as it
	    // interferes with the sync solve.

	    if (gravity->NoSync() == 0) {

		phi_new.minus(comp_minus_level_phi, 0, 1, 0);

		for (int n = 0; n < BL_SPACEDIM; ++n)
		    gravity->get_grad_phi_curr(level)[n]->minus(*comp_minus_level_grad_phi[n], 0, 1, 0);

	    }

	    // In any event we can now clear this memory, as we no longer need it.

	    comp_minus_level_phi.clear();
	    comp_minus_level_grad_phi.clear();

	}

    }

}
#endif

void Castro::construct_old_gravity_source(MultiFab& source, MultiFab& state, Real time, Real dt)
{

#ifdef SELF_GRAVITY
    const MultiFab& phi_old = get_old_data(PhiGrav_Type);
    const MultiFab& grav_old = get_old_data(Gravity_Type);
#endif

    if (!do_grav) return;

    // Gravitational source term for the time-level n data.

    const Real* dx = geom.CellSize();
    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state, true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();

	ca_gsrc(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		ARLIM_3D(domlo), ARLIM_3D(domhi),
		BL_TO_FORTRAN_3D(state[mfi]),
#ifdef SELF_GRAVITY
		BL_TO_FORTRAN_3D(phi_old[mfi]),
		BL_TO_FORTRAN_3D(grav_old[mfi]),
#endif
		BL_TO_FORTRAN_3D(source[mfi]),
		ZFILL(dx),dt,&time);

    }

}

void Castro::construct_new_gravity_source(MultiFab& source, MultiFab& state_old, MultiFab& state_new, Real time, Real dt)
{

#ifdef SELF_GRAVITY
    MultiFab& phi_old = get_old_data(PhiGrav_Type);
    MultiFab& phi_new = get_new_data(PhiGrav_Type);

    MultiFab& grav_old = get_old_data(Gravity_Type);
    MultiFab& grav_new = get_new_data(Gravity_Type);
#endif

    if (!do_grav) return;

    const Real *dx = geom.CellSize();
    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	for (MFIter mfi(state_new, true); mfi.isValid(); ++mfi)
	{
	    const Box& bx = mfi.tilebox();

	    ca_corrgsrc(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			ARLIM_3D(domlo), ARLIM_3D(domhi),
			BL_TO_FORTRAN_3D(state_old[mfi]),
			BL_TO_FORTRAN_3D(state_new[mfi]),
#ifdef SELF_GRAVITY
			BL_TO_FORTRAN_3D(phi_old[mfi]),
			BL_TO_FORTRAN_3D(phi_new[mfi]),
			BL_TO_FORTRAN_3D(grav_old[mfi]),
			BL_TO_FORTRAN_3D(grav_new[mfi]),
#endif
			BL_TO_FORTRAN_3D(volume[mfi]),
			BL_TO_FORTRAN_3D((*mass_fluxes[0])[mfi]),
			BL_TO_FORTRAN_3D((*mass_fluxes[1])[mfi]),
			BL_TO_FORTRAN_3D((*mass_fluxes[2])[mfi]),
			BL_TO_FORTRAN_3D(source[mfi]),
			ZFILL(dx),dt,&time);

	}
    }

}
