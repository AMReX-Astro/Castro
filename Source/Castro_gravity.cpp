#include "Castro.H"
#include "Castro_F.H"

#ifdef SELF_GRAVITY
#include "Gravity.H"

void
Castro::construct_old_gravity(int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle, Real time)
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
	PArray<MultiFab> comp_gphi(BL_SPACEDIM, PArrayManage);

        if (gravity->NoComposite() != 1 && level < parent->finestLevel()) {

	    comp_phi.define(phi_old.boxArray(), phi_old.nComp(), phi_old.nGrow(), Fab_allocate);
	    MultiFab::Copy(comp_phi, phi_old, 0, 0, phi_old.nComp(), phi_old.nGrow());

	    for (int n = 0; n < BL_SPACEDIM; ++n) {
		comp_gphi.set(n, new MultiFab(getEdgeBoxArray(n), 1, 0));
		comp_gphi[n].copy(gravity->get_grad_phi_prev(level)[n], 0, 0, 1);
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
			       gravity->get_grad_phi_prev(level),
			       is_new);

        if (gravity->NoComposite() != 1 && level < parent->finestLevel()) {

	    // Subtract the level solve from the composite solution.

	    gravity->create_comp_minus_level_grad_phi(level,
						      comp_phi,
						      comp_gphi,
						      comp_minus_level_phi,
						      comp_minus_level_grad_phi);

	    // Copy the composite data back. This way the forcing
	    // uses the most accurate data we have.

	    MultiFab::Copy(phi_old, comp_phi, 0, 0, phi_old.nComp(), phi_old.nGrow());

	    for (int n = 0; n < BL_SPACEDIM; ++n)
		gravity->get_grad_phi_prev(level)[n].copy(comp_gphi[n], 0, 0, 1);

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
Castro::construct_new_gravity(int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle, Real time)
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

	if (level < parent->finestLevel() && gravity->NoComposite() != 1)
	    phi_new.minus(comp_minus_level_phi, 0, 1, 0);

	if (verbose && ParallelDescriptor::IOProcessor()) {
	    std::cout << " " << '\n';
	    std::cout << "... new-time level solve at level " << level << '\n';
	}

	int is_new = 1;

	gravity->solve_for_phi(level,
			       phi_new,
			       gravity->get_grad_phi_curr(level),
			       is_new);

	if (level < parent->finestLevel() && gravity->NoComposite() != 1) {

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
	    // If we are doing the sync solve, it should make that faster.

	    phi_new.plus(comp_minus_level_phi, 0, 1, 0);
	    gravity->plus_grad_phi_curr(level, comp_minus_level_grad_phi);

	    // We can clear this memory, we no longer need it.

	    comp_minus_level_phi.clear();
	    AMReX::FillNull(comp_minus_level_grad_phi);

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

}

void Castro::construct_old_gravity_source(Real time, Real dt)
{

    MultiFab& phi_old = get_old_data(PhiGrav_Type);
    MultiFab& grav_old = get_old_data(Gravity_Type);

    old_sources[grav_src]->setVal(0.0);

    if (!do_grav) return;

    // Gravitational source term for the time-level n data.

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
    for (MFIter mfi(Sborder,true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.growntilebox();

	Real mom_added[3] = { 0.0 };

	ca_gsrc(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		ARLIM_3D(domlo), ARLIM_3D(domhi),
		BL_TO_FORTRAN_3D(phi_old[mfi]),
		BL_TO_FORTRAN_3D(grav_old[mfi]),
		BL_TO_FORTRAN_3D(Sborder[mfi]),
		BL_TO_FORTRAN_3D((*old_sources[grav_src])[mfi]),
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
}

void Castro::construct_new_gravity_source(Real time, Real dt)
{
    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    MultiFab& phi_old = get_old_data(PhiGrav_Type);
    MultiFab& phi_new = get_new_data(PhiGrav_Type);

    MultiFab& grav_old = get_old_data(Gravity_Type);
    MultiFab& grav_new = get_new_data(Gravity_Type);

    new_sources[grav_src]->setVal(0.0);

    if (!do_grav) return;

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
			BL_TO_FORTRAN_3D((*new_sources[grav_src])[mfi]),
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
}
#else
// This is the constant gravity version
void Castro::construct_old_gravity_source(Real time, Real dt)
{
    old_sources[grav_src]->setVal(0.0);

    if (!do_grav) return;

    // Gravitational source term for the time-level n data.

    const Real* dx = geom.CellSize();
    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();

    for (MFIter mfi(Sborder,true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.growntilebox();

	Real mom_added[3] = { 0.0 };

	ca_gsrc(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		ARLIM_3D(domlo), ARLIM_3D(domhi),
		BL_TO_FORTRAN_3D(Sborder[mfi]),
		BL_TO_FORTRAN_3D((*old_sources[grav_src])[mfi]),
		BL_TO_FORTRAN_3D(volume[mfi]),
		ZFILL(dx),dt,&time);
    }
}

// This is the constant gravity version
void Castro::construct_new_gravity_source(Real time, Real dt)
{
    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    new_sources[grav_src]->setVal(0.0);

    if (!do_grav) return;

    const Real *dx = geom.CellSize();
    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();

    {
	for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
	{
	    const Box& bx = mfi.tilebox();

	    Real mom_added[3] = { 0.0 };

	    ca_corrgsrc(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			ARLIM_3D(domlo), ARLIM_3D(domhi),
			BL_TO_FORTRAN_3D(S_old[mfi]),
			BL_TO_FORTRAN_3D(S_new[mfi]),
			BL_TO_FORTRAN_3D((*new_sources[grav_src])[mfi]),
			BL_TO_FORTRAN_3D((*fluxes[0])[mfi]),
			BL_TO_FORTRAN_3D((*fluxes[1])[mfi]),
			BL_TO_FORTRAN_3D((*fluxes[2])[mfi]),
			ZFILL(dx),dt,&time,
			BL_TO_FORTRAN_3D(volume[mfi]));
	}
    }
}
#endif
