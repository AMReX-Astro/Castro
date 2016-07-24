#include <winstd.H>

#include "Castro.H"
#include "Castro_F.H"

void Castro::construct_old_rotation(int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle, Real time)
{

    MultiFab& phirot_old = get_old_data(PhiRot_Type);
    MultiFab& rot_old = get_old_data(Rotation_Type);

    // Fill the old rotation data.

    if (do_rotation) {

      fill_rotation_field(phirot_old, rot_old, (*Sborder), time);

    } else {

      phirot_old.setVal(0.);
      rot_old.setVal(0.);

    }

}



void Castro::construct_new_rotation(int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle, Real time)
{

    MultiFab& phirot_new = get_new_data(PhiRot_Type);
    MultiFab& rot_new = get_new_data(Rotation_Type);

    MultiFab& S_new = get_new_data(State_Type);

    // Fill the old rotation data.

    if (do_rotation) {

      fill_rotation_field(phirot_new, rot_new, S_new, time);

    } else {

      phirot_new.setVal(0.);
      rot_new.setVal(0.);

    }

}



void Castro::construct_old_rotation_source(Real time, Real dt)
{
    MultiFab& phirot_old = get_old_data(PhiRot_Type);
    MultiFab& rot_old = get_old_data(Rotation_Type);

    if (do_rotation) {

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
	for (MFIter mfi((*Sborder),true); mfi.isValid(); ++mfi)
	{
	    const Box& bx = mfi.tilebox();

	    Real mom_added[3] = { 0.0 };

	    ca_rsrc(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		    ARLIM_3D(domlo), ARLIM_3D(domhi),
		    BL_TO_FORTRAN_3D(phirot_old[mfi]),
		    BL_TO_FORTRAN_3D(rot_old[mfi]),
		    BL_TO_FORTRAN_3D((*Sborder)[mfi]),
		    BL_TO_FORTRAN_3D(old_sources[rot_src][mfi]),
		    BL_TO_FORTRAN_3D(volume[mfi]),
		    ZFILL(dx),dt,&time,
		    E_added,mom_added);

	    xmom_added += mom_added[0];
	    ymom_added += mom_added[1];
	    zmom_added += mom_added[2];

	}

        if (print_energy_diagnostics)
        {
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

		std::cout << "(rho E) added from rot. source  terms          : " << E_added << std::endl;
		std::cout << "xmom added from rot. source terms              : " << xmom_added << std::endl;
		std::cout << "ymom added from rot. source terms              : " << ymom_added << std::endl;
		std::cout << "zmom added from rot. source terms              : " << zmom_added << std::endl;
	    }
#ifdef BL_LAZY
	    });
#endif
        }

	add_force_to_sources(rot_old, *sources_for_hydro, (*Sborder));

    }

}



void Castro::construct_new_rotation_source(MultiFab fluxes[], Real time, Real dt)
{
    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    MultiFab& phirot_old = get_old_data(PhiRot_Type);
    MultiFab& rot_old = get_old_data(Rotation_Type);

    MultiFab& phirot_new = get_new_data(PhiRot_Type);
    MultiFab& rot_new = get_new_data(Rotation_Type);

    if (do_rotation) {

	// Now do corrector part of rotation source term update

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

		ca_corrrsrc(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			    ARLIM_3D(domlo), ARLIM_3D(domhi),
			    BL_TO_FORTRAN_3D(phirot_old[mfi]),
			    BL_TO_FORTRAN_3D(phirot_new[mfi]),
			    BL_TO_FORTRAN_3D(rot_old[mfi]),
			    BL_TO_FORTRAN_3D(rot_new[mfi]),
			    BL_TO_FORTRAN_3D(S_old[mfi]),
			    BL_TO_FORTRAN_3D(S_new[mfi]),
			    BL_TO_FORTRAN_3D(new_sources[rot_src][mfi]),
			    BL_TO_FORTRAN_3D(fluxes[0][mfi]),
			    BL_TO_FORTRAN_3D(fluxes[1][mfi]),
			    BL_TO_FORTRAN_3D(fluxes[2][mfi]),
			    ZFILL(dx),dt,&time,
			    BL_TO_FORTRAN_3D(volume[mfi]),
			    E_added,mom_added);

		xmom_added += mom_added[0];
		ymom_added += mom_added[1];
		zmom_added += mom_added[2];
	    }
	}

        if (print_energy_diagnostics)
        {
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

		std::cout << "(rho E) added from rot. corr.  terms          : " << E_added << std::endl;
		std::cout << "xmom added from rot. corr. terms              : " << xmom_added << std::endl;
		std::cout << "ymom added from rot. corr. terms              : " << ymom_added << std::endl;
		std::cout << "zmom added from rot. corr. terms              : " << zmom_added << std::endl;
	    }
#ifdef BL_LAZY
	    });
#endif
        }

	// Add this to the source term array if we're using the source term predictor.
	// If not, don't bother because sources isn't actually used in the update after this point.

	if (source_term_predictor == 1)
	  add_force_to_sources(rot_new, *sources_for_hydro, S_new);

    }

}



void Castro::fill_rotation_field(MultiFab& phi, MultiFab& rot, MultiFab& state, Real time)
{
    const Real* dx = geom.CellSize();

    phi.setVal(0.0);    

    int ng = phi.nGrow();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(phi, true); mfi.isValid(); ++mfi)
    {

      const Box& bx = mfi.growntilebox(ng);

      ca_fill_rotational_potential(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()), 
				   BL_TO_FORTRAN_3D(phi[mfi]),
				   ZFILL(dx),time);

    }

    rot.setVal(0.0);

    ng = state.nGrow();

    if (ng > rot.nGrow())
      BoxLib::Error("State MF has more ghost cells than rotation MF.");
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state, true); mfi.isValid(); ++mfi)
    {

      const Box& bx = mfi.growntilebox(ng);

      ca_fill_rotational_acceleration(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()), 
				      BL_TO_FORTRAN_3D(rot[mfi]),
				      BL_TO_FORTRAN_3D(state[mfi]),
				      ZFILL(dx),time);

    }

}


