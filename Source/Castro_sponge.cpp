#ifdef SPONGE
#include "Castro.H"
#include "Castro_F.H"

void
Castro::construct_old_sponge_source(Real time, Real dt)
{
    int ng = Sborder.nGrow();

    old_sources[sponge_src].setVal(0.0);

    if (!time_center_sponge || !do_sponge) return;

    update_sponge_params(&time);

    Real E_added    = 0.;
    Real xmom_added = 0.;
    Real ymom_added = 0.;
    Real zmom_added = 0.;

    const Real *dx = geom.CellSize();

#ifdef _OPENMP
#pragma omp parallel reduction(+:E_added,xmom_added,ymom_added,zmom_added)
#endif
    for (MFIter mfi(Sborder,true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();

	Real mom_added[3] = { 0.0 };

	ca_sponge(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		  BL_TO_FORTRAN_3D(Sborder[mfi]),
		  BL_TO_FORTRAN_3D(old_sources[sponge_src][mfi]),
		  BL_TO_FORTRAN_3D(volume[mfi]),
		  ZFILL(dx), dt, &time,
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
		    xmom_added = foo[1];
		    ymom_added = foo[2];
		    zmom_added = foo[3];

		    std::cout << "(rho E) added from sponge                      : " << E_added << std::endl;
		    std::cout << "xmom added from sponge                         : " << xmom_added << std::endl;
		    std::cout << "ymom added from sponge                         : " << ymom_added << std::endl;
		    std::cout << "zmom added from sponge                         : " << zmom_added << std::endl;
		}
#ifdef BL_LAZY
	    });
#endif
    }

}

void
Castro::construct_new_sponge_source(Real time, Real dt)
{
    MultiFab& S_new = get_new_data(State_Type);

    int ng = 0;

    new_sources[sponge_src].setVal(0.0);

    if (!do_sponge) return;

    update_sponge_params(&time);

    Real E_added    = 0.;
    Real xmom_added = 0.;
    Real ymom_added = 0.;
    Real zmom_added = 0.;

    const Real *dx = geom.CellSize();

#ifdef _OPENMP
#pragma omp parallel reduction(+:E_added,xmom_added,ymom_added,zmom_added)
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();

	Real mom_added[3] = { 0.0 };

	ca_sponge(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		  BL_TO_FORTRAN_3D(S_new[mfi]),
		  BL_TO_FORTRAN_3D(new_sources[sponge_src][mfi]),
		  BL_TO_FORTRAN_3D(volume[mfi]),
		  ZFILL(dx), dt, &time,
		  E_added,mom_added);

	xmom_added += mom_added[0];
	ymom_added += mom_added[1];
	zmom_added += mom_added[2];

    }

    // Time center the source term.

    if (time_center_sponge) {
	new_sources[sponge_src].mult(0.5);

	MultiFab::Saxpy(new_sources[sponge_src],-0.5,old_sources[sponge_src],0,0,NUM_STATE,ng);
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
		    xmom_added = foo[1];
		    ymom_added = foo[2];
		    zmom_added = foo[3];

		    std::cout << "(rho E) added from sponge                      : " << E_added << std::endl;
		    std::cout << "xmom added from sponge                         : " << xmom_added << std::endl;
		    std::cout << "ymom added from sponge                         : " << ymom_added << std::endl;
		    std::cout << "zmom added from sponge                         : " << zmom_added << std::endl;
		}
#ifdef BL_LAZY
	    });
#endif
    }
}
#endif
