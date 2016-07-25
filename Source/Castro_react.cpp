#include <winstd.H>

#include "Castro.H"
#include "Castro_F.H"

using std::string;


void
Castro::strang_react_first_half(Real time, Real dt)
{

    // Reactions are expensive and we would usually rather do a communication
    // step than burn on the ghost zones. So what we will do here is create a mask
    // that indicates that we want to turn on the valid interior zones but NOT
    // on the ghost zones that are interior to the level. However, we DO want to
    // burn on the ghost zones that are on the coarse-fine interfaces, since that
    // is going to be more accurate than interpolating from coarse zones. So we will
    // not mask out those zones, and the subsequent FillBoundary call will not
    // interfere with it.

    MultiFab& reactions = get_old_data(Reactions_Type);

    MultiFab& state = (*Sborder);

    const int ng = state.nGrow();

    const iMultiFab& interior_mask = build_interior_boundary_mask(ng);

    react_state(state, reactions, interior_mask, time, dt, ng);

    reset_internal_energy(state);

    BoxLib::fill_boundary(state, geom);

}



void
Castro::strang_react_second_half(Real time, Real dt)
{

    MultiFab& reactions = get_new_data(Reactions_Type);

    MultiFab& state = get_new_data(State_Type);

    const int ng = state.nGrow();

    const iMultiFab& interior_mask = build_interior_boundary_mask(ng);

    react_state(state, reactions, interior_mask, time, dt, ng);

    reset_internal_energy(state);

    BoxLib::fill_boundary(state, geom);

}



void
Castro::react_state(MultiFab& s, MultiFab& r, const iMultiFab& mask, Real time, Real dt_react, int ngrow)
{

    BL_PROFILE("Castro::react_state()");

    const Real strt_time = ParallelDescriptor::second();

    r.setVal(0.0);

    if (do_react == 1)
    {

      if (verbose && ParallelDescriptor::IOProcessor())
	std::cout << "\n" << "... Entering burner and doing half-timestep of burning." << "\n";

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(s, true); mfi.isValid(); ++mfi)
	{

	  const Box& bx = mfi.growntilebox(ngrow);

	  // Note that box is *not* necessarily just the valid region!
	  ca_react_state(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			 BL_TO_FORTRAN_3D(s[mfi]),
			 BL_TO_FORTRAN_3D(r[mfi]),
#ifdef TAU
			 BL_TO_FORTRAN_3D((*tau_diff)[mfi]),
#endif
			 BL_TO_FORTRAN_3D(mask[mfi]),
			 time, dt_react);

	}

	if (verbose) {

	  Real e_added = r.sum(NumSpec + 1);

	  if (ParallelDescriptor::IOProcessor() && e_added != 0.0)
	    std::cout << "... (rho e) added from burning: " << e_added << std::endl;

	}

	if (verbose && ParallelDescriptor::IOProcessor())
	  std::cout << "... Leaving burner after completing half-timestep of burning." << "\n";

    }

    if (verbose > 0 && do_react)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

	if (ParallelDescriptor::IOProcessor())
	  std::cout << "Castro::react_state() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
	});
#endif
    }

}
