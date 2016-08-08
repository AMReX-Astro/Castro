#include <winstd.H>

#include "Castro.H"
#include "Castro_F.H"

using std::string;

#ifndef SDC

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

    MultiFab& state = Sborder;

    const int ng = state.nGrow();

    const iMultiFab& interior_mask = build_interior_boundary_mask(ng);

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "\n" << "... Entering burner and doing half-timestep of burning." << "\n";

    react_state(state, reactions, interior_mask, time, dt, ng);

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... Leaving burner after completing half-timestep of burning." << "\n";

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

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "\n" << "... Entering burner and doing half-timestep of burning." << "\n";

    react_state(state, reactions, interior_mask, time, dt, ng);

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... Leaving burner after completing half-timestep of burning." << "\n";

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
			 BL_TO_FORTRAN_3D(tau_diff[mfi]),
#endif
			 BL_TO_FORTRAN_3D(mask[mfi]),
			 time, dt_react);

	}

	if (verbose) {

	  Real e_added = r.sum(NumSpec + 1);

	  if (ParallelDescriptor::IOProcessor() && e_added != 0.0)
	    std::cout << "... (rho e) added from burning: " << e_added << std::endl;

	}

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

#else

void
Castro::react_state(Real time, Real dt)
{
    BL_PROFILE("Castro::react_state()");

    const Real strt_time = ParallelDescriptor::second();

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "\n" << "... Entering burner and doing full timestep of burning." << "\n";

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    // Build the burning mask, in case the state has ghost zones.

    const int ng = S_new.nGrow();
    const iMultiFab& interior_mask = build_interior_boundary_mask(ng);

    // Create a MultiFab with all of the non-reacting source terms.

    MultiFab A_src(grids, NUM_STATE, ng, Fab_allocate);
    sum_of_sources(A_src);

    MultiFab& reactions = get_old_data(Reactions_Type);

    reactions.setVal(0.0);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
    {

	const Box& bx = mfi.growntilebox(ng);

	FArrayBox& uold    = S_old[mfi];
	FArrayBox& unew    = S_new[mfi];
	FArrayBox& a       = A_src[mfi];
	FArrayBox& r       = reactions[mfi];
	const IArrayBox& m = interior_mask[mfi];

	ca_react_state(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		       uold.dataPtr(), ARLIM_3D(uold.loVect()), ARLIM_3D(uold.hiVect()),
		       unew.dataPtr(), ARLIM_3D(unew.loVect()), ARLIM_3D(unew.hiVect()),
		       a.dataPtr(), ARLIM_3D(a.loVect()), ARLIM_3D(a.hiVect()),
		       r.dataPtr(), ARLIM_3D(r.loVect()), ARLIM_3D(r.hiVect()),
		       m.dataPtr(), ARLIM_3D(m.loVect()), ARLIM_3D(m.hiVect()),
		       time, dt);

    }

    if (ng > 0)
        BoxLib::fill_boundary(S_new, geom);

    if (verbose) {

        Real e_added = reactions.sum(NumSpec + 1);

	if (ParallelDescriptor::IOProcessor() && e_added != 0.0)
	    std::cout << "... (rho e) added from burning: " << e_added << std::endl;

	if (ParallelDescriptor::IOProcessor())
	    std::cout << "... Leaving burner after completing full timestep of burning." << "\n";

        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time, IOProc);

	if (ParallelDescriptor::IOProcessor())
	  std::cout << "Castro::react_state() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
	});
#endif

    }

}

#endif
