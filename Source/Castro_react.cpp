#include <winstd.H>

#include "Castro.H"
#include "Castro_F.H"

using std::string;

#ifdef REACTIONS
void
#ifdef TAU
Castro::react_half_dt(MultiFab& state, MultiFab& tau_diff, Real time, Real dt) 
#else
Castro::react_half_dt(MultiFab& state, Real time, Real dt) 
#endif
{
    BL_PROFILE("Castro::react_second_half_dt()");

    // Note that here we only react on the valid region of the MultiFab (S_new has no ghost cells).
    if (do_react == 1) 
    {

      if (verbose && ParallelDescriptor::IOProcessor())
	std::cout << "\n" << "... Entering burner and doing half-timestep of burning." << "\n";

        MultiFab& ReactMF = get_new_data(Reactions_Type);
	ReactMF.setVal(0.);

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(state, true); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
#ifdef TAU
	  reactState(state[mfi], state[mfi], ReactMF[mfi], tau_diff[mfi], bx, time, 0.5*dt);
#else
	  reactState(state[mfi], state[mfi], ReactMF[mfi], bx, time, 0.5*dt);
#endif
	}

        reset_internal_energy(state);

	if (verbose && ParallelDescriptor::IOProcessor())
	  std::cout << "... Leaving burner after completing half-timestep of burning." << "\n\n";

    }

}

void
Castro::reactState(FArrayBox&        Snew,
                   FArrayBox&        Sold,
                   FArrayBox&        ReactionTerms,
#ifdef TAU
                   FArrayBox&        tau,
#endif
                   const Box&        box,
                   Real              time,
                   Real              dt_react)
{
    BL_PROFILE("Castro::reactState()");

    const Real strt_time = ParallelDescriptor::second();

    const Real cur_time = state[State_Type].curTime();

    // Note that box is *not* necessarily just the valid region!
    BL_FORT_PROC_CALL(CA_REACT_STATE,ca_react_state)
                    (ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()), 
 	             BL_TO_FORTRAN_3D(Sold),
                     BL_TO_FORTRAN_3D(Snew),
                     BL_TO_FORTRAN_3D(ReactionTerms),
#ifdef TAU
                     BL_TO_FORTRAN_3d(tau),
#endif
                     time,dt_react);

    if (verbose > 1)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);
	if (ParallelDescriptor::IOProcessor()) 
	    std::cout << "reactState time = " << run_time << '\n';
#ifdef BL_LAZY
	});
#endif
    }

}
#endif
