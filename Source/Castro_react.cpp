#include <winstd.H>

#include "Castro.H"
#include "Castro_F.H"

using std::string;

#ifdef REACTIONS
void
#ifdef TAU
Castro::react_first_half_dt(MultiFab& S_old, MultiFab& tau_diff, Real time, Real dt) 
#else
Castro::react_first_half_dt(MultiFab& S_old, Real time, Real dt) 
#endif
{
    BL_PROFILE("Castro::react_first_half_dt()");

    if (do_react == 1)
    {

      if (verbose && ParallelDescriptor::IOProcessor())
	std::cout << "\n" << "... Entering burner and doing first half-timestep of burning." << "\n";

      MultiFab& ReactMF_old = get_old_data(Reactions_Type);
      MultiFab& ReactMF     = get_new_data(Reactions_Type);
      ReactMF_old.setVal(0.);
      ReactMF.setVal(0.);

#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(S_old,true); mfi.isValid(); ++mfi) 
      {
	 // Note that here we react on the valid region *and* the ghost cells (i.e. the whole FAB)
	 const Box& bx   = mfi.tilebox();

#ifdef TAU
	 reactState(S_old[mfi], S_old[mfi], ReactMF[mfi], tau_diff[mfi], bx, time, 0.5*dt);
#else
	 reactState(S_old[mfi], S_old[mfi], ReactMF[mfi], bx, time, 0.5*dt);
#endif

	 // Synchronize (rho e) and (rho E)
	 BL_FORT_PROC_CALL(RESET_INTERNAL_E,reset_internal_e)
             (ARLIM_3D(bx.loVect()),ARLIM_3D(bx.hiVect()),BL_TO_FORTRAN_3D(S_old[mfi]), verbose);
      }
      if (verbose && ParallelDescriptor::IOProcessor())
	std::cout << "... Leaving burner after completing first half-timestep of burning." << "\n\n";

    }

}

void
#ifdef TAU
Castro::react_second_half_dt(MultiFab& S_new, MultiFab& tau_diff, Real time, Real dt, int ngrow) 
#else
Castro::react_second_half_dt(MultiFab& S_new, Real time, Real dt, int ngrow) 
#endif
{
    BL_PROFILE("Castro::react_second_half_dt()");

    const Real strt_time = ParallelDescriptor::second();

    const Real cur_time = state[State_Type].curTime();

    // Note that here we only react on the valid region of the MultiFab but we may need ghost
    // cells in order to compute things.
    if (do_react == 1) 
    {

      if (verbose && ParallelDescriptor::IOProcessor())
	std::cout << "\n" << "... Entering burner and doing second half-timestep of burning." << "\n";

        MultiFab& ReactMF = get_new_data(Reactions_Type);
        if (ngrow > 0) 
        {
            for (FillPatchIterator Sfpi(*this, S_new, 1, cur_time, State_Type, 0, NUM_STATE);
                 Sfpi.isValid(); ++Sfpi)
            {
	        const Box& bx = Sfpi.tilebox();
#ifdef TAU
                reactState(Sfpi(), Sfpi(), ReactMF[Sfpi], tau_diff[Sfpi], bx, time, 0.5*dt);
#else
                reactState(Sfpi(), Sfpi(), ReactMF[Sfpi], bx, time, 0.5*dt);
#endif
                S_new[Sfpi].copy(Sfpi());
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel
#endif
  	    for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
#ifdef TAU
                reactState(S_new[mfi], S_new[mfi], ReactMF[mfi], tau_diff[mfi], bx, time, 0.5*dt);
#else
                reactState(S_new[mfi], S_new[mfi], ReactMF[mfi], bx, time, 0.5*dt);
#endif
            }
        }
        ReactMF.mult(1.0/dt);
        reset_internal_energy(S_new);
	if (verbose && ParallelDescriptor::IOProcessor())
	  std::cout << "... Leaving burner after completing second half-timestep of burning." << "\n\n";

    }

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
}
#endif
