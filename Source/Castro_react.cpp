#include <winstd.H>

#include "Castro.H"
#include "Castro_F.H"

using std::string;

#ifdef REACTIONS
void
#ifdef TAU
Castro::react_half_dt(MultiFab& s, MultiFab& r, MultiFab& tau_diff, Real time, Real dt, int ngrow) 
#else
Castro::react_half_dt(MultiFab& s, MultiFab& r, Real time, Real dt, int ngrow) 
#endif
{
    BL_PROFILE("Castro::react_half_dt()");

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
			 BL_TO_FORTRAN_3D(tau_diff[mfi]),
#endif
			 time, 0.5 * dt);

	}

        reset_internal_energy(s);

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
	  std::cout << "Castro::react_half_dt() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
	});
#endif
    }
}
#endif
