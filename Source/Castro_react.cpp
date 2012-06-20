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
    // Make sure to zero these even if do_react == 0.
    MultiFab& ReactMF_old = get_old_data(Reactions_Type);
    ReactMF_old.setVal(0.);
    MultiFab& ReactMF = get_new_data(Reactions_Type);
    ReactMF.setVal(0.);
    if (do_react == 1)
    {
#ifdef TAU
        strang_chem(S_old,ReactMF,tau_diff,time,dt);
#else
        strang_chem(S_old,ReactMF,time,dt);
#endif
        reset_internal_energy(S_old);
    }
}

void
#ifdef TAU
Castro::react_second_half_dt(MultiFab& S_new, MultiFab& tau_diff, Real cur_time, Real dt) 
#else
Castro::react_second_half_dt(MultiFab& S_new, Real cur_time, Real dt) 
#endif
{
    if (do_react == 1) 
    {
        MultiFab& ReactMF = get_new_data(Reactions_Type);

#ifdef TAU
        strang_chem(S_new,ReactMF,tau_diff,cur_time,dt);
#else
        strang_chem(S_new,ReactMF,cur_time,dt);
#endif
        ReactMF.mult(1.0/dt);
        reset_internal_energy(S_new);
    }
}

void
Castro::strang_chem (MultiFab&  state,
                     MultiFab&  React_mf,
#ifdef TAU
                     MultiFab&  tau,
#endif
                     Real       time,
                     Real       dt)
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::strang_chem(MultiFab&,...");
    const Real strt_time = ParallelDescriptor::second();

    for (MFIter Smfi(state); Smfi.isValid(); ++Smfi)
    {
        FArrayBox& fb   = state[Smfi];
        const Box& bx   = Smfi.validbox();
#ifdef TAU
        FArrayBox& tfab = tau[Smfi];
        reactState(fb, fb, React_mf[Smfi], tfab, bx, time, 0.5*dt);
#else
        reactState(fb, fb, React_mf[Smfi], bx, time, 0.5*dt);
#endif
    }

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

//      if (ParallelDescriptor::IOProcessor()) 
//          std::cout << "strang_chem time = " << run_time << '\n';
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
    BL_FORT_PROC_CALL(CA_REACT_STATE,ca_react_state)
                     (box.loVect(), box.hiVect(), 
                     BL_TO_FORTRAN(Sold),
                     BL_TO_FORTRAN(Snew),
                     BL_TO_FORTRAN(ReactionTerms),
#ifdef TAU
                     BL_TO_FORTRAN(tau),
#endif
                     time,dt_react);
}
#endif
