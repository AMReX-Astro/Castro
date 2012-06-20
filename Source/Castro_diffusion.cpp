#include <winstd.H>

#include "Castro.H"
#include "Castro_F.H"

using std::string;

#ifdef DIFFUSION
#include "Diffusion.H"


void
#ifdef TAU
Castro::add_diffusion_to_old_source (MultiFab& ext_src_old, MultiFab& OldTempDiffTerm, Real prev_time, MultiFab& tau_diff)
#else
Castro::add_diffusion_to_old_source (MultiFab& ext_src_old, MultiFab& OldTempDiffTerm, Real prev_time)
#endif
{
    // Define an explicit temperature update.
    OldTempDiffTerm.setVal(0.);
    if (diffuse_temp == 1) {
#ifdef TAU
       getTempDiffusionTerm(prev_time,OldTempDiffTerm,&tau_diff);
#else
       getTempDiffusionTerm(prev_time,OldTempDiffTerm);
#endif
      MultiFab::Add(ext_src_old,OldTempDiffTerm,0,Eden,1,0);
      MultiFab::Add(ext_src_old,OldTempDiffTerm,0,Eint,1,0);
    }
    geom.FillPeriodicBoundary(ext_src_old,0,NUM_STATE);
}

void
#ifdef TAU
Castro::time_center_diffusion(MultiFab& S_new, MultiFab& OldTempDiffTerm, Real cur_time, Real dt, MultiFab& tau_diff)
#else
Castro::time_center_diffusion(MultiFab& S_new, MultiFab& OldTempDiffTerm, Real cur_time, Real dt)
#endif
{
        // Correct the temperature update so that it will be time-centered.
        MultiFab NewTempDiffTerm(grids,1,1);
        NewTempDiffTerm.setVal(0.);
        if (diffuse_temp == 1) {
#ifdef TAU
           getTempDiffusionTerm(cur_time,NewTempDiffTerm,&tau_diff);
#else
           getTempDiffusionTerm(cur_time,NewTempDiffTerm);
#endif
           NewTempDiffTerm.mult( 0.5*dt);
           OldTempDiffTerm.mult(-0.5*dt);
           // Subtract off half of the old source term, and add half of the new.
           MultiFab::Add(S_new,OldTempDiffTerm,0,Eden,1,0);
           MultiFab::Add(S_new,OldTempDiffTerm,0,Eint,1,0);
           MultiFab::Add(S_new,NewTempDiffTerm,0,Eden,1,0);
           MultiFab::Add(S_new,NewTempDiffTerm,0,Eint,1,0);
           computeTemp(S_new);
        }
}

void
#ifdef TAU
Castro::full_diffusion_update (MultiFab& S_new, Real prev_time, Real cur_time, Real dt, MultiFab& tau_diff)
#else
Castro::full_diffusion_update (MultiFab& S_new, Real prev_time, Real cur_time, Real dt)
#endif
{
        if (diffuse_temp == 1) {
           // Define an explicit temperature update.
           MultiFab OldTempDiffTerm(grids,1,1);
           OldTempDiffTerm.setVal(0.);
#ifdef TAU
           getTempDiffusionTerm(prev_time,OldTempDiffTerm,&tau_diff);
#else
           getTempDiffusionTerm(prev_time,OldTempDiffTerm);
#endif
           OldTempDiffTerm.mult(dt);
           MultiFab::Add(S_new,OldTempDiffTerm,0,Eden,1,0);
           MultiFab::Add(S_new,OldTempDiffTerm,0,Eint,1,0);
           computeTemp(S_new);

           // Correct the temperature update so that it will be time-centered.
           MultiFab NewTempDiffTerm(grids,1,1);
           NewTempDiffTerm.setVal(0.);
#ifdef TAU
           getTempDiffusionTerm(cur_time,NewTempDiffTerm,&tau_diff);
#else
           getTempDiffusionTerm(cur_time,NewTempDiffTerm);
#endif
           NewTempDiffTerm.mult( 0.5*dt);
           OldTempDiffTerm.mult(-0.5*dt);
           // Subtract off half of the old source term, and add half of the new.
           MultiFab::Add(S_new,OldTempDiffTerm,0,Eden,1,0);
           MultiFab::Add(S_new,OldTempDiffTerm,0,Eint,1,0);
           MultiFab::Add(S_new,NewTempDiffTerm,0,Eden,1,0);
           MultiFab::Add(S_new,NewTempDiffTerm,0,Eint,1,0);
           computeTemp(S_new);
        }
}
#endif
