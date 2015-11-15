#include <winstd.H>

#include "Castro.H"
#include "Castro_F.H"

using std::string;

#ifdef DIFFUSION
#include "Diffusion.H"

// **********************************************************************************************

void
#ifdef TAU
Castro::add_temp_diffusion_to_source (MultiFab& ext_src, MultiFab& TempDiffTerm, Real t, MultiFab& tau_diff)
#else
Castro::add_temp_diffusion_to_source (MultiFab& ext_src, MultiFab& TempDiffTerm, Real t)
#endif
{
    // Define an explicit temperature update.
    TempDiffTerm.setVal(0.);
    if (diffuse_temp == 1) {
#ifdef TAU
       getTempDiffusionTerm(t,TempDiffTerm,&tau_diff);
#else
       getTempDiffusionTerm(t,TempDiffTerm);
#endif
       int ng = std::min(ext_src.nGrow(),TempDiffTerm.nGrow());
       MultiFab::Add(ext_src,TempDiffTerm,0,Eden,1,ng);
       MultiFab::Add(ext_src,TempDiffTerm,0,Eint,1,ng);
    }
}

// **********************************************************************************************

void
#ifdef TAU
Castro::time_center_temp_diffusion(MultiFab& S_new, MultiFab& OldTempDiffTerm, Real cur_time, Real dt, MultiFab& tau_diff)
#else
Castro::time_center_temp_diffusion(MultiFab& S_new, MultiFab& OldTempDiffTerm, Real cur_time, Real dt)
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

// **********************************************************************************************

void
#ifdef TAU
Castro::full_temp_diffusion_update (MultiFab& S_new, Real prev_time, Real cur_time, Real dt, MultiFab& tau_diff)
#else
Castro::full_temp_diffusion_update (MultiFab& S_new, Real prev_time, Real cur_time, Real dt)
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

// **********************************************************************************************

void
Castro::add_spec_diffusion_to_source (MultiFab& ext_src, MultiFab& SpecDiffTerm, Real t)
{
    // Define an explicit species update.
    SpecDiffTerm.setVal(0.);
    if (diffuse_spec == 1) {
       getSpecDiffusionTerm(t,SpecDiffTerm);
       int ng = std::min(ext_src.nGrow(),SpecDiffTerm.nGrow());
       MultiFab::Add(ext_src,SpecDiffTerm,0,FirstSpec,NumSpec,ng);
    }
}

// **********************************************************************************************

void
Castro::time_center_spec_diffusion(MultiFab& S_new, MultiFab& OldSpecDiffTerm, Real cur_time, Real dt)
{
        // Correct the species update so that it will be time-centered.
        MultiFab NewSpecDiffTerm(grids,NumSpec,1);
        NewSpecDiffTerm.setVal(0.);
        if (diffuse_spec == 1) {
           getSpecDiffusionTerm(cur_time,NewSpecDiffTerm);
           NewSpecDiffTerm.mult( 0.5*dt);
           OldSpecDiffTerm.mult(-0.5*dt);
           // Subtract off half of the old source term, and add half of the new.
           MultiFab::Add(S_new,OldSpecDiffTerm,0,FirstSpec,NumSpec,0);
           MultiFab::Add(S_new,NewSpecDiffTerm,0,FirstSpec,NumSpec,0);
           computeTemp(S_new);
        }
}

// **********************************************************************************************

void
Castro::full_spec_diffusion_update (MultiFab& S_new, Real prev_time, Real cur_time, Real dt)
{
        if (diffuse_spec == 1) {
           // Define an explicit species update.
           MultiFab OldSpecDiffTerm(grids,NumSpec,1);
           OldSpecDiffTerm.setVal(0.);
           getSpecDiffusionTerm(prev_time,OldSpecDiffTerm);
           OldSpecDiffTerm.mult(dt);
           MultiFab::Add(S_new,OldSpecDiffTerm,0,FirstSpec,NumSpec,0);

           // Correct the species update so that it will be time-centered.
           MultiFab NewSpecDiffTerm(grids,NumSpec,1);
           NewSpecDiffTerm.setVal(0.);
           getSpecDiffusionTerm(cur_time,NewSpecDiffTerm);
           NewSpecDiffTerm.mult( 0.5*dt);
           OldSpecDiffTerm.mult(-0.5*dt);
           // Subtract off half of the old source term, and add half of the new.
           MultiFab::Add(S_new,OldSpecDiffTerm,0,FirstSpec,NumSpec,0);
           MultiFab::Add(S_new,NewSpecDiffTerm,0,FirstSpec,NumSpec,0);

           computeTemp(S_new);
        }
}

#if (BL_SPACEDIM == 1)
// **********************************************************************************************

void
Castro::add_viscous_term_to_source(MultiFab& ext_src, MultiFab& ViscousTermforMomentum, 
                                   MultiFab& ViscousTermforEnergy, Real t)
{
    // Define an explicit viscous term
    ViscousTermforMomentum.setVal(0.);
    ViscousTermforEnergy.setVal(0.);
    if (diffuse_vel == 1) {
       getViscousTerm(t,ViscousTermforMomentum,ViscousTermforEnergy);
       int ng = std::min(ext_src.nGrow(),ViscousTermforMomentum.nGrow());
       MultiFab::Add(ext_src,ViscousTermforMomentum,0,Xmom,1,ng);

       ng = std::min(ext_src.nGrow(),ViscousTermforEnergy.nGrow());
       MultiFab::Add(ext_src,ViscousTermforEnergy  ,0,Eden,1,ng);
    }
}

// **********************************************************************************************

void
Castro::time_center_viscous_term(MultiFab& S_new, MultiFab& OldViscousTermforMomentum,
                                 MultiFab& OldViscousTermforEnergy, Real cur_time, Real dt)
{
        // Correct the species update so that it will be time-centered.
        MultiFab NewViscousTermforMomentum(grids,BL_SPACEDIM,1);
        MultiFab NewViscousTermforEnergy  (grids,1,1);
        NewViscousTermforMomentum.setVal(0.);
        NewViscousTermforEnergy.setVal(0.);
        if (diffuse_vel == 1) {

           getViscousTerm(cur_time,NewViscousTermforMomentum,NewViscousTermforEnergy);

           NewViscousTermforMomentum.mult( 0.5*dt);
           OldViscousTermforMomentum.mult(-0.5*dt);
           MultiFab::Add(S_new,OldViscousTermforMomentum,0,Xmom,BL_SPACEDIM,0);
           MultiFab::Add(S_new,NewViscousTermforMomentum,0,Xmom,BL_SPACEDIM,0);

           NewViscousTermforEnergy.mult( 0.5*dt);
           OldViscousTermforEnergy.mult(-0.5*dt);
           MultiFab::Add(S_new,OldViscousTermforEnergy  ,0,Eden,1          ,0);
           MultiFab::Add(S_new,NewViscousTermforEnergy  ,0,Eden,1          ,0);

           computeTemp(S_new);
        }
}
#endif
#endif
