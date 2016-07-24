#include <winstd.H>

#include "Castro.H"
#include "Castro_F.H"

using std::string;

#ifdef DIFFUSION
#include "Diffusion.H"

void
Castro::construct_old_diff_source(PArray<MultiFab>& old_sources,
				  MultiFab& sources_for_hydro,
				  Real time, Real dt)
{
    add_temp_diffusion_to_source(old_sources[diff_src],*OldTempDiffTerm,time);

#if (BL_SPACEDIM == 1)
    add_spec_diffusion_to_source(old_sources[diff_src],*OldSpecDiffTerm,time);
    add_viscous_term_to_source(old_sources[diff_src],*OldViscousTermforMomentum,OldViscousTermforEnergy,time);
#endif

    BoxLib::fill_boundary(old_sources[diff_src], geom);

    // Add to the hydro source terms.

    MultiFab::Add(sources_for_hydro,old_sources[diff_src],0,0,NUM_STATE,NUM_GROW);
}

void
Castro::construct_new_diff_source(PArray<MultiFab>& old_sources,
				  PArray<MultiFab>& new_sources,
				  MultiFab& sources_for_hydro,
				  Real time, Real dt)
{
    add_temp_diffusion_to_source(new_sources[diff_src],*NewTempDiffTerm,time);

#if (BL_SPACEDIM == 1)
    add_spec_diffusion_to_source(new_sources[diff_src],*NewSpecDiffTerm,time);
    add_viscous_term_to_source(new_sources[diff_src],*NewViscousTermforMomentum,NewViscousTermforEnergy,time);
#endif

    // Time center the source term.

    old_sources[diff_src].mult(-0.5);
    new_sources[diff_src].mult( 0.5);

    MultiFab::Add(new_sources[diff_src],old_sources[diff_src],0,0,NUM_STATE,0);

    old_sources[diff_src].mult(-2.0);

    // Add to the hydro source terms.

    MultiFab::Add(sources_for_hydro,new_sources[diff_src],0,0,NUM_STATE,0);

}

// **********************************************************************************************

void
Castro::add_temp_diffusion_to_source (MultiFab& ext_src, MultiFab& DiffTerm, Real t)
{
    // Define an explicit temperature update.
    DiffTerm.setVal(0.);
    if (diffuse_temp == 1) {
#ifdef TAU
       getTempDiffusionTerm(t,DiffTerm,tau_diff);
#else
       getTempDiffusionTerm(t,DiffTerm);
#endif
    } else if (diffuse_enth == 1) {
       getEnthDiffusionTerm(t,DiffTerm);
    }

    if (diffuse_temp == 1 or diffuse_enth == 1) {
       int ng = std::min(ext_src.nGrow(),DiffTerm.nGrow());
       MultiFab::Add(ext_src,DiffTerm,0,Eden,1,ng);
       MultiFab::Add(ext_src,DiffTerm,0,Eint,1,ng);
    }
}

// **********************************************************************************************

void
#ifdef TAU
Castro::time_center_temp_diffusion(MultiFab& S_new, MultiFab& OldDiffTerm, Real cur_time, Real dt, MultiFab& tau_diff)
#else
Castro::time_center_temp_diffusion(MultiFab& S_new, MultiFab& OldDiffTerm, Real cur_time, Real dt)
#endif
{
        // Correct the temperature update so that it will be time-centered.
        MultiFab NewDiffTerm(grids,1,1);
        NewDiffTerm.setVal(0.);
        if (diffuse_temp == 1) {
#ifdef TAU
           getTempDiffusionTerm(cur_time,NewDiffTerm,&tau_diff);
#else
           getTempDiffusionTerm(cur_time,NewDiffTerm);
#endif
        } else if (diffuse_enth == 1) {
           getEnthDiffusionTerm(cur_time,NewDiffTerm);
        }
        if (diffuse_temp == 1 or diffuse_enth == 1) {
           NewDiffTerm.mult( 0.5*dt);
           OldDiffTerm.mult(-0.5*dt);
           // Subtract off half of the old source term, and add half of the new.
           MultiFab::Add(S_new,OldDiffTerm,0,Eden,1,0);
           MultiFab::Add(S_new,OldDiffTerm,0,Eint,1,0);
           MultiFab::Add(S_new,NewDiffTerm,0,Eden,1,0);
           MultiFab::Add(S_new,NewDiffTerm,0,Eint,1,0);
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
        MultiFab OldDiffTerm, NewDiffTerm;
        if (diffuse_temp == 1 or diffuse_enth == 1) {
           // Define an explicit temperature/enthalpy update.
           OldDiffTerm.define(grids,1,1,Fab_allocate);
           OldDiffTerm.setVal(0.);
        }

        if (diffuse_temp == 1) {
#ifdef TAU
           getTempDiffusionTerm(prev_time,OldDiffTerm,&tau_diff);
#else
           getTempDiffusionTerm(prev_time,OldDiffTerm);
#endif
        } else if (diffuse_enth == 1) {
           getEnthDiffusionTerm(prev_time,OldDiffTerm);
        }

        if (diffuse_temp == 1 or diffuse_enth == 1) {
           OldDiffTerm.mult(dt);
           MultiFab::Add(S_new,OldDiffTerm,0,Eden,1,0);
           MultiFab::Add(S_new,OldDiffTerm,0,Eint,1,0);
           computeTemp(S_new);

           // Correct the temperature update so that it will be time-centered.
           NewDiffTerm.define(grids,1,1,Fab_allocate);
           NewDiffTerm.setVal(0.);
        }

        if (diffuse_temp == 1) {
#ifdef TAU
           getTempDiffusionTerm(cur_time,NewDiffTerm,&tau_diff);
#else
           getTempDiffusionTerm(cur_time,NewDiffTerm);
#endif
        } else if (diffuse_enth == 1) {
           getEnthDiffusionTerm(cur_time,NewDiffTerm);
        }

        if (diffuse_temp == 1 or diffuse_enth == 1) {
           NewDiffTerm.mult( 0.5*dt);
           OldDiffTerm.mult(-0.5*dt);
           // Subtract off half of the old source term, and add half of the new.
           MultiFab::Add(S_new,OldDiffTerm,0,Eden,1,0);
           MultiFab::Add(S_new,OldDiffTerm,0,Eint,1,0);
           MultiFab::Add(S_new,NewDiffTerm,0,Eden,1,0);
           MultiFab::Add(S_new,NewDiffTerm,0,Eint,1,0);
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
