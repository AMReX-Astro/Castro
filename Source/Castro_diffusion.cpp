#include <winstd.H>

#include "Castro.H"
#include "Castro_F.H"

using std::string;

#include "Diffusion.H"

void
Castro::construct_old_diff_source(Real time, Real dt)
{
    int ng = (*Sborder).nGrow();

    old_sources[diff_src].setVal(0.0, ng);

    add_temp_diffusion_to_source(old_sources[diff_src],*OldTempDiffTerm,time);

#if (BL_SPACEDIM == 1)
    add_spec_diffusion_to_source(old_sources[diff_src],*OldSpecDiffTerm,time);
    add_viscous_term_to_source(old_sources[diff_src],*OldViscousTermforMomentum,OldViscousTermforEnergy,time);
#endif

    BoxLib::fill_boundary(old_sources[diff_src], geom);

    // Add to the hydro source terms.

    MultiFab::Add(*sources_for_hydro,old_sources[diff_src],0,0,NUM_STATE,NUM_GROW);
}

void
Castro::construct_new_diff_source(Real time, Real dt)
{
    MultiFab& S_new = get_new_data(State_Type);

    int ng = S_new.nGrow();

    new_sources[diff_src].setVal(0.0, ng);

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

    if (source_term_predictor == 1)
        MultiFab::Add(*sources_for_hydro,new_sources[diff_src],0,0,NUM_STATE,0);

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

void
Castro::getTempDiffusionTerm (Real time, MultiFab& TempDiffTerm)
{
    BL_PROFILE("Castro::getTempDiffusionTerm()");

   MultiFab& S_old = get_old_data(State_Type);
   if (verbose && ParallelDescriptor::IOProcessor())
      std::cout << "Calculating diffusion term at time " << time << std::endl;

#ifdef TAU
   if (tau_diff == 0)
     std::cerr << "ERROR:tau must be defined if USE_TAU = TRUE " << std::endl;
#endif

   // Fill coefficients at this level.
   PArray<MultiFab> coeffs(BL_SPACEDIM,PArrayManage);
   PArray<MultiFab> coeffs_temporary(3,PArrayManage); // This is what we pass to the dimension-agnostic Fortran
   for (int dir = 0; dir < 3; dir++) {
       if (dir < BL_SPACEDIM) {
	 coeffs.set(dir,new MultiFab(getEdgeBoxArray(dir), 1, 0, Fab_allocate));
	 coeffs_temporary.set(dir,new MultiFab(getEdgeBoxArray(dir), 1, 0, Fab_allocate));
       } else {
	 coeffs_temporary.set(dir,new MultiFab(grids, 1, 0, Fab_allocate));
       }
   }

   // Fill temperature at this level.
   MultiFab Temperature(grids,1,1,Fab_allocate);

   {
       FillPatchIterator fpi(*this,S_old,1,time,State_Type,0,NUM_STATE);
       MultiFab& state_old = fpi.get_mf();

       MultiFab::Copy(Temperature, state_old, Temp, 0, 1, 1);

       for (MFIter mfi(state_old); mfi.isValid(); ++mfi)
       {
	   const Box& bx = grids[mfi.index()];

	   ca_fill_temp_cond(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			     BL_TO_FORTRAN_3D(state_old[mfi]),
#ifdef TAU
			     BL_TO_FORTRAN_3D((*tau_diff)[mfi]),
#endif
			     BL_TO_FORTRAN_3D(coeffs_temporary[0][mfi]),
			     BL_TO_FORTRAN_3D(coeffs_temporary[1][mfi]),
			     BL_TO_FORTRAN_3D(coeffs_temporary[2][mfi]));
       }
   }

   // Now copy the temporary array results back to the
   // correctly dimensioned coeffs array.

   for (int dir = 0; dir < BL_SPACEDIM; dir++)
     MultiFab::Copy(coeffs[dir], coeffs_temporary[dir], 0, 0, 1, 0);

   if (Geometry::isAnyPeriodic())
     for (int d = 0; d < BL_SPACEDIM; d++)
       geom.FillPeriodicBoundary(coeffs[d]);

   MultiFab CrseTemp;
   if (level > 0) {
       // Fill temperature at next coarser level, if it exists.
       const BoxArray& crse_grids = getLevel(level-1).boxArray();
       CrseTemp.define(crse_grids,1,1,Fab_allocate);
       FillPatch(getLevel(level-1),CrseTemp,1,time,State_Type,Temp,1);
   }

   diffusion->applyop(level,Temperature,CrseTemp,TempDiffTerm,coeffs);

   // Extrapolate to ghost cells
   if (TempDiffTerm.nGrow() > 0) {
       for (MFIter mfi(TempDiffTerm); mfi.isValid(); ++mfi)
       {
	   const Box& bx = mfi.validbox();
	   ca_tempdiffextrap(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			     BL_TO_FORTRAN_3D(TempDiffTerm[mfi]));
       }
   }
}

void
Castro::getEnthDiffusionTerm (Real time, MultiFab& DiffTerm)
{
    BL_PROFILE("Castro::getEnthDiffusionTerm()");

   MultiFab& S_old = get_old_data(State_Type);
   if (verbose && ParallelDescriptor::IOProcessor())
      std::cout << "Calculating diffusion term at time " << time << std::endl;

   // Fill coefficients at this level.
   PArray<MultiFab> coeffs(BL_SPACEDIM,PArrayManage);
   PArray<MultiFab> coeffs_temporary(3,PArrayManage); // This is what we pass to the dimension-agnostic Fortran
   for (int dir = 0; dir < 3; dir++) {
       if (dir < BL_SPACEDIM) {
	 coeffs.set(dir,new MultiFab(getEdgeBoxArray(dir), 1, 0, Fab_allocate));
	 coeffs_temporary.set(dir,new MultiFab(getEdgeBoxArray(dir), 1, 0, Fab_allocate));
       } else {
	 coeffs_temporary.set(dir,new MultiFab(grids, 1, 0, Fab_allocate));
       }
   }

   // Define enthalpy at this level.
   MultiFab Enthalpy(grids,1,1,Fab_allocate);
   {
       FillPatchIterator fpi(*this,S_old,1,time,State_Type,0,NUM_STATE);
       MultiFab& state_old = fpi.get_mf();

       for (MFIter mfi(state_old); mfi.isValid(); ++mfi)
       {
	   const Box& bx = grids[mfi.index()];
	   make_enthalpy(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
	                 BL_TO_FORTRAN_3D(state_old[mfi]),
	                 BL_TO_FORTRAN_3D(Enthalpy[mfi]));

	   ca_fill_enth_cond(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			     BL_TO_FORTRAN_3D(state_old[mfi]),
			     BL_TO_FORTRAN_3D(coeffs_temporary[0][mfi]),
			     BL_TO_FORTRAN_3D(coeffs_temporary[1][mfi]),
			     BL_TO_FORTRAN_3D(coeffs_temporary[2][mfi]));
       }
   }

   // Now copy the temporary array results back to the
   // correctly dimensioned coeffs array.

   for (int dir = 0; dir < BL_SPACEDIM; dir++)
     MultiFab::Copy(coeffs[dir], coeffs_temporary[dir], 0, 0, 1, 0);

   if (Geometry::isAnyPeriodic())
     for (int d = 0; d < BL_SPACEDIM; d++)
       geom.FillPeriodicBoundary(coeffs[d]);

   MultiFab CrseEnth, CrseState;
   if (level > 0) {
       // Fill temperature at next coarser level, if it exists.
       const BoxArray& crse_grids = getLevel(level-1).boxArray();
       CrseEnth.define (crse_grids,1,1,Fab_allocate);
       CrseState.define(crse_grids,NUM_STATE,1,Fab_allocate);
       FillPatch(getLevel(level-1),CrseState,1,time,State_Type,Density,NUM_STATE);

       for (MFIter mfi(CrseState); mfi.isValid(); ++mfi)
       {
	   const Box& bx = crse_grids[mfi.index()];
	   make_enthalpy(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
	                 BL_TO_FORTRAN_3D(CrseState[mfi]),
	                 BL_TO_FORTRAN_3D( CrseEnth[mfi]));
       }
   }

   diffusion->applyop(level,Enthalpy,CrseEnth,DiffTerm,coeffs);

   // Extrapolate to ghost cells
   if (DiffTerm.nGrow() > 0) {
       for (MFIter mfi(DiffTerm); mfi.isValid(); ++mfi)
       {
	   const Box& bx = mfi.validbox();
	   ca_tempdiffextrap(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			     BL_TO_FORTRAN_3D(DiffTerm[mfi]));
       }
   }
}

void
Castro::getSpecDiffusionTerm (Real time, MultiFab& SpecDiffTerm)
{
    BL_PROFILE("Castro::getSpecDiffusionTerm()");

   MultiFab& S_old = get_old_data(State_Type);
   if (verbose && ParallelDescriptor::IOProcessor())
      std::cout << "Calculating species diffusion term at time " << time << std::endl;

   // Fill coefficients at this level.
   PArray<MultiFab> coeffs(BL_SPACEDIM,PArrayManage);
   PArray<MultiFab> coeffs_temporary(3,PArrayManage); // This is what we pass to the dimension-agnostic Fortran
   for (int dir = 0; dir < 3; dir++) {
       if (dir < BL_SPACEDIM) {
	 coeffs.set(dir,new MultiFab(getEdgeBoxArray(dir), 1, 0, Fab_allocate));
	 coeffs_temporary.set(dir,new MultiFab(getEdgeBoxArray(dir), 1, 0, Fab_allocate));
       } else {
	 coeffs_temporary.set(dir,new MultiFab(grids, 1, 0, Fab_allocate));
       }
   }

   FillPatchIterator fpi(*this,S_old,1,time,State_Type,0,NUM_STATE);
   MultiFab& state_old = fpi.get_mf();

   for (MFIter mfi(state_old); mfi.isValid(); ++mfi)
   {
       const Box& bx = grids[mfi.index()];

       ca_fill_spec_coeff(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			  BL_TO_FORTRAN_3D(state_old[mfi]),
			  BL_TO_FORTRAN_3D(coeffs_temporary[0][mfi]),
			  BL_TO_FORTRAN_3D(coeffs_temporary[1][mfi]),
			  BL_TO_FORTRAN_3D(coeffs_temporary[2][mfi]));
   }

   // Now copy the temporary array results back to the
   // correctly dimensioned coeffs array.
   for (int dir = 0; dir < BL_SPACEDIM; dir++)
     MultiFab::Copy(coeffs[dir], coeffs_temporary[dir], 0, 0, 1, 0);

   if (Geometry::isAnyPeriodic())
     for (int d = 0; d < BL_SPACEDIM; d++)
       geom.FillPeriodicBoundary(coeffs[d]);

   // Create MultiFabs that only hold the data for one species at a time.
   MultiFab Species(grids,1,1,Fab_allocate);
   MultiFab     SDT(grids,SpecDiffTerm.nComp(),SpecDiffTerm.nGrow(),Fab_allocate);
   SDT.setVal(0.);
   MultiFab CrseSpec, CrseDen;
   if (level > 0) {
       const BoxArray& crse_grids = getLevel(level-1).boxArray();
       CrseSpec.define(crse_grids,1,1,Fab_allocate);
        CrseDen.define(crse_grids,1,1,Fab_allocate);
   }

   // Fill one species at a time at this level.
   for (int ispec = 0; ispec < NumSpec; ispec++)
   {
       MultiFab::Copy  (Species, state_old, FirstSpec+ispec, 0, 1, 1);
       MultiFab::Divide(Species, state_old, Density        , 0, 1, 1);

       // Fill temperature at next coarser level, if it exists.
       if (level > 0)
       {
           FillPatch(getLevel(level-1),CrseSpec,1,time,State_Type,FirstSpec+ispec,1);
           FillPatch(getLevel(level-1),CrseDen ,1,time,State_Type,Density        ,1);
           MultiFab::Divide(CrseSpec, CrseDen, 0, 0, 1, 1);
       }

       diffusion->applyop(level,Species,CrseSpec,SDT,coeffs);

       // Extrapolate to ghost cells
       if (SDT.nGrow() > 0) {
           for (MFIter mfi(SDT); mfi.isValid(); ++mfi)
           {
	       const Box& bx = mfi.validbox();
	       ca_tempdiffextrap(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
				 BL_TO_FORTRAN_3D(SDT[mfi]));
           }
       }
       // Copy back into SpecDiffTerm from the temporary SDT
       MultiFab::Copy(SpecDiffTerm, SDT, 0, ispec, 1, 1);
   }
}

#if (BL_SPACEDIM == 1)
// **********************************************
// Note: this currently just gets the term that looks like div(2 mu grad(u)) which is
//       only part of the viscous term.  We assume that the coefficient that is filled is "2 mu"
// **********************************************
void
Castro::getViscousTerm (Real time, MultiFab& ViscousTermforMomentum, MultiFab& ViscousTermforEnergy)
{
    BL_PROFILE("Castro::getViscousTerm()");

   if (verbose && ParallelDescriptor::IOProcessor())
      std::cout << "Calculating viscous term at time " << time << std::endl;

   getFirstViscousTerm(time,ViscousTermforMomentum);

   MultiFab SecndTerm(grids,ViscousTermforMomentum.nComp(),ViscousTermforMomentum.nGrow(),Fab_allocate);
   getSecndViscousTerm(time,SecndTerm);
   MultiFab::Add(ViscousTermforMomentum, SecndTerm, 0, 0, ViscousTermforMomentum.nComp(), 0);

   getViscousTermForEnergy(time,ViscousTermforEnergy);
}

void
Castro::getFirstViscousTerm (Real time, MultiFab& ViscousTerm)
{
   MultiFab& S_old = get_old_data(State_Type);

   // Fill coefficients at this level.
   PArray<MultiFab> coeffs(BL_SPACEDIM,PArrayManage);
   PArray<MultiFab> coeffs_temporary(3,PArrayManage); // This is what we pass to the dimension-agnostic Fortran
   for (int dir = 0; dir < 3; dir++) {
       if (dir < BL_SPACEDIM) {
	 coeffs.set(dir,new MultiFab(getEdgeBoxArray(dir), 1, 0, Fab_allocate));
	 coeffs_temporary.set(dir,new MultiFab(getEdgeBoxArray(dir), 1, 0, Fab_allocate));
       } else {
	 coeffs_temporary.set(dir,new MultiFab(grids, 1, 0, Fab_allocate));
       }
   }

   // Fill velocity at this level.
   MultiFab Vel(grids,1,1,Fab_allocate);

   FillPatchIterator fpi(*this,S_old,1,time,State_Type,0,NUM_STATE);
   MultiFab& state_old = fpi.get_mf();

   // Remember this is just 1-d
   MultiFab::Copy  (Vel, state_old, Xmom   , 0, 1, 1);
   MultiFab::Divide(Vel, state_old, Density, 0, 1, 1);

   for (MFIter mfi(state_old); mfi.isValid(); ++mfi)
   {
       const Box& bx = grids[mfi.index()];

       ca_fill_first_visc_coeff(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
				BL_TO_FORTRAN_3D(state_old[mfi]),
				BL_TO_FORTRAN_3D(coeffs_temporary[0][mfi]),
				BL_TO_FORTRAN_3D(coeffs_temporary[1][mfi]),
				BL_TO_FORTRAN_3D(coeffs_temporary[2][mfi]));
   }

   // Now copy the temporary array results back to the
   // correctly dimensioned coeffs array.
   for (int dir = 0; dir < BL_SPACEDIM; dir++)
     MultiFab::Copy(coeffs[dir], coeffs_temporary[dir], 0, 0, 1, 0);

   if (Geometry::isAnyPeriodic())
     for (int d = 0; d < BL_SPACEDIM; d++)
       geom.FillPeriodicBoundary(coeffs[d]);

   MultiFab CrseVel, CrseDen;
   if (level > 0) {
       // Fill temperature at next coarser level, if it exists.
       const BoxArray& crse_grids = getLevel(level-1).boxArray();
       CrseVel.define(crse_grids,1,1,Fab_allocate);
       CrseDen.define(crse_grids,1,1,Fab_allocate);
       FillPatch(getLevel(level-1),CrseVel ,1,time,State_Type,Xmom   ,1);
       FillPatch(getLevel(level-1),CrseDen ,1,time,State_Type,Density,1);
       MultiFab::Divide(CrseVel, CrseDen, 0, 0, 1, 1);
   }
   diffusion->applyop(level,Vel,CrseVel,ViscousTerm,coeffs);

   // Extrapolate to ghost cells
   if (ViscousTerm.nGrow() > 0) {
       for (MFIter mfi(ViscousTerm); mfi.isValid(); ++mfi)
       {
	   const Box& bx = mfi.validbox();
	   ca_tempdiffextrap(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			     BL_TO_FORTRAN_3D(ViscousTerm[mfi]));
       }
   }
}

void
Castro::getSecndViscousTerm (Real time, MultiFab& ViscousTerm)
{
   MultiFab& S_old = get_old_data(State_Type);

   // Fill coefficients at this level.
   PArray<MultiFab> coeffs(BL_SPACEDIM,PArrayManage);
   PArray<MultiFab> coeffs_temporary(3,PArrayManage); // This is what we pass to the dimension-agnostic Fortran
   for (int dir = 0; dir < 3; dir++) {
       if (dir < BL_SPACEDIM) {
	 coeffs.set(dir,new MultiFab(getEdgeBoxArray(dir), 1, 0, Fab_allocate));
	 coeffs_temporary.set(dir,new MultiFab(getEdgeBoxArray(dir), 1, 0, Fab_allocate));
       } else {
	 coeffs_temporary.set(dir,new MultiFab(grids, 1, 0, Fab_allocate));
       }
   }

   // Fill velocity at this level.
   MultiFab Vel(grids,1,1,Fab_allocate);

   FillPatchIterator fpi(*this,S_old,1,time,State_Type,0,NUM_STATE);
   MultiFab& state_old = fpi.get_mf();

   // Remember this is just 1-d
   MultiFab::Copy  (Vel, state_old, Xmom   , 0, 1, 1);
   MultiFab::Divide(Vel, state_old, Density, 0, 1, 1);

   for (MFIter mfi(state_old); mfi.isValid(); ++mfi)
   {
       const Box& bx = grids[mfi.index()];

       ca_fill_secnd_visc_coeff(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
				BL_TO_FORTRAN_3D(state_old[mfi]),
				BL_TO_FORTRAN_3D(coeffs_temporary[0][mfi]),
				BL_TO_FORTRAN_3D(coeffs_temporary[1][mfi]),
				BL_TO_FORTRAN_3D(coeffs_temporary[2][mfi]));
   }

   // Now copy the temporary array results back to the
   // correctly dimensioned coeffs array.
   for (int dir = 0; dir < BL_SPACEDIM; dir++)
     MultiFab::Copy(coeffs[dir], coeffs_temporary[dir], 0, 0, 1, 0);

   if (Geometry::isAnyPeriodic())
     for (int d = 0; d < BL_SPACEDIM; d++)
       geom.FillPeriodicBoundary(coeffs[d]);

   MultiFab CrseVel, CrseDen;
   if (level > 0) {
       // Fill temperature at next coarser level, if it exists.
       const BoxArray& crse_grids = getLevel(level-1).boxArray();
       CrseVel.define(crse_grids,1,1,Fab_allocate);
       CrseDen.define(crse_grids,1,1,Fab_allocate);
       FillPatch(getLevel(level-1),CrseVel ,1,time,State_Type,Xmom   ,1);
       FillPatch(getLevel(level-1),CrseDen ,1,time,State_Type,Density,1);
       MultiFab::Divide(CrseVel, CrseDen, 0, 0, 1, 1);
   }
   diffusion->applyViscOp(level,Vel,CrseVel,ViscousTerm,coeffs);

   // Extrapolate to ghost cells
   if (ViscousTerm.nGrow() > 0) {
       for (MFIter mfi(ViscousTerm); mfi.isValid(); ++mfi)
       {
	   const Box& bx = mfi.validbox();
	   ca_tempdiffextrap(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			     BL_TO_FORTRAN_3D(ViscousTerm[mfi]));
       }
   }
}

void
Castro::getViscousTermForEnergy (Real time, MultiFab& ViscousTerm)
{
   MultiFab& S_old = get_old_data(State_Type);

   // Fill coefficients at this level.
   PArray<MultiFab> coeffs(BL_SPACEDIM,PArrayManage);
   PArray<MultiFab> coeffs_temporary(3,PArrayManage); // This is what we pass to the dimension-agnostic Fortran
   for (int dir = 0; dir < 3; dir++) {
       if (dir < BL_SPACEDIM) {
	 coeffs.set(dir,new MultiFab(getEdgeBoxArray(dir), 1, 0, Fab_allocate));
	 coeffs_temporary.set(dir,new MultiFab(getEdgeBoxArray(dir), 1, 0, Fab_allocate));
       } else {
	 coeffs_temporary.set(dir,new MultiFab(grids, 1, 0, Fab_allocate));
       }
   }

   FillPatchIterator fpi(*this,S_old,2,time,State_Type,0,NUM_STATE);
   MultiFab& state_old = fpi.get_mf();

   const Geometry& fine_geom = parent->Geom(parent->finestLevel());
   const Real*       dx_fine = fine_geom.CellSize();

   // Remember this is just 1-d
   int coord_type = Geometry::Coord();
   for (MFIter mfi(state_old); mfi.isValid(); ++mfi)
   {
       const Box& bx = grids[mfi.index()];

       ca_compute_div_tau_u(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			    BL_TO_FORTRAN_3D(ViscousTerm[mfi]),
			    BL_TO_FORTRAN_3D(state_old[mfi]),
			    ZFILL(dx_fine),&coord_type);
   }

   // Extrapolate to ghost cells
   if (ViscousTerm.nGrow() > 0) {
       for (MFIter mfi(ViscousTerm); mfi.isValid(); ++mfi)
       {
	   const Box& bx = mfi.validbox();
	   ca_tempdiffextrap(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			     BL_TO_FORTRAN_3D(ViscousTerm[mfi]));
       }
   }
}
#endif
