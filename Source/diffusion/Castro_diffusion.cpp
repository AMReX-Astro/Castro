
#include "Castro.H"
#include "Castro_F.H"

using std::string;

#include "Diffusion.H"

using namespace amrex;

void
Castro::construct_old_diff_source(MultiFab& source, MultiFab& state, Real time, Real dt)
{
    MultiFab TempDiffTerm(grids, dmap, 1, 1);
    MultiFab SpecDiffTerm(grids, dmap, NumSpec, 1);
    MultiFab ViscousTermforMomentum(grids, dmap, BL_SPACEDIM, 1);
    MultiFab ViscousTermforEnergy(grids, dmap, 1, 1);

    add_temp_diffusion_to_source(source, state, TempDiffTerm, time);

#if (BL_SPACEDIM == 1)
    add_spec_diffusion_to_source(source, state, SpecDiffTerm, time);
    add_viscous_term_to_source(source, state, ViscousTermforMomentum, ViscousTermforEnergy, time);
#endif
}

void
Castro::construct_new_diff_source(MultiFab& source, MultiFab& state_old, MultiFab& state_new, Real time, Real dt)
{
    MultiFab TempDiffTerm(grids, dmap, 1, 1);
    MultiFab SpecDiffTerm(grids, dmap, NumSpec, 1);
    MultiFab ViscousTermforMomentum(grids, dmap, BL_SPACEDIM, 1);
    MultiFab ViscousTermforEnergy(grids, dmap, 1, 1);

    Real mult_factor = 0.5;

    add_temp_diffusion_to_source(source, state_new, TempDiffTerm, time, mult_factor);

#if (BL_SPACEDIM == 1)
    add_spec_diffusion_to_source(source, state_new, SpecDiffTerm, time, mult_factor);
    add_viscous_term_to_source(source, state_new, ViscousTermforMomentum, ViscousTermforEnergy, time, mult_factor);
#endif

    // Time center the source term.

    mult_factor = -0.5;
    Real old_time = time - dt;

    add_temp_diffusion_to_source(source, state_old, TempDiffTerm, old_time, mult_factor);

#if (BL_SPACEDIM == 1)
    add_spec_diffusion_to_source(source, state_old, SpecDiffTerm, old_time, mult_factor);
    add_viscous_term_to_source(source, state_old, ViscousTermforMomentum, ViscousTermforEnergy, old_time, mult_factor);
#endif


}

// **********************************************************************************************

void
Castro::add_temp_diffusion_to_source (MultiFab& ext_src, MultiFab& state, MultiFab& DiffTerm, Real t, Real mult_factor)
{
    // Define an explicit temperature update.
    DiffTerm.setVal(0.);
    if (diffuse_temp == 1) {
        getTempDiffusionTerm(t, state, DiffTerm);
    } else if (diffuse_enth == 1) {
        getEnthDiffusionTerm(t, state, DiffTerm);
    }

    if (diffuse_temp == 1 or diffuse_enth == 1) {
       MultiFab::Saxpy(ext_src,mult_factor,DiffTerm,0,Eden,1,0);
       MultiFab::Saxpy(ext_src,mult_factor,DiffTerm,0,Eint,1,0);
    }
}

// **********************************************************************************************

#if (BL_SPACEDIM == 1)
void
Castro::add_spec_diffusion_to_source (MultiFab& ext_src, MultiFab& state, MultiFab& SpecDiffTerm, Real t, Real mult_factor)
{
    // Define an explicit species update.
    SpecDiffTerm.setVal(0.);
    if (diffuse_spec == 1) {
       getSpecDiffusionTerm(t, state, SpecDiffTerm);
       MultiFab::Saxpy(ext_src,mult_factor,SpecDiffTerm,0,FirstSpec,NumSpec,0);
    }
}
#endif

// **********************************************************************************************

#if (BL_SPACEDIM == 1)
void
Castro::add_viscous_term_to_source(MultiFab& ext_src, MultiFab& state, MultiFab& ViscousTermforMomentum,
                                   MultiFab& ViscousTermforEnergy, Real t, Real mult_factor)
{
    // Define an explicit viscous term
    ViscousTermforMomentum.setVal(0.);
    ViscousTermforEnergy.setVal(0.);
    if (diffuse_vel == 1) {
       getViscousTerm(t,state,ViscousTermforMomentum,ViscousTermforEnergy);
       MultiFab::Saxpy(ext_src,mult_factor,ViscousTermforMomentum,0,Xmom,1,0);
       MultiFab::Saxpy(ext_src,mult_factor,ViscousTermforEnergy  ,0,Eden,1,0);
    }
}
#endif

// **********************************************************************************************

void
Castro::getTempDiffusionTerm (Real time, MultiFab& state, MultiFab& TempDiffTerm)
{
    BL_PROFILE("Castro::getTempDiffusionTerm()");

   if (verbose && ParallelDescriptor::IOProcessor())
      std::cout << "Calculating diffusion term at time " << time << std::endl;

   // Fill coefficients at this level.
   Vector<std::unique_ptr<MultiFab> >coeffs(BL_SPACEDIM);
   Vector<std::unique_ptr<MultiFab> >coeffs_temporary(3); // This is what we pass to the dimension-agnostic Fortran
   for (int dir = 0; dir < 3; dir++) {
       if (dir < BL_SPACEDIM) {
	   coeffs[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, 1, 0));
	   coeffs_temporary[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, 1, 0));
       } else {
	   coeffs_temporary[dir].reset(new MultiFab(grids, dmap, 1, 0));
       }
   }

   // Fill temperature at this level.
   MultiFab Temperature(grids,dmap,1,1);

   {
       FillPatchIterator fpi(*this, state, 1, time, State_Type, 0, NUM_STATE);
       MultiFab& grown_state = fpi.get_mf();

       MultiFab::Copy(Temperature, grown_state, Temp, 0, 1, 1);

       for (MFIter mfi(grown_state); mfi.isValid(); ++mfi)
       {
	   const Box& bx = grids[mfi.index()];

	   ca_fill_temp_cond(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			     BL_TO_FORTRAN_ANYD(grown_state[mfi]),
			     BL_TO_FORTRAN_ANYD((*coeffs_temporary[0])[mfi]),
			     BL_TO_FORTRAN_ANYD((*coeffs_temporary[1])[mfi]),
			     BL_TO_FORTRAN_ANYD((*coeffs_temporary[2])[mfi]));
       }
   }

   // Now copy the temporary array results back to the
   // correctly dimensioned coeffs array.

   for (int dir = 0; dir < BL_SPACEDIM; dir++)
     MultiFab::Copy(*coeffs[dir], *coeffs_temporary[dir], 0, 0, 1, 0);

   MultiFab CrseTemp;
   if (level > 0) {
       // Fill temperature at next coarser level, if it exists.
       const BoxArray& crse_grids = getLevel(level-1).boxArray();
       const DistributionMapping& crse_dmap = getLevel(level-1).DistributionMap();
       CrseTemp.define(crse_grids,crse_dmap,1,1);
       FillPatch(getLevel(level-1),CrseTemp,1,time,State_Type,Temp,1);
   }

   diffusion->applyop(level,Temperature,CrseTemp,TempDiffTerm,coeffs);

   // Extrapolate to ghost cells
   if (TempDiffTerm.nGrow() > 0) {
       for (MFIter mfi(TempDiffTerm); mfi.isValid(); ++mfi)
       {
	   const Box& bx = mfi.validbox();
	   ca_tempdiffextrap(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			     BL_TO_FORTRAN_ANYD(TempDiffTerm[mfi]));
       }
   }
}

void
Castro::getEnthDiffusionTerm (Real time, MultiFab& state, MultiFab& DiffTerm)
{
    BL_PROFILE("Castro::getEnthDiffusionTerm()");

   if (verbose && ParallelDescriptor::IOProcessor())
      std::cout << "Calculating diffusion term at time " << time << std::endl;

   // Fill coefficients at this level.
   Vector<std::unique_ptr<MultiFab> > coeffs(BL_SPACEDIM);
   Vector<std::unique_ptr<MultiFab> > coeffs_temporary(3); // This is what we pass to the dimension-agnostic Fortran
   for (int dir = 0; dir < 3; dir++) {
       if (dir < BL_SPACEDIM) {
	   coeffs[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, 1, 0));
	   coeffs_temporary[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, 1, 0));
       } else {
	   coeffs_temporary[dir].reset(new MultiFab(grids, dmap, 1, 0));
       }
   }

   // Define enthalpy at this level.
   MultiFab Enthalpy(grids,dmap,1,1);
   {
       FillPatchIterator fpi(*this, state, 1, time, State_Type, 0, NUM_STATE);
       const MultiFab& grown_state = fpi.get_mf();

       for (MFIter mfi(grown_state); mfi.isValid(); ++mfi)
       {
	   const Box& bx = grids[mfi.index()];
	   make_enthalpy(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
	                 BL_TO_FORTRAN_ANYD(grown_state[mfi]),
	                 BL_TO_FORTRAN_ANYD(Enthalpy[mfi]));

	   ca_fill_enth_cond(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			     BL_TO_FORTRAN_ANYD(grown_state[mfi]),
			     BL_TO_FORTRAN_ANYD((*coeffs_temporary[0])[mfi]),
			     BL_TO_FORTRAN_ANYD((*coeffs_temporary[1])[mfi]),
			     BL_TO_FORTRAN_ANYD((*coeffs_temporary[2])[mfi]));
       }
   }

   // Now copy the temporary array results back to the
   // correctly dimensioned coeffs array.

   for (int dir = 0; dir < BL_SPACEDIM; dir++)
     MultiFab::Copy(*coeffs[dir], *coeffs_temporary[dir], 0, 0, 1, 0);

   MultiFab CrseEnth, CrseState;
   if (level > 0) {
       // Fill temperature at next coarser level, if it exists.
       const BoxArray& crse_grids = getLevel(level-1).boxArray();
       const DistributionMapping& crse_dmap = getLevel(level-1).DistributionMap();
       CrseEnth.define (crse_grids,crse_dmap,1,1);
       CrseState.define(crse_grids,crse_dmap,NUM_STATE,1);
       FillPatch(getLevel(level-1),CrseState,1,time,State_Type,Density,NUM_STATE);

       for (MFIter mfi(CrseState); mfi.isValid(); ++mfi)
       {
	   const Box& bx = crse_grids[mfi.index()];
	   make_enthalpy(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
	                 BL_TO_FORTRAN_ANYD(CrseState[mfi]),
	                 BL_TO_FORTRAN_ANYD( CrseEnth[mfi]));
       }
   }

   diffusion->applyop(level,Enthalpy,CrseEnth,DiffTerm,coeffs);

   // Extrapolate to ghost cells
   if (DiffTerm.nGrow() > 0) {
       for (MFIter mfi(DiffTerm); mfi.isValid(); ++mfi)
       {
	   const Box& bx = mfi.validbox();
	   ca_tempdiffextrap(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			     BL_TO_FORTRAN_ANYD(DiffTerm[mfi]));
       }
   }
}

#if (BL_SPACEDIM == 1)
void
Castro::getSpecDiffusionTerm (Real time, MultiFab& state, MultiFab& SpecDiffTerm)
{
  BL_PROFILE("Castro::getSpecDiffusionTerm()");

  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "Calculating species diffusion term at time " << time << std::endl;

  // Fill coefficients at this level.
  Vector<std::unique_ptr<MultiFab> > coeffs(BL_SPACEDIM);
  Vector<std::unique_ptr<MultiFab> > coeffs_temporary(3); // This is what we pass to the dimension-agnostic Fortran
   for (int dir = 0; dir < 3; dir++) {
       if (dir < BL_SPACEDIM) {
	   coeffs[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, 1, 0));
	   coeffs_temporary[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, 1, 0));
       } else {
	   coeffs_temporary[dir].reset(new MultiFab(grids, dmap, 1, 0));
       }
   }

   FillPatchIterator fpi(*this, state, 1, time, State_Type, 0, NUM_STATE);
   MultiFab& grown_state = fpi.get_mf();

   for (MFIter mfi(grown_state); mfi.isValid(); ++mfi)
   {
       const Box& bx = grids[mfi.index()];

       ca_fill_spec_coeff(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			  BL_TO_FORTRAN_ANYD(grown_state[mfi]),
			  BL_TO_FORTRAN_ANYD((*coeffs_temporary[0])[mfi]),
			  BL_TO_FORTRAN_ANYD((*coeffs_temporary[1])[mfi]),
			  BL_TO_FORTRAN_ANYD((*coeffs_temporary[2])[mfi]));
   }

   // Now copy the temporary array results back to the
   // correctly dimensioned coeffs array.
   for (int dir = 0; dir < BL_SPACEDIM; dir++)
     MultiFab::Copy(*coeffs[dir], *coeffs_temporary[dir], 0, 0, 1, 0);

   // Create MultiFabs that only hold the data for one species at a time.
   MultiFab Species(grids,dmap,1,1);
   MultiFab     SDT(grids,dmap,SpecDiffTerm.nComp(),SpecDiffTerm.nGrow());
   SDT.setVal(0.);
   MultiFab CrseSpec, CrseDen;
   if (level > 0) {
       const BoxArray& crse_grids = getLevel(level-1).boxArray();
       const DistributionMapping& crse_dmap = getLevel(level-1).DistributionMap();
       CrseSpec.define(crse_grids,crse_dmap,1,1);
       CrseDen.define (crse_grids,crse_dmap,1,1);
   }

   // Fill one species at a time at this level.
   for (int ispec = 0; ispec < NumSpec; ispec++)
   {
       MultiFab::Copy  (Species, grown_state, FirstSpec+ispec, 0, 1, 1);
       MultiFab::Divide(Species, grown_state, Density        , 0, 1, 1);

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
				 BL_TO_FORTRAN_ANYD(SDT[mfi]));
           }
       }
       // Copy back into SpecDiffTerm from the temporary SDT
       MultiFab::Copy(SpecDiffTerm, SDT, 0, ispec, 1, 1);
   }
}
#endif

#if (BL_SPACEDIM == 1)
// **********************************************
/// @note this currently just gets the term that looks like div(2 mu grad(u)) which is
///       only part of the viscous term.  We assume that the coefficient that is filled is "2 mu"
// **********************************************
void
Castro::getViscousTerm (Real time, MultiFab& state, MultiFab& ViscousTermforMomentum, MultiFab& ViscousTermforEnergy)
{
    BL_PROFILE("Castro::getViscousTerm()");

   if (verbose && ParallelDescriptor::IOProcessor())
      std::cout << "Calculating viscous term at time " << time << std::endl;

   getFirstViscousTerm(time,state,ViscousTermforMomentum);

   MultiFab SecndTerm(grids,dmap,ViscousTermforMomentum.nComp(),ViscousTermforMomentum.nGrow());
   getSecndViscousTerm(time,state,SecndTerm);
   MultiFab::Add(ViscousTermforMomentum, SecndTerm, 0, 0, ViscousTermforMomentum.nComp(), 0);

   getViscousTermForEnergy(time,state,ViscousTermforEnergy);
}

void
Castro::getFirstViscousTerm (Real time, MultiFab& state, MultiFab& ViscousTerm)
{

   // Fill coefficients at this level.
   Vector<std::unique_ptr<MultiFab> > coeffs(BL_SPACEDIM);
   Vector<std::unique_ptr<MultiFab> > coeffs_temporary(3); // This is what we pass to the dimension-agnostic Fortran
   for (int dir = 0; dir < 3; dir++) {
       if (dir < BL_SPACEDIM) {
	   coeffs[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, 1, 0));
	   coeffs_temporary[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, 1, 0));
       } else {
	   coeffs_temporary[dir].reset(new MultiFab(grids, dmap, 1, 0));
       }
   }

   // Fill velocity at this level.
   MultiFab Vel(grids,dmap,1,1);

   FillPatchIterator fpi(*this, state, 1, time, State_Type, 0, NUM_STATE);
   MultiFab& grown_state = fpi.get_mf();

   // Remember this is just 1-d
   MultiFab::Copy  (Vel, grown_state, Xmom   , 0, 1, 1);
   MultiFab::Divide(Vel, grown_state, Density, 0, 1, 1);

   for (MFIter mfi(grown_state); mfi.isValid(); ++mfi)
   {
       const Box& bx = grids[mfi.index()];

       ca_fill_first_visc_coeff(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
				BL_TO_FORTRAN_ANYD(grown_state[mfi]),
				BL_TO_FORTRAN_ANYD((*coeffs_temporary[0])[mfi]),
				BL_TO_FORTRAN_ANYD((*coeffs_temporary[1])[mfi]),
				BL_TO_FORTRAN_ANYD((*coeffs_temporary[2])[mfi]));
   }

   // Now copy the temporary array results back to the
   // correctly dimensioned coeffs array.
   for (int dir = 0; dir < BL_SPACEDIM; dir++)
     MultiFab::Copy(*coeffs[dir], *coeffs_temporary[dir], 0, 0, 1, 0);

   MultiFab CrseVel, CrseDen;
   if (level > 0) {
       // Fill temperature at next coarser level, if it exists.
       const BoxArray& crse_grids = getLevel(level-1).boxArray();
       const DistributionMapping& crse_dmap = getLevel(level-1).DistributionMap();
       CrseVel.define(crse_grids,crse_dmap,1,1);
       CrseDen.define(crse_grids,crse_dmap,1,1);
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
			     BL_TO_FORTRAN_ANYD(ViscousTerm[mfi]));
       }
   }
}

void
Castro::getSecndViscousTerm (Real time, MultiFab& state, MultiFab& ViscousTerm)
{

   // Fill coefficients at this level.
   Vector<std::unique_ptr<MultiFab> > coeffs(BL_SPACEDIM);
   Vector<std::unique_ptr<MultiFab> > coeffs_temporary(3); // This is what we pass to the dimension-agnostic Fortran
   for (int dir = 0; dir < 3; dir++) {
       if (dir < BL_SPACEDIM) {
	   coeffs[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, 1, 0));
	   coeffs_temporary[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, 1, 0));
       } else {
	   coeffs_temporary[dir].reset(new MultiFab(grids, dmap, 1, 0));
       }
   }

   // Fill velocity at this level.
   MultiFab Vel(grids,dmap,1,1);

   FillPatchIterator fpi(*this, state, 1, time, State_Type, 0, NUM_STATE);
   MultiFab& grown_state = fpi.get_mf();

   // Remember this is just 1-d
   MultiFab::Copy  (Vel, grown_state, Xmom   , 0, 1, 1);
   MultiFab::Divide(Vel, grown_state, Density, 0, 1, 1);

   for (MFIter mfi(grown_state); mfi.isValid(); ++mfi)
   {
       const Box& bx = grids[mfi.index()];

       ca_fill_secnd_visc_coeff(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
				BL_TO_FORTRAN_ANYD(grown_state[mfi]),
				BL_TO_FORTRAN_ANYD((*coeffs_temporary[0])[mfi]),
				BL_TO_FORTRAN_ANYD((*coeffs_temporary[1])[mfi]),
				BL_TO_FORTRAN_ANYD((*coeffs_temporary[2])[mfi]));
   }

   // Now copy the temporary array results back to the
   // correctly dimensioned coeffs array.
   for (int dir = 0; dir < BL_SPACEDIM; dir++)
     MultiFab::Copy(*coeffs[dir], *coeffs_temporary[dir], 0, 0, 1, 0);

   MultiFab CrseVel, CrseDen;
   if (level > 0) {
       // Fill temperature at next coarser level, if it exists.
       const BoxArray& crse_grids = getLevel(level-1).boxArray();
       const DistributionMapping& crse_dmap = getLevel(level-1).DistributionMap();
       CrseVel.define(crse_grids,crse_dmap,1,1);
       CrseDen.define(crse_grids,crse_dmap,1,1);
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
			     BL_TO_FORTRAN_ANYD(ViscousTerm[mfi]));
       }
   }
}

void
Castro::getViscousTermForEnergy (Real time, MultiFab& state, MultiFab& ViscousTerm)
{

   // Fill coefficients at this level.
   Vector<std::unique_ptr<MultiFab> > coeffs(BL_SPACEDIM);
   Vector<std::unique_ptr<MultiFab> > coeffs_temporary(3); // This is what we pass to the dimension-agnostic Fortran
   for (int dir = 0; dir < 3; dir++) {
       if (dir < BL_SPACEDIM) {
	   coeffs[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, 1, 0));
	   coeffs_temporary[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, 1, 0));
       } else {
	   coeffs_temporary[dir].reset(new MultiFab(grids, dmap, 1, 0));
       }
   }

   FillPatchIterator fpi(*this, state, 2, time, State_Type, 0, NUM_STATE);
   MultiFab& grown_state = fpi.get_mf();

   const Geometry& fine_geom = parent->Geom(parent->finestLevel());
   const Real*       dx_fine = fine_geom.CellSize();

   // Remember this is just 1-d
   int coord_type = Geometry::Coord();
   for (MFIter mfi(grown_state); mfi.isValid(); ++mfi)
   {
       const Box& bx = grids[mfi.index()];

       ca_compute_div_tau_u(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			    BL_TO_FORTRAN_ANYD(ViscousTerm[mfi]),
			    BL_TO_FORTRAN_ANYD(grown_state[mfi]),
			    ZFILL(dx_fine),&coord_type);
   }

   // Extrapolate to ghost cells
   if (ViscousTerm.nGrow() > 0) {
       for (MFIter mfi(ViscousTerm); mfi.isValid(); ++mfi)
       {
	   const Box& bx = mfi.validbox();
	   ca_tempdiffextrap(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			     BL_TO_FORTRAN_ANYD(ViscousTerm[mfi]));
       }
   }
}
#endif
