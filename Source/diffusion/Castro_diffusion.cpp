
#include "Castro.H"
#include "Castro_F.H"

using std::string;

#include "Diffusion.H"

using namespace amrex;

void
Castro::construct_old_diff_source(MultiFab& source, MultiFab& state, Real time, Real dt)
{
    MultiFab TempDiffTerm(grids, dmap, 1, 0);

    add_temp_diffusion_to_source(source, state, TempDiffTerm, time);

}

void
Castro::construct_new_diff_source(MultiFab& source, MultiFab& state_old, MultiFab& state_new, Real time, Real dt)
{
    MultiFab TempDiffTerm(grids, dmap, 1, 0);

    Real mult_factor = 0.5;

    add_temp_diffusion_to_source(source, state_new, TempDiffTerm, time, mult_factor);

    // Time center the source term.

    mult_factor = -0.5;
    Real old_time = time - dt;

    add_temp_diffusion_to_source(source, state_old, TempDiffTerm, old_time, mult_factor);

}

// **********************************************************************************************

void
Castro::add_temp_diffusion_to_source (MultiFab& ext_src, MultiFab& state, MultiFab& DiffTerm, Real t, Real mult_factor)
{
    // Define an explicit temperature update.
    DiffTerm.setVal(0.);
    if (diffuse_temp == 1) {
        getTempDiffusionTerm(t, state, DiffTerm);
    }

    if (diffuse_temp == 1) {
       MultiFab::Saxpy(ext_src,mult_factor,DiffTerm,0,Eden,1,0);
       MultiFab::Saxpy(ext_src,mult_factor,DiffTerm,0,Eint,1,0);
    }
}


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

}
