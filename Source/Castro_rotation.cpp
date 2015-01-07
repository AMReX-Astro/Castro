#include <winstd.H>

#include "Castro.H"
#include "Castro_F.H"

using std::string;

#ifdef ROTATION

void Castro::add_rotation_to_old_source(MultiFab& ext_src_old, MultiFab& OldRotationTerms, Real old_time)
{
  const Real* dx = geom.CellSize();
  const Real* prob_lo = geom.ProbLo();
  
  MultiFab& S_old = get_old_data(State_Type);
  const int ncomp = S_old.nComp();

  OldRotationTerms.setVal(0.);

  if (do_rotation == 1) {
    {
	FillPatchIterator old_fpi(*this,S_old,NUM_GROW,old_time,State_Type,Density,ncomp);
	MultiFab& old_state = old_fpi.get_mf();
#ifdef _OPENMP
#pragma omp parallel	  
#endif
	for (MFIter mfi(OldRotationTerms,true); mfi.isValid(); ++mfi)
	{
	    const Box& bx = mfi.growntilebox();
	    
	    BL_FORT_PROC_CALL(CA_ROTATE,ca_rotate)
		(bx.loVect(), bx.hiVect(),
		 BL_TO_FORTRAN(old_state[mfi]),
		 BL_TO_FORTRAN(OldRotationTerms[mfi]),
		 prob_lo, dx);
	}
    }

    // Add the source terms to ext_src_old
    MultiFab::Add(ext_src_old,OldRotationTerms,0,Xmom,BL_SPACEDIM,1);
    MultiFab::Add(ext_src_old,OldRotationTerms,BL_SPACEDIM,Eden,1,1);
  }
}

void Castro::time_center_rotation(MultiFab& S_new, MultiFab& OldRotationTerms, Real cur_time, Real dt)
{
  // Correct the momentum and energy update with time-centered rotation terms
  MultiFab NewRotationTerms(grids,BL_SPACEDIM+1,1);
  NewRotationTerms.setVal(0.);

  const Real* dx = geom.CellSize();
  const Real* prob_lo = geom.ProbLo();
  const int ncomp = S_new.nComp();

  if (do_rotation == 1) {
    {
	FillPatchIterator New_fpi(*this,S_new,NUM_GROW,cur_time,State_Type,Density,ncomp);
	MultiFab& new_state = New_fpi.get_mf();
#ifdef _OPENMP
#pragma omp parallel	  
#endif
	for (MFIter mfi(NewRotationTerms,true); mfi.isValid(); ++mfi)
	{
	    const Box& bx = mfi.tilebox();
	
	    BL_FORT_PROC_CALL(CA_ROTATE,ca_rotate)
		(bx.loVect(), bx.hiVect(),
		 BL_TO_FORTRAN(new_state[mfi]),
		 BL_TO_FORTRAN(NewRotationTerms[mfi]),
		 prob_lo, dx);
	}
    }

    NewRotationTerms.mult( 0.5*dt);
    OldRotationTerms.mult(-0.5*dt);

    MultiFab::Add(S_new,OldRotationTerms,0,Xmom,BL_SPACEDIM,0);
    MultiFab::Add(S_new,OldRotationTerms,BL_SPACEDIM,Eden,1,0);
    MultiFab::Add(S_new,NewRotationTerms,0,Xmom,BL_SPACEDIM,0);
    MultiFab::Add(S_new,NewRotationTerms,BL_SPACEDIM,Eden,1,0);
  }
}


#endif
