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
    for (FillPatchIterator Old_fpi(*this,S_old,4,old_time,State_Type,Density,ncomp);
	 Old_fpi.isValid();++Old_fpi)
      {
	const Box& bx = grids[Old_fpi.index()];

	BL_FORT_PROC_CALL(CA_ROTATE,ca_rotate)
	  (bx.loVect(), bx.hiVect(),
	   BL_TO_FORTRAN(Old_fpi()),
	   BL_TO_FORTRAN(OldRotationTerms[Old_fpi]),
	   prob_lo, dx);
      }

    // Add the source terms to ext_src_old
    MultiFab::Add(ext_src_old,OldRotationTerms,0,Xmom,BL_SPACEDIM,0);
    MultiFab::Add(ext_src_old,OldRotationTerms,BL_SPACEDIM,Eden,1,0);

    geom.FillPeriodicBoundary(ext_src_old,0,NUM_STATE);
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
    for (FillPatchIterator New_fpi(*this,S_new,4,cur_time,State_Type,Density,ncomp);
	 New_fpi.isValid();++New_fpi)
      {
	const Box& bx = grids[New_fpi.index()];
	
	BL_FORT_PROC_CALL(CA_ROTATE,ca_rotate)
	  (bx.loVect(), bx.hiVect(),
	   BL_TO_FORTRAN(New_fpi()),
	   BL_TO_FORTRAN(NewRotationTerms[New_fpi]),
	   prob_lo, dx);
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
