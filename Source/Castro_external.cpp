#include "Castro.H"
#include "Castro_F.H"

void
Castro::construct_old_ext_source(PArray<MultiFab>& old_sources, MultiFab& sources_for_hydro, MultiFab& S_old, Real time, Real dt)
{
    int ng = S_old.nGrow();

    old_sources.set(ext_src, new MultiFab(grids,NUM_STATE,ng));
    old_sources[ext_src].setVal(0.0,ng);

    if (add_ext_src) {
      fill_ext_source(time, dt, S_old, S_old, old_sources[ext_src], ng);
      BoxLib::fill_boundary(old_sources[ext_src], geom);
      MultiFab::Add(sources_for_hydro,old_sources[ext_src],0,0,NUM_STATE,ng);
    }
}



void
Castro::construct_new_ext_source(PArray<MultiFab>& old_sources, PArray<MultiFab>& new_sources, MultiFab& sources_for_hydro,
				 MultiFab& S_old, MultiFab& S_new, Real time, Real dt)
{
    int ng = S_new.nGrow();

    new_sources.set(ext_src, new MultiFab(grids,NUM_STATE,ng));
    new_sources[ext_src].setVal(0.0,ng);

    fill_ext_source(time, dt, S_old, S_new, new_sources[ext_src], ng);

    // Time center the source term.

    old_sources[ext_src].mult(-0.5);
    new_sources[ext_src].mult( 0.5);

    MultiFab::Add(new_sources[ext_src],old_sources[ext_src],0,0,NUM_STATE,0);

    old_sources[ext_src].mult(-2.0);

}



void
Castro::fill_ext_source (Real time, Real dt, MultiFab& state_old, MultiFab& state_new, MultiFab& ext_src, int ng)
{
    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(ext_src,true); mfi.isValid(); ++mfi)
    {

        const Box& bx = mfi.growntilebox(ng);

#ifdef DIMENSION_AGNOSTIC
        BL_FORT_PROC_CALL(CA_EXT_SRC,ca_ext_src)
	  (ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
	   BL_TO_FORTRAN_3D(state_old[mfi]),
	   BL_TO_FORTRAN_3D(state_new[mfi]),
	   BL_TO_FORTRAN_3D(ext_src[mfi]),
	   ZFILL(prob_lo),ZFILL(dx),&time,&dt);
#else
	BL_FORT_PROC_CALL(CA_EXT_SRC,ca_ext_src)
	  (bx.loVect(), bx.hiVect(),
	   BL_TO_FORTRAN(state_old[mfi]),
	   BL_TO_FORTRAN(state_new[mfi]),
	   BL_TO_FORTRAN(ext_src[mfi]),
	   prob_lo,dx,&time,&dt);
#endif
    }
}
