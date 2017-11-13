#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::construct_old_ext_source(Real time, Real dt)
{
    int ng = Sborder.nGrow();

    if (!add_ext_src) return;

    MultiFab src(grids, dmap, NUM_STATE, ng);

    src.setVal(0.0);

    fill_ext_source(time, dt, Sborder, Sborder, src, ng);

    Real mult_factor = 1.0;

    MultiFab::Saxpy(old_sources, mult_factor, src, 0, 0, NUM_STATE, ng);
}



void
Castro::construct_new_ext_source(Real time, Real dt)
{
    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    int ng = 0;

    if (!add_ext_src) return;

    MultiFab src(grids, dmap, NUM_STATE, ng);

    src.setVal(0.0);

    // Subtract off the old-time value first.

    Real old_time = time - dt;
    Real mult_factor = -0.5;

    fill_ext_source(old_time, dt, S_old, S_old, src, ng);

    MultiFab::Saxpy(new_sources, mult_factor, src, 0, 0, NUM_STATE, ng);

    // Time center with the new data.

    src.setVal(0.0);

    mult_factor = 0.5;

    fill_ext_source(time, dt, S_old, S_new, src, ng);

    MultiFab::Saxpy(new_sources, mult_factor, src, 0, 0, NUM_STATE, ng);

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
