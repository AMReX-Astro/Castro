#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::construct_old_ext_source(MultiFab& source, MultiFab& state, Real time, Real dt)
{
    if (!add_ext_src) return;

    MultiFab ext_src(grids, dmap, NSRC, 0);

    ext_src.setVal(0.0);

    fill_ext_source(time, dt, state, state, ext_src);

    Real mult_factor = 1.0;

    MultiFab::Saxpy(source, mult_factor, ext_src, 0, 0, NSRC, 0);
}



void
Castro::construct_new_ext_source(MultiFab& source, MultiFab& state_old, MultiFab& state_new, Real time, Real dt)
{
    if (!add_ext_src) return;

    MultiFab ext_src(grids, dmap, NSRC, 0);

    ext_src.setVal(0.0);

    // Subtract off the old-time value first.

    Real old_time = time - dt;
    Real mult_factor = -0.5;

    fill_ext_source(old_time, dt, state_old, state_old, ext_src);

    MultiFab::Saxpy(source, mult_factor, ext_src, 0, 0, NSRC, 0);

    // Time center with the new data.

    ext_src.setVal(0.0);

    mult_factor = 0.5;

    fill_ext_source(time, dt, state_old, state_new, ext_src);

    MultiFab::Saxpy(source, mult_factor, ext_src, 0, 0, NSRC, 0);

}



void
Castro::fill_ext_source (Real time, Real dt, MultiFab& state_old, MultiFab& state_new, MultiFab& ext_src)
{
    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(ext_src,true); mfi.isValid(); ++mfi)
    {

        const Box& bx = mfi.tilebox();

#ifdef AMREX_DIMENSION_AGNOSTIC
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
