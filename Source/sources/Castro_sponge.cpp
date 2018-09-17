#ifdef SPONGE
#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::construct_old_sponge_source(MultiFab& source, MultiFab& state, Real time, Real dt)
{

	if (!do_sponge) return;

	update_sponge_params(&time);

	const Real *dx = geom.CellSize();

	const Real mult_factor = 1.0;
#ifndef AMREX_USE_CUDA
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(state, true); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();

		ca_sponge(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		          BL_TO_FORTRAN_ANYD(state[mfi]),
		          BL_TO_FORTRAN_ANYD(source[mfi]),
		          BL_TO_FORTRAN_ANYD(volume[mfi]),
		          ZFILL(dx), dt, time, mult_factor);

	}
#else

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(state, true); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();

#pragma gpu
		ca_sponge(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
		          BL_TO_FORTRAN_ANYD(state[mfi]),
		          BL_TO_FORTRAN_ANYD(source[mfi]),
		          BL_TO_FORTRAN_ANYD(volume[mfi]),
		          AMREX_REAL_ANYD(dx), dt, time, mult_factor);

	}
#endif

}

void
Castro::construct_new_sponge_source(MultiFab& source, MultiFab& state_old, MultiFab& state_new, Real time, Real dt)
{

	if (!do_sponge) return;

	const Real *dx = geom.CellSize();

	const Real mult_factor_old = -0.5;
	const Real mult_factor_new =  0.5;

	// First, subtract half of the old-time source.
	// Note that the sponge parameters are still current
	// at this point from their evaluation at the old time.
#ifndef AMREX_USE_CUDA
#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(state_old, true); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();

		ca_sponge(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		          BL_TO_FORTRAN_ANYD(state_old[mfi]),
		          BL_TO_FORTRAN_ANYD(source[mfi]),
		          BL_TO_FORTRAN_ANYD(volume[mfi]),
		          ZFILL(dx), dt, time, mult_factor_old);

	}

	// Now update to the new-time sponge parameter values
	// and then evaluate the new-time part of the corrector.

	update_sponge_params(&time);

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(state_new, true); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();

		ca_sponge(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		          BL_TO_FORTRAN_ANYD(state_new[mfi]),
		          BL_TO_FORTRAN_ANYD(source[mfi]),
		          BL_TO_FORTRAN_ANYD(volume[mfi]),
		          ZFILL(dx), dt, time, mult_factor_new);

	}

#else

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(state_old, true); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
#pragma gpu
		ca_sponge(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
		          BL_TO_FORTRAN_ANYD(state_old[mfi]),
		          BL_TO_FORTRAN_ANYD(source[mfi]),
		          BL_TO_FORTRAN_ANYD(volume[mfi]),
		          AMREX_REAL_ANYD(dx), dt, time, mult_factor_old);

	}

// Now update to the new-time sponge parameter values
// and then evaluate the new-time part of the corrector.

	update_sponge_params(&time);

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(state_new, true); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.tilebox();
#pragma gpu
		ca_sponge(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
		          BL_TO_FORTRAN_ANYD(state_new[mfi]),
		          BL_TO_FORTRAN_ANYD(source[mfi]),
		          BL_TO_FORTRAN_ANYD(volume[mfi]),
		          AMREX_REAL_ANYD(dx), dt, time, mult_factor_new);

	}

#endif

}

void Castro::sponge_finalize()
{
    ca_deallocate_sponge_params();
}


#endif
