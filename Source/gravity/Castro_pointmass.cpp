#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::pointmass_update(Real time, Real dt)
{
    int finest_level = parent->finestLevel();

    if (level == finest_level && point_mass_fix_solution)
    {

        MultiFab& S_old = get_old_data(State_Type);
	MultiFab& S_new = get_new_data(State_Type);

        Real mass_change_at_center = 0.0;

	const Real* dx = geom.CellSize();

#ifdef _OPENMP
#pragma omp parallel reduction(+:mass_change_at_center)
#endif
        for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {

	    const Box& bx = mfi.tilebox();

	    pm_compute_delta_mass(&mass_change_at_center,
				  ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
				  BL_TO_FORTRAN_3D(S_old[mfi]),
				  BL_TO_FORTRAN_3D(S_new[mfi]),
				  BL_TO_FORTRAN_3D(volume[mfi]),
				  ZFILL(geom.ProbLo()), ZFILL(dx),
				  &time, &dt);

	}

	ParallelDescriptor::ReduceRealSum(mass_change_at_center);

	if (mass_change_at_center > 0.0)
        {
	    point_mass += mass_change_at_center;

	    set_pointmass(&point_mass);

#ifdef _OPENMP
#pragma omp parallel
#endif
	    for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();

		pm_fix_solution(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
				BL_TO_FORTRAN_3D(S_old[mfi]), BL_TO_FORTRAN_3D(S_new[mfi]),
				ZFILL(geom.ProbLo()), ZFILL(dx), &time, &dt);
             }
          }
    }
}
