/* Implementations of functions in Problem.H go here */

#include "Castro.H"
#include "Castro_F.H"
#include "Problem_F.H"

void
Castro::flame_width_properties (Real time, Real& T_max, Real& T_min, Real& grad_T_max)
{
    BL_PROFILE("Castro::flame_width_properties()");

    const Real* dx = geom.CellSize();

    MultiFab* mf = derive("Temp",time,1);

    BL_ASSERT(mf != 0);

#ifdef _OPENMP
#pragma omp parallel reduction(max:T_max,grad_T_max) reduction(min:T_min)
#endif    
    for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

	BL_FORT_PROC_CALL(FLAME_WIDTH_TEMP,flame_width_temp)
            (BL_TO_FORTRAN_3D(fab),
	     ARLIM_3D(lo),ARLIM_3D(hi),
	     ZFILL(dx),&time,
	     &T_max, &T_min, &grad_T_max);
    }

    delete mf;

}



void
Castro::flame_speed_properties (Real old_time, Real new_time, Real& rho_fuel_dot_old, Real& rho_fuel_dot_new)
{
    BL_PROFILE("Castro::flame_speed_properties()");

    const Real* dx = geom.CellSize();

    MultiFab* mf_old = derive("omegadot_fuel",old_time,0);
    MultiFab* mf_new = derive("omegadot_fuel",new_time,0);

    BL_ASSERT(mf_old != 0);
    BL_ASSERT(mf_new != 0);

    Real rho_fuel_dot_old_temp = 0.0;
    Real rho_fuel_dot_new_temp = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:rho_fuel_dot_old_temp,rho_fuel_dot_new_temp)
#endif    
    for (MFIter mfi(*mf_new,true); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab_old = (*mf_old)[mfi];
	FArrayBox& fab_new = (*mf_new)[mfi];

        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

	BL_FORT_PROC_CALL(FLAME_SPEED_DATA,flame_speed_data)
            (BL_TO_FORTRAN_3D(fab_old),
	     BL_TO_FORTRAN_3D(fab_new),
	     ARLIM_3D(lo),ARLIM_3D(hi),
	     ZFILL(dx),
	     &rho_fuel_dot_old_temp,&rho_fuel_dot_new_temp);
    }

    rho_fuel_dot_old += rho_fuel_dot_old_temp;
    rho_fuel_dot_new += rho_fuel_dot_new_temp;
    
    delete mf_old;
    delete mf_new;

}
