#include <winstd.H>

#include "Castro.H"
#include "Castro_F.H"

void Castro::fill_rotation_field(MultiFab& phi, MultiFab& rot, MultiFab& state, Real time)
{
    const Real* dx = geom.CellSize();

    phi.setVal(0.0);    

    int ng = phi.nGrow();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(phi, true); mfi.isValid(); ++mfi)
    {

      const Box& bx = mfi.growntilebox(ng);

      BL_FORT_PROC_CALL(CA_FILL_ROTATIONAL_POTENTIAL,ca_fill_rotational_potential)
		(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()), 
		 BL_TO_FORTRAN_3D(phi[mfi]),
		 ZFILL(dx),time);

    }

    rot.setVal(0.0);

    ng = 0;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(rot, true); mfi.isValid(); ++mfi)
    {

      const Box& bx = mfi.growntilebox(ng);

      BL_FORT_PROC_CALL(CA_FILL_ROTATIONAL_ACCELERATION,ca_fill_rotational_acceleration)
		(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()), 
		 BL_TO_FORTRAN_3D(rot[mfi]),
		 BL_TO_FORTRAN_3D(state[mfi]),
		 ZFILL(dx),time);

    }

}
