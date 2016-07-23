#include <winstd.H>

#include "Castro.H"
#include "Castro_F.H"

void Castro::construct_old_rotation(int amr_iteration, int amr_ncycle,
				    int sub_iteration, int sub_ncycle,
				    Real time, MultiFab& state)
{

    MultiFab& phirot_old = get_old_data(PhiRot_Type);
    MultiFab& rot_old = get_old_data(Rotation_Type);

    // Fill the old rotation data.

    if (do_rotation) {

      fill_rotation_field(phirot_old, rot_old, state, time);

    } else {

      phirot_old.setVal(0.);
      rot_old.setVal(0.);

    }

}



void Castro::construct_new_rotation(int amr_iteration, int amr_ncycle,
				    int sub_iteration, int sub_ncycle,
				    Real time, MultiFab& state)
{

    MultiFab& phirot_new = get_new_data(PhiRot_Type);
    MultiFab& rot_new = get_new_data(Rotation_Type);

    // Fill the old rotation data.

    if (do_rotation) {

      fill_rotation_field(phirot_new, rot_new, state, time);

    } else {

      phirot_new.setVal(0.);
      rot_new.setVal(0.);

    }

}



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

      ca_fill_rotational_potential(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()), 
				   BL_TO_FORTRAN_3D(phi[mfi]),
				   ZFILL(dx),time);

    }

    rot.setVal(0.0);

    ng = state.nGrow();

    if (ng > rot.nGrow())
      BoxLib::Error("State MF has more ghost cells than rotation MF.");
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state, true); mfi.isValid(); ++mfi)
    {

      const Box& bx = mfi.growntilebox(ng);

      ca_fill_rotational_acceleration(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()), 
				      BL_TO_FORTRAN_3D(rot[mfi]),
				      BL_TO_FORTRAN_3D(state[mfi]),
				      ZFILL(dx),time);

    }

}
