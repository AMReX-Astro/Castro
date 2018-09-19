
#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void Castro::construct_old_rotation_source(MultiFab& source, MultiFab& state, Real time, Real dt)
{
    MultiFab& phirot_old = get_old_data(PhiRot_Type);
    MultiFab& rot_old = get_old_data(Rotation_Type);

    // Fill the rotation data.

    if (!do_rotation) {

	phirot_old.setVal(0.0);
	rot_old.setVal(0.0);

	return;

    }

    fill_rotation_field(phirot_old, rot_old, state, time);

    const Real *dx = geom.CellSize();
    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state, true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();

	ca_rsrc(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		ARLIM_3D(domlo), ARLIM_3D(domhi),
		BL_TO_FORTRAN_ANYD(phirot_old[mfi]),
		BL_TO_FORTRAN_ANYD(rot_old[mfi]),
		BL_TO_FORTRAN_ANYD(state[mfi]),
		BL_TO_FORTRAN_ANYD(source[mfi]),
		BL_TO_FORTRAN_ANYD(volume[mfi]),
		ZFILL(dx),dt,&time);

    }

}



void Castro::construct_new_rotation_source(MultiFab& source, MultiFab& state_old, MultiFab& state_new, Real time, Real dt)
{

    MultiFab& phirot_old = get_old_data(PhiRot_Type);
    MultiFab& rot_old = get_old_data(Rotation_Type);

    MultiFab& phirot_new = get_new_data(PhiRot_Type);
    MultiFab& rot_new = get_new_data(Rotation_Type);

    // Fill the rotation data.

    if (!do_rotation) {

      phirot_new.setVal(0.);
      rot_new.setVal(0.);

      return;

    }

    fill_rotation_field(phirot_new, rot_new, state_new, time);

    // Now do corrector part of rotation source term update

    const Real *dx = geom.CellSize();
    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	for (MFIter mfi(state_new, true); mfi.isValid(); ++mfi)
	{
	    const Box& bx = mfi.tilebox();

	    ca_corrrsrc(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			ARLIM_3D(domlo), ARLIM_3D(domhi),
			BL_TO_FORTRAN_ANYD(phirot_old[mfi]),
			BL_TO_FORTRAN_ANYD(phirot_new[mfi]),
			BL_TO_FORTRAN_ANYD(rot_old[mfi]),
			BL_TO_FORTRAN_ANYD(rot_new[mfi]),
			BL_TO_FORTRAN_ANYD(state_old[mfi]),
			BL_TO_FORTRAN_ANYD(state_new[mfi]),
			BL_TO_FORTRAN_ANYD(source[mfi]),
			BL_TO_FORTRAN_ANYD((*mass_fluxes[0])[mfi]),
			BL_TO_FORTRAN_ANYD((*mass_fluxes[1])[mfi]),
			BL_TO_FORTRAN_ANYD((*mass_fluxes[2])[mfi]),
			ZFILL(dx),dt,&time,
			BL_TO_FORTRAN_ANYD(volume[mfi]));
	}
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
				   BL_TO_FORTRAN_ANYD(phi[mfi]),
				   ZFILL(dx),time);

    }

    rot.setVal(0.0);

    ng = state.nGrow();

    if (ng > rot.nGrow())
      amrex::Error("State MF has more ghost cells than rotation MF.");
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state, true); mfi.isValid(); ++mfi)
    {

      const Box& bx = mfi.growntilebox(ng);

      ca_fill_rotational_acceleration(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()), 
				      BL_TO_FORTRAN_ANYD(rot[mfi]),
				      BL_TO_FORTRAN_ANYD(state[mfi]),
				      ZFILL(dx),time);

    }

}


