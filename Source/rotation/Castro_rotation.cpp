
#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void Castro::construct_old_rotation_source(Real time, Real dt)
{
    MultiFab& phirot_old = get_old_data(PhiRot_Type);
    MultiFab& rot_old = get_old_data(Rotation_Type);

    MultiFab& old_sources = get_old_data(Source_Type);

    // Fill the rotation data.

    if (!do_rotation) {

	phirot_old.setVal(0.0);
	rot_old.setVal(0.0);

	return;

    }

    fill_rotation_field(phirot_old, rot_old, Sborder, time);

    const Real *dx = geom.CellSize();
    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(Sborder,true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();

	ca_rsrc(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		ARLIM_3D(domlo), ARLIM_3D(domhi),
		BL_TO_FORTRAN_3D(phirot_old[mfi]),
		BL_TO_FORTRAN_3D(rot_old[mfi]),
		BL_TO_FORTRAN_3D(Sborder[mfi]),
		BL_TO_FORTRAN_3D(old_sources[mfi]),
		BL_TO_FORTRAN_3D(volume[mfi]),
		ZFILL(dx),dt,&time);

    }

}



void Castro::construct_new_rotation_source(Real time, Real dt)
{
    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    MultiFab& phirot_old = get_old_data(PhiRot_Type);
    MultiFab& rot_old = get_old_data(Rotation_Type);

    MultiFab& phirot_new = get_new_data(PhiRot_Type);
    MultiFab& rot_new = get_new_data(Rotation_Type);

    MultiFab& new_sources = get_new_data(Source_Type);

    // Fill the rotation data.

    if (!do_rotation) {

      phirot_new.setVal(0.);
      rot_new.setVal(0.);

      return;

    }

    fill_rotation_field(phirot_new, rot_new, S_new, time);

    // Now do corrector part of rotation source term update

    const Real *dx = geom.CellSize();
    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
	{
	    const Box& bx = mfi.tilebox();

	    ca_corrrsrc(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			ARLIM_3D(domlo), ARLIM_3D(domhi),
			BL_TO_FORTRAN_3D(phirot_old[mfi]),
			BL_TO_FORTRAN_3D(phirot_new[mfi]),
			BL_TO_FORTRAN_3D(rot_old[mfi]),
			BL_TO_FORTRAN_3D(rot_new[mfi]),
			BL_TO_FORTRAN_3D(S_old[mfi]),
			BL_TO_FORTRAN_3D(S_new[mfi]),
			BL_TO_FORTRAN_3D(new_sources[mfi]),
			BL_TO_FORTRAN_3D((*fluxes[0])[mfi]),
			BL_TO_FORTRAN_3D((*fluxes[1])[mfi]),
			BL_TO_FORTRAN_3D((*fluxes[2])[mfi]),
			ZFILL(dx),dt,&time,
			BL_TO_FORTRAN_3D(volume[mfi]));
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
				   BL_TO_FORTRAN_3D(phi[mfi]),
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
				      BL_TO_FORTRAN_3D(rot[mfi]),
				      BL_TO_FORTRAN_3D(state[mfi]),
				      ZFILL(dx),time);

    }

}


