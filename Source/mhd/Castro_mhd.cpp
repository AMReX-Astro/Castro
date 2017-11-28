#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::just_the_mhd(Real time, Real dt)
{
      if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... mhd ...!!! " << std::endl << std::endl;

      hydro_source.setVal(0.0);

      AmrLevel::FillPatchAdd(*this, sources_for_hydro, NUM_GROW, time, Source_Type, 0, NUM_STATE);

      const int finest_level = parent->finestLevel();

      const Real *dx = geom.CellSize();
      Real courno = -1.0e+200;

      MultiFab& S_new = get_new_data(State_Type);
      MultiFab& Bx_old= get_old_data(Mag_Type_x);
      MultiFab& By_old= get_old_data(Mag_Type_y);
      MultiFab& Bz_old= get_old_data(Mag_Type_z);
      MultiFab& Bx_new= get_new_data(Mag_Type_x);
      MultiFab& By_new= get_new_data(Mag_Type_y);
      MultiFab& Bz_new= get_new_data(Mag_Type_z);


      MultiFab grav_vector(grids, dmap, BL_SPACEDIM, 3);
      grav_vector.setVal(0.);
        

#ifdef _OPENMP
#pragma omp parallel reduction(+:mass:courno)
#endif
    { 
        FArrayBox flux[BL_SPACEDIM], u_gdnv[BL_SPACEDIM], E[BL_SPACEDIM];

        FArrayBox q, qaux, src_q;

        int priv_nstep_fsp = -1;

        Real cflLoc = -1.0e+200;
        int is_finest_level = (level == finest_level) ? 1 : 0;
        const int*  domain_lo = geom.Domain().loVect();
        const int*  domain_hi = geom.Domain().hiVect();

        for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi)
	{
	   const Box& bx = mfi.tilebox();
	   const Box& qbx = amrex::grow(bx, NUM_GROW);

           const int* lo = bx.loVect();
           const int* hi = bx.hiVect();

           FArrayBox &statein  = Sborder[mfi];
           FArrayBox &stateout = S_new[mfi];

           FArrayBox &source_in  = sources_for_hydro[mfi];
           FArrayBox &source_out = hydro_source[mfi];
           
	   q.resize(qbx, QVAR);
	   
	   qaux.resize(qbx, NQAUX);
	   src_q.resize(qbx, QVAR);

           FArrayBox& Bx  = Bx_old[mfi];
           FArrayBox& By  = By_old[mfi]; 
	   FArrayBox& Bz  = Bz_old[mfi];

	   FArrayBox& Bxout = Bx_new[mfi];
	   FArrayBox& Byout = By_new[mfi];
	   FArrayBox& Bzout = Bz_new[mfi];

           Real se  = 0;
	   Real ske = 0; 

           ca_advance_mhd
            (&time, bx.loVect(), bx.hiVect(),
             BL_TO_FORTRAN(statein),
             BL_TO_FORTRAN(stateout),
             BL_TO_FORTRAN(Bx),
             BL_TO_FORTRAN(By),
             BL_TO_FORTRAN(Bz),
             BL_TO_FORTRAN(Bxout),
             BL_TO_FORTRAN(Byout),
             BL_TO_FORTRAN(Bzout),
             BL_TO_FORTRAN(u_gdnv[0]),
             BL_TO_FORTRAN(u_gdnv[1]),
             BL_TO_FORTRAN(u_gdnv[2]),
             BL_TO_FORTRAN(sources_for_hydro[mfi]),
             BL_TO_FORTRAN(grav_vector[mfi]),
             dx, &dt,
             BL_TO_FORTRAN(flux[0]),
             BL_TO_FORTRAN(flux[1]),
             BL_TO_FORTRAN(flux[2]),
             BL_TO_FORTRAN(E[0]),
             BL_TO_FORTRAN(E[1]),
             BL_TO_FORTRAN(E[2]),
             &cflLoc, &se, &ske, &print_fortran_warnings);


	}

    } 



}
