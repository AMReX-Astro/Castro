#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::just_the_mhd(Real time, Real dt)
{
      if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... mhd ...!!! " << std::endl << std::endl;

      hydro_source.setVal(0.0);

      const int finest_level = parent->finestLevel();

      const Real *dx = geom.CellSize();
      Real courno = -1.0e+200;

      const int*  domain_lo = geom.Domain().loVect();
      const int*  domain_hi = geom.Domain().hiVect();


      MultiFab& S_new = get_new_data(State_Type);
      MultiFab& Bx_new= get_new_data(Mag_Type_x);
      MultiFab& By_new= get_new_data(Mag_Type_y);
      MultiFab& Bz_new= get_new_data(Mag_Type_z);


      MultiFab electric[BL_SPACEDIM];
      for (int j = 0; j < BL_SPACEDIM; j++)
      {
        electric[j].define(getEdgeBoxArray(j), dmap, 1, 0);
        electric[j].setVal(0.0);
      }


      BL_ASSERT(NUM_GROW == 4);


#ifdef _OPENMP
#pragma omp parallel reduction(+:mass:courno)
#endif
    { 
        FArrayBox flux[BL_SPACEDIM], E[BL_SPACEDIM];


        int priv_nstep_fsp = -1;

        Real cflLoc = -1.0e+200;
        int is_finest_level = (level == finest_level) ? 1 : 0;
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
           const Box& bx = mfi.tilebox();

           const int* lo = bx.loVect();
           const int* hi = bx.hiVect();

           FArrayBox &statein  = Sborder[mfi];
           FArrayBox &stateout = S_new[mfi];

           FArrayBox &source_in  = sources_for_hydro[mfi];
           FArrayBox &source_out = hydro_source[mfi]; 

           FArrayBox& Bx  = Bx_old_tmp[mfi];
           FArrayBox& By  = By_old_tmp[mfi]; 
           FArrayBox& Bz  = Bz_old_tmp[mfi];

           FArrayBox& Bxout = Bx_new[mfi];
           FArrayBox& Byout = By_new[mfi];
           FArrayBox& Bzout = Bz_new[mfi];

           Real se  = 0;
           Real ske = 0;

           //Allocate fabs for fluxes
           for (int i = 0; i < BL_SPACEDIM; i++){
              const Box& bxtmp = amrex::surroundingNodes(bx,i);
              flux[i].resize(bxtmp,NUM_STATE);
              E[i].resize(bxtmp,NUM_STATE);
              
           }


           ca_advance_mhd
            (&time, bx.loVect(), bx.hiVect(),
             BL_TO_FORTRAN_3D(statein),
             BL_TO_FORTRAN_3D(source_out),
             BL_TO_FORTRAN_3D(Bx),
             BL_TO_FORTRAN_3D(By),
             BL_TO_FORTRAN_3D(Bz),
             BL_TO_FORTRAN_3D(Bxout),
             BL_TO_FORTRAN_3D(Byout),
             BL_TO_FORTRAN_3D(Bzout),
             BL_TO_FORTRAN_3D(sources_for_hydro[mfi]),
             dx, &dt,
             D_DECL(BL_TO_FORTRAN_3D(flux[0]),
             BL_TO_FORTRAN_3D(flux[1]),
             BL_TO_FORTRAN_3D(flux[2])),
             BL_TO_FORTRAN_3D(E[0]),
             BL_TO_FORTRAN_3D(E[1]),
             BL_TO_FORTRAN_3D(E[2]),
             &cflLoc, &se, &ske, &print_fortran_warnings);

           for (int i = 0; i < BL_SPACEDIM; i++){
             (*fluxes[i])[mfi].plus(flux[i], mfi.nodaltilebox(i),0,0,NUM_STATE);
             
             (*mass_fluxes[i])[mfi].copy(flux[i],mfi.nodaltilebox(i),Density,mfi.nodaltilebox(i),0,1);
             electric[i][mfi].copy(E[i], mfi.nodaltilebox(i));
           }


        }

    } 



}
