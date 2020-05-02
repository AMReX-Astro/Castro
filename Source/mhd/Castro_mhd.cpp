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


      //MultiFab electric[BL_SPACEDIM];
      //for (int j = 0; j < BL_SPACEDIM; j++)
      //{
      //  electric[j].define(getEdgeBoxArray(j), dmap, 1, 0);
      //  electric[j].setVal(0.0);
      //}


      BL_ASSERT(NUM_GROW == 4);


#ifdef _OPENMP
#pragma omp parallel reduction(+:mass:courno)
#endif
    {

      FArrayBox flux[AMREX_SPACEDIM], E[AMREX_SPACEDIM];
      FArrayBox cs[AMREX_SPACEDIM];

      FArrayBox bcc;
      FArrayBox q;
      FArrayBox srcQ;

      FArrayBox qm;
      FArrayBox qp;

      FArrayBox flxx;
      FArrayBox flxy;
      FArrayBox flxz;

      FArrayBox Extmp;
      FArrayBox Eytmp;
      FArrayBox Eztmp;

      for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {

          const Box& bx = mfi.tilebox();

          const Box& obx = amrex::grow(bx, 1);

          // box with NUM_GROW ghost cells for PPM stuff
          const Box& bx_gc = amrex::grow(bx, NUM_GROW);

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


          // allocate fabs for fluxes
          for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
            const Box& bxtmp = amrex::surroundingNodes(bx, i);
            flux[i].resize(bxtmp, NUM_STATE);
            E[i].resize(bxtmp, NUM_STATE);
          }

          // Calculate primitives based on conservatives
          bcc.resize(bc_gc, 3);
          q.resize(bc_gc, NQ);
          srcQ.resize(obx, NQSRC);

          for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
            cs[idir].resize(bc_gc, 1);
          }

          const int* lo_gc = bx_gc.loVect();
          const int* hi_gc = bx_gc.hiVect();

          ctoprim_mhd(lo_gc, hi_gc,
                      BL_TO_FORTRAN_ANYD(statein),
                      BL_TO_FORTRAN_ANYD(bcc),
                      BL_TO_FORTRAN_ANYD(Bx),
                      BL_TO_FORTRAN_ANYD(By),
                      BL_TO_FORTRAN_ANYD(Bz),
                      BL_TO_FORTRAN_ANYD(q),
                      BL_TO_FORTRAN_ANYD(cs[0]),
                      BL_TO_FORTRAN_ANYD(cs[1]),
                      BL_TO_FORTRAN_ANYD(cs[2]),
                      BL_TO_FORTRAN_ANYD(source_in),
                      BL_TO_FORTRAN_ANYD(srcQ));

          const int* lo1 = obx.loVect();
          const int* hi1 = obx.hiVect();

          srctoprim_mhd(lo1, hi1,
                        BL_TO_FORTRAN_ANYD(q),
                        BL_TO_FORTRAN_ANYD(source_in),
                        BL_TO_FORTRAN_ANYD(srcQ));

          check_for_mhd_cfl_violation(lo, hi,
                                      BL_TO_FORTRAN_ANYD(q),
                                      BL_TO_FORTRAN_ANYD(cs[0]),
                                      BL_TO_FORTRAN_ANYD(cs[1]),
                                      BL_TO_FORTRAN_ANYD(cs[2]));

          flatn.resize(obx, 1);
          auto flatn_arr = flatn.array();
          amrex::ParallelFor(obx,
          [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            flatn_arr(i,j,k) = 0.0;
          });

          mhd_flatten(lo1, hi1,
                      BL_TO_FORTRAN_ANYD(q),
                      BL_TO_FORTRAN_ANYD(flatn));


          // Interpolate Cell centered values to faces
          qp.resize(bx_gc, NQ * AMREX_SPACEDIM);
          qm.resize(bx_gc, NQ * AMREX_SPACEDIM);

          plm(lo, hi,
              BL_TO_FORTRAN_ANYD(q)
              BL_TO_FORTRAN_ANYD(flatn),
              BL_TO_FORTRAN_ANYD(Bx),
              BL_TO_FORTRAN_ANYD(By),
              BL_TO_FORTRAN_ANYD(Bz),
              BL_TO_FORTRAN_ANYD(qp),
              BL_TO_FORTRAN_ANYD(qm),
              BL_TO_FORTRAN_ANYD(srcQ),
              dx, dt);


          // Corner Couple and find the correct fluxes + electric fields

          // need to revisit these box sizes
          nbx = amrex::surroundingNodes(bx, 0);
          nbxf = amrex::grow(nbx, IntVect(2, 3, 3));

          nby = amrex::surroundingNodes(bx, 1);
          nbyf = amrex::grow(nby, IntVect(3, 2, 3));

          nbz = amrex::surroundingNodes(bx, 2);
          nbzf = amrex::grow(nbz, IntVect(3, 3, 2));

          flxx.resize(nbxf, NVAR+3);
          Extmp.resize(nbxf);

          flxy.resize(nbyf, NVAR+3);
          Eytmp.resize(nbyf);

          flxz.resize(nbzf, NVAR+3);
          Eztmp.resize(nbzf);

          amrex::ParallelFor(nbxf, NVAR+3,
          [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            flxx(i,j,k,n) = 0.0_rt;
          });

          amrex::ParallelFor(nbyf, NVAR+3,
          [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            flxy(i,j,k,n) = 0.0_rt;
          });

          amrex::ParallelFor(nbzf, NVAR+3,
          [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
          {
            flxz(i,j,k,n) = 0.0_rt;
          });

          corner_transport(BL_TO_FORTRAN_ANYD(q),
                           BL_TO_FORTRAN_ANYD(qm),
                           BL_TO_FORTRAN_ANYD(qp),
                           BL_TO_FORTRAN_ANYD(flxx),
                           BL_TO_FORTRAN_ANYD(flxy),
                           BL_TO_FORTRAN_ANYD(flxz),
                           BL_TO_FORTRAN_ANYD(Extmp),
                           BL_TO_FORTRAN_ANYD(Eytmp),
                           BL_TO_FORTRAN_ANYD(Eztmp),
                           dx, dt);

          // Conservative update
          consup(lo, hi,
                 BL_TO_FORTRAN_ANYD(statein),
                 BL_TO_FORTRAN_ANYD(source_out),
                 BL_TO_FORTRAN(bcc),
                 BL_TO_FORTRAN(flxx),
                 BL_TO_FORTRAN(flxy),
                 BL_TO_FORTRAN(flxz),
                 dx, dt);


          // magnetic update
          amrex::ParallelFor(nbx,
          [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            Bxout(i,j,k) = 0.0_rt;
          });

          amrex::ParallelFor(nby,
          [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            Byout(i,j,k) = 0.0_rt;
          });

          amrex::ParallelFor(nbz,
          [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            Bzout(i,j,k) = 0.0_rt;
          });

          magup(lo, hi,
                BL_TO_FORTRAN_ANYD(bxin),
                BL_TO_FORTRAN_ANYD(byin),
                BL_TO_FORTRAN_ANYD(bzin),
                BL_TO_FORTRAN_ANYD(bxout),
                BL_TO_FORTRAN_ANYD(byout),
                BL_TO_FORTRAN_ANYD(bzout),
                BL_TO_FORTRAN_ANYD(Extmp),
                BL_TO_FORTRAN_ANYD(Eytmp),
                BL_TO_FORTRAN_ANYD(Eztmp),
                dx, dt);


          // store the fluxes
          amrex::ParallelFor(nbx,
          [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            flux1(i,j,k,URHO) = flxx(i,j,k,URHO);
            flux1(i,j,k,UMX) = flxx(i,j,k,UMX);
            flux1(i,j,k,UMY) = flxx(i,j,k,UMY);
            flux1(i,j,k,UMZ) = flxx(i,j,k,UMZ);
            flux1(i,j,k,UEDEN) = flxx(i,j,k,UEDEN);
            for (int n = 0; n < NumSpec; n++) {
              flux1(i,j,k,UFS+n) = flxx(i,j,k,UFS+n);
            }
          });

          amrex::ParallelFor(nby,
          [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            flux2(i,j,k,URHO) = flxy(i,j,k,URHO);
            flux2(i,j,k,UMX) = flxy(i,j,k,UMX);
            flux2(i,j,k,UMY) = flxy(i,j,k,UMY);
            flux2(i,j,k,UMZ) = flxy(i,j,k,UMZ);
            flux2(i,j,k,UEDEN) = flxy(i,j,k,UEDEN);
            for (int n = 0; n < NumSpec; n++) {
              flux2(i,j,k,UFS+n) = flxy(i,j,k,UFS+n);
            }
          });

          amrex::ParallelFor(nbz,
          [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
            flux3(i,j,k,URHO) = flxz(i,j,k,URHO);
            flux3(i,j,k,UMX) = flxz(i,j,k,UMX);
            flux3(i,j,k,UMY) = flxz(i,j,k,UMY);
            flux3(i,j,k,UMZ) = flxz(i,j,k,UMZ);
            flux3(i,j,k,UEDEN) = flxz(i,j,k,UEDEN);
            for (int n = 0; n < NumSpec; n++) {
              flux3(i,j,k,UFS+n) = flxz(i,j,k,UFS+n);
            }
          });

          // not sure if this is needed

          //Ex(ex_lo(1):ex_hi(1),ex_lo(2):ex_hi(2), ex_lo(3):ex_hi(3)) = Extemp(ex_lo(1):ex_hi(1),ex_lo(2):ex_hi(2),ex_lo(3):ex_hi(3))
          //Ey(ey_lo(1):ey_hi(1),ey_lo(2):ey_hi(2), ey_lo(3):ey_hi(3)) = Eytemp(ey_lo(1):ey_hi(1),ey_lo(2):ey_hi(2),ey_lo(3):ey_hi(3))
          //Ez(ez_lo(1):ez_hi(1),ez_lo(2):ez_hi(2), ez_lo(3):ez_hi(3)) = Eztemp(ez_lo(1):ez_hi(1),ez_lo(2):ez_hi(2),ez_lo(3):ez_hi(3))

          for (int i = 0; i < BL_SPACEDIM; i++){
            (*fluxes[i])[mfi].plus(flux[i], mfi.nodaltilebox(i),0,0,NUM_STATE);

            (*mass_fluxes[i])[mfi].copy(flux[i],mfi.nodaltilebox(i),Density,mfi.nodaltilebox(i),0,1);
            //electric[i][mfi].copy(E[i], mfi.nodaltilebox(i));
          }

        }

    }

}
