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

      const auto dx = geom.CellSizeArray();
      const Real* dx_f = geom.CellSize();

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
#pragma omp parallel
#endif
    {

      FArrayBox flux[AMREX_SPACEDIM], E[AMREX_SPACEDIM];

      FArrayBox q;
      FArrayBox qaux;
      FArrayBox srcQ;

      FArrayBox flatn;

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
          auto u_arr = statein.array();

          FArrayBox &stateout = S_new[mfi];

          FArrayBox &source_in  = sources_for_hydro[mfi];
          auto src_arr = source_in.array();

          FArrayBox &hydro_update = hydro_source[mfi];
          auto update_arr = hydro_update.array();

          FArrayBox& Bx  = Bx_old_tmp[mfi];
          auto Bx_arr = Bx.array();
          auto elix_Bx = Bx.elixir();

          FArrayBox& By  = By_old_tmp[mfi];
          auto By_arr = By.array();
          auto elix_By = By.elixir();

          FArrayBox& Bz  = Bz_old_tmp[mfi];
          auto Bz_arr = Bz.array();
          auto elix_Bz = Bz.elixir();

          FArrayBox& Bxout = Bx_new[mfi];
          auto Bxo_arr = Bxout.array();

          FArrayBox& Byout = By_new[mfi];
          auto Byo_arr = Byout.array();

          FArrayBox& Bzout = Bz_new[mfi];
          auto Bzo_arr = Bzout.array();


          // allocate fabs for fluxes
          for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
            const Box& bxtmp = amrex::surroundingNodes(bx, idir);
            flux[idir].resize(bxtmp, NUM_STATE);
            E[idir].resize(bxtmp, NUM_STATE);
          }

          // Calculate primitives based on conservatives
          q.resize(bx_gc, NQ);
          auto q_arr = q.array();
          auto elix_q = q.elixir();

          qaux.resize(bx_gc, NQAUX);
          auto qaux_arr = qaux.array();
          auto elix_qaux = qaux.elixir();

          srcQ.resize(bx_gc, NQSRC);
          auto src_q_arr = srcQ.array();
          auto elix_src_q = srcQ.elixir();

          const int* lo_gc = bx_gc.loVect();
          const int* hi_gc = bx_gc.hiVect();


          ctoprim(bx_gc, time,
                  u_arr,
                  Bx_arr, By_arr, Bz_arr,
                  q_arr, qaux_arr);

          const int* lo1 = obx.loVect();
          const int* hi1 = obx.hiVect();


          src_to_prim(bx_gc, q_arr, src_arr, src_q_arr);

          check_for_mhd_cfl_violation(bx, dt, q_arr, qaux_arr);

          flatn.resize(bx_gc, 1);
          auto flatn_arr = flatn.array();
          auto elix_flatn = flatn.elixir();

          amrex::ParallelFor(obx,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
          {
            flatn_arr(i,j,k) = 0.0;
          });

          mhd_flatten(lo1, hi1,
                      BL_TO_FORTRAN_ANYD(q),
                      BL_TO_FORTRAN_ANYD(flatn));


          // Interpolate Cell centered values to faces
          qp.resize(bx_gc, NQ * AMREX_SPACEDIM);
          qm.resize(bx_gc, NQ * AMREX_SPACEDIM);

          const Box& nbx = amrex::surroundingNodes(bx, 0);
          const Box& nby = amrex::surroundingNodes(bx, 1);
          const Box& nbz = amrex::surroundingNodes(bx, 2);

          const Box& nbxi = amrex::grow(bx, IntVect(3, 3, 3));

          plm(nbxi.loVect(), nbxi.hiVect(), 1,
              BL_TO_FORTRAN_ANYD(q),
              BL_TO_FORTRAN_ANYD(qaux),
              BL_TO_FORTRAN_ANYD(flatn),
              BL_TO_FORTRAN_ANYD(Bx),
              BL_TO_FORTRAN_ANYD(By),
              BL_TO_FORTRAN_ANYD(Bz),
              BL_TO_FORTRAN_ANYD(qp),
              BL_TO_FORTRAN_ANYD(qm),
              BL_TO_FORTRAN_ANYD(srcQ),
              dx_f, dt);

          const Box& nbyi = amrex::grow(bx, IntVect(3, 3, 3));

          plm(nbyi.loVect(), nbyi.hiVect(), 2,
              BL_TO_FORTRAN_ANYD(q),
              BL_TO_FORTRAN_ANYD(qaux),
              BL_TO_FORTRAN_ANYD(flatn),
              BL_TO_FORTRAN_ANYD(Bx),
              BL_TO_FORTRAN_ANYD(By),
              BL_TO_FORTRAN_ANYD(Bz),
              BL_TO_FORTRAN_ANYD(qp),
              BL_TO_FORTRAN_ANYD(qm),
              BL_TO_FORTRAN_ANYD(srcQ),
              dx_f, dt);

          const Box& nbzi = amrex::grow(bx, IntVect(3, 3, 3));

          plm(nbzi.loVect(), nbzi.hiVect(), 3,
              BL_TO_FORTRAN_ANYD(q),
              BL_TO_FORTRAN_ANYD(qaux),
              BL_TO_FORTRAN_ANYD(flatn),
              BL_TO_FORTRAN_ANYD(Bx),
              BL_TO_FORTRAN_ANYD(By),
              BL_TO_FORTRAN_ANYD(Bz),
              BL_TO_FORTRAN_ANYD(qp),
              BL_TO_FORTRAN_ANYD(qm),
              BL_TO_FORTRAN_ANYD(srcQ),
              dx_f, dt);


          // Corner Couple and find the correct fluxes + electric fields

          // need to revisit these box sizes
          const Box& nbxf = amrex::grow(nbx, IntVect(2, 3, 3));
          const Box& nbyf = amrex::grow(nby, IntVect(3, 2, 3));
          const Box& nbzf = amrex::grow(nbz, IntVect(3, 3, 2));

          flxx.resize(nbxf, NUM_STATE+3);
          auto flxx_arr = flxx.array();
          auto elix_flxx = flxx.elixir();

          Extmp.resize(nbxf);
          auto Ex_arr = Extmp.array();
          auto elix_Ex = Extmp.elixir();

          flxy.resize(nbyf, NUM_STATE+3);
          auto flxy_arr = flxy.array();
          auto elix_flxy = flxy.elixir();

          Eytmp.resize(nbyf);
          auto Ey_arr = Eytmp.array();
          auto elix_Ey = Eytmp.elixir();

          flxz.resize(nbzf, NUM_STATE+3);
          auto flxz_arr = flxz.array();
          auto elix_flxz = flxz.elixir();

          Eztmp.resize(nbzf);
          auto Ez_arr = Eztmp.array();
          auto elix_Ez = Eztmp.elixir();

          corner_transport(lo, hi,
                           BL_TO_FORTRAN_ANYD(q),
                           BL_TO_FORTRAN_ANYD(qm),
                           BL_TO_FORTRAN_ANYD(qp),
                           BL_TO_FORTRAN_ANYD(flxx),
                           BL_TO_FORTRAN_ANYD(flxy),
                           BL_TO_FORTRAN_ANYD(flxz),
                           BL_TO_FORTRAN_ANYD(Extmp),
                           BL_TO_FORTRAN_ANYD(Eytmp),
                           BL_TO_FORTRAN_ANYD(Eztmp),
                           dx_f, dt);

          // Conservative update

          consup_mhd(bx, update_arr, flxx_arr, flxy_arr, flxz_arr);

          // magnetic update

          amrex::ParallelFor(nbx,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
          {
            Bxo_arr(i,j,k) = Bx_arr(i,j,k) + dt/dx[0] *
              ((Ey_arr(i,j,k+1) - Ey_arr(i,j,k)) - (Ez_arr(i,j+1,k) - Ez_arr(i,j,k)));
          });

          amrex::ParallelFor(nby,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
          {
            Byo_arr(i,j,k) = By_arr(i,j,k) + dt/dx[1] *
              ((Ez_arr(i+1,j,k) - Ez_arr(i,j,k)) - (Ex_arr(i,j,k+1) - Ex_arr(i,j,k)));
          });

          amrex::ParallelFor(nbz,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
          {
            Bzo_arr(i,j,k) = Bz_arr(i,j,k) + dt/dx[2] *
              ((Ex_arr(i,j+1,k) - Ex_arr(i,j,k)) - (Ey_arr(i+1,j,k) - Ey_arr(i,j,k)));
          });


          // store the fluxes -- it looks like we don't need these temporary fluxes?
          auto flux0_arr = flux[0].array();
          amrex::ParallelFor(nbx,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
          {
            flux0_arr(i,j,k,URHO) = flxx_arr(i,j,k,URHO);
            flux0_arr(i,j,k,UMX) = flxx_arr(i,j,k,UMX);
            flux0_arr(i,j,k,UMY) = flxx_arr(i,j,k,UMY);
            flux0_arr(i,j,k,UMZ) = flxx_arr(i,j,k,UMZ);
            flux0_arr(i,j,k,UEDEN) = flxx_arr(i,j,k,UEDEN);
            flux0_arr(i,j,k,UEINT) = flxx_arr(i,j,k,UEINT);
            flux0_arr(i,j,k,UTEMP) = 0.0;
#ifdef SHOCK_VAR
            flux0_arr(i,j,k,USHK) = 0.0;
#endif
            for (int n = 0; n < NumSpec; n++) {
              flux0_arr(i,j,k,UFS+n) = flxx_arr(i,j,k,UFS+n);
            }
          });

          auto flux1_arr = flux[1].array();
          amrex::ParallelFor(nby,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
          {
            flux1_arr(i,j,k,URHO) = flxy_arr(i,j,k,URHO);
            flux1_arr(i,j,k,UMX) = flxy_arr(i,j,k,UMX);
            flux1_arr(i,j,k,UMY) = flxy_arr(i,j,k,UMY);
            flux1_arr(i,j,k,UMZ) = flxy_arr(i,j,k,UMZ);
            flux1_arr(i,j,k,UEDEN) = flxy_arr(i,j,k,UEDEN);
            flux1_arr(i,j,k,UEINT) = flxy_arr(i,j,k,UEINT);
            flux1_arr(i,j,k,UTEMP) = 0.0;
#ifdef SHOCK_VAR
            flux1_arr(i,j,k,USHK) = 0.0;
#endif
            for (int n = 0; n < NumSpec; n++) {
              flux1_arr(i,j,k,UFS+n) = flxy_arr(i,j,k,UFS+n);
            }
          });

          auto flux2_arr = flux[2].array();
          amrex::ParallelFor(nbz,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
          {
            flux2_arr(i,j,k,URHO) = flxz_arr(i,j,k,URHO);
            flux2_arr(i,j,k,UMX) = flxz_arr(i,j,k,UMX);
            flux2_arr(i,j,k,UMY) = flxz_arr(i,j,k,UMY);
            flux2_arr(i,j,k,UMZ) = flxz_arr(i,j,k,UMZ);
            flux2_arr(i,j,k,UEDEN) = flxz_arr(i,j,k,UEDEN);
            flux2_arr(i,j,k,UEINT) = flxz_arr(i,j,k,UEINT);
            flux2_arr(i,j,k,UTEMP) = 0.0;
#ifdef SHOCK_VAR
            flux2_arr(i,j,k,USHK) = 0.0;
#endif
            for (int n = 0; n < NumSpec; n++) {
              flux2_arr(i,j,k,UFS+n) = flxz_arr(i,j,k,UFS+n);
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
