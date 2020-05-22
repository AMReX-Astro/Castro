#include "Castro.H"
#include "Castro_F.H"
#include <iostream>
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


      BL_ASSERT(NUM_GROW == 6);


#ifdef _OPENMP
#pragma omp parallel
#endif
    {

      FArrayBox flux[AMREX_SPACEDIM], E[AMREX_SPACEDIM];

      FArrayBox q;
      FArrayBox qaux;
      FArrayBox srcQ;

      FArrayBox flatn;

      FArrayBox qx_left;
      FArrayBox qx_right;
      FArrayBox qy_left;
      FArrayBox qy_right;
      FArrayBox qz_left;
      FArrayBox qz_right;

      FArrayBox flxx1D;
      FArrayBox flxy1D;
      FArrayBox flxz1D;

      FArrayBox ux_left;
      FArrayBox ux_right;
      FArrayBox uy_left;
      FArrayBox uy_right;
      FArrayBox uz_left;
      FArrayBox uz_right;

      FArrayBox qtmp_left;
      FArrayBox qtmp_right;

      FArrayBox flx_xy;
      FArrayBox flx_xz;

      FArrayBox flx_yx;
      FArrayBox flx_yz;

      FArrayBox flx_zx;
      FArrayBox flx_zy;

      FArrayBox q2D;

      for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {

          const Box& bx = mfi.tilebox();
          const Box& obx = amrex::grow(bx, 1);
          const Box& gbx = amrex::grow(bx, 2);

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


          // allocate fabs for fluxes and electric field

          const Box& nbx = amrex::surroundingNodes(bx, 0);
          const Box& nby = amrex::surroundingNodes(bx, 1);
          const Box& nbz = amrex::surroundingNodes(bx, 2);

          const Box& nbxf = amrex::grow(nbx, IntVect(0, 1, 1));
          const Box& nbyf = amrex::grow(nby, IntVect(1, 0, 1));
          const Box& nbzf = amrex::grow(nbz, IntVect(1, 1, 0));

          // need to revisit these box sizes
          const Box& nbxe = amrex::grow(nbx, IntVect(2, 3, 3));
          const Box& nbye = amrex::grow(nby, IntVect(3, 2, 3));
          const Box& nbze = amrex::grow(nbz, IntVect(3, 3, 2));

          flux[0].resize(nbxf, NUM_STATE+3);
          auto flxx_arr = flux[0].array();
          auto elix_flxx = flux[0].elixir();

          E[0].resize(nbxe);
          auto Ex_arr = E[0].array();
          auto elix_Ex = E[0].elixir();

          flux[1].resize(nbyf, NUM_STATE+3);
          auto flxy_arr = flux[1].array();
          auto elix_flxy = flux[1].elixir();

          E[1].resize(nbye);
          auto Ey_arr = E[1].array();
          auto elix_Ey = E[1].elixir();

          flux[2].resize(nbzf, NUM_STATE+3);
          auto flxz_arr = flux[2].array();
          auto elix_flxz = flux[2].elixir();

          E[2].resize(nbze);
          auto Ez_arr = E[2].array();
          auto elix_Ez = E[2].elixir();


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
          qx_left.resize(bx_gc, NQ);
          auto qx_left_arr = qx_left.array();
          auto elix_qx_left = qx_left.elixir();

          qx_right.resize(bx_gc, NQ);
          auto qx_right_arr = qx_right.array();
          auto elix_qx_right = qx_right.elixir();

          qy_left.resize(bx_gc, NQ);
          auto qy_left_arr = qy_left.array();
          auto elix_qy_left = qy_left.elixir();

          qy_right.resize(bx_gc, NQ);
          auto qy_right_arr = qy_right.array();
          auto elix_qy_right = qy_right.elixir();

          qz_left.resize(bx_gc, NQ);
          auto qz_left_arr = qz_left.array();
          auto elix_qz_left = qz_left.elixir();

          qz_right.resize(bx_gc, NQ);
          auto qz_right_arr = qz_right.array();
          auto elix_qz_right = qz_right.elixir();

          const Box& nbxi = amrex::grow(bx, IntVect(3, 3, 3));

          plm(nbxi.loVect(), nbxi.hiVect(), 1,
              BL_TO_FORTRAN_ANYD(q),
              BL_TO_FORTRAN_ANYD(qaux),
              BL_TO_FORTRAN_ANYD(flatn),
              BL_TO_FORTRAN_ANYD(Bx),
              BL_TO_FORTRAN_ANYD(By),
              BL_TO_FORTRAN_ANYD(Bz),
              BL_TO_FORTRAN_ANYD(qx_left),
              BL_TO_FORTRAN_ANYD(qx_right),
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
              BL_TO_FORTRAN_ANYD(qy_left),
              BL_TO_FORTRAN_ANYD(qy_right),
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
              BL_TO_FORTRAN_ANYD(qz_left),
              BL_TO_FORTRAN_ANYD(qz_right),
              BL_TO_FORTRAN_ANYD(srcQ),
              dx_f, dt);


          // Corner Couple and find the correct fluxes + electric fields

          // Do the corner coupling and the CT updates

          // MM CTU Step 1
          // Calculate Flux 1D, eq.35

          // x-dir
          // [lo(1)-2, lo(2)-3, lo(3)-3] [hi(1)+3, hi(2)+3, hi(3)+3]
          const Box& bfx = amrex::grow(nbx, IntVect(2, 3, 3));

          flxx1D.resize(bfx, NUM_STATE+3);
          auto flxx1D_arr = flxx1D.array();
          auto elix_flxx1D = flxx1D.elixir();

          hlld(bfx.loVect(), bfx.hiVect(),
               BL_TO_FORTRAN_ANYD(qx_left),
               BL_TO_FORTRAN_ANYD(qx_right),
               BL_TO_FORTRAN_ANYD(flxx1D), 1);

          // y-dir
          // [lo(1)-3, lo(2)-2, lo(3)-3] [hi(1)+3, hi(2)+3, hi(3)+3]
          const Box& bfy = amrex::grow(nby, IntVect(3, 2, 3));

          flxy1D.resize(bfy, NUM_STATE+3);
          auto flxy1D_arr = flxy1D.array();
          auto elix_flxy1D = flxy1D.elixir();

          hlld(bfy.loVect(), bfy.hiVect(),
               BL_TO_FORTRAN_ANYD(qy_left),
               BL_TO_FORTRAN_ANYD(qy_right),
               BL_TO_FORTRAN_ANYD(flxy1D), 2);

          // z-dir
          // [lo(1)-3, lo(2)-3, lo(3)-2] [hi(1)+3, hi(2)+3, hi(3)+3]
          const Box& bfz = amrex::grow(nbz, IntVect(3, 3, 2));

          flxz1D.resize(bfz, NUM_STATE+3);
          auto flxz1D_arr = flxz1D.array();
          auto elix_flxz1D = flxz1D.elixir();

          hlld(bfz.loVect(), bfz.hiVect(),
               BL_TO_FORTRAN_ANYD(qz_left),
               BL_TO_FORTRAN_ANYD(qz_right),
               BL_TO_FORTRAN_ANYD(flxz1D), 3);


          // Prim to Cons

          ux_left.resize(gbx, NUM_STATE+3);
          auto ux_left_arr = ux_left.array();
          auto elix_ux_left = ux_left.elixir();

          ux_right.resize(gbx, NUM_STATE+3);
          auto ux_right_arr = ux_right.array();
          auto elix_ux_right = ux_right.elixir();

          PrimToCons(gbx, qx_left_arr, ux_left_arr);
          PrimToCons(gbx, qx_right_arr, ux_right_arr);

          uy_left.resize(gbx, NUM_STATE+3);
          auto uy_left_arr = uy_left.array();
          auto elix_uy_left = uy_left.elixir();

          uy_right.resize(gbx, NUM_STATE+3);
          auto uy_right_arr = uy_right.array();
          auto elix_uy_right = uy_right.elixir();

          PrimToCons(gbx, qy_left_arr, uy_left_arr);
          PrimToCons(gbx, qy_right_arr, uy_right_arr);

          uz_left.resize(gbx, NUM_STATE+3);
          auto uz_left_arr = uz_left.array();
          auto elix_uz_left = uz_left.elixir();

          uz_right.resize(gbx, NUM_STATE+3);
          auto uz_right_arr = uz_right.array();
          auto elix_uz_right = uz_right.elixir();

          PrimToCons(gbx, qz_left_arr, uz_left_arr);
          PrimToCons(gbx, qz_right_arr, uz_right_arr);

          // MM CTU Step 2
          // Use "1D" fluxes To interpolate Temporary Edge Centered Electric Fields, eq.36

          // [lo(1)-2, lo(2)-2, lo(3)-2][hi(1)+2, hi(2)+3, hi(3)+3]
          Box eebx = amrex::grow(bx, 2);
          eebx.growHi(1);
          eebx.growHi(2);

          electric_edge_x(eebx, q_arr, Ex_arr, flxy1D_arr, flxz1D_arr);


          // [lo(1)-2, lo(2)-2, lo(3)-2][hi(1)+3, hi(2)+2, hi(3)+3]
          Box eeby = amrex::grow(bx, 2);
          eeby.growHi(0);
          eeby.growHi(2);

          electric_edge_y(eeby, q_arr, Ey_arr, flxx1D_arr, flxz1D_arr);


          // [lo(1)-2, lo(2)-2, lo(3)-2][hi(1)+3, hi(2)+3, hi(3)+2]
          Box eebz = amrex::grow(bx, 2);
          eebz.growHi(0);
          eebz.growHi(1);

          electric_edge_z(eebz, q_arr, Ez_arr, flxx1D_arr, flxy1D_arr);


          // MM CTU Steps 3, 4, and 5
          // Corner Couple, eq. 37, 38 and 39 Correct Conservative vars using Transverse Fluxes

          // X direction

          // affected by Y Flux
          // [lo(1)-1, lo(2)-2, lo(3)-2] [hi(1)+2, hi(2)+2, hi(2)+2]
          const Box& ccbx = amrex::grow(nbx, IntVect(1, 2, 2));

          qtmp_left.resize(gbx, NQ);
          auto qtmp_left_arr = qtmp_left.array();
          auto elix_qtmp_left = qtmp_left.elixir();

          qtmp_right.resize(gbx, NQ);
          auto qtmp_right_arr = qtmp_right.array();
          auto elix_qtmp_right = qtmp_right.elixir();

          corner_couple(ccbx,
                        qtmp_right_arr, qtmp_left_arr,
                        ux_right_arr, ux_left_arr,
                        flxy1D_arr, Ex_arr, Ez_arr,
                        0, 1, 2, dt);

          // Calculate Flux 2D eq. 40
          // F^{x|y}
          flx_xy.resize(ccbx, NUM_STATE+3);
          auto flx_xy_arr = flx_xy.array();
          auto elix_flx_xy = flx_xy.elixir();

          hlld(ccbx.loVect(), ccbx.hiVect(),
               BL_TO_FORTRAN_ANYD(qtmp_left),
               BL_TO_FORTRAN_ANYD(qtmp_right),
               BL_TO_FORTRAN_ANYD(flx_xy), 1);

          // affected by Z Flux
          corner_couple(ccbx,
                        qtmp_right_arr, qtmp_left_arr,
                        ux_right_arr, ux_left_arr,
                        flxz1D_arr, Ex_arr, Ey_arr,
                        0, 2, 1, dt);

          // F^{x|z}
          flx_xz.resize(ccbx, NUM_STATE+3);
          auto flx_xz_arr = flx_xz.array();
          auto elix_flx_xz = flx_xz.elixir();

          hlld(ccbx.loVect(), ccbx.hiVect(),
               BL_TO_FORTRAN_ANYD(qtmp_left),
               BL_TO_FORTRAN_ANYD(qtmp_right),
               BL_TO_FORTRAN_ANYD(flx_xz), 1);


          // Y direction

          // affected by X Flux
          // [lo(1)-2, lo(2)-1, lo(3)-2] [hi(1)+2, hi(2)+2, hi(3)+2]
          const Box& ccby = amrex::grow(nby, IntVect(2, 1, 2));

          corner_couple(ccby,
                        qtmp_right_arr, qtmp_left_arr,
                        uy_right_arr, uy_left_arr,
                        flxx1D_arr, Ey_arr, Ez_arr,
                        1, 0, 2, dt);

          // F^{y|x}
          flx_yx.resize(ccby, NUM_STATE+3);
          auto flx_yx_arr = flx_yx.array();
          auto elix_flx_yx = flx_yx.elixir();

          hlld(ccby.loVect(), ccby.hiVect(),
               BL_TO_FORTRAN_ANYD(qtmp_left),
               BL_TO_FORTRAN_ANYD(qtmp_right),
               BL_TO_FORTRAN_ANYD(flx_yx), 2);

          // affected by Z Flux

          corner_couple(ccby,
                        qtmp_right_arr, qtmp_left_arr,
                        uy_right_arr, uy_left_arr,
                        flxz1D_arr, Ey_arr, Ex_arr,
                        1, 2, 0, dt);

          // F^{y|z}
          flx_yz.resize(ccby, NUM_STATE+3);
          auto flx_yz_arr = flx_yz.array();
          auto elix_flx_yz = flx_yz.elixir();

          hlld(ccby.loVect(), ccby.hiVect(),
               BL_TO_FORTRAN_ANYD(qtmp_left),
               BL_TO_FORTRAN_ANYD(qtmp_right),
               BL_TO_FORTRAN_ANYD(flx_yz), 2);

          // Z direction

          // affected by X Flux
          // [lo(1)-2, lo(2)-2, lo(3)-1] [hi(1)+2, hi(2)+2, hi(3)+2]
          const Box& ccbz = amrex::grow(nbz, IntVect(2, 2, 1));

          corner_couple(ccbz,
                        qtmp_right_arr, qtmp_left_arr,
                        uz_right_arr, uz_left_arr,
                        flxx1D_arr, Ez_arr, Ey_arr,
                        2, 0, 1, dt);

          // F^{z|x}
          flx_zx.resize(ccbz, NUM_STATE+3);
          auto flx_zx_arr = flx_zx.array();
          auto elix_flx_zx = flx_zx.elixir();

          hlld(ccbz.loVect(), ccbz.hiVect(),
               BL_TO_FORTRAN_ANYD(qtmp_left),
               BL_TO_FORTRAN_ANYD(qtmp_right),
               BL_TO_FORTRAN_ANYD(flx_zx), 3);

          // affected by Y Flux

          corner_couple(ccbz,
                        qtmp_right_arr, qtmp_left_arr,
                        uz_right_arr, uz_left_arr,
                        flxy1D_arr, Ez_arr, Ex_arr,
                        2, 1, 0, dt);

          // F^{z|y}
          flx_zy.resize(ccbz, NUM_STATE+3);
          auto flx_zy_arr = flx_zy.array();
          auto elix_flx_zy = flx_zy.elixir();

          hlld(ccbz.loVect(), ccbz.hiVect(),
               BL_TO_FORTRAN_ANYD(qtmp_left),
               BL_TO_FORTRAN_ANYD(qtmp_right),
               BL_TO_FORTRAN_ANYD(flx_zy), 3);


          // MM CTU Step 6
          // Use Averaged 2D fluxes to interpolate temporary Edge Centered Electric Fields, reuse "flx1D"
          // eq. 42 and 43

          amrex::ParallelFor(ccbx, NUM_STATE+3,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int n) noexcept
          {
            flxx1D_arr(i,j,k,n) = 0.5_rt * (flx_xy_arr(i,j,k,n) + flx_xz_arr(i,j,k,n));
          });

          amrex::ParallelFor(ccby, NUM_STATE+3,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int n) noexcept
          {
            flxy1D_arr(i,j,k,n) = 0.5_rt * (flx_yx_arr(i,j,k,n) + flx_yz_arr(i,j,k,n));
          });

          amrex::ParallelFor(ccbz, NUM_STATE+3,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int n) noexcept
          {
            flxz1D_arr(i,j,k,n) = 0.5_rt * (flx_zx_arr(i,j,k,n) + flx_zy_arr(i,j,k,n));
          });


          // eq. 41
          // [lo(1)-1, lo(2)-1, lo(3)-1][hi(1)+1, hi(2)+2, hi(3)+2]
          Box eebx2 = amrex::grow(bx, 1);
          eebx2.growHi(1);
          eebx2.growHi(2);

          electric_edge_x(eebx2, q_arr, Ex_arr, flxy1D_arr, flxz1D_arr);


          // [lo(1)-1, lo(2)-1, lo(3)-1][hi(1)+2, hi(2)+1, hi(3)+2]
          Box eeby2 = amrex::grow(bx, 1);
          eeby2.growHi(0);
          eeby2.growHi(2);

          electric_edge_y(eeby2, q_arr, Ey_arr, flxx1D_arr, flxz1D_arr);


          // [lo(1)-1, lo(2)-1, lo(3)-1][hi(1)+2, hi(2)+2, hi(3)+1]
          Box eebz2 = amrex::grow(bx, 1);
          eebz2.growHi(0);
          eebz2.growHi(1);

          electric_edge_z(eebz2, q_arr, Ez_arr, flxx1D_arr, flxy1D_arr);


          // MM CTU Step 7, 8, and 9
          // Half Step conservative vars eq.44, eq.45, eq.46
          // Here we reuse qtmp_left/right to denote the half-time conservative state

          // for x direction
          // [lo(1), lo(2)-1, lo(3)-1][hi(1)+1, hi(2)+1, hi(3)+1]
          const Box& nbx1 = amrex::grow(nbx, IntVect(0, 1, 1));

          half_step(nbx1,
                    qtmp_right_arr, qtmp_left_arr,
                    ux_right_arr, ux_left_arr,
                    flx_yz_arr, flx_zy_arr,
                    Ex_arr, Ey_arr, Ez_arr,
                    0, 1, 2, dt);

          // Final Fluxes eq.47

          // We need to compute these on a box 1 larger in the transverse directions
          // than we'd need for hydro alone due to the electric update

          hlld(nbx1.loVect(), nbx1.hiVect(),
               BL_TO_FORTRAN_ANYD(qtmp_left),
               BL_TO_FORTRAN_ANYD(qtmp_right),
               BL_TO_FORTRAN_ANYD(flux[0]), 1);

          // for y direction
          const Box& nby1 = amrex::grow(nby, IntVect(1, 0, 1));

          half_step(nby1,
                    qtmp_right_arr, qtmp_left_arr,
                    uy_right_arr, uy_left_arr,
                    flx_xz_arr, flx_zx_arr,
                    Ey_arr, Ex_arr, Ez_arr,
                    1, 0, 2, dt);

          hlld(nby1.loVect(), nby1.hiVect(),
               BL_TO_FORTRAN_ANYD(qtmp_left),
               BL_TO_FORTRAN_ANYD(qtmp_right),
               BL_TO_FORTRAN_ANYD(flux[1]), 2);

          // for z direction
          const Box& nbz1 = amrex::grow(nbz, IntVect(1, 1, 0));

          half_step(nbz1,
                    qtmp_right_arr, qtmp_left_arr,
                    uz_right_arr, uz_left_arr,
                    flx_xy_arr, flx_yx_arr,
                    Ez_arr, Ex_arr, Ey_arr,
                    2, 0, 1, dt);

          hlld(nbz1.loVect(), nbz1.hiVect(),
               BL_TO_FORTRAN_ANYD(qtmp_left),
               BL_TO_FORTRAN_ANYD(qtmp_right),
               BL_TO_FORTRAN_ANYD(flux[2]), 3);


          // MM CTU Step 10
          // Primitive update eq. 48
          q2D.resize(obx, NQ);
          auto q2D_arr = q2D.array();
          auto elix_q2D = q2D.elixir();

          prim_half(obx, q2D_arr, q_arr,
                    flxx1D_arr, flxy1D_arr, flxz1D_arr, dt);

          // Final Electric Field Update eq.48

          // [lo(1), lo(2), lo(3)][hi(1), hi(2)+1, hi(3)+1]
          Box eebxf = mfi.tilebox();
          eebxf.growHi(1, 1);
          eebxf.growHi(2, 1);

          electric_edge_x(eebxf, q2D_arr, Ex_arr, flxy_arr, flxz_arr);

          // [lo(1), lo(2), lo(3)][hi(1)+1, hi(2), hi(3)+1]
          Box eebyf = mfi.tilebox();
          eebyf.growHi(0, 1);
          eebyf.growHi(2, 1);

          electric_edge_y(eebyf, q2D_arr, Ey_arr, flxx_arr, flxz_arr);

          // [lo(1), lo(2), lo(3)][hi(1)+1, hi(2)+1 ,hi(3)]
          Box eebzf = mfi.tilebox();
          eebzf.growHi(0, 1);
          eebzf.growHi(1, 1);

          electric_edge_z(eebzf, q2D_arr, Ez_arr, flxx_arr, flxy_arr);

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

