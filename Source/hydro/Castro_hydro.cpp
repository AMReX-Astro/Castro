#include "Castro.H"
#include "Castro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

using namespace amrex;

void
Castro::construct_hydro_source(Real time, Real dt)
{

  BL_PROFILE("Castro::construct_hydro_source()");

  const Real strt_time = ParallelDescriptor::second();

  // this constructs the hydrodynamic source (essentially the flux
  // divergence) using the CTU framework for unsplit hydrodynamics

  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "... Entering hydro advance" << std::endl << std::endl;

  hydro_source.setVal(0.0);

  int finest_level = parent->finestLevel();

  const Real *dx = geom.CellSize();

  const int* domain_lo = geom.Domain().loVect();
  const int* domain_hi = geom.Domain().hiVect();

  MultiFab& S_new = get_new_data(State_Type);

#ifdef RADIATION
  MultiFab& Er_new = get_new_data(Rad_Type);

  if (!Radiation::rad_hydro_combined) {
    amrex::Abort("Castro::construct_hydro_source -- we don't implement a mode where we have radiation, but it is not coupled to hydro");
  }

  int nstep_fsp = -1;
#endif

  // note: the radiation consup currently does not fill these
  Real mass_lost       = 0.;
  Real xmom_lost       = 0.;
  Real ymom_lost       = 0.;
  Real zmom_lost       = 0.;
  Real eden_lost       = 0.;
  Real xang_lost       = 0.;
  Real yang_lost       = 0.;
  Real zang_lost       = 0.;

  BL_PROFILE_VAR("Castro::advance_hydro_ca_umdrv()", CA_UMDRV);

#ifdef _OPENMP
#ifdef RADIATION
#pragma omp parallel reduction(max:nstep_fsp) \
                     reduction(+:mass_lost,xmom_lost,ymom_lost,zmom_lost) \
                     reduction(+:eden_lost,xang_lost,yang_lost,zang_lost)
#else
#pragma omp parallel reduction(+:mass_lost,xmom_lost,ymom_lost,zmom_lost) \
                     reduction(+:eden_lost,xang_lost,yang_lost,zang_lost)
#endif
#endif
  {

#ifdef RADIATION
    int priv_nstep_fsp = -1;
#endif

    // Declare local storage now. This should be done outside the MFIter loop,
    // and then we will resize the Fabs in each MFIter loop iteration. Then,
    // we apply an Elixir to ensure that their memory is saved until it is no
    // longer needed (only relevant for the asynchronous case, usually on GPUs).
    
    FArrayBox flatn;
#ifdef RADIATION
    FArrayBox flatg;
#endif
    FArrayBox dq;
    FArrayBox Ip, Im, Ip_src, Im_src, Ip_gc, Im_gc;
    FArrayBox sm, sp;
    FArrayBox shk;
    FArrayBox qxm, qxp;
#if AMREX_SPACEDIM >= 2
    FArrayBox qym, qyp;
#endif
#if AMREX_SPACEDIM == 3
    FArrayBox qzm, qzp;
#endif
    FArrayBox div;
    FArrayBox q_int;
#ifdef RADIATION
    FArrayBox lambda_int;
#endif
#if AMREX_SPACEDIM >= 2
    FArrayBox ftmp1, ftmp2;
#ifdef RADIATION
    FArrayBox rftmp1, rftmp2;
#endif
    FArrayBox qgdnvtmp1, qgdnvtmp2;
    FArrayBox ql, qr;
#endif
    FArrayBox flux[AMREX_SPACEDIM], qe[AMREX_SPACEDIM];
#ifdef RADIATION
    FArrayBox rad_flux[AMREX_SPACEDIM];
#endif
#if AMREX_SPACEDIM <= 2
    FArrayBox pradial;
#endif
#if AMREX_SPACEDIM == 3
    FArrayBox qmyx, qpyx;
    FArrayBox qmzx, qpzx;
    FArrayBox qmxy, qpxy;
    FArrayBox qmzy, qpzy;
    FArrayBox qmxz, qpxz;
    FArrayBox qmyz, qpyz;
#endif
    FArrayBox pdivu;
    
    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      // the valid region box
      const Box& bx = mfi.tilebox();

      const Box& obx = amrex::grow(bx, 1);

      flatn.resize(obx, 1);
      Elixir elix_flatn = flatn.elixir();

#ifdef RADIATION
      flatg.resize(obx, 1);
      Elixir elix_flatg = flatg.elixir();
#endif

      // compute the flattening coefficient

      if (first_order_hydro == 1) {
        flatn.setVal(0.0, obx);
      } else if (use_flattening == 1) {
#ifdef RADIATION
        ca_rad_flatten(ARLIM_3D(obx.loVect()), ARLIM_3D(obx.hiVect()),
                       BL_TO_FORTRAN_ANYD(q[mfi]),
                       BL_TO_FORTRAN_ANYD(flatn),
                       BL_TO_FORTRAN_ANYD(flatg));
#else
#pragma gpu
        ca_uflatten(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
                    BL_TO_FORTRAN_ANYD(q[mfi]),
                    BL_TO_FORTRAN_ANYD(flatn), QPRES+1);
#endif
      } else {
        flatn.setVal(1.0, obx);
      }

      const Box& xbx = amrex::surroundingNodes(bx, 0);
      const Box& gxbx = amrex::grow(xbx, 1);
#if AMREX_SPACEDIM >= 2
      const Box& ybx = amrex::surroundingNodes(bx, 1);
      const Box& gybx = amrex::grow(ybx, 1);
#endif
#if AMREX_SPACEDIM == 3
      const Box& zbx = amrex::surroundingNodes(bx, 2);
      const Box& gzbx = amrex::grow(zbx, 1);
#endif

      dq.resize(obx, AMREX_SPACEDIM*NQ);
      Elixir elix_dq = dq.elixir();

      Ip.resize(obx, AMREX_SPACEDIM*3*NQ);
      Elixir elix_Ip = Ip.elixir();

      Im.resize(obx, AMREX_SPACEDIM*3*NQ);
      Elixir elix_Im = Im.elixir();

      Ip_src.resize(obx, AMREX_SPACEDIM*3*QVAR);
      Elixir elix_Ip_src = Ip_src.elixir();

      Im_src.resize(obx, AMREX_SPACEDIM*3*QVAR);
      Elixir elix_Im_src = Im_src.elixir();

      Ip_gc.resize(obx, AMREX_SPACEDIM*3);
      Elixir elix_Ip_gc = Ip_gc.elixir();

      Im_gc.resize(obx,AMREX_SPACEDIM*3);
      Elixir elix_Im_gc = Im_gc.elixir();

      sm.resize(obx, AMREX_SPACEDIM);
      Elixir elix_sm = sm.elixir();

      sp.resize(obx, AMREX_SPACEDIM);
      Elixir elix_sp = sp.elixir();

      shk.resize(obx, 1);
      Elixir elix_shk = shk.elixir();

      qxm.resize(obx, NQ);
      Elixir elix_qxm = qxm.elixir();

      qxp.resize(obx, NQ);
      Elixir elix_qxp = qxp.elixir();

#if AMREX_SPACEDIM >= 2
      qym.resize(obx, NQ);
      Elixir elix_qym = qym.elixir();

      qyp.resize(obx, NQ);
      Elixir elix_qyp = qyp.elixir();
#endif

#if AMREX_SPACEDIM == 3
      qzm.resize(obx, NQ);
      Elixir elix_qzm = qzm.elixir();

      qzp.resize(obx, NQ);
      Elixir elix_qzp = qzp.elixir();
#endif

#pragma gpu
      ctu_normal_states(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
                        AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                        BL_TO_FORTRAN_ANYD(q[mfi]),
                        BL_TO_FORTRAN_ANYD(flatn),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(src_q[mfi]),
                        BL_TO_FORTRAN_ANYD(shk),
                        BL_TO_FORTRAN_ANYD(Ip),
                        BL_TO_FORTRAN_ANYD(Im),
                        BL_TO_FORTRAN_ANYD(Ip_src),
                        BL_TO_FORTRAN_ANYD(Im_src),
                        BL_TO_FORTRAN_ANYD(Ip_gc),
                        BL_TO_FORTRAN_ANYD(Im_gc),
                        BL_TO_FORTRAN_ANYD(dq),
                        BL_TO_FORTRAN_ANYD(sm),
                        BL_TO_FORTRAN_ANYD(sp),
                        BL_TO_FORTRAN_ANYD(qxm),
                        BL_TO_FORTRAN_ANYD(qxp),
#if AMREX_SPACEDIM >= 2
                        BL_TO_FORTRAN_ANYD(qym),
                        BL_TO_FORTRAN_ANYD(qyp),
#endif
#if AMREX_SPACEDIM == 3
                        BL_TO_FORTRAN_ANYD(qzm),
                        BL_TO_FORTRAN_ANYD(qzp),
#endif
                        AMREX_REAL_ANYD(dx), dt,
#if (AMREX_SPACEDIM < 3)
                        BL_TO_FORTRAN_ANYD(dLogArea[0][mfi]),
#endif
                        AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      div.resize(obx, 1);
      Elixir elix_div = div.elixir();

      // compute divu -- we'll use this later when doing the artifical viscosity
#pragma gpu
      divu(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
           BL_TO_FORTRAN_ANYD(q[mfi]),
           AMREX_REAL_ANYD(dx),
           BL_TO_FORTRAN_ANYD(div));

#ifdef AMREX_USE_CUDA
      Gpu::Device::synchronize();
#endif

      q_int.resize(obx, NQ);
      Elixir elix_q_int = q_int.elixir();

#ifdef RADIATION
      lambda_int.resize(obx, Radiation::nGroups);
      Elixir elix_lambda_int = lambda_int.elixir();
#endif

      flux[0].resize(gxbx, NUM_STATE);
      Elixir elix_flux_x = flux[0].elixir();

      qe[0].resize(gxbx, NGDNV);
      Elixir elix_qe_x = qe[0].elixir();

#ifdef RADIATION
      rad_flux[0].resize(gxbx, Radiation::nGroups);
      Elixir elix_rad_flux_x = rad_flux[0].elixir();
#endif

#if AMREX_SPACEDIM >= 2
      flux[1].resize(gybx, NUM_STATE);
      Elixir elix_flux_y = flux[1].elixir();

      qe[1].resize(gybx, NGDNV);
      Elixir elix_qe_y = qe[1].elixir();

#ifdef RADIATION
      rad_flux[1].resize(gybx, Radiation::nGroups);
      Elixir elix_rad_flux_y = rad_flux[1].elixir();
#endif
#endif

#if AMREX_SPACEDIM == 3
      flux[2].resize(gzbx, NUM_STATE);
      Elixir elix_flux_z = flux[2].elixir();

      qe[2].resize(gzbx, NGDNV);
      Elixir elix_qe_z = qe[2].elixir();

#ifdef RADIATION
      rad_flux[2].resize(gzbx, Radiation::nGroups);
      Elixir elix_rad_flux_z = rad_flux[2].elixir();
#endif
#endif

#if AMREX_SPACEDIM <= 2
      if (!Geometry::IsCartesian()) {
          pradial.resize(xbx, 1);
          Elixir elix_pradial = pradial.elixir();
      }
#endif

#if AMREX_SPACEDIM == 1
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(xbx.loVect()), AMREX_INT_ANYD(xbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qxm),
                          BL_TO_FORTRAN_ANYD(qxp), 1, 1,
                          BL_TO_FORTRAN_ANYD(flux[0]),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rad_flux[0]),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qe[0]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          1, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

#endif // 1-d



#if AMREX_SPACEDIM >= 2
      ftmp1.resize(obx, NUM_STATE);
      Elixir elix_ftmp1 = ftmp1.elixir();

      ftmp2.resize(obx, NUM_STATE);
      Elixir elix_ftmp2 = ftmp2.elixir();

#ifdef RADIATION
      rftmp1.resize(obx, Radiation::nGroups);
      Elixir elix_rftmp1 = rftmp1.elixir();

      rftmp2.resize(obx, Radiation::nGroups);
      Elixir elix_rftmp2 = rftmp2.elixir();
#endif

      qgdnvtmp1.resize(obx, NGDNV);
      Elixir elix_qgdnvtmp1 = qgdnvtmp1.elixir();

      qgdnvtmp2.resize(obx, NGDNV);
      Elixir elix_qgdnvtmp2 = qgdnvtmp2.elixir();

      ql.resize(obx, NQ);
      Elixir elix_ql = ql.elixir();

      qr.resize(obx, NQ);
      Elixir elix_qr = qr.elixir();
#endif



#if AMREX_SPACEDIM == 2

      const amrex::Real hdt = 0.5*dt;
      const amrex::Real hdtdx = 0.5*dt/dx[0];
      const amrex::Real hdtdy = 0.5*dt/dx[1];

      // compute F^x
      // [lo(1), lo(2)-1, 0], [hi(1)+1, hi(2)+1, 0]
      const Box& cxbx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,1,0)));

      // ftmp1 = fx
      // rftmp1 = rfx
      // qgdnvtmp1 = qgdnxv
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(cxbx.loVect()), AMREX_INT_ANYD(cxbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qxm),
                          BL_TO_FORTRAN_ANYD(qxp), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp1),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          1, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // compute F^y
      // [lo(1)-1, lo(2), 0], [hi(1)+1, hi(2)+1, 0]
      const Box& cybx = amrex::grow(ybx, IntVect(AMREX_D_DECL(1,0,0)));

      // ftmp2 = fy
      // rftmp2 = rfy
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(cybx.loVect()), AMREX_INT_ANYD(cybx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qym),
                          BL_TO_FORTRAN_ANYD(qyp), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp2),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp2),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qe[1]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          2, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // add the transverse flux difference in y to the x states
      // [lo(1), lo(2), 0], [hi(1)+1, hi(2), 0]

      // ftmp2 = fy
      // rftmp2 = rfy
#pragma gpu
      transy_on_xstates(AMREX_INT_ANYD(xbx.loVect()), AMREX_INT_ANYD(xbx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qxm),
                        BL_TO_FORTRAN_ANYD(ql),
                        BL_TO_FORTRAN_ANYD(qxp),
                        BL_TO_FORTRAN_ANYD(qr),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp2),
#ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp2),
#endif
                        BL_TO_FORTRAN_ANYD(qe[1]),
                        hdtdy);

      // solve the final Riemann problem axross the x-interfaces

#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(xbx.loVect()), AMREX_INT_ANYD(xbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(ql),
                          BL_TO_FORTRAN_ANYD(qr), 1, 1,
                          BL_TO_FORTRAN_ANYD(flux[0]),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rad_flux[0]),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qe[0]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          1, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // add the transverse flux difference in x to the y states
      // [lo(1), lo(2), 0], [hi(1), hi(2)+1, 0]

      // ftmp1 = fx
      // rftmp1 = rfx
      // qgdnvtmp1 = qgdnvx

#pragma gpu
      transx_on_ystates(AMREX_INT_ANYD(ybx.loVect()), AMREX_INT_ANYD(ybx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qym),
                        BL_TO_FORTRAN_ANYD(ql),
                        BL_TO_FORTRAN_ANYD(qyp),
                        BL_TO_FORTRAN_ANYD(qr),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1),
#ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp1),
#endif
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                        BL_TO_FORTRAN_ANYD(area[0][mfi]),
                        BL_TO_FORTRAN_ANYD(volume[mfi]),
                        hdt, hdtdx);

      // solve the final Riemann problem axross the y-interfaces

#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(ybx.loVect()), AMREX_INT_ANYD(ybx.hiVect()),
                          BL_TO_FORTRAN_ANYD(ql),
                          BL_TO_FORTRAN_ANYD(qr), 1, 1,
                          BL_TO_FORTRAN_ANYD(flux[1]),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rad_flux[1]),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qe[1]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          2, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));
#endif // 2-d



#if AMREX_SPACEDIM == 3

      const amrex::Real hdt = 0.5*dt;

      const amrex::Real hdtdx = 0.5*dt/dx[0];
      const amrex::Real hdtdy = 0.5*dt/dx[1];
      const amrex::Real hdtdz = 0.5*dt/dx[2];

      const amrex::Real cdtdx = dt/dx[0]/3.0;
      const amrex::Real cdtdy = dt/dx[1]/3.0;
      const amrex::Real cdtdz = dt/dx[2]/3.0;

      // compute F^x
      // [lo(1), lo(2)-1, lo(3)-1], [hi(1)+1, hi(2)+1, hi(3)+1]
      const Box& cxbx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,1,1)));

      // ftmp1 = fx
      // rftmp1 = rfx
      // qgdnvtmp1 = qgdnxv
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(cxbx.loVect()), AMREX_INT_ANYD(cxbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qxm),
                          BL_TO_FORTRAN_ANYD(qxp), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp1),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          1, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // [lo(1), lo(2), lo(3)-1], [hi(1), hi(2)+1, hi(3)+1]
      const Box& tyxbx = amrex::grow(ybx, IntVect(AMREX_D_DECL(0,0,1)));

      qmyx.resize(tyxbx, NQ);
      Elixir elix_qmyx = qmyx.elixir();

      qpyx.resize(tyxbx, NQ);
      Elixir elix_qpyx = qpyx.elixir();

      // ftmp1 = fx
      // rftmp1 = rfx
      // qgdnvtmp1 = qgdnvx
#pragma gpu
      transx_on_ystates(AMREX_INT_ANYD(tyxbx.loVect()), AMREX_INT_ANYD(tyxbx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qym),
                        BL_TO_FORTRAN_ANYD(qmyx),
                        BL_TO_FORTRAN_ANYD(qyp),
                        BL_TO_FORTRAN_ANYD(qpyx),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1),
#ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp1),
#endif
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                        hdt, cdtdx);

      // [lo(1), lo(2)-1, lo(3)], [hi(1), hi(2)+1, hi(3)+1]
      const Box& tzxbx = amrex::grow(zbx, IntVect(AMREX_D_DECL(0,1,0)));

      qmzx.resize(tzxbx, NQ);
      Elixir elix_qmzx = qmzx.elixir();

      qpzx.resize(tzxbx, NQ);
      Elixir elix_qpzx = qpzx.elixir();

#pragma gpu
      transx_on_zstates(AMREX_INT_ANYD(tzxbx.loVect()), AMREX_INT_ANYD(tzxbx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qzm),
                        BL_TO_FORTRAN_ANYD(qmzx),
                        BL_TO_FORTRAN_ANYD(qzp),
                        BL_TO_FORTRAN_ANYD(qpzx),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1),
#ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp1),
#endif
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                        hdt, cdtdx);

      // compute F^y
      // [lo(1)-1, lo(2), lo(3)-1], [hi(1)+1, hi(2)+1, hi(3)+1]
      const Box& cybx = amrex::grow(ybx, IntVect(AMREX_D_DECL(1,0,1)));

      // ftmp1 = fy
      // rftmp1 = rfy
      // qgdnvtmp1 = qgdnvy
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(cybx.loVect()), AMREX_INT_ANYD(cybx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qym),
                          BL_TO_FORTRAN_ANYD(qyp), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp1),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          2, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // [lo(1), lo(2), lo(3)-1], [hi(1)+1, hi(2), lo(3)+1]
      const Box& txybx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,0,1)));

      qmxy.resize(txybx, NQ);
      Elixir elix_qmxy = qmxy.elixir();

      qpxy.resize(txybx, NQ);
      Elixir elix_qpxy = qpxy.elixir();

      // ftmp1 = fy
      // rftmp1 = rfy
      // qgdnvtmp1 = qgdnvy
#pragma gpu
      transy_on_xstates(AMREX_INT_ANYD(txybx.loVect()), AMREX_INT_ANYD(txybx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qxm),
                        BL_TO_FORTRAN_ANYD(qmxy),
                        BL_TO_FORTRAN_ANYD(qxp),
                        BL_TO_FORTRAN_ANYD(qpxy),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1),
#ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp1),
#endif
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                        cdtdy);

      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2), lo(3)+1]
      const Box& tzybx = amrex::grow(zbx, IntVect(AMREX_D_DECL(1,0,0)));

      qmzy.resize(tzybx, NQ);
      Elixir elix_qmzy = qmzy.elixir();

      qpzy.resize(tzybx, NQ);
      Elixir elix_qpzy = qpzy.elixir();

      // ftmp1 = fy
      // rftmp1 = rfy
      // qgdnvtmp1 = qgdnvy
#pragma gpu
      transy_on_zstates(AMREX_INT_ANYD(tzybx.loVect()), AMREX_INT_ANYD(tzybx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qzm),
                        BL_TO_FORTRAN_ANYD(qmzy),
                        BL_TO_FORTRAN_ANYD(qzp),
                        BL_TO_FORTRAN_ANYD(qpzy),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1),
#ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp1),
#endif
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                        cdtdy);

      // compute F^z
      // [lo(1)-1, lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, hi(3)+1]
      const Box& czbx = amrex::grow(zbx, IntVect(AMREX_D_DECL(1,1,0)));

      // ftmp1 = fz
      // rftmp1 = rfz
      // qgdnvtmp1 = qgdnvz
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(czbx.loVect()), AMREX_INT_ANYD(czbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qzm),
                          BL_TO_FORTRAN_ANYD(qzp), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp1),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          3, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // [lo(1)-1, lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, lo(3)]
      const Box& txzbx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,1,0)));

      qmxz.resize(txzbx, NQ);
      Elixir elix_qmxz = qmxz.elixir();

      qpxz.resize(txzbx, NQ);
      Elixir elix_qpxz = qpxz.elixir();

      // ftmp1 = fz
      // rftmp1 = rfz
      // qgdnvtmp1 = qgdnvz
#pragma gpu
      transz_on_xstates(AMREX_INT_ANYD(txzbx.loVect()), AMREX_INT_ANYD(txzbx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qxm),
                        BL_TO_FORTRAN_ANYD(qmxz),
                        BL_TO_FORTRAN_ANYD(qxp),
                        BL_TO_FORTRAN_ANYD(qpxz),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1),
#ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp1),
#endif
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                        cdtdz);

      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2)+1, lo(3)]
      const Box& tyzbx = amrex::grow(ybx, IntVect(AMREX_D_DECL(1,0,0)));

      qmyz.resize(tyzbx, NQ);
      Elixir elix_qmyz = qmyz.elixir();

      qpyz.resize(tyzbx, NQ);
      Elixir elix_qpyz = qpyz.elixir();

      // ftmp1 = fz
      // rftmp1 = rfz
      // qgdnvtmp1 = qgdnvz
#pragma gpu
      transz_on_ystates(AMREX_INT_ANYD(tyzbx.loVect()), AMREX_INT_ANYD(tyzbx.hiVect()),
                        BL_TO_FORTRAN_ANYD(qym),
                        BL_TO_FORTRAN_ANYD(qmyz),
                        BL_TO_FORTRAN_ANYD(qyp),
                        BL_TO_FORTRAN_ANYD(qpyz),
                        BL_TO_FORTRAN_ANYD(qaux[mfi]),
                        BL_TO_FORTRAN_ANYD(ftmp1),
#ifdef RADIATION
                        BL_TO_FORTRAN_ANYD(rftmp1),
#endif
                        BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                        cdtdz);

      // we now have q?zx, q?yx, q?zy, q?xy, q?yz, q?xz

      //
      // Use qx?, q?yz, q?zy to compute final x-flux
      //

      // compute F^{y|z}
      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2)+1, hi(3)]
      const Box& cyzbx = amrex::grow(ybx, IntVect(AMREX_D_DECL(1,0,0)));

      // ftmp1 = fyz
      // rftmp1 = rfyz
      // qgdnvtmp1 = qgdnvyz
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(cyzbx.loVect()), AMREX_INT_ANYD(cyzbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qmyz),
                          BL_TO_FORTRAN_ANYD(qpyz), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp1),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          2, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // compute F^{z|y}
      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)+1]
      const Box& czybx = amrex::grow(zbx, IntVect(AMREX_D_DECL(1,0,0)));

      // ftmp2 = fzy
      // rftmp2 = rfzy
      // qgdnvtmp2 = qgdnvzy
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(czybx.loVect()), AMREX_INT_ANYD(czybx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qmzy),
                          BL_TO_FORTRAN_ANYD(qpzy), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp2),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp2),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp2),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          3, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // compute the corrected x interface states and fluxes
      // [lo(1), lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)]

#pragma gpu
      transyz(AMREX_INT_ANYD(xbx.loVect()), AMREX_INT_ANYD(xbx.hiVect()),
              BL_TO_FORTRAN_ANYD(qxm),
              BL_TO_FORTRAN_ANYD(ql),
              BL_TO_FORTRAN_ANYD(qxp),
              BL_TO_FORTRAN_ANYD(qr),
              BL_TO_FORTRAN_ANYD(qaux[mfi]),
              BL_TO_FORTRAN_ANYD(ftmp1),
#ifdef RADIATION
              BL_TO_FORTRAN_ANYD(rftmp1),
#endif
              BL_TO_FORTRAN_ANYD(ftmp2),
#ifdef RADIATION
              BL_TO_FORTRAN_ANYD(rftmp2),
#endif
              BL_TO_FORTRAN_ANYD(qgdnvtmp1),
              BL_TO_FORTRAN_ANYD(qgdnvtmp2),
              hdt, hdtdy, hdtdz);

#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(cxbx.loVect()), AMREX_INT_ANYD(cxbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(ql),
                          BL_TO_FORTRAN_ANYD(qr), 1, 1,
                          BL_TO_FORTRAN_ANYD(flux[0]),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rad_flux[0]),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qe[0]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          1, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      //
      // Use qy?, q?zx, q?xz to compute final y-flux
      //

      // compute F^{z|x}
      // [lo(1), lo(2)-1, lo(3)], [hi(1), hi(2)+1, hi(3)+1]
      const Box& czxbx = amrex::grow(zbx, IntVect(AMREX_D_DECL(0,1,0)));

      // ftmp1 = fzx
      // rftmp1 = rfzx
      // qgdnvtmp1 = qgdnvzx
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(czxbx.loVect()), AMREX_INT_ANYD(czxbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qmzx),
                          BL_TO_FORTRAN_ANYD(qpzx), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp1),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          3, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // compute F^{x|z}
      // [lo(1), lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, hi(3)]
      const Box& cxzbx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,1,0)));

      // ftmp2 = fxz
      // rftmp2 = rfxz
      // qgdnvtmp2 = qgdnvxz
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(cxzbx.loVect()), AMREX_INT_ANYD(cxzbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qmxz),
                          BL_TO_FORTRAN_ANYD(qpxz), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp2),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp2),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp2),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          1, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // Compute the corrected y interface states and fluxes
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)]

#pragma gpu
      transxz(AMREX_INT_ANYD(ybx.loVect()), AMREX_INT_ANYD(ybx.hiVect()),
              BL_TO_FORTRAN_ANYD(qym),
              BL_TO_FORTRAN_ANYD(ql),
              BL_TO_FORTRAN_ANYD(qyp),
              BL_TO_FORTRAN_ANYD(qr),
              BL_TO_FORTRAN_ANYD(qaux[mfi]),
              BL_TO_FORTRAN_ANYD(ftmp2),
#ifdef RADIATION
              BL_TO_FORTRAN_ANYD(rftmp2),
#endif
              BL_TO_FORTRAN_ANYD(ftmp1),
#ifdef RADIATION
              BL_TO_FORTRAN_ANYD(rftmp1),
#endif
              BL_TO_FORTRAN_ANYD(qgdnvtmp2),
              BL_TO_FORTRAN_ANYD(qgdnvtmp1),
              hdt, hdtdx, hdtdz);

      // Compute the final F^y
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)]
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(ybx.loVect()), AMREX_INT_ANYD(ybx.hiVect()),
                          BL_TO_FORTRAN_ANYD(ql),
                          BL_TO_FORTRAN_ANYD(qr), 1, 1,
                          BL_TO_FORTRAN_ANYD(flux[1]),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rad_flux[1]),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qe[1]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          2, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      //
      // Use qz?, q?xy, q?yx to compute final z-flux
      //

      // compute F^{x|y}
      // [lo(1), lo(2), lo(3)-1], [hi(1)+1, hi(2), hi(3)+1]
      const Box& cxybx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,0,1)));

      // ftmp1 = fxy
      // rftmp1 = rfxy
      // qgdnvtmp1 = qgdnvxy
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(cxybx.loVect()), AMREX_INT_ANYD(cxybx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qmxy),
                          BL_TO_FORTRAN_ANYD(qpxy), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp1),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp1),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp1),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          1, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // compute F^{y|x}
      // [lo(1), lo(2), lo(3)-1], [hi(1), hi(2)+dg(2), hi(3)+1]
      const Box& cyxbx = amrex::grow(ybx, IntVect(AMREX_D_DECL(0,0,1)));

      // ftmp2 = fyx
      // rftmp2 = rfyx
      // qgdnvtmp2 = qgdnvyx
#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(cyxbx.loVect()), AMREX_INT_ANYD(cyxbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qmyx),
                          BL_TO_FORTRAN_ANYD(qpyx), 1, 1,
                          BL_TO_FORTRAN_ANYD(ftmp2),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rftmp2),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qgdnvtmp2),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          2, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      // compute the corrected z interface states and fluxes
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1]

#pragma gpu
      transxy(AMREX_INT_ANYD(zbx.loVect()), AMREX_INT_ANYD(zbx.hiVect()),
              BL_TO_FORTRAN_ANYD(qzm),
              BL_TO_FORTRAN_ANYD(ql),
              BL_TO_FORTRAN_ANYD(qzp),
              BL_TO_FORTRAN_ANYD(qr),
              BL_TO_FORTRAN_ANYD(qaux[mfi]),
              BL_TO_FORTRAN_ANYD(ftmp1),
#ifdef RADIATION
              BL_TO_FORTRAN_ANYD(rftmp1),
#endif
              BL_TO_FORTRAN_ANYD(ftmp2),
#ifdef RADIATION
              BL_TO_FORTRAN_ANYD(rftmp2),
#endif
              BL_TO_FORTRAN_ANYD(qgdnvtmp1),
              BL_TO_FORTRAN_ANYD(qgdnvtmp2),
              hdt, hdtdx, hdtdy);

      // compute the final z fluxes F^z
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1]

#pragma gpu
      cmpflx_plus_godunov(AMREX_INT_ANYD(zbx.loVect()), AMREX_INT_ANYD(zbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(ql),
                          BL_TO_FORTRAN_ANYD(qr), 1, 1,
                          BL_TO_FORTRAN_ANYD(flux[2]),
                          BL_TO_FORTRAN_ANYD(q_int),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(rad_flux[2]),
                          BL_TO_FORTRAN_ANYD(lambda_int),
#endif
                          BL_TO_FORTRAN_ANYD(qe[2]),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          BL_TO_FORTRAN_ANYD(shk),
                          3, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

#endif // 3-d



      // clean the fluxes

      for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

          const Box& nbx = amrex::surroundingNodes(bx, idir);

          int idir_f = idir + 1;

#pragma gpu
          ctu_clean_fluxes(AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
                           idir_f,
                           BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                           BL_TO_FORTRAN_ANYD(q[mfi]),
                           BL_TO_FORTRAN_ANYD(flux[idir]),
#ifdef RADIATION
                           BL_TO_FORTRAN_ANYD(Erborder[mfi]),
                           BL_TO_FORTRAN_ANYD(rad_flux[idir]),
#endif
                           BL_TO_FORTRAN_ANYD(area[idir][mfi]),
                           BL_TO_FORTRAN_ANYD(volume[mfi]),
                           BL_TO_FORTRAN_ANYD(div),
                           AMREX_REAL_ANYD(dx), dt);

      }



      pdivu.resize(bx, 1);
      Elixir elix_pdivu = pdivu.elixir();

      // conservative update

#pragma gpu
      ctu_consup(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                 BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                 BL_TO_FORTRAN_ANYD(q[mfi]),
                 BL_TO_FORTRAN_ANYD(shk),
                 BL_TO_FORTRAN_ANYD(S_new[mfi]),
                 BL_TO_FORTRAN_ANYD(hydro_source[mfi]),
                 BL_TO_FORTRAN_ANYD(flux[0]),
#if AMREX_SPACEDIM >= 2
                 BL_TO_FORTRAN_ANYD(flux[1]),
#endif
#if AMREX_SPACEDIM == 3
                 BL_TO_FORTRAN_ANYD(flux[2]),
#endif
#ifdef RADIATION
                 BL_TO_FORTRAN_ANYD(Erborder[mfi]),
                 BL_TO_FORTRAN_ANYD(Er_new[mfi]),
                 BL_TO_FORTRAN_ANYD(rad_flux[0]),
#if AMREX_SPACEDIM >= 2
                 BL_TO_FORTRAN_ANYD(rad_flux[1]),
#endif
#if AMREX_SPACEDIM == 3
                 BL_TO_FORTRAN_ANYD(rad_flux[2]),
#endif
                 &priv_nstep_fsp,
#endif
                 BL_TO_FORTRAN_ANYD(qe[0]),
#if AMREX_SPACEDIM >= 2
                 BL_TO_FORTRAN_ANYD(qe[1]),
#endif
#if AMREX_SPACEDIM == 3
                 BL_TO_FORTRAN_ANYD(qe[2]),
#endif
                 BL_TO_FORTRAN_ANYD(area[0][mfi]),
#if AMREX_SPACEDIM >= 2
                 BL_TO_FORTRAN_ANYD(area[1][mfi]),
#endif
#if AMREX_SPACEDIM == 3
                 BL_TO_FORTRAN_ANYD(area[2][mfi]),
#endif
                 BL_TO_FORTRAN_ANYD(volume[mfi]),
                 BL_TO_FORTRAN_ANYD(pdivu),
                 AMREX_REAL_ANYD(dx), dt);

#ifdef RADIATION
      nstep_fsp = std::max(nstep_fsp, priv_nstep_fsp);
#endif

#if AMREX_SPACEDIM <= 2
      Array4<Real> pradial_fab = pradial.array();
#endif

      for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

        const Box& nbx = amrex::surroundingNodes(bx, idir);

#pragma gpu
        scale_flux(AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
#if AMREX_SPACEDIM == 1
                   BL_TO_FORTRAN_ANYD(qe[idir]),
#endif
                   BL_TO_FORTRAN_ANYD(flux[idir]),
                   BL_TO_FORTRAN_ANYD(area[idir][mfi]), dt);

#ifdef RADIATION
#pragma gpu
        scale_rad_flux(AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
                       BL_TO_FORTRAN_ANYD(rad_flux[idir]),
                       BL_TO_FORTRAN_ANYD(area[idir][mfi]), dt);
#endif

        if (idir == 0) {
            // get the scaled radial pressure -- we need to treat this specially
            Array4<Real> const qex_fab = qe[idir].array();
            const int prescomp = GDPRES;

#if AMREX_SPACEDIM == 1
            if (!Geometry::IsCartesian()) {
                AMREX_PARALLEL_FOR_3D(nbx, i, j, k,
                {
                    pradial_fab(i,j,k) = qex_fab(i,j,k,prescomp) * dt;
                });
            }
#endif

#if AMREX_SPACEDIM == 2
            if (!mom_flux_has_p[0][0]) {
                AMREX_PARALLEL_FOR_3D(nbx, i, j, k,
                {
                    pradial_fab(i,j,k) = qex_fab(i,j,k,prescomp) * dt;
                });
            }
#endif
        }

        // Store the fluxes from this advance.

        // For normal integration we want to add the fluxes from this advance
        // since we may be subcycling the timestep. But for SDC integration
        // we want to copy the fluxes since we expect that there will not be
        // subcycling and we only want the last iteration's fluxes.

        Array4<Real> const flux_fab = (flux[idir]).array();
        Array4<Real> fluxes_fab = (*fluxes[idir]).array(mfi);
        const int numcomp = NUM_STATE;

        AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(idir), numcomp, i, j, k, n,
        {
#ifndef SDC
            fluxes_fab(i,j,k,n) += flux_fab(i,j,k,n);
#else
            fluxes_fab(i,j,k,n) = flux_fab(i,j,k,n);
#endif
        });

#ifdef RADIATION
        Array4<Real> const rad_flux_fab = (rad_flux[idir]).array();
        Array4<Real> rad_fluxes_fab = (*rad_fluxes[idir]).array(mfi);
        const int radcomp = Radiation::nGroups;

        AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(idir), radcomp, i, j, k, n,
        {
#ifndef SDC
            rad_fluxes_fab(i,j,k,n) += rad_flux_fab(i,j,k,n);
#else
            rad_fluxes_fab(i,j,k,n) = rad_flux_fab(i,j,k,n);
#endif
        });
#endif

        Array4<Real> mass_fluxes_fab = (*mass_fluxes[idir]).array(mfi);
        const int dens_comp = Density;

        AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(idir), 1, i, j, k, n,
        {
            mass_fluxes_fab(i,j,k,0) = flux_fab(i,j,k,dens_comp);
        });

      } // idir loop

#if AMREX_SPACEDIM <= 2
      if (!Geometry::IsCartesian()) {

          Array4<Real> P_radial_fab = P_radial.array(mfi);

          AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(0), 1, i, j, k, n,
          {
#ifndef SDC
              P_radial_fab(i,j,k,0) += pradial_fab(i,j,k,0);
#else
              P_radial_fab(i,j,k,0) = pradial_fab(i,j,k,0);
#endif
          });

      }
#endif

      if (track_grid_losses == 1) {

#ifdef AMREX_USE_CUDA
          Gpu::Device::streamSynchronize();
#endif

          ca_track_grid_losses(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                               D_DECL(BL_TO_FORTRAN_ANYD(flux[0]),
                                      BL_TO_FORTRAN_ANYD(flux[1]),
                                      BL_TO_FORTRAN_ANYD(flux[2])),
                               mass_lost, xmom_lost, ymom_lost, zmom_lost,
                               eden_lost, xang_lost, yang_lost, zang_lost);
      }

    } // MFIter loop

  } // OMP loop


  BL_PROFILE_VAR_STOP(CA_UMDRV);

#ifdef RADIATION
  if (radiation->verbose>=1) {
#ifdef BL_LAZY
    Lazy::QueueReduction( [=] () mutable {
#endif
       ParallelDescriptor::ReduceIntMax(nstep_fsp, ParallelDescriptor::IOProcessorNumber());
       if (ParallelDescriptor::IOProcessor() && nstep_fsp > 0) {
         std::cout << "Radiation f-space advection on level " << level
                   << " takes as many as " << nstep_fsp;
         if (nstep_fsp == 1) {
           std::cout<< " substep.\n";
         }
         else {
           std::cout<< " substeps.\n";
         }
       }
#ifdef BL_LAZY
     });
#endif
  }
#else
  // Flush Fortran output

  if (verbose)
    flush_output();

  if (track_grid_losses)
    {
      material_lost_through_boundary_temp[0] += mass_lost;
      material_lost_through_boundary_temp[1] += xmom_lost;
      material_lost_through_boundary_temp[2] += ymom_lost;
      material_lost_through_boundary_temp[3] += zmom_lost;
      material_lost_through_boundary_temp[4] += eden_lost;
      material_lost_through_boundary_temp[5] += xang_lost;
      material_lost_through_boundary_temp[6] += yang_lost;
      material_lost_through_boundary_temp[7] += zang_lost;
    }

  if (print_update_diagnostics)
    {

      bool local = true;
      Vector<Real> hydro_update = evaluate_source_change(hydro_source, dt, local);

#ifdef BL_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
         ParallelDescriptor::ReduceRealSum(hydro_update.dataPtr(), hydro_update.size(), ParallelDescriptor::IOProcessorNumber());

         if (ParallelDescriptor::IOProcessor())
           std::cout << std::endl << "  Contributions to the state from the hydro source:" << std::endl;

         print_source_change(hydro_update);

#ifdef BL_LAZY
      });
#endif
    }
#endif

  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "... Leaving hydro advance" << std::endl << std::endl;

  if (verbose > 0)
    {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

	if (ParallelDescriptor::IOProcessor())
	  std::cout << "Castro::construct_hydro_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
	});
#endif
    }

}



void
Castro::construct_mol_hydro_source(Real time, Real dt)
{

  BL_PROFILE("Castro::construct_mol_hydro_source()");

  // this constructs the hydrodynamic source (essentially the flux
  // divergence) using method of lines integration.  The output, as a
  // update to the state, is stored in the k_mol array of multifabs.

  const Real strt_time = ParallelDescriptor::second();

  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "... hydro MOL stage " << mol_iteration << std::endl;


  // we'll add each stage's contribution to -div{F(U)} as we compute them
  if (mol_iteration == 0) {
    hydro_source.setVal(0.0);
  }


  const Real *dx = geom.CellSize();

  MultiFab& S_new = get_new_data(State_Type);

  MultiFab& k_stage = *k_mol[mol_iteration];

#ifdef RADIATION
  MultiFab& Er_new = get_new_data(Rad_Type);

  if (!Radiation::rad_hydro_combined) {
    amrex::Abort("Castro::construct_mol_hydro_source -- we don't implement a mode where we have radiation, but it is not coupled to hydro");
  }

  int nstep_fsp = -1;
#endif

  BL_PROFILE_VAR("Castro::advance_hydro_ca_umdrv()", CA_UMDRV);

  const int* domain_lo = geom.Domain().loVect();
  const int* domain_hi = geom.Domain().hiVect();

#ifndef AMREX_USE_CUDA

#ifdef _OPENMP
#ifdef RADIATION
#pragma omp parallel reduction(max:nstep_fsp)
#endif
#endif
  {

    FArrayBox flux[AMREX_SPACEDIM];
#if (AMREX_SPACEDIM <= 2)
    FArrayBox pradial(Box::TheUnitBox(),1);
#endif
#ifdef RADIATION
    FArrayBox rad_flux[AMREX_SPACEDIM];
#endif

#ifdef RADIATION
    int priv_nstep_fsp = -1;
#endif
    // The fourth order stuff cannot do tiling because of the Laplacian corrections
    for (MFIter mfi(S_new, (fourth_order) ? no_tile_size : hydro_tile_size); mfi.isValid(); ++mfi)
      {
	const Box& bx  = mfi.tilebox();

	const int* lo = bx.loVect();
	const int* hi = bx.hiVect();

	FArrayBox &statein  = Sborder[mfi];
	FArrayBox &stateout = S_new[mfi];

	FArrayBox &source_in  = sources_for_hydro[mfi];

	// the output of this will be stored in the correct stage MF
	FArrayBox &source_out = k_stage[mfi];
	FArrayBox &source_hydro_only = hydro_source[mfi];

#ifdef RADIATION
	FArrayBox &Er = Erborder[mfi];
	FArrayBox &lam = lamborder[mfi];
	FArrayBox &Erout = Er_new[mfi];
#endif

	// All cate fabs for fluxes
	for (int i = 0; i < AMREX_SPACEDIM ; i++)  {
	  const Box& bxtmp = amrex::surroundingNodes(bx,i);
	  flux[i].resize(bxtmp,NUM_STATE);
#ifdef RADIATION
	  rad_flux[i].resize(bxtmp,Radiation::nGroups);
#endif
	}

#if (AMREX_SPACEDIM <= 2)
	if (!Geometry::IsCartesian()) {
	  pradial.resize(amrex::surroundingNodes(bx,0),1);
	}
#endif
        if (fourth_order) {
          ca_fourth_single_stage
            (ARLIM_3D(lo), ARLIM_3D(hi), &time, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi),
             &(b_mol[mol_iteration]),
             BL_TO_FORTRAN_ANYD(statein),
             BL_TO_FORTRAN_ANYD(stateout),
             BL_TO_FORTRAN_ANYD(q[mfi]),
             BL_TO_FORTRAN_ANYD(q_bar[mfi]),
             BL_TO_FORTRAN_ANYD(qaux[mfi]),
             BL_TO_FORTRAN_ANYD(qaux_bar[mfi]),
             BL_TO_FORTRAN_ANYD(source_in),
             BL_TO_FORTRAN_ANYD(source_out),
             BL_TO_FORTRAN_ANYD(source_hydro_only),
             ZFILL(dx), &dt,
             D_DECL(BL_TO_FORTRAN_ANYD(flux[0]),
                    BL_TO_FORTRAN_ANYD(flux[1]),
                    BL_TO_FORTRAN_ANYD(flux[2])),
             D_DECL(BL_TO_FORTRAN_ANYD(area[0][mfi]),
                    BL_TO_FORTRAN_ANYD(area[1][mfi]),
                    BL_TO_FORTRAN_ANYD(area[2][mfi])),
#if (AMREX_SPACEDIM < 3)
             BL_TO_FORTRAN_ANYD(pradial),
             BL_TO_FORTRAN_ANYD(dLogArea[0][mfi]),
#endif
             BL_TO_FORTRAN_ANYD(volume[mfi]),
             verbose);

        } else {
          ca_mol_single_stage
            (ARLIM_3D(lo), ARLIM_3D(hi), &time, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi),
             &(b_mol[mol_iteration]),
             BL_TO_FORTRAN_ANYD(statein),
             BL_TO_FORTRAN_ANYD(stateout),
             BL_TO_FORTRAN_ANYD(q[mfi]),
             BL_TO_FORTRAN_ANYD(qaux[mfi]),
             BL_TO_FORTRAN_ANYD(source_in),
             BL_TO_FORTRAN_ANYD(source_out),
             BL_TO_FORTRAN_ANYD(source_hydro_only),
             ZFILL(dx), &dt,
             D_DECL(BL_TO_FORTRAN_ANYD(flux[0]),
                    BL_TO_FORTRAN_ANYD(flux[1]),
                    BL_TO_FORTRAN_ANYD(flux[2])),
             D_DECL(BL_TO_FORTRAN_ANYD(area[0][mfi]),
                    BL_TO_FORTRAN_ANYD(area[1][mfi]),
                    BL_TO_FORTRAN_ANYD(area[2][mfi])),
#if (AMREX_SPACEDIM < 3)
             BL_TO_FORTRAN_ANYD(pradial),
             BL_TO_FORTRAN_ANYD(dLogArea[0][mfi]),
#endif
             BL_TO_FORTRAN_ANYD(volume[mfi]),
             verbose);
        }

	// Store the fluxes from this advance -- we weight them by the
	// integrator weight for this stage
	for (int i = 0; i < AMREX_SPACEDIM ; i++) {
	  (*fluxes    [i])[mfi].saxpy(b_mol[mol_iteration], flux[i],
				      mfi.nodaltilebox(i), mfi.nodaltilebox(i), 0, 0, NUM_STATE);
#ifdef RADIATION
	  (*rad_fluxes[i])[mfi].saxpy(b_mol[mol_iteration], rad_flux[i],
				      mfi.nodaltilebox(i), mfi.nodaltilebox(i), 0, 0, Radiation::nGroups);
#endif
	}

#if (AMREX_SPACEDIM <= 2)
	if (!Geometry::IsCartesian()) {
	  P_radial[mfi].saxpy(b_mol[mol_iteration], pradial,
                              mfi.nodaltilebox(0), mfi.nodaltilebox(0), 0, 0, 1);
	}
#endif
      } // MFIter loop

#ifdef RADIATION
    nstep_fsp = std::max(nstep_fsp, priv_nstep_fsp);
#endif
  }  // end of omp parallel region

#else
  // CUDA version

#ifndef RADIATION

  MultiFab flatn;
  flatn.define(grids, dmap, 1, 1);

  MultiFab div;
  div.define(grids, dmap, 1, 1);

  MultiFab qm;
  qm.define(grids, dmap, AMREX_SPACEDIM*NQ, 2);

  MultiFab qp;
  qp.define(grids, dmap, AMREX_SPACEDIM*NQ, 2);

  MultiFab shk;
  shk.define(grids, dmap, 1, 1);


  MultiFab flux[AMREX_SPACEDIM];
  MultiFab qe[AMREX_SPACEDIM];
  MultiFab qi[AMREX_SPACEDIM];

  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      flux[i].define(getEdgeBoxArray(i), dmap, NUM_STATE, 0);
      qe[i].define(getEdgeBoxArray(i), dmap, NGDNV, 0);
      qi[i].define(getEdgeBoxArray(i), dmap, NQ, 0);
  }


#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      const Box& obx = mfi.growntilebox(1);

      // Compute divergence of velocity field.
#pragma gpu
      divu(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
           BL_TO_FORTRAN_ANYD(q[mfi]),
           AMREX_REAL_ANYD(dx),
           BL_TO_FORTRAN_ANYD(div[mfi]));

      // Compute flattening coefficient for slope calculations.
#pragma gpu
      ca_uflatten
          (AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
           BL_TO_FORTRAN_ANYD(q[mfi]),
           BL_TO_FORTRAN_ANYD(flatn[mfi]), QPRES+1);

      // Do PPM reconstruction to the zone edges.
      int put_on_edges = 1;

#pragma gpu
      ca_ppm_reconstruct
          (AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()), put_on_edges,
           BL_TO_FORTRAN_ANYD(q[mfi]), NQ, 1, NQ,
           BL_TO_FORTRAN_ANYD(flatn[mfi]),
           BL_TO_FORTRAN_ANYD(qm[mfi]),
           BL_TO_FORTRAN_ANYD(qp[mfi]), NQ, 1, NQ);

      // Compute the shk variable
#pragma gpu
      ca_shock
        (AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
         BL_TO_FORTRAN_ANYD(q[mfi]),
         BL_TO_FORTRAN_ANYD(shk[mfi]),
         AMREX_REAL_ANYD(dx));

  } // MFIter loop



#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

          const Box& ebx = mfi.nodaltilebox(idir);

          int idir_f = idir + 1;

#pragma gpu
          ca_construct_flux_cuda
              (AMREX_INT_ANYD(ebx.loVect()), AMREX_INT_ANYD(ebx.hiVect()),
               AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
               AMREX_REAL_ANYD(dx), dt,
               idir_f,
               BL_TO_FORTRAN_ANYD(Sborder[mfi]),
               BL_TO_FORTRAN_ANYD(div[mfi]),
               BL_TO_FORTRAN_ANYD(qaux[mfi]),
               BL_TO_FORTRAN_ANYD(shk[mfi]),
               BL_TO_FORTRAN_ANYD(qm[mfi]),
               BL_TO_FORTRAN_ANYD(qp[mfi]),
               BL_TO_FORTRAN_ANYD(qi[idir][mfi]),
               BL_TO_FORTRAN_ANYD(flux[idir][mfi]),
               BL_TO_FORTRAN_ANYD(area[idir][mfi]));

#pragma gpu
          ca_store_godunov_state
            (AMREX_INT_ANYD(ebx.loVect()), AMREX_INT_ANYD(ebx.hiVect()),
             BL_TO_FORTRAN_ANYD(qi[idir][mfi]),
             BL_TO_FORTRAN_ANYD(qe[idir][mfi]));

          // Store the fluxes from this advance -- we weight them by the
          // integrator weight for this stage
#ifdef AMREX_USE_CUDA
          Cuda::Device::synchronize();  // because saxpy below is run on cpu
#endif
          (*fluxes[idir])[mfi].saxpy(b_mol[mol_iteration], flux[idir][mfi], ebx, ebx, 0, 0, NUM_STATE);

      }

  } // MFIter loop


#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      const Box& bx = mfi.tilebox();

#pragma gpu
      ca_construct_hydro_update_cuda
          (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
           AMREX_REAL_ANYD(dx), dt,
           BL_TO_FORTRAN_ANYD(qe[0][mfi]),
           BL_TO_FORTRAN_ANYD(qe[1][mfi]),
           BL_TO_FORTRAN_ANYD(qe[2][mfi]),
           BL_TO_FORTRAN_ANYD(flux[0][mfi]),
           BL_TO_FORTRAN_ANYD(flux[1][mfi]),
           BL_TO_FORTRAN_ANYD(flux[2][mfi]),
           BL_TO_FORTRAN_ANYD(area[0][mfi]),
           BL_TO_FORTRAN_ANYD(area[1][mfi]),
           BL_TO_FORTRAN_ANYD(area[2][mfi]),
           BL_TO_FORTRAN_ANYD(volume[mfi]),
           BL_TO_FORTRAN_ANYD(sources_for_hydro[mfi]),
           BL_TO_FORTRAN_ANYD(k_stage[mfi]));

  } // MFIter loop

#endif // RADIATION

#endif // CUDA check

  BL_PROFILE_VAR_STOP(CA_UMDRV);

  // Flush Fortran output

  if (verbose)
    flush_output();


  if (print_update_diagnostics)
    {

      bool local = true;
      Vector<Real> hydro_update = evaluate_source_change(k_stage, dt, local);

#ifdef BL_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
	  ParallelDescriptor::ReduceRealSum(hydro_update.dataPtr(), hydro_update.size(), ParallelDescriptor::IOProcessorNumber());

	  if (ParallelDescriptor::IOProcessor())
	    std::cout << std::endl << "  Contributions to the state from the hydro source:" << std::endl;

	  print_source_change(hydro_update);

#ifdef BL_LAZY
	});
#endif
    }

    if (verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

	if (ParallelDescriptor::IOProcessor())
	  std::cout << "Castro::construct_mol_hydro_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
	});
#endif
    }


}



void
Castro::cons_to_prim(const Real time)
{

#ifdef RADIATION
    AmrLevel::FillPatch(*this, Erborder, NUM_GROW, time, Rad_Type, 0, Radiation::nGroups);

    MultiFab lamborder(grids, dmap, Radiation::nGroups, NUM_GROW);
    if (radiation->pure_hydro) {
      lamborder.setVal(0.0, NUM_GROW);
    }
    else {
      radiation->compute_limiter(level, grids, Sborder, Erborder, lamborder);
    }
#endif

    MultiFab& S_new = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

        const Box& qbx = mfi.growntilebox(NUM_GROW);

        // Convert the conservative state to the primitive variable state.
        // This fills both q and qaux.

#pragma gpu
        ca_ctoprim(AMREX_INT_ANYD(qbx.loVect()), AMREX_INT_ANYD(qbx.hiVect()),
                   BL_TO_FORTRAN_ANYD(Sborder[mfi]),
#ifdef RADIATION
                   BL_TO_FORTRAN_ANYD(Erborder[mfi]),
                   BL_TO_FORTRAN_ANYD(lamborder[mfi]),
#endif
                   BL_TO_FORTRAN_ANYD(q[mfi]),
                   BL_TO_FORTRAN_ANYD(qaux[mfi]));

        // Convert the source terms expressed as sources to the conserved state to those
        // expressed as sources for the primitive state.
        if (time_integration_method == CornerTransportUpwind) {
#pragma gpu
            ca_srctoprim(BL_TO_FORTRAN_BOX(qbx),
                         BL_TO_FORTRAN_ANYD(q[mfi]),
                         BL_TO_FORTRAN_ANYD(qaux[mfi]),
                         BL_TO_FORTRAN_ANYD(sources_for_hydro[mfi]),
                         BL_TO_FORTRAN_ANYD(src_q[mfi]));
        }

#ifndef RADIATION

        // Add in the reactions source term; only done in SDC.

#ifdef SDC
#ifdef REACTIONS
        MultiFab& SDC_react_source = get_new_data(SDC_React_Type);

        if (do_react)
	    src_q[mfi].plus(SDC_react_source[mfi],qbx,qbx,0,0,QVAR);
#endif
#endif
#endif

    }

}


#ifndef AMREX_USE_CUDA
void
Castro::cons_to_prim_fourth(const Real time)
{
  // convert the conservative state cell averages to primitive cell
  // averages with 4th order accuracy

    MultiFab& S_new = get_new_data(State_Type);

    // we don't support radiation here
#ifdef RADIATION
    amrex::Abort("radiation not supported to fourth order");
#else
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      const Box& qbx = mfi.growntilebox(NUM_GROW);
      const Box& qbxm1 = mfi.growntilebox(NUM_GROW-1);

      // note: these conversions are using a growntilebox, so it
      // will include ghost cells

      // convert U_avg to U_cc -- this will use a Laplacian
      // operation and will result in U_cc defined only on
      // NUM_GROW-1 ghost cells at the end.
      FArrayBox U_cc;
      U_cc.resize(qbx, NUM_STATE);

      ca_make_cell_center(BL_TO_FORTRAN_BOX(qbxm1),
                          BL_TO_FORTRAN_FAB(Sborder[mfi]),
                          BL_TO_FORTRAN_FAB(U_cc));

      // convert U_avg to q_bar -- this will be done on all NUM_GROW
      // ghost cells.
      ca_ctoprim(BL_TO_FORTRAN_BOX(qbx),
                 BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                 BL_TO_FORTRAN_ANYD(q_bar[mfi]),
                 BL_TO_FORTRAN_ANYD(qaux_bar[mfi]));

      // this is what we should construct the flattening coefficient
      // from

      // convert U_cc to q_cc (we'll store this temporarily in q,
      // qaux).  This will remain valid only on the NUM_GROW-1 ghost
      // cells.
      ca_ctoprim(BL_TO_FORTRAN_BOX(qbxm1),
                 BL_TO_FORTRAN_ANYD(U_cc),
                 BL_TO_FORTRAN_ANYD(q[mfi]),
                 BL_TO_FORTRAN_ANYD(qaux[mfi]));
    }


    // check for NaNs
    check_for_nan(q);
    check_for_nan(q_bar);


#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      const Box& qbxm1 = mfi.growntilebox(NUM_GROW-1);

      // now convert q, qaux into 4th order accurate averages
      // this will create q, qaux in NUM_GROW-1 ghost cells, but that's
      // we need here

      ca_make_fourth_average(BL_TO_FORTRAN_BOX(qbxm1),
                             BL_TO_FORTRAN_FAB(q[mfi]),
                             BL_TO_FORTRAN_FAB(q_bar[mfi]));

      // not sure if we need to convert qaux this way, or if we can
      // just evaluate it (we may not need qaux at all actually)
      ca_make_fourth_average(BL_TO_FORTRAN_BOX(qbxm1),
                             BL_TO_FORTRAN_FAB(qaux[mfi]),
                             BL_TO_FORTRAN_FAB(qaux_bar[mfi]));

    }

    check_for_nan(q_bar);
#endif // RADIATION
}
#endif



void
Castro::check_for_cfl_violation(const Real dt)
{

    Real courno = -1.0e+200;

    const Real *dx = geom.CellSize();

    MultiFab& S_new = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel reduction(max:courno)
#endif
    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

#pragma gpu
        ca_compute_cfl(BL_TO_FORTRAN_BOX(bx),
                       BL_TO_FORTRAN_ANYD(q[mfi]),
                       BL_TO_FORTRAN_ANYD(qaux[mfi]),
                       dt, AMREX_REAL_ANYD(dx), AMREX_MFITER_REDUCE_MAX(&courno), print_fortran_warnings);

    }

    ParallelDescriptor::ReduceRealMax(courno);

    if (courno > 1.0) {
        amrex::Print() << "WARNING -- EFFECTIVE CFL AT LEVEL " << level << " IS " << courno << std::endl << std::endl;

        cfl_violation = 1;
    }

}
