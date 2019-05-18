#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

using namespace amrex;

void
Castro::construct_ctu_hydro_source(Real time, Real dt)
{

  BL_PROFILE("Castro::construct_ctu_hydro_source()");

  const Real strt_time = ParallelDescriptor::second();

  // this constructs the hydrodynamic source (essentially the flux
  // divergence) using the CTU framework for unsplit hydrodynamics

  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "... Entering construct_ctu_hydro_source" << std::endl << std::endl;

  hydro_source.setVal(0.0);

  int finest_level = parent->finestLevel();

  const Real *dx = geom.CellSize();

  const int* domain_lo = geom.Domain().loVect();
  const int* domain_hi = geom.Domain().hiVect();

  MultiFab& S_new = get_new_data(State_Type);

#ifdef RADIATION
  MultiFab& Er_new = get_new_data(Rad_Type);

  if (!Radiation::rad_hydro_combined) {
    amrex::Abort("Castro::construct_ctu_hydro_source -- we don't implement a mode where we have radiation, but it is not coupled to hydro");
  }

  int nstep_fsp = -1;
#endif

  Real mass_lost = 0.;
  Real xmom_lost = 0.;
  Real ymom_lost = 0.;
  Real zmom_lost = 0.;
  Real eden_lost = 0.;
  Real xang_lost = 0.;
  Real yang_lost = 0.;
  Real zang_lost = 0.;

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
    FArrayBox dq_core, dq_pass;
    FArrayBox Ip_core, Im_core, Ip_core_src, Im_core_src, Ip_gc, Im_gc;
    FArrayBox Ip_pass, Im_pass, Ip_pass_src, Im_pass_src;
#ifdef RADIATION
    FArrayBox Ip_rad, Im_rad;
#endif
    FArrayBox sm, sp;
    FArrayBox shk;

    FArrayBox qxm_core, qxp_core;
    FArrayBox qxm_pass, qxp_pass;
#ifdef RADIATION
    FArrayBox qxm_rad, qxp_rad;
#endif

#if AMREX_SPACEDIM >= 2
    FArrayBox qym_core, qyp_core;
    FArrayBox qym_pass, qyp_pass;
#ifdef RADIATION
    FArrayBox qym_rad, qyp_rad;
#endif
#endif

#if AMREX_SPACEDIM == 3
    FArrayBox qzm_core, qzp_core;
    FArrayBox qzm_pass, qzp_pass;
#ifdef RADIATION
    FArrayBox qzm_rad, qzp_rad;
#endif
#endif

    FArrayBox div;
    FArrayBox q_int_core;
    FArrayBox q_int_pass;
#ifdef RADIATION
    FArrayBox q_int_rad;
    FArrayBox lambda_int;
#endif

#if AMREX_SPACEDIM >= 2
    FArrayBox ftmp1, ftmp2;
#ifdef RADIATION
    FArrayBox rftmp1, rftmp2;
#endif
    FArrayBox qgdnvtmp1, qgdnvtmp2;

    FArrayBox ql_core, qr_core;
    FArrayBox ql_pass, qr_pass;
#ifdef RADIATION
    FArrayBox ql_rad, qr_rad;
#endif
#endif

    FArrayBox flux[AMREX_SPACEDIM], qe[AMREX_SPACEDIM];
#ifdef RADIATION
    FArrayBox rad_flux[AMREX_SPACEDIM];
#endif
#if AMREX_SPACEDIM <= 2
    FArrayBox pradial;
#endif
#if AMREX_SPACEDIM == 3
    FArrayBox qmyx_core, qpyx_core;
    FArrayBox qmyx_pass, qpyx_pass;
#ifdef RADIATION
    FArrayBox qmyx_rad, qpyx_rad;
#endif

    FArrayBox qmzx_core, qpzx_core;
    FArrayBox qmzx_pass, qpzx_pass;
#ifdef RADIATION
    FArrayBox qmzx_rad, qpzx_rad;
#endif

    FArrayBox qmxy_core, qpxy_core;
    FArrayBox qmxy_pass, qpxy_pass;
#ifdef RADIATION
    FArrayBox qmxy_rad, qpxy_rad;
#endif

    FArrayBox qmzy_core, qpzy_core;
    FArrayBox qmzy_pass, qpzy_pass;
#ifdef RADIATION
    FArrayBox qmzy_rad, qpzy_rad;
#endif

    FArrayBox qmxz_core, qpxz_core;
    FArrayBox qmxz_pass, qpxz_pass;
#ifdef RADIATION
    FArrayBox qmxz_rad, qpxz_rad;
#endif

    FArrayBox qmyz_core, qpyz_core;
    FArrayBox qmyz_pass, qpyz_pass;
#ifdef RADIATION
    FArrayBox qmyz_rad, qpyz_rad;
#endif

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

      Array4<Real> const flatn_arr = flatn.array();

      if (first_order_hydro == 1) {
        AMREX_PARALLEL_FOR_3D(obx, i, j, k, { flatn_arr(i,j,k) = 0.0; });
      } else if (use_flattening == 1) {
#ifdef RADIATION
        ca_rad_flatten(ARLIM_3D(obx.loVect()), ARLIM_3D(obx.hiVect()),
                       BL_TO_FORTRAN_ANYD(q_core[mfi]),
                       BL_TO_FORTRAN_ANYD(q_rad[mfi]),
                       BL_TO_FORTRAN_ANYD(flatn),
                       BL_TO_FORTRAN_ANYD(flatg));
#else
#pragma gpu box(obx)
        ca_uflatten(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
                    BL_TO_FORTRAN_ANYD(q_core[mfi]),
                    BL_TO_FORTRAN_ANYD(q_core[mfi]), NQC, QPRES+1,
                    BL_TO_FORTRAN_ANYD(flatn));
#endif
      } else {
        AMREX_PARALLEL_FOR_3D(obx, i, j, k, { flatn_arr(i,j,k) = 1.0; });
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

      shk.resize(obx, 1);
      Elixir elix_shk = shk.elixir();

      // qxm

      qxm_core.resize(obx, NQC);
      Elixir elix_qxm_core = qxm_core.elixir();

      qxm_pass.resize(obx, NQP);
      Elixir elix_qxm_pass = qxm_pass.elixir();

#ifdef RADIATION
      qxm_rad.resize(obx, NQR);
      Elixir elix_qxm_rad = qxm_rad.elixir();
#endif

      // qxp

      qxp_core.resize(obx, NQC);
      Elixir elix_qxp_core = qxp_core.elixir();

      qxp_pass.resize(obx, NQP);
      Elixir elix_qxp_pass = qxp_pass.elixir();

#ifdef RADIATION
      qxp_rad.resize(obx, NQR);
      Elixir elix_qxp_rad = qxp_rad.elixir();
#endif

#if AMREX_SPACEDIM >= 2
      // qym

      qym_core.resize(obx, NQC);
      Elixir elix_qym_core = qym_core.elixir();

      qym_pass.resize(obx, NQP);
      Elixir elix_qym_pass = qym_pass.elixir();

#ifdef RADIATION
      qym_rad.resize(obx, NQR);
      Elixir elix_qym_rad = qym_rad.elixir();
#endif

      // qyp

      qyp_core.resize(obx, NQC);
      Elixir elix_qyp_core = qyp_core.elixir();

      qyp_pass.resize(obx, NQP);
      Elixir elix_qyp_pass = qyp_pass.elixir();

#ifdef RADIATION
      qyp_rad.resize(obx, NQR);
      Elixir elix_qyp_rad = qyp_rad.elixir();
#endif
#endif

#if AMREX_SPACEDIM == 3
      // qzm

      qzm_core.resize(obx, NQC);
      Elixir elix_qzm_core = qzm_core.elixir();

      qzm_pass.resize(obx, NQP);
      Elixir elix_qzm_pass = qzm_pass.elixir();

#ifdef RADIATION
      qzm_rad.resize(obx, NQR);
      Elixir elix_qzm_rad = qzm_rad.elixir();
#endif

      // qzp

      qzp_core.resize(obx, NQC);
      Elixir elix_qzp_core = qzp_core.elixir();

      qzp_pass.resize(obx, NQP);
      Elixir elix_qzp_pass = qzp_pass.elixir();

#ifdef RADIATION
      qzp_rad.resize(obx, NQR);
      Elixir elix_qzp_rad = qzp_rad.elixir();
#endif
#endif

      if (ppm_type == 0) {

        dq_core.resize(obx, NQC);
        Elixir elix_dq_core = dq_core.elixir();

        dq_pass.resize(obx, NQP);
        Elixir elix_dq_pass = dq_pass.elixir();

#pragma gpu box(obx)
        ctu_plm_states(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
                       AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                       BL_TO_FORTRAN_ANYD(q_core[mfi]),
                       BL_TO_FORTRAN_ANYD(q_pass[mfi]),
                       BL_TO_FORTRAN_ANYD(flatn),
                       BL_TO_FORTRAN_ANYD(qaux[mfi]),
                       BL_TO_FORTRAN_ANYD(q_core_src[mfi]),
#ifdef PRIM_SPECIES_HAVE_SOURCES
                       BL_TO_FORTRAN_ANYD(q_pass_src[mfi]),
#endif
                       BL_TO_FORTRAN_ANYD(shk),
                       BL_TO_FORTRAN_ANYD(dq_core),
                       BL_TO_FORTRAN_ANYD(dq_pass),
                       BL_TO_FORTRAN_ANYD(qxm_core),
                       BL_TO_FORTRAN_ANYD(qxp_core),
                       BL_TO_FORTRAN_ANYD(qxm_pass),
                       BL_TO_FORTRAN_ANYD(qxp_pass),
#if AMREX_SPACEDIM >= 2
                       BL_TO_FORTRAN_ANYD(qym_core),
                       BL_TO_FORTRAN_ANYD(qyp_core),
                       BL_TO_FORTRAN_ANYD(qym_pass),
                       BL_TO_FORTRAN_ANYD(qyp_pass),
#endif
#if AMREX_SPACEDIM == 3
                       BL_TO_FORTRAN_ANYD(qzm_core),
                       BL_TO_FORTRAN_ANYD(qzp_core),
                       BL_TO_FORTRAN_ANYD(qzm_pass),
                       BL_TO_FORTRAN_ANYD(qzp_pass),
#endif
                       AMREX_REAL_ANYD(dx), dt,
#if (AMREX_SPACEDIM < 3)
                       BL_TO_FORTRAN_ANYD(dLogArea[0][mfi]),
#endif
                       AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      } else {

        // first the core and radiation variables
        Ip_core.resize(obx, 3*NQC);
        Elixir elix_Ip_core = Ip_core.elixir();

        Ip_pass.resize(obx, 3*NQP);
        Elixir elix_Ip_pass = Ip_pass.elixir();

        Im_core.resize(obx, 3*NQC);
        Elixir elix_Im_core = Im_core.elixir();

        Im_pass.resize(obx, 3*NQP);
        Elixir elix_Im_pass = Im_pass.elixir();

        Ip_core_src.resize(obx, 3*NQC_SRC);
        Elixir elix_Ip_core_src = Ip_core_src.elixir();

        Im_core_src.resize(obx, 3*NQC_SRC);
        Elixir elix_Im_core_src = Im_core_src.elixir();

#ifdef PRIM_SPECIES_HAVE_SOURCES
        Ip_pass_src.resize(obx, 3*NQP_SRC);
        Elixir elix_Ip_pass_src = Ip_pass_src.elixir();

        Im_pass_src.resize(obx, 3*NQP_SRC);
        Elixir elix_Im_pass_src = Im_pass_src.elixir();
#endif

        Ip_gc.resize(obx, 3);
        Elixir elix_Ip_gc = Ip_gc.elixir();

        Im_gc.resize(obx, 3);
        Elixir elix_Im_gc = Im_gc.elixir();

        sm.resize(obx, AMREX_SPACEDIM);
        Elixir elix_sm = sm.elixir();

        sp.resize(obx, AMREX_SPACEDIM);
        Elixir elix_sp = sp.elixir();

#pragma gpu box(obx)
        ctu_ppm_states(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
                       AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                       BL_TO_FORTRAN_ANYD(q_core[mfi]),
                       BL_TO_FORTRAN_ANYD(q_pass[mfi]),
#ifdef RADIATION
                       BL_TO_FORTRAN_ANYD(q_rad[mfi]),
#endif
                       BL_TO_FORTRAN_ANYD(flatn),
                       BL_TO_FORTRAN_ANYD(qaux[mfi]),
                       BL_TO_FORTRAN_ANYD(q_core_src[mfi]),
#ifdef PRIM_SPECIES_HAVE_SOURCES
                       BL_TO_FORTRAN_ANYD(q_pass_src[mfi]),
#endif
                       BL_TO_FORTRAN_ANYD(shk),
                       BL_TO_FORTRAN_ANYD(Ip_core),
                       BL_TO_FORTRAN_ANYD(Im_core),
                       BL_TO_FORTRAN_ANYD(Ip_pass),
                       BL_TO_FORTRAN_ANYD(Im_pass),
#ifdef RADIATION
                       BL_TO_FORTRAN_ANYD(Ip_rad),
                       BL_TO_FORTRAN_ANYD(Im_rad),
#endif
                       BL_TO_FORTRAN_ANYD(Ip_core_src),
                       BL_TO_FORTRAN_ANYD(Im_core_src),
#ifdef PRIM_SPECIES_HAVE_SOURCES
                       BL_TO_FORTRAN_ANYD(Ip_pass_src),
                       BL_TO_FORTRAN_ANYD(Im_pass_src),
#endif
                       BL_TO_FORTRAN_ANYD(Ip_gc),
                       BL_TO_FORTRAN_ANYD(Im_gc),
                       BL_TO_FORTRAN_ANYD(sm),
                       BL_TO_FORTRAN_ANYD(sp),
                       BL_TO_FORTRAN_ANYD(qxm_core),
                       BL_TO_FORTRAN_ANYD(qxp_core),
                       BL_TO_FORTRAN_ANYD(qxm_pass),
                       BL_TO_FORTRAN_ANYD(qxp_pass),
#ifdef RADIATION
                       BL_TO_FORTRAN_ANYD(qxm_rad),
                       BL_TO_FORTRAN_ANYD(qxp_rad),
#endif
#if AMREX_SPACEDIM >= 2
                       BL_TO_FORTRAN_ANYD(qym_core),
                       BL_TO_FORTRAN_ANYD(qyp_core),
                       BL_TO_FORTRAN_ANYD(qym_pass),
                       BL_TO_FORTRAN_ANYD(qyp_pass),
#ifdef RADIATION
                       BL_TO_FORTRAN_ANYD(qym_rad),
                       BL_TO_FORTRAN_ANYD(qyp_rad),
#endif
#endif
#if AMREX_SPACEDIM == 3
                       BL_TO_FORTRAN_ANYD(qzm_core),
                       BL_TO_FORTRAN_ANYD(qzp_core),
                       BL_TO_FORTRAN_ANYD(qzm_pass),
                       BL_TO_FORTRAN_ANYD(qzp_pass),
#ifdef RADIATION
                       BL_TO_FORTRAN_ANYD(qzp_rad),
                       BL_TO_FORTRAN_ANYD(qzp_rad),
#endif
#endif
                       AMREX_REAL_ANYD(dx), dt,
#if (AMREX_SPACEDIM < 3)
                       BL_TO_FORTRAN_ANYD(dLogArea[0][mfi]),
#endif
                       AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));
      }

      div.resize(obx, 1);
      Elixir elix_div = div.elixir();

      // compute divu -- we'll use this later when doing the artifical viscosity
#pragma gpu box(obx)
      divu(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
           BL_TO_FORTRAN_ANYD(q_core[mfi]),
           AMREX_REAL_ANYD(dx),
           BL_TO_FORTRAN_ANYD(div));

      q_int_core.resize(obx, NQC);
      Elixir elix_q_int_core = q_int_core.elixir();

      q_int_pass.resize(obx, NQP);
      Elixir elix_q_int_pass = q_int_pass.elixir();

#ifdef RADIATION
      q_int_rad.resize(obx, NQR);
      Elixir elix_q_int_rad = q_int_rad.elixir();

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
      }
      Elixir elix_pradial = pradial.elixir();
#endif

#if AMREX_SPACEDIM == 1
#pragma gpu box(xbx)
      cmpflx_plus_godunov(AMREX_INT_ANYD(xbx.loVect()), AMREX_INT_ANYD(xbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qxm_core),
                          BL_TO_FORTRAN_ANYD(qxp_core),
                          BL_TO_FORTRAN_ANYD(qxm_pass),
                          BL_TO_FORTRAN_ANYD(qxp_pass),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(qxm_rad),
                          BL_TO_FORTRAN_ANYD(qxp_rad),
#endif
                          1, 1,
                          BL_TO_FORTRAN_ANYD(flux[0]),
                          BL_TO_FORTRAN_ANYD(q_int_core),
                          BL_TO_FORTRAN_ANYD(q_int_pass),
#ifdef RADIATION
                          BL_TO_FORTRAN_ANYD(q_int_rad),
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

      ql_core.resize(obx, NQC);
      Elixir elix_ql_core = ql_core.elixir();

      ql_pass.resize(obx, NQP);
      Elixir elix_ql_pass = ql_pass.elixir();

#ifdef RADIATION
      ql_rad.resize(obx, NQR);
      Elixir elix_ql_rad = ql_rad.elixir();
#endif

      qr_core.resize(obx, NQC);
      Elixir elix_qr_core = qr_core.elixir();

      qr_pass.resize(obx, NQP);
      Elixir elix_qr_pass = qr_pass.elixir();

#ifdef RADIATION
      qr_rad.resize(obx, NQR);
      Elixir elix_qr_rad = qr_rad.elixir();
#endif

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
#pragma gpu box(cxbx)
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
#pragma gpu box(cybx)
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
#pragma gpu box(xbx)
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

#pragma gpu box(xbx)
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

#pragma gpu box(ybx)
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

#pragma gpu box(ybx)
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
#pragma gpu box(cxbx)
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
#pragma gpu box(tyxbx)
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

#pragma gpu box(tzxbx)
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
#pragma gpu box(cybx)
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
#pragma gpu box(txybx)
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
#pragma gpu box(tzybx)
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
#pragma gpu box(czbx)
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
#pragma gpu box(txzbx)
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
#pragma gpu box(tyzbx)
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
#pragma gpu box(cyzbx)
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
#pragma gpu box(czybx)
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

#pragma gpu box(xbx)
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

#pragma gpu box(xbx)
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

      //
      // Use qy?, q?zx, q?xz to compute final y-flux
      //

      // compute F^{z|x}
      // [lo(1), lo(2)-1, lo(3)], [hi(1), hi(2)+1, hi(3)+1]
      const Box& czxbx = amrex::grow(zbx, IntVect(AMREX_D_DECL(0,1,0)));

      // ftmp1 = fzx
      // rftmp1 = rfzx
      // qgdnvtmp1 = qgdnvzx
#pragma gpu box(czxbx)
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
#pragma gpu box(cxzbx)
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

#pragma gpu box(ybx)
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
#pragma gpu box(ybx)
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
#pragma gpu box(cxybx)
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
#pragma gpu box(cyxbx)
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

#pragma gpu box(zbx)
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

#pragma gpu box(zbx)
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

          Array4<Real> const flux_arr = (flux[idir]).array();
          const int temp_comp = Temp;
#ifdef SHOCK_VAR
          const int shk_comp = Shock;
#endif

          // Zero out shock and temp fluxes -- these are physically meaningless here
          AMREX_PARALLEL_FOR_3D(nbx, i, j, k,
          {
              flux_arr(i,j,k,temp_comp) = 0.e0;
#ifdef SHOCK_VAR
              flux_arr(i,j,k,shk_comp) = 0.e0;
#endif
          });

#pragma gpu box(nbx)
          apply_av(AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
                   idir_f, AMREX_REAL_ANYD(dx),
                   BL_TO_FORTRAN_ANYD(div),
                   BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                   BL_TO_FORTRAN_ANYD(flux[idir]));

#ifdef RADIATION
#pragma gpu box(nbx)
          apply_av_rad(AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
                       idir_f, AMREX_REAL_ANYD(dx),
                       BL_TO_FORTRAN_ANYD(div),
                       BL_TO_FORTRAN_ANYD(Erborder[mfi]),
                       BL_TO_FORTRAN_ANYD(rad_flux[idir]));
#endif

          if (limit_fluxes_on_small_dens == 1) {
#pragma gpu box(nbx)
              limit_hydro_fluxes_on_small_dens
                  (AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
                   idir_f,
                   BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                   BL_TO_FORTRAN_ANYD(q_core[mfi]),
                   BL_TO_FORTRAN_ANYD(volume[mfi]),
                   BL_TO_FORTRAN_ANYD(flux[idir]),
                   BL_TO_FORTRAN_ANYD(area[idir][mfi]),
                   dt, AMREX_REAL_ANYD(dx));
          }

#pragma gpu box(nbx)
          normalize_species_fluxes(AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
                                   BL_TO_FORTRAN_ANYD(flux[idir]));

      }



      pdivu.resize(bx, 1);
      Elixir elix_pdivu = pdivu.elixir();

      // conservative update

#pragma gpu box(bx)
      ctu_consup(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                 BL_TO_FORTRAN_ANYD(Sborder[mfi]),
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

#pragma gpu box(nbx)
        scale_flux(AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
#if AMREX_SPACEDIM == 1
                   BL_TO_FORTRAN_ANYD(qe[idir]),
#endif
                   BL_TO_FORTRAN_ANYD(flux[idir]),
                   BL_TO_FORTRAN_ANYD(area[idir][mfi]), dt);

#ifdef RADIATION
#pragma gpu box(nbx)
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
        // since we may be subcycling the timestep. But for simplified SDC integration
        // we want to copy the fluxes since we expect that there will not be
        // subcycling and we only want the last iteration's fluxes.

        Array4<Real> const flux_fab = (flux[idir]).array();
        Array4<Real> fluxes_fab = (*fluxes[idir]).array(mfi);
        const int numcomp = NUM_STATE;

        if (time_integration_method == SimplifiedSpectralDeferredCorrections) {

            AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(idir), numcomp, i, j, k, n,
            {
                fluxes_fab(i,j,k,n) = flux_fab(i,j,k,n);
            });

        } else {

            AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(idir), numcomp, i, j, k, n,
            {
                fluxes_fab(i,j,k,n) += flux_fab(i,j,k,n);
            });

        }

#ifdef RADIATION
        Array4<Real> const rad_flux_fab = (rad_flux[idir]).array();
        Array4<Real> rad_fluxes_fab = (*rad_fluxes[idir]).array(mfi);
        const int radcomp = Radiation::nGroups;

        if (time_integration_method == SimplifiedSpectralDeferredCorrections) {

            AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(idir), radcomp, i, j, k, n,
            {
                rad_fluxes_fab(i,j,k,n) = rad_flux_fab(i,j,k,n);
            });

        } else {

            AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(idir), radcomp, i, j, k, n,
            {
                rad_fluxes_fab(i,j,k,n) += rad_flux_fab(i,j,k,n);
            });

        }
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

          if (time_integration_method == SimplifiedSpectralDeferredCorrections) {

              AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(0), 1, i, j, k, n,
              {
                  P_radial_fab(i,j,k,0) = pradial_fab(i,j,k,0);
              });

          } else {

              AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(0), 1, i, j, k, n,
              {
                  P_radial_fab(i,j,k,0) += pradial_fab(i,j,k,0);
              });

          }

      }
#endif

      if (track_grid_losses == 1) {

#pragma gpu box(bx)
          ca_track_grid_losses(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                               BL_TO_FORTRAN_ANYD(flux[0]),
#if AMREX_SPACEDIM >= 2
                               BL_TO_FORTRAN_ANYD(flux[1]),
#endif
#if AMREX_SPACEDIM == 3
                               BL_TO_FORTRAN_ANYD(flux[2]),
#endif
                               AMREX_MFITER_REDUCE_SUM(&mass_lost),
                               AMREX_MFITER_REDUCE_SUM(&xmom_lost),
                               AMREX_MFITER_REDUCE_SUM(&ymom_lost),
                               AMREX_MFITER_REDUCE_SUM(&zmom_lost),
                               AMREX_MFITER_REDUCE_SUM(&eden_lost),
                               AMREX_MFITER_REDUCE_SUM(&xang_lost),
                               AMREX_MFITER_REDUCE_SUM(&yang_lost),
                               AMREX_MFITER_REDUCE_SUM(&zang_lost));
      }

    } // MFIter loop

  } // OMP loop

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
    std::cout << "... Leaving construct_ctu_hydro_sources()" << std::endl << std::endl;

  if (verbose > 0)
    {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

	if (ParallelDescriptor::IOProcessor())
	  std::cout << "Castro::construct_ctu_hydro_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
	});
#endif
    }

}
