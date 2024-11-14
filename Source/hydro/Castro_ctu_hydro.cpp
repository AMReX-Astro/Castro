#include <Castro.H>
#include <Castro_util.H>

#ifdef RADIATION
#include <Radiation.H>
#endif

#ifdef HYBRID_MOMENTUM
#include <hybrid.H>
#endif

#include <advection_util.H>

using namespace amrex;

advance_status
Castro::construct_ctu_hydro_source(Real time, Real dt)  // NOLINT(readability-convert-member-functions-to-static)
{
  amrex::ignore_unused(time);
  amrex::ignore_unused(dt);

  advance_status status {};

  if (!do_hydro) {
      return status;
  }

#ifndef TRUE_SDC

  BL_PROFILE("Castro::construct_ctu_hydro_source()");

  const Real strt_time = ParallelDescriptor::second();

  // this constructs the hydrodynamic source (essentially the flux
  // divergence) using the CTU framework for unsplit hydrodynamics

  if (verbose) {
      amrex::Print() << "... Entering construct_ctu_hydro_source() on level " << level << std::endl << std::endl;
  }

#ifdef HYBRID_MOMENTUM
  GeometryData geomdata = geom.data();
#endif

#if AMREX_SPACEDIM <= 2
  int coord = geom.Coord();
#endif

#if AMREX_SPACEDIM >= 2
  const Real *dx = geom.CellSize();
#endif

  MultiFab& S_new = get_new_data(State_Type);

#ifdef SIMPLIFIED_SDC
#ifdef REACTIONS
  MultiFab& SDC_react_source = get_new_data(Simplified_SDC_React_Type);
#endif
#endif

  // we will treat the hydro source as any other source term

#ifdef RADIATION
  MultiFab& Er_new = get_new_data(Rad_Type);

  if (!Radiation::rad_hydro_combined) {
    amrex::Abort("Castro::construct_ctu_hydro_source -- we don't implement a mode where we have radiation, but it is not coupled to hydro");
  }

  int nstep_fsp = -1;

  AmrLevel::FillPatch(*this, Erborder, NUM_GROW, time, Rad_Type, 0, Radiation::nGroups);

  MultiFab lamborder(grids, dmap, Radiation::nGroups, NUM_GROW);
  if (radiation->pure_hydro) {
      lamborder.setVal(0.0, NUM_GROW);
  }
  else {
      radiation->compute_limiter(level, grids, Sborder, Erborder, lamborder);
  }
#endif

  // Record a running total of the number of bytes allocated as temporary Fab data.

  size_t fab_size = 0;
#ifdef AMREX_USE_GPU
  size_t mf_size = 0;
#endif
  IntVect maximum_tile_size{0};

  // Our strategy for launching work on GPUs in the hydro is incompatible with OpenMP,
  // throw an error if the user is trying this. If this were to ever change, it would
  // require at minimum doing a safe atomic update on the running fab_size total.

#if defined(AMREX_USE_OMP) && defined(AMREX_USE_GPU)
  amrex::Error("USE_OMP=TRUE and USE_GPU=TRUE are not concurrently supported in Castro");
#endif

#ifdef AMREX_USE_GPU
   if (castro::hydro_memory_footprint_ratio > 0.0) {
       // If we haven't done any tuning yet, set the tile size to an arbitrary
       // small value to start with.

       if (hydro_tile_size_has_been_tuned == 0) {
           hydro_tile_size[0] = 16;
#if AMREX_SPACEDIM >= 2
           hydro_tile_size[1] = 16;
#endif
#if AMREX_SPACEDIM == 3
           hydro_tile_size[2] = 16;
#endif
       }

       // Run through boxes on this level and see if any of them are
       // bigger than the biggest box from our previous tuning. If so,
       // we need to re-compute the tile size.
       for (MFIter mfi(S_new, false); mfi.isValid(); ++mfi) {
           if (mfi.validbox().numPts() > largest_box_from_hydro_tile_size_tuning) {
               hydro_tile_size_has_been_tuned = 0;
           }

           // Also, sum up the number of bytes in the state data.
           mf_size += S_new[mfi].nBytes();
       }
   }
#endif

#ifdef _OPENMP
#ifdef RADIATION
#pragma omp parallel reduction(max:nstep_fsp)
#else
#pragma omp parallel
#endif
#endif
  {

#ifdef RADIATION
    int priv_nstep_fsp = -1;
#endif

    // Declare local storage now. This should be done outside the
    // MFIter loop, and then we will resize the Fabs in each MFIter
    // loop iteration. We use the async arenato ensure that their
    // memory is saved until it is no longer needed (only relevant for
    // the asynchronous case, usually on GPUs).

    FArrayBox shk(The_Async_Arena());
    FArrayBox q(The_Async_Arena()), qaux(The_Async_Arena());
    FArrayBox rho_inv(The_Async_Arena());
    FArrayBox src_q(The_Async_Arena());
    FArrayBox qxm(The_Async_Arena()), qxp(The_Async_Arena());
#if AMREX_SPACEDIM >= 2
    FArrayBox qym(The_Async_Arena()), qyp(The_Async_Arena());
#endif
#if AMREX_SPACEDIM == 3
    FArrayBox qzm(The_Async_Arena()), qzp(The_Async_Arena());
#endif
    FArrayBox div(The_Async_Arena());
#if AMREX_SPACEDIM >= 2
    FArrayBox ftmp1(The_Async_Arena()), ftmp2(The_Async_Arena());
#ifdef RADIATION
    FArrayBox rftmp1(The_Async_Arena()), rftmp2(The_Async_Arena());
#endif
    FArrayBox qgdnvtmp1(The_Async_Arena()), qgdnvtmp2(The_Async_Arena());
    FArrayBox ql(The_Async_Arena()), qr(The_Async_Arena());
#endif
    Vector<FArrayBox> flux, qe;
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        flux.push_back(FArrayBox(The_Async_Arena()));
        qe.push_back(FArrayBox(The_Async_Arena()));
    }

#ifdef RADIATION
    Vector<FArrayBox> rad_flux;
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        rad_flux.push_back(FArrayBox(The_Async_Arena()));
    }
#endif
#if AMREX_SPACEDIM <= 2
    FArrayBox pradial(The_Async_Arena());
#endif
#if AMREX_SPACEDIM == 3
    FArrayBox qmyx(The_Async_Arena()), qpyx(The_Async_Arena());
    FArrayBox qmzx(The_Async_Arena()), qpzx(The_Async_Arena());
    FArrayBox qmxy(The_Async_Arena()), qpxy(The_Async_Arena());
    FArrayBox qmzy(The_Async_Arena()), qpzy(The_Async_Arena());
    FArrayBox qmxz(The_Async_Arena()), qpxz(The_Async_Arena());
    FArrayBox qmyz(The_Async_Arena()), qpyz(The_Async_Arena());
#endif

    MultiFab& old_source = get_old_data(Source_Type);

    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      // the valid region box
      const Box& bx = mfi.tilebox();

      const Box& obx = amrex::grow(bx, 1);

      // Compute the primitive variables (both q and qaux) from
      // the conserved variables.

      const Box& qbx = amrex::grow(bx, NUM_GROW);
      const Box& qbx3 = amrex::grow(bx, 3);

#ifdef RADIATION
      q.resize(qbx, NQ);
#else
      // note: we won't store the passives in q, so we'll compute their
      // primitive versions on demand as needed
      q.resize(qbx, NQTHERM);
#endif
      fab_size += q.nBytes();
      Array4<Real> const q_arr = q.array();

      qaux.resize(qbx, NQAUX);
      fab_size += qaux.nBytes();
      Array4<Real> const qaux_arr = qaux.array();

      Array4<Real const> const U_old_arr = Sborder.array(mfi);

      rho_inv.resize(qbx3, 1);
      fab_size += rho_inv.nBytes();
      Array4<Real> const rho_inv_arr = rho_inv.array();

      amrex::ParallelFor(qbx3,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
          rho_inv_arr(i,j,k) = 1.0 / U_old_arr(i,j,k,URHO);
      });

      ctoprim(qbx, time, U_old_arr,
#ifdef RADIATION
              Erborder.array(mfi), lamborder.array(mfi),
#endif
              q_arr, qaux_arr);

#if AMREX_SPACEDIM == 2
      Array4<Real const> const areax_arr = area[0].array(mfi);
      Array4<Real const> const areay_arr = area[1].array(mfi);
      Array4<Real> const vol_arr = volume.array(mfi);
#endif

#if AMREX_SPACEDIM < 3
      Array4<Real const> const dLogAreaX_arr = (dLogArea[0]).array(mfi);
#endif
#if AMREX_SPACEDIM == 2
      Array4<Real const> const dLogAreaY_arr = (dLogArea[1]).array(mfi);
#endif

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
      fab_size += shk.nBytes();

      // Multidimensional shock detection

      // this is a local shock variable used only for the Riemann
      // solver -- this will never be used to update the value in the
      // conserved state.

      Array4<Real> const shk_arr = shk.array();

      // get the primitive variable hydro sources

      src_q.resize(qbx3, NQSRC);
      fab_size += src_q.nBytes();
      Array4<Real> const src_q_arr = src_q.array();

      Array4<Real> const old_src_arr = old_source.array(mfi);
      Array4<Real> const src_corr_arr = source_corrector.array(mfi);

      amrex::ParallelFor(qbx3,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
          hydro::src_to_prim(i, j, k, dt, U_old_arr, q_arr, old_src_arr, src_corr_arr, src_q_arr);
      });

      if (hybrid_riemann == 1) {
        shock(obx, q_arr, old_src_arr, shk_arr);
      }
      else {
        amrex::ParallelFor(obx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
          shk_arr(i,j,k) = 0.0;
        });
      }

      // work on the interface states

      qxm.resize(obx, NQ);
      fab_size += qxm.nBytes();

      qxp.resize(obx, NQ);
      fab_size += qxp.nBytes();

      Array4<Real> const qxm_arr = qxm.array();
      Array4<Real> const qxp_arr = qxp.array();

#if AMREX_SPACEDIM >= 2
      qym.resize(obx, NQ);
      fab_size += qym.nBytes();

      qyp.resize(obx, NQ);
      fab_size += qyp.nBytes();

      Array4<Real> const qym_arr = qym.array();
      Array4<Real> const qyp_arr = qyp.array();

#endif

#if AMREX_SPACEDIM == 3
      qzm.resize(obx, NQ);
      fab_size += qzm.nBytes();

      qzp.resize(obx, NQ);
      fab_size += qzp.nBytes();

      Array4<Real> const qzm_arr = qzm.array();
      Array4<Real> const qzp_arr = qzp.array();

#endif

      if (ppm_type == 0) {

        ctu_plm_states(obx, bx,
                       U_old_arr, rho_inv_arr,
                       q_arr,
                       qaux_arr,
                       src_q_arr,
                       qxm_arr, qxp_arr,
#if AMREX_SPACEDIM >= 2
                       qym_arr, qyp_arr,
#endif
#if AMREX_SPACEDIM == 3
                       qzm_arr, qzp_arr,
#endif
#if AMREX_SPACEDIM < 3
                       dLogAreaX_arr,
#endif
#if AMREX_SPACEDIM == 2
                       dLogAreaY_arr,
#endif
                       dt);

      } else {

#ifdef RADIATION
        ctu_ppm_rad_states(obx, bx,
                           U_old_arr, rho_inv_arr,
                           q_arr, qaux_arr, src_q_arr,
                           qxm_arr, qxp_arr,
#if AMREX_SPACEDIM >= 2
                           qym_arr, qyp_arr,
#endif
#if AMREX_SPACEDIM == 3
                           qzm_arr, qzp_arr,
#endif
#if AMREX_SPACEDIM < 3
                           dLogAreaX_arr,
#endif
#if AMREX_SPACEDIM == 2
                           dLogAreaY_arr,
#endif
                           dt);
#else

        ctu_ppm_states(obx, bx,
                       U_old_arr, rho_inv_arr,
                       q_arr, qaux_arr, src_q_arr,
                       qxm_arr, qxp_arr,
#if AMREX_SPACEDIM >= 2
                       qym_arr, qyp_arr,
#endif
#if AMREX_SPACEDIM == 3
                       qzm_arr, qzp_arr,
#endif
#if AMREX_SPACEDIM < 3
                           dLogAreaX_arr,
#endif
#if AMREX_SPACEDIM == 2
                           dLogAreaY_arr,
#endif
                       dt);
#endif

      }

      div.resize(obx, 1);
      fab_size += div.nBytes();
      auto div_arr = div.array();

      // compute divu -- we'll use this later when doing the artificial viscosity
      divu(obx, q_arr, div_arr);

      flux[0].resize(gxbx, NUM_STATE);
      fab_size += flux[0].nBytes();
      Array4<Real> const flux0_arr = (flux[0]).array();

      qe[0].resize(gxbx, NGDNV);
      auto qex_arr = qe[0].array();
      fab_size += qe[0].nBytes();

#ifdef RADIATION
      rad_flux[0].resize(gxbx, Radiation::nGroups);
      fab_size += rad_flux[0].nBytes();
      auto rad_flux0_arr = (rad_flux[0]).array();
#endif

#if AMREX_SPACEDIM >= 2
      flux[1].resize(gybx, NUM_STATE);
      fab_size += flux[1].nBytes();
      Array4<Real> const flux1_arr = (flux[1]).array();

      qe[1].resize(gybx, NGDNV);
      auto qey_arr = qe[1].array();
      fab_size += qe[1].nBytes();

#ifdef RADIATION
      rad_flux[1].resize(gybx, Radiation::nGroups);
      fab_size += rad_flux[1].nBytes();
      auto const rad_flux1_arr = (rad_flux[1]).array();
#endif
#endif

#if AMREX_SPACEDIM == 3
      flux[2].resize(gzbx, NUM_STATE);
      fab_size += flux[2].nBytes();
      Array4<Real> const flux2_arr = (flux[2]).array();

      qe[2].resize(gzbx, NGDNV);
      auto qez_arr = qe[2].array();
      fab_size += qe[2].nBytes();

#ifdef RADIATION
      rad_flux[2].resize(gzbx, Radiation::nGroups);
      fab_size += rad_flux[2].nBytes();
      auto const rad_flux2_arr = (rad_flux[2]).array();
#endif
#endif

#if AMREX_SPACEDIM <= 2
      if (!Geom().IsCartesian()) {
          pradial.resize(xbx, 1);
      }
      fab_size += pradial.nBytes();
#endif

#ifdef SIMPLIFIED_SDC
#ifdef REACTIONS
      Array4<Real> const sdc_src_arr = SDC_react_source.array(mfi);
#endif
#endif

#if AMREX_SPACEDIM == 1

#ifdef SIMPLIFIED_SDC
#ifdef REACTIONS
      add_sdc_source_to_states(xbx, 0, dt,
                               qxm_arr, qxp_arr, sdc_src_arr);
#endif
#endif

      // compute the fluxes through the x-interface

      cmpflx_plus_godunov(xbx,
                          qxm_arr, qxp_arr,
                          flux0_arr,
#ifdef RADIATION
                          rad_flux0_arr,
#endif
                          qex_arr,
                          qaux_arr,
                          shk_arr,
                          0, false);

#endif // 1-d



#if AMREX_SPACEDIM >= 2
      ftmp1.resize(obx, NUM_STATE);
      auto ftmp1_arr = ftmp1.array();
      fab_size += ftmp1.nBytes();

      ftmp2.resize(obx, NUM_STATE);
      auto ftmp2_arr = ftmp2.array();
      fab_size += ftmp2.nBytes();

#ifdef RADIATION
      rftmp1.resize(obx, Radiation::nGroups);
      auto rftmp1_arr = rftmp1.array();
      fab_size += rftmp1.nBytes();

      rftmp2.resize(obx, Radiation::nGroups);
      auto rftmp2_arr = rftmp2.array();
      fab_size += rftmp2.nBytes();
#endif

      qgdnvtmp1.resize(obx, NGDNV);
      auto qgdnvtmp1_arr = qgdnvtmp1.array();
      fab_size += qgdnvtmp1.nBytes();

#if AMREX_SPACEDIM == 3
      qgdnvtmp2.resize(obx, NGDNV);
      auto qgdnvtmp2_arr = qgdnvtmp2.array();
      fab_size += qgdnvtmp2.nBytes();
#endif

      ql.resize(obx, NQ);
      auto ql_arr = ql.array();
      fab_size += ql.nBytes();

      qr.resize(obx, NQ);
      auto qr_arr = qr.array();
      fab_size += qr.nBytes();
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
      cmpflx_plus_godunov(cxbx,
                          qxm_arr, qxp_arr,
                          ftmp1_arr,
#ifdef RADIATION
                          rftmp1_arr,
#endif
                          qgdnvtmp1_arr,
                          qaux_arr, shk_arr,
                          0, false);

      // compute F^y
      // [lo(1)-1, lo(2), 0], [hi(1)+1, hi(2)+1, 0]
      const Box& cybx = amrex::grow(ybx, IntVect(AMREX_D_DECL(1,0,0)));

      // ftmp2 = fy
      // rftmp2 = rfy
      cmpflx_plus_godunov(cybx,
                          qym_arr, qyp_arr,
                          ftmp2_arr,
#ifdef RADIATION
                          rftmp2_arr,
#endif
                          qey_arr,
                          qaux_arr, shk_arr,
                          1, false);

      // add the transverse flux difference in y to the x states
      // [lo(1), lo(2), 0], [hi(1)+1, hi(2), 0]

      // ftmp2 = fy
      // rftmp2 = rfy
      trans_single(xbx, 1, 0,
                   qxm_arr, ql_arr,
                   qxp_arr, qr_arr,
                   qaux_arr,
                   ftmp2_arr,
#ifdef RADIATION
                   rftmp2_arr,
#endif
                   qey_arr,
                   areay_arr,
                   vol_arr,
                   hdt, hdtdy);

      reset_edge_state_thermo(xbx, ql_arr);

      reset_edge_state_thermo(xbx, qr_arr);

      // solve the final Riemann problem axross the x-interfaces

#ifdef SIMPLIFIED_SDC
#ifdef REACTIONS
      add_sdc_source_to_states(xbx, 0, dt,
                               ql_arr, qr_arr, sdc_src_arr);
#endif

#endif

      cmpflx_plus_godunov(xbx,
                          ql_arr, qr_arr,
                          flux0_arr,
#ifdef RADIATION
                          rad_flux0_arr,
#endif
                          qex_arr,
                          qaux_arr, shk_arr,
                          0, false);

      // add the transverse flux difference in x to the y states
      // [lo(1), lo(2), 0], [hi(1), hi(2)+1, 0]

      // ftmp1 = fx
      // rftmp1 = rfx
      // qgdnvtmp1 = qgdnvx

      trans_single(ybx, 0, 1,
                   qym_arr, ql_arr,
                   qyp_arr, qr_arr,
                   qaux_arr,
                   ftmp1_arr,
#ifdef RADIATION
                   rftmp1_arr,
#endif
                   qgdnvtmp1_arr,
                   areax_arr,
                   vol_arr,
                   hdt, hdtdx);

      reset_edge_state_thermo(ybx, ql_arr);

      reset_edge_state_thermo(ybx, qr_arr);


      // solve the final Riemann problem axross the y-interfaces

#ifdef SIMPLIFIED_SDC
#ifdef REACTIONS
      add_sdc_source_to_states(ybx, 1, dt,
                               ql_arr, qr_arr, sdc_src_arr);
#endif
#endif

      cmpflx_plus_godunov(ybx,
                          ql_arr, qr_arr,
                          flux1_arr,
#ifdef RADIATION
                          rad_flux1_arr,
#endif
                          qey_arr,
                          qaux_arr, shk_arr,
                          1, false);
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
      cmpflx_plus_godunov(cxbx,
                          qxm_arr, qxp_arr,
                          ftmp1_arr,
#ifdef RADIATION
                          rftmp1_arr,
#endif
                          qgdnvtmp1_arr,
                          qaux_arr, shk_arr,
                          0, false);


      // [lo(1), lo(2), lo(3)-1], [hi(1), hi(2)+1, hi(3)+1]
      const Box& tyxbx = amrex::grow(ybx, IntVect(AMREX_D_DECL(0,0,1)));

      qmyx.resize(tyxbx, NQ);
      auto qmyx_arr = qmyx.array();
      fab_size += qmyx.nBytes();

      qpyx.resize(tyxbx, NQ);
      auto qpyx_arr = qpyx.array();
      fab_size += qpyx.nBytes();

      // ftmp1 = fx
      // rftmp1 = rfx
      // qgdnvtmp1 = qgdnvx
      trans_single(tyxbx, 0, 1,
                   qym_arr, qmyx_arr,
                   qyp_arr, qpyx_arr,
                   qaux_arr,
                   ftmp1_arr,
#ifdef RADIATION
                   rftmp1_arr,
#endif
                   qgdnvtmp1_arr,
                   hdt, cdtdx);

      reset_edge_state_thermo(tyxbx, qmyx.array());

      reset_edge_state_thermo(tyxbx, qpyx.array());

      // [lo(1), lo(2)-1, lo(3)], [hi(1), hi(2)+1, hi(3)+1]
      const Box& tzxbx = amrex::grow(zbx, IntVect(AMREX_D_DECL(0,1,0)));

      qmzx.resize(tzxbx, NQ);
      auto qmzx_arr = qmzx.array();
      fab_size += qmzx.nBytes();

      qpzx.resize(tzxbx, NQ);
      auto qpzx_arr = qpzx.array();
      fab_size += qpzx.nBytes();

      trans_single(tzxbx, 0, 2,
                   qzm_arr, qmzx_arr,
                   qzp_arr, qpzx_arr,
                   qaux_arr,
                   ftmp1_arr,
#ifdef RADIATION
                   rftmp1_arr,
#endif
                   qgdnvtmp1_arr,
                   hdt, cdtdx);

      reset_edge_state_thermo(tzxbx, qmzx.array());

      reset_edge_state_thermo(tzxbx, qpzx.array());

      // compute F^y
      // [lo(1)-1, lo(2), lo(3)-1], [hi(1)+1, hi(2)+1, hi(3)+1]
      const Box& cybx = amrex::grow(ybx, IntVect(AMREX_D_DECL(1,0,1)));

      // ftmp1 = fy
      // rftmp1 = rfy
      // qgdnvtmp1 = qgdnvy
      cmpflx_plus_godunov(cybx,
                          qym_arr, qyp_arr,
                          ftmp1_arr,
#ifdef RADIATION
                          rftmp1_arr,
#endif
                          qgdnvtmp1_arr,
                          qaux_arr, shk_arr,
                          1, false);

      // [lo(1), lo(2), lo(3)-1], [hi(1)+1, hi(2), lo(3)+1]
      const Box& txybx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,0,1)));

      qmxy.resize(txybx, NQ);
      auto qmxy_arr = qmxy.array();
      fab_size += qmxy.nBytes();

      qpxy.resize(txybx, NQ);
      auto qpxy_arr = qpxy.array();
      fab_size += qpxy.nBytes();

      // ftmp1 = fy
      // rftmp1 = rfy
      // qgdnvtmp1 = qgdnvy
      trans_single(txybx, 1, 0,
                   qxm_arr, qmxy_arr,
                   qxp_arr, qpxy_arr,
                   qaux_arr,
                   ftmp1_arr,
#ifdef RADIATION
                   rftmp1_arr,
#endif
                   qgdnvtmp1_arr,
                   hdt, cdtdy);

      reset_edge_state_thermo(txybx, qmxy.array());

      reset_edge_state_thermo(txybx, qpxy.array());

      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2), lo(3)+1]
      const Box& tzybx = amrex::grow(zbx, IntVect(AMREX_D_DECL(1,0,0)));

      qmzy.resize(tzybx, NQ);
      auto qmzy_arr = qmzy.array();
      fab_size += qmzy.nBytes();

      qpzy.resize(tzybx, NQ);
      auto qpzy_arr = qpzy.array();
      fab_size += qpzy.nBytes();

      // ftmp1 = fy
      // rftmp1 = rfy
      // qgdnvtmp1 = qgdnvy
      trans_single(tzybx, 1, 2,
                   qzm_arr, qmzy_arr,
                   qzp_arr, qpzy_arr,
                   qaux_arr,
                   ftmp1_arr,
#ifdef RADIATION
                   rftmp1_arr,
#endif
                   qgdnvtmp1_arr,
                   hdt, cdtdy);

      reset_edge_state_thermo(tzybx, qmzy.array());

      reset_edge_state_thermo(tzybx, qpzy.array());

      // compute F^z
      // [lo(1)-1, lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, hi(3)+1]
      const Box& czbx = amrex::grow(zbx, IntVect(AMREX_D_DECL(1,1,0)));

      // ftmp1 = fz
      // rftmp1 = rfz
      // qgdnvtmp1 = qgdnvz
      cmpflx_plus_godunov(czbx,
                          qzm_arr, qzp_arr,
                          ftmp1_arr,
#ifdef RADIATION
                          rftmp1_arr,
#endif
                          qgdnvtmp1_arr,
                          qaux_arr, shk_arr,
                          2, false);

      // [lo(1)-1, lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, lo(3)]
      const Box& txzbx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,1,0)));

      qmxz.resize(txzbx, NQ);
      auto qmxz_arr = qmxz.array();
      fab_size += qmxz.nBytes();

      qpxz.resize(txzbx, NQ);
      auto qpxz_arr = qpxz.array();
      fab_size += qpxz.nBytes();

      // ftmp1 = fz
      // rftmp1 = rfz
      // qgdnvtmp1 = qgdnvz
      trans_single(txzbx, 2, 0,
                   qxm_arr, qmxz_arr,
                   qxp_arr, qpxz_arr,
                   qaux_arr,
                   ftmp1_arr,
#ifdef RADIATION
                   rftmp1_arr,
#endif
                   qgdnvtmp1_arr,
                   hdt, cdtdz);

      reset_edge_state_thermo(txzbx, qmxz.array());

      reset_edge_state_thermo(txzbx, qpxz.array());

      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2)+1, lo(3)]
      const Box& tyzbx = amrex::grow(ybx, IntVect(AMREX_D_DECL(1,0,0)));

      qmyz.resize(tyzbx, NQ);
      auto qmyz_arr = qmyz.array();
      fab_size += qmyz.nBytes();

      qpyz.resize(tyzbx, NQ);
      auto qpyz_arr = qpyz.array();
      fab_size += qpyz.nBytes();

      // ftmp1 = fz
      // rftmp1 = rfz
      // qgdnvtmp1 = qgdnvz
      trans_single(tyzbx, 2, 1,
                   qym_arr, qmyz_arr,
                   qyp_arr, qpyz_arr,
                   qaux_arr,
                   ftmp1_arr,
#ifdef RADIATION
                   rftmp1_arr,
#endif
                   qgdnvtmp1_arr,
                   hdt, cdtdz);

      reset_edge_state_thermo(tyzbx, qmyz.array());

      reset_edge_state_thermo(tyzbx, qpyz.array());

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
      cmpflx_plus_godunov(cyzbx,
                          qmyz_arr, qpyz_arr,
                          ftmp1_arr,
#ifdef RADIATION
                          rftmp1_arr,
#endif
                          qgdnvtmp1_arr,
                          qaux_arr, shk_arr,
                          1, false);

      // compute F^{z|y}
      // [lo(1)-1, lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)+1]
      const Box& czybx = amrex::grow(zbx, IntVect(AMREX_D_DECL(1,0,0)));

      // ftmp2 = fzy
      // rftmp2 = rfzy
      // qgdnvtmp2 = qgdnvzy
      cmpflx_plus_godunov(czybx,
                          qmzy_arr, qpzy_arr,
                          ftmp2_arr,
#ifdef RADIATION
                          rftmp2_arr,
#endif
                          qgdnvtmp2_arr,
                          qaux_arr, shk_arr,
                          2, false);

      // compute the corrected x interface states and fluxes
      // [lo(1), lo(2), lo(3)], [hi(1)+1, hi(2), hi(3)]

      trans_final(xbx, 0, 1, 2,
                  qxm_arr, ql_arr,
                  qxp_arr, qr_arr,
                  qaux_arr,
                  ftmp1_arr,
#ifdef RADIATION
                  rftmp1_arr,
#endif
                  ftmp2_arr,
#ifdef RADIATION
                  rftmp2_arr,
#endif
                  qgdnvtmp1_arr,
                  qgdnvtmp2_arr,
                  hdtdy, hdtdz);

      reset_edge_state_thermo(xbx, ql_arr);

      reset_edge_state_thermo(xbx, qr_arr);

#ifdef SIMPLIFIED_SDC
#ifdef REACTIONS
      add_sdc_source_to_states(xbx, 0, dt,
                               ql_arr, qr_arr, sdc_src_arr);
#endif
#endif


      cmpflx_plus_godunov(xbx,
                          ql_arr, qr_arr,
                          flux0_arr,
#ifdef RADIATION
                          rad_flux0_arr,
#endif
                          qex_arr,
                          qaux_arr, shk_arr,
                          0, false);

      //
      // Use qy?, q?zx, q?xz to compute final y-flux
      //

      // compute F^{z|x}
      // [lo(1), lo(2)-1, lo(3)], [hi(1), hi(2)+1, hi(3)+1]
      const Box& czxbx = amrex::grow(zbx, IntVect(AMREX_D_DECL(0,1,0)));

      // ftmp1 = fzx
      // rftmp1 = rfzx
      // qgdnvtmp1 = qgdnvzx
      cmpflx_plus_godunov(czxbx,
                          qmzx_arr, qpzx_arr,
                          ftmp1_arr,
#ifdef RADIATION
                          rftmp1_arr,
#endif
                          qgdnvtmp1_arr,
                          qaux_arr, shk_arr,
                          2, false);

      // compute F^{x|z}
      // [lo(1), lo(2)-1, lo(3)], [hi(1)+1, hi(2)+1, hi(3)]
      const Box& cxzbx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,1,0)));

      // ftmp2 = fxz
      // rftmp2 = rfxz
      // qgdnvtmp2 = qgdnvxz
      cmpflx_plus_godunov(cxzbx,
                          qmxz_arr, qpxz_arr,
                          ftmp2_arr,
#ifdef RADIATION
                          rftmp2_arr,
#endif
                          qgdnvtmp2_arr,
                          qaux_arr, shk_arr,
                          0, false);

      // Compute the corrected y interface states and fluxes
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)]

      trans_final(ybx, 1, 0, 2,
                  qym_arr, ql_arr,
                  qyp_arr, qr_arr,
                  qaux_arr,
                  ftmp2_arr,
#ifdef RADIATION
                  rftmp2_arr,
#endif
                  ftmp1_arr,
#ifdef RADIATION
                  rftmp1_arr,
#endif
                  qgdnvtmp2_arr,
                  qgdnvtmp1_arr,
                  hdtdx, hdtdz);

      reset_edge_state_thermo(ybx, ql_arr);

      reset_edge_state_thermo(ybx, qr_arr);

#ifdef SIMPLIFIED_SDC
#ifdef REACTIONS
      add_sdc_source_to_states(ybx, 1, dt,
                               ql_arr, qr_arr, sdc_src_arr);
#endif
#endif


      // Compute the final F^y
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2)+1, hi(3)]
      cmpflx_plus_godunov(ybx,
                          ql_arr, qr_arr,
                          flux1_arr,
#ifdef RADIATION
                          rad_flux1_arr,
#endif
                          qey_arr,
                          qaux_arr, shk_arr,
                          1, false);

      //
      // Use qz?, q?xy, q?yx to compute final z-flux
      //

      // compute F^{x|y}
      // [lo(1), lo(2), lo(3)-1], [hi(1)+1, hi(2), hi(3)+1]
      const Box& cxybx = amrex::grow(xbx, IntVect(AMREX_D_DECL(0,0,1)));

      // ftmp1 = fxy
      // rftmp1 = rfxy
      // qgdnvtmp1 = qgdnvxy
      cmpflx_plus_godunov(cxybx,
                          qmxy_arr, qpxy_arr,
                          ftmp1_arr,
#ifdef RADIATION
                          rftmp1_arr,
#endif
                          qgdnvtmp1_arr,
                          qaux_arr, shk_arr,
                          0, false);

      // compute F^{y|x}
      // [lo(1), lo(2), lo(3)-1], [hi(1), hi(2)+dg(2), hi(3)+1]
      const Box& cyxbx = amrex::grow(ybx, IntVect(AMREX_D_DECL(0,0,1)));

      // ftmp2 = fyx
      // rftmp2 = rfyx
      // qgdnvtmp2 = qgdnvyx
      cmpflx_plus_godunov(cyxbx,
                          qmyx_arr, qpyx_arr,
                          ftmp2_arr,
#ifdef RADIATION
                          rftmp2_arr,
#endif
                          qgdnvtmp2_arr,
                          qaux_arr, shk_arr,
                          1, false);

      // compute the corrected z interface states and fluxes
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1]

      trans_final(zbx, 2, 0, 1,
                  qzm_arr, ql_arr,
                  qzp_arr, qr_arr,
                  qaux_arr,
                  ftmp1_arr,
#ifdef RADIATION
                  rftmp1_arr,
#endif
                  ftmp2_arr,
#ifdef RADIATION
                  rftmp2_arr,
#endif
                  qgdnvtmp1_arr,
                  qgdnvtmp2_arr,
                  hdtdx, hdtdy);

      reset_edge_state_thermo(zbx, ql_arr);

      reset_edge_state_thermo(zbx, qr_arr);

#ifdef SIMPLIFIED_SDC
#ifdef REACTIONS
      add_sdc_source_to_states(zbx, 2, dt,
                               ql_arr, qr_arr, sdc_src_arr);
#endif
#endif

      // compute the final z fluxes F^z
      // [lo(1), lo(2), lo(3)], [hi(1), hi(2), hi(3)+1]

      cmpflx_plus_godunov(zbx,
                          ql_arr, qr_arr,
                          flux2_arr,
#ifdef RADIATION
                          rad_flux2_arr,
#endif
                          qez_arr,
                          qaux_arr, shk_arr,
                          2, false);

#endif // 3-d



      // clean the fluxes

      for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

          const Box& nbx = amrex::surroundingNodes(bx, idir);

          Array4<Real> const flux_arr = (flux[idir]).array();

          // Zero out shock and temp fluxes -- these are physically meaningless here
          amrex::ParallelFor(nbx,
          [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
          {
              flux_arr(i,j,k,UTEMP) = 0.e0;
#ifdef SHOCK_VAR
              flux_arr(i,j,k,USHK) = 0.e0;
#endif
#ifdef NSE_NET
              flux_arr(i,j,k,UMUP) = 0.e0;
              flux_arr(i,j,k,UMUN) = 0.e0;
#endif
          });

          apply_av(nbx, idir, div_arr, U_old_arr, flux_arr);

#ifdef RADIATION
          Array4<Real> const rad_flux_arr = (rad_flux[idir]).array();
          Array4<Real const> const Erin_arr = Erborder.array(mfi);

          apply_av_rad(nbx, idir, div_arr, Erin_arr, rad_flux_arr);
#endif

          if (limit_fluxes_on_small_dens == 1) {
              limit_hydro_fluxes_on_small_dens
                  (nbx, idir,
                   Sborder.array(mfi),
                   volume.array(mfi),
                   flux[idir].array(),
                   area[idir].array(mfi),
                   dt);
          }

          normalize_species_fluxes(nbx, flux_arr);

      }



      // conservative update
      Array4<Real> const update_arr = S_new.array(mfi);

      Array4<Real> const flx_arr = (flux[0]).array();
      Array4<Real> const qx_arr = (qe[0]).array();

#if AMREX_SPACEDIM >= 2
      Array4<Real> const fly_arr = (flux[1]).array();
      Array4<Real> const qy_arr = (qe[1]).array();
#endif

#if AMREX_SPACEDIM == 3
      Array4<Real> const flz_arr = (flux[2]).array();
      Array4<Real> const qz_arr = (qe[2]).array();
#endif

      consup_hydro(bx,
                   update_arr,
                   flx_arr, qx_arr,
#if AMREX_SPACEDIM >= 2
                   fly_arr, qy_arr,
#endif
#if AMREX_SPACEDIM == 3
                   flz_arr, qz_arr,
#endif
                   dt);


#ifdef HYBRID_MOMENTUM
      auto dx_arr = geom.CellSizeArray();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

          GpuArray<Real, 3> loc;

          position(i, j, k, geomdata, loc);

          for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
              loc[dir] -= problem::center[dir];

          Real R = amrex::max(std::sqrt(loc[0] * loc[0] + loc[1] * loc[1]),
                              std::numeric_limits<Real>::min());
          Real RInv = 1.0_rt / R;

          update_arr(i,j,k,UMR) -= dt * (loc[0] * RInv) * (qx_arr(i+1,j,k,GDPRES) - qx_arr(i,j,k,GDPRES)) / dx_arr[0];
#if AMREX_SPACEDIM >= 2
          update_arr(i,j,k,UMR) -= dt * (loc[1] * RInv) * (qy_arr(i,j+1,k,GDPRES) - qy_arr(i,j,k,GDPRES)) / dx_arr[1];
#endif
      });
#endif

#ifdef RADIATION
      ctu_rad_consup(bx,
                     update_arr,
                     Erborder.array(mfi),
                     Er_new.array(mfi),
                     (rad_flux[0]).array(),
                     (qe[0]).array(),
                     (area[0]).array(mfi),
#if AMREX_SPACEDIM >= 2
                     (rad_flux[1]).array(),
                     (qe[1]).array(),
                     (area[1]).array(mfi),
#endif
#if AMREX_SPACEDIM == 3
                     (rad_flux[2]).array(),
                     (qe[2]).array(),
                     (area[2]).array(mfi),
#endif
                     priv_nstep_fsp,
                     volume.array(mfi),
                     dt);

      nstep_fsp = std::max(nstep_fsp, priv_nstep_fsp);
#endif


      for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

        const Box& nbx = amrex::surroundingNodes(bx, idir);

        Array4<Real> const flux_arr = (flux[idir]).array();
        Array4<Real const> const area_arr = (area[idir]).array(mfi);

        scale_flux(nbx,
#if AMREX_SPACEDIM == 1
                   qex_arr,
#endif
                   flux_arr, area_arr, dt);

#ifdef RADIATION
        Array4<Real> const rad_flux_arr = (rad_flux[idir]).array();
        scale_rad_flux(nbx, rad_flux_arr, area_arr, dt);
#endif

        if (idir == 0) {
#if AMREX_SPACEDIM <= 2
            Array4<Real> pradial_fab = pradial.array();
#endif

            // get the scaled radial pressure -- we need to treat this specially
#if AMREX_SPACEDIM <= 2
            if (!mom_flux_has_p(0, 0, coord)) {
                amrex::ParallelFor(nbx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    pradial_fab(i,j,k) = qex_arr(i,j,k,GDPRES) * dt;
                });
            }

#endif
        }

        // Store the fluxes from this advance. For simplified SDC integration we
        // only need to do this on the last iteration.

        bool add_fluxes = true;

        if (time_integration_method == SimplifiedSpectralDeferredCorrections &&
            sdc_iteration != sdc_iters - 1) {
            add_fluxes = false;
        }

        if (add_fluxes) {

            Array4<Real> const flux_fab = (flux[idir]).array();
            Array4<Real> fluxes_fab = (*fluxes[idir]).array(mfi);

            amrex::ParallelFor(mfi.nodaltilebox(idir), NUM_STATE,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                fluxes_fab(i,j,k,n) += flux_fab(i,j,k,n);
            });

#ifdef RADIATION
            Array4<Real> const rad_flux_fab = (rad_flux[idir]).array();
            Array4<Real> rad_fluxes_fab = (*rad_fluxes[idir]).array(mfi);

            amrex::ParallelFor(mfi.nodaltilebox(idir), Radiation::nGroups,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                rad_fluxes_fab(i,j,k,n) += rad_flux_fab(i,j,k,n);
            });
#endif

#if AMREX_SPACEDIM <= 2
            if (idir == 0 && !mom_flux_has_p(0, 0, coord)) {
                Array4<Real> pradial_fab = pradial.array();
                Array4<Real> P_radial_fab = P_radial.array(mfi);

                amrex::ParallelFor(mfi.nodaltilebox(0),
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    P_radial_fab(i,j,k,0) += pradial_fab(i,j,k,0);
                });
            }

#endif

        } // add_fluxes

        Array4<Real> const flux_fab = (flux[idir]).array();
        Array4<Real> mass_fluxes_fab = (*mass_fluxes[idir]).array(mfi);

        amrex::ParallelFor(mfi.nodaltilebox(idir),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // This is a copy, not an add, since we need mass_fluxes to be
            // only this subcycle's data when we evaluate the gravitational
            // forces.

            mass_fluxes_fab(i,j,k,0) = flux_fab(i,j,k,URHO);
        });

      } // idir loop

#ifdef AMREX_USE_GPU
      if (castro::hydro_memory_footprint_ratio > 0.0) {

          if (hydro_tile_size_has_been_tuned == 0) {

              // Keep a running record of the largest box we've encountered.

              largest_box_from_hydro_tile_size_tuning = amrex::max(largest_box_from_hydro_tile_size_tuning,
                                                                   mfi.validbox().numPts());

              // If we're tuning the hydro tile size during this timestep, we will record
              // the total amount of additional memory allocated, relative to the size of S_new.
              // Then we will reset the tile size so that it is no larger than the requested
              // memory footprint.

              // This could be generalized in the future to operate with more granularity
              // than the MFIter loop boundary. We could have potential synchronization
              // points prior to each of the above kernel launches.

              maximum_tile_size[0] = amrex::max(maximum_tile_size[0], bx.bigEnd(0) - mfi.validbox().smallEnd(0) + 1);
#if AMREX_SPACEDIM >= 2
              maximum_tile_size[1] = amrex::max(maximum_tile_size[1], bx.bigEnd(1) - mfi.validbox().smallEnd(1) + 1);
#endif
#if AMREX_SPACEDIM ==3
              maximum_tile_size[2] = amrex::max(maximum_tile_size[2], bx.bigEnd(2) - mfi.validbox().smallEnd(2) + 1);
#endif

              if (fab_size >= castro::hydro_memory_footprint_ratio * mf_size) {
                  // If we reached the memory limit, set the tile size to the current
                  // maximum tile size.
                  Gpu::synchronize();
                  hydro_tile_size = maximum_tile_size;
                  hydro_tile_size_has_been_tuned = 1;
              }
              else if (mfi.tileIndex() == mfi.length() - 1) {
                  // If we reached the last tile and we haven't gone over the
                  // memory limit, effectively disable tiling.
                  hydro_tile_size[0] = 1024;
#if AMREX_SPACEDIM >= 2
                  hydro_tile_size[1] = 1024;
#endif
#if AMREX_SPACEDIM == 3
                  hydro_tile_size[2] = 1024;
#endif
                  hydro_tile_size_has_been_tuned = 1;
              }

          }
          else {
              // If we have already tuned the parameter, then synchronize each time the
              // outstanding number of active bytes is larger than our ratio.

              if (fab_size >= castro::hydro_memory_footprint_ratio * mf_size) {
                  Gpu::synchronize();

                  // Reset the counter for the next sequence of tiles.
                  fab_size = 0;
              }
          }
      }
#endif

    } // MFIter loop

  } // OMP loop

#ifdef RADIATION
  if (radiation->verbose>=1) {
      amrex::Real llevel = level;
#ifdef BL_LAZY
    Lazy::QueueReduction( [=] () mutable {
#endif
       ParallelDescriptor::ReduceIntMax(nstep_fsp, ParallelDescriptor::IOProcessorNumber());
       if (ParallelDescriptor::IOProcessor() && nstep_fsp > 0) {
         std::cout << "Radiation f-space advection on level " << llevel
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
#endif

  // Check for small/negative densities and X > 1 or X < 0.

  status = check_for_negative_density();

  if (status.success == false) {
      return status;
  }

  // Sync up state after hydro source.

  clean_state(
#ifdef MHD
               Bx_new, By_new, Bz_new,
#endif
               S_new, time + dt, 0);

  // Check for NaN's.

  check_for_nan(S_new);

#ifdef GRAVITY
  // Must define new value of "center" after advecting on the grid

  if (moving_center == 1) {
      define_new_center(S_new, time);
  }
#endif

  // Perform reflux (for non-subcycling advances).

  if (parent->subcyclingMode() == "None") {
      if (do_reflux == 1) {
          FluxRegCrseInit();
          FluxRegFineAdd();
      }
  }

  if (verbose) {
      amrex::Print() << "... Leaving construct_ctu_hydro_source() on level " << level << std::endl << std::endl;
  }

  if (verbose > 0)
    {
      const int IOProc = ParallelDescriptor::IOProcessorNumber();
      amrex::Real run_time = ParallelDescriptor::second() - strt_time;
      amrex::Real llevel = level;

#ifdef BL_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        amrex::Print() << "Castro::construct_ctu_hydro_source() time = " << run_time
                       << " on level " << llevel << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }

#endif

  return status;
}
