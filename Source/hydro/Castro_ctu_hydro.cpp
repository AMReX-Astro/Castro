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
    std::cout << "... Entering construct_ctu_hydro_source()" << std::endl << std::endl;

  hydro_source.setVal(0.0);

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
    FArrayBox dq;
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

    size_t starting_size = MultiFab::queryMemUsage("AmrLevel_Level_" + std::to_string(level));
    size_t current_size = starting_size;

    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      size_t fab_size = 0;

      // the valid region box
      const Box& bx = mfi.tilebox();

      const Box& obx = amrex::grow(bx, 1);

      flatn.resize(obx, 1);
      Elixir elix_flatn = flatn.elixir();
      fab_size += flatn.nBytes();

#ifdef RADIATION
      flatg.resize(obx, 1);
      Elixir elix_flatg = flatg.elixir();
      fab_size += flatg.nBytes();
#endif

      // If we are oversubscribing the GPU, performance of the hydro will be constrained
      // due to its heavy memory requirements. We can help the situation by prefetching in
      // all the data we will need, and then prefetching it out at the end. This at least
      // improves performance by mitigating the number of unified memory page faults.

      // Unfortunately in CUDA there is no easy way to see actual current memory usage when
      // using unified memory; querying CUDA for free memory usage will only tell us whether
      // we've oversubscribed at any point, not whether we're currently oversubscribing, but
      // this is still a good heuristic in most cases.

      bool oversubscribed = false;

#ifdef AMREX_USE_CUDA
      if (Gpu::Device::freeMemAvailable() < 0.05 * Gpu::Device::totalGlobalMem()) {
          oversubscribed = true;
      }
#endif

      if (oversubscribed) {
          q[mfi].prefetchToDevice();
          qaux[mfi].prefetchToDevice();
          src_q[mfi].prefetchToDevice();
          volume[mfi].prefetchToDevice();
          Sborder[mfi].prefetchToDevice();
          hydro_source[mfi].prefetchToDevice();
          for (int i = 0; i < AMREX_SPACEDIM; ++i) {
              area[i][mfi].prefetchToDevice();
              (*fluxes[i])[mfi].prefetchToDevice();
          }
#if AMREX_SPACEDIM < 3
          dLogArea[0][mfi].prefetchToDevice();
          P_radial[mfi].prefetchToDevice();
#endif
#ifdef RADIATION
          Erborder[mfi].prefetchToDevice();
          Er_new[mfi].prefetchToDevice();
#endif
      }

      // compute the flattening coefficient

      Array4<Real const> const q_arr = q.array(mfi);
      Array4<Real> const flatn_arr = flatn.array();
#ifdef RADIATION
      Array4<Real> const flatg_arr = flatg.array();
#endif

      if (first_order_hydro == 1) {
        AMREX_PARALLEL_FOR_3D(obx, i, j, k, { flatn_arr(i,j,k) = 0.0; });
      } else if (use_flattening == 1) {

        uflatten(obx, q_arr, flatn_arr, QPRES);

#ifdef RADIATION
        uflatten(obx, q_arr, flatg_arr, QPTOT);

        AMREX_PARALLEL_FOR_3D(obx, i, j, k,
        {
          flatn_arr(i,j,k) = flatn_arr(i,j,k) * flatg_arr(i,j,k);

          if (Radiation::flatten_pp_threshold > 0.0) {
            if ( q_arr(i-1,j,k,QU) + q_arr(i,j-1,k,QV) + q_arr(i,j,k-1,QW) >
                 q_arr(i+1,j,k,QU) + q_arr(i,j+1,k,QV) + q_arr(i,j,k+1,QW) ) {

              if (q_arr(i,j,k,QPRES) < Radiation::flatten_pp_threshold * q_arr(i,j,k,QPTOT)) {
                flatn_arr(i,j,k) = 0.0;
              }
            }
          }
        });
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
      fab_size += shk.nBytes();

      Array4<Real> const shk_arr = shk.array();

      // Multidimensional shock detection
      // Used for the hybrid Riemann solver

#ifdef SHOCK_VAR
      bool compute_shock = true;
#else
      bool compute_shock = false;
#endif

      if (hybrid_riemann == 1 || compute_shock) {
#pragma gpu box(obx)
          ca_shock(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
                   BL_TO_FORTRAN_ANYD(q[mfi]),
                   BL_TO_FORTRAN_ANYD(shk),
                   AMREX_REAL_ANYD(dx));
      }
      else {
        AMREX_PARALLEL_FOR_3D(obx, i, j, k, { shk_arr(i,j,k) = 0.0; });
      }

      qxm.resize(obx, NQ);
      Elixir elix_qxm = qxm.elixir();
      fab_size += shk.nBytes();

      qxp.resize(obx, NQ);
      Elixir elix_qxp = qxp.elixir();
      fab_size += qxp.nBytes();

#if AMREX_SPACEDIM >= 2
      qym.resize(obx, NQ);
      Elixir elix_qym = qym.elixir();
      fab_size += qym.nBytes();

      qyp.resize(obx, NQ);
      Elixir elix_qyp = qyp.elixir();
      fab_size += qyp.nBytes();
#endif

#if AMREX_SPACEDIM == 3
      qzm.resize(obx, NQ);
      Elixir elix_qzm = qzm.elixir();
      fab_size += qzm.nBytes();

      qzp.resize(obx, NQ);
      Elixir elix_qzp = qzp.elixir();
      fab_size += qzp.nBytes();
#endif

      if (ppm_type == 0) {

        dq.resize(obx, NQ);
        Elixir elix_dq = dq.elixir();
        fab_size += dq.nBytes();

#pragma gpu box(obx)
        ctu_plm_states(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
                       AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                       BL_TO_FORTRAN_ANYD(q[mfi]),
                       BL_TO_FORTRAN_ANYD(flatn),
                       BL_TO_FORTRAN_ANYD(qaux[mfi]),
                       BL_TO_FORTRAN_ANYD(src_q[mfi]),
                       BL_TO_FORTRAN_ANYD(dq),
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

      } else {

#pragma gpu box(obx)
        ctu_ppm_states(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
                       AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                       BL_TO_FORTRAN_ANYD(q[mfi]),
                       BL_TO_FORTRAN_ANYD(flatn),
                       BL_TO_FORTRAN_ANYD(qaux[mfi]),
                       BL_TO_FORTRAN_ANYD(src_q[mfi]),
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
      }

      div.resize(obx, 1);
      Elixir elix_div = div.elixir();
      fab_size += div.nBytes();

      // compute divu -- we'll use this later when doing the artifical viscosity
#pragma gpu box(obx)
      divu(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
           BL_TO_FORTRAN_ANYD(q[mfi]),
           AMREX_REAL_ANYD(dx),
           BL_TO_FORTRAN_ANYD(div));

      q_int.resize(obx, NQ);
      Elixir elix_q_int = q_int.elixir();
      fab_size += q_int.nBytes();

#ifdef RADIATION
      lambda_int.resize(obx, Radiation::nGroups);
      Elixir elix_lambda_int = lambda_int.elixir();
      fab_size += lambda_int.nBytes();
#endif

      flux[0].resize(gxbx, NUM_STATE);
      Elixir elix_flux_x = flux[0].elixir();
      fab_size += flux[0].nBytes();

      qe[0].resize(gxbx, NGDNV);
      Elixir elix_qe_x = qe[0].elixir();
      fab_size += qe[0].nBytes();

#ifdef RADIATION
      rad_flux[0].resize(gxbx, Radiation::nGroups);
      Elixir elix_rad_flux_x = rad_flux[0].elixir();
      fab_size += rad_flux[0].nBytes();
#endif

#if AMREX_SPACEDIM >= 2
      flux[1].resize(gybx, NUM_STATE);
      Elixir elix_flux_y = flux[1].elixir();
      fab_size += flux[1].nBytes();

      qe[1].resize(gybx, NGDNV);
      Elixir elix_qe_y = qe[1].elixir();
      fab_size += qe[1].nBytes();

#ifdef RADIATION
      rad_flux[1].resize(gybx, Radiation::nGroups);
      Elixir elix_rad_flux_y = rad_flux[1].elixir();
      fab_size += rad_flux[1].nBytes();
#endif
#endif

#if AMREX_SPACEDIM == 3
      flux[2].resize(gzbx, NUM_STATE);
      Elixir elix_flux_z = flux[2].elixir();
      fab_size += flux[2].nBytes();

      qe[2].resize(gzbx, NGDNV);
      Elixir elix_qe_z = qe[2].elixir();
      fab_size += qe[2].nBytes();

#ifdef RADIATION
      rad_flux[2].resize(gzbx, Radiation::nGroups);
      Elixir elix_rad_flux_z = rad_flux[2].elixir();
      fab_size += rad_flux[2].nBytes();
#endif
#endif

#if AMREX_SPACEDIM <= 2
      if (!Geom().IsCartesian()) {
          pradial.resize(xbx, 1);
      }
      Elixir elix_pradial = pradial.elixir();
      fab_size += pradial.nBytes();
#endif

#if AMREX_SPACEDIM == 1
#pragma gpu box(xbx)
      cmpflx_plus_godunov(AMREX_INT_ANYD(xbx.loVect()), AMREX_INT_ANYD(xbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qxm),
                          BL_TO_FORTRAN_ANYD(qxp),
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
      fab_size += ftmp1.nBytes();

      ftmp2.resize(obx, NUM_STATE);
      Elixir elix_ftmp2 = ftmp2.elixir();
      fab_size += ftmp2.nBytes();

#ifdef RADIATION
      rftmp1.resize(obx, Radiation::nGroups);
      Elixir elix_rftmp1 = rftmp1.elixir();
      fab_size += rftmp1.nBytes();

      rftmp2.resize(obx, Radiation::nGroups);
      Elixir elix_rftmp2 = rftmp2.elixir();
      fab_size += rftmp2.nBytes();
#endif

      qgdnvtmp1.resize(obx, NGDNV);
      Elixir elix_qgdnvtmp1 = qgdnvtmp1.elixir();
      fab_size += qgdnvtmp1.nBytes();

      qgdnvtmp2.resize(obx, NGDNV);
      Elixir elix_qgdnvtmp2 = qgdnvtmp2.elixir();
      fab_size += qgdnvtmp2.nBytes();

      ql.resize(obx, NQ);
      Elixir elix_ql = ql.elixir();
      fab_size += ql.nBytes();

      qr.resize(obx, NQ);
      Elixir elix_qr = qr.elixir();
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
#pragma gpu box(cxbx)
      cmpflx_plus_godunov(AMREX_INT_ANYD(cxbx.loVect()), AMREX_INT_ANYD(cxbx.hiVect()),
                          BL_TO_FORTRAN_ANYD(qxm),
                          BL_TO_FORTRAN_ANYD(qxp),
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
                          BL_TO_FORTRAN_ANYD(qyp),
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
                          BL_TO_FORTRAN_ANYD(qr),
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
                          BL_TO_FORTRAN_ANYD(qr),
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
                          BL_TO_FORTRAN_ANYD(qxp),
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
      fab_size += qmyx.nBytes();

      qpyx.resize(tyxbx, NQ);
      Elixir elix_qpyx = qpyx.elixir();
      fab_size += qpyx.nBytes();

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
      fab_size += qmzx.nBytes();

      qpzx.resize(tzxbx, NQ);
      Elixir elix_qpzx = qpzx.elixir();
      fab_size += qpzx.nBytes();

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
                          BL_TO_FORTRAN_ANYD(qyp),
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
      fab_size += qmxy.nBytes();

      qpxy.resize(txybx, NQ);
      Elixir elix_qpxy = qpxy.elixir();
      fab_size += qpxy.nBytes();

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
      fab_size += qmzy.nBytes();

      qpzy.resize(tzybx, NQ);
      Elixir elix_qpzy = qpzy.elixir();
      fab_size += qpzy.nBytes();

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
                          BL_TO_FORTRAN_ANYD(qzp),
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
      fab_size += qmxz.nBytes();

      qpxz.resize(txzbx, NQ);
      Elixir elix_qpxz = qpxz.elixir();
      fab_size += qpxz.nBytes();

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
      fab_size += qmyz.nBytes();

      qpyz.resize(tyzbx, NQ);
      Elixir elix_qpyz = qpyz.elixir();
      fab_size += qpyz.nBytes();

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
                          BL_TO_FORTRAN_ANYD(qpyz),
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
                          BL_TO_FORTRAN_ANYD(qpzy),
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
                          BL_TO_FORTRAN_ANYD(qr),
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
                          BL_TO_FORTRAN_ANYD(qpzx),
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
                          BL_TO_FORTRAN_ANYD(qpxz),
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
                          BL_TO_FORTRAN_ANYD(qr),
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
                          BL_TO_FORTRAN_ANYD(qpxy),
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
                          BL_TO_FORTRAN_ANYD(qpyx),
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
                          BL_TO_FORTRAN_ANYD(qr),
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

          // Zero out shock and temp fluxes -- these are physically meaningless here
          AMREX_PARALLEL_FOR_3D(nbx, i, j, k,
          {
              flux_arr(i,j,k,UTEMP) = 0.e0;
#ifdef SHOCK_VAR
              flux_arr(i,j,k,USHK) = 0.e0;
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
                   BL_TO_FORTRAN_ANYD(q[mfi]),
                   BL_TO_FORTRAN_ANYD(volume[mfi]),
                   BL_TO_FORTRAN_ANYD(flux[idir]),
                   BL_TO_FORTRAN_ANYD(area[idir][mfi]),
                   dt, AMREX_REAL_ANYD(dx));
          }

          if (limit_fluxes_on_large_vel == 1) {
#pragma gpu box(nbx)
              limit_hydro_fluxes_on_large_vel
                  (AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
                   idir_f,
                   BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                   BL_TO_FORTRAN_ANYD(q[mfi]),
                   BL_TO_FORTRAN_ANYD(volume[mfi]),
                   BL_TO_FORTRAN_ANYD(flux[idir]),
                   BL_TO_FORTRAN_ANYD(area[idir][mfi]),
                   dt, AMREX_REAL_ANYD(dx));
          }

#pragma gpu box(nbx)
          normalize_species_fluxes(AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
                                   BL_TO_FORTRAN_ANYD(flux[idir]));

      }



      // conservative update
      Array4<Real> const update_arr = hydro_source.array(mfi);

      Array4<Real> const flx_arr = (flux[0]).array();
      Array4<Real> const qx_arr = (qe[0]).array();
      Array4<Real> const areax_arr = (area[0]).array(mfi);

#if AMREX_SPACEDIM >= 2
      Array4<Real> const fly_arr = (flux[1]).array();
      Array4<Real> const qy_arr = (qe[1]).array();
      Array4<Real> const areay_arr = (area[1]).array(mfi);
#endif

#if AMREX_SPACEDIM == 3
      Array4<Real> const flz_arr = (flux[2]).array();
      Array4<Real> const qz_arr = (qe[2]).array();
      Array4<Real> const areaz_arr = (area[2]).array(mfi);
#endif

      Array4<Real> const vol_arr = volume.array(mfi);

      consup_hydro(bx,
                   shk_arr,
                   update_arr,
                   flx_arr, qx_arr, areax_arr,
#if AMREX_SPACEDIM >= 2
                   fly_arr, qy_arr, areay_arr,
#endif
#if AMREX_SPACEDIM == 3
                   flz_arr, qz_arr, areaz_arr,
#endif
                   vol_arr,
                   dt);


#ifdef HYBRID_MOMENTUM
#pragma gpu box(bx)
    add_hybrid_advection_source(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                                dt,
                                BL_TO_FORTRAN_ANYD(hydro_source[mfi]),
                                BL_TO_FORTRAN_ANYD(qe[0]),
                                BL_TO_FORTRAN_ANYD(qe[1]),
                                BL_TO_FORTRAN_ANYD(qe[2]));
#endif

#ifdef RADIATION
#pragma gpu box(bx)
      ctu_rad_consup(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                     BL_TO_FORTRAN_ANYD(hydro_source[mfi]),
                     BL_TO_FORTRAN_ANYD(Erborder[mfi]),
                     BL_TO_FORTRAN_ANYD(S_new[mfi]),
                     BL_TO_FORTRAN_ANYD(Er_new[mfi]),
                     BL_TO_FORTRAN_ANYD(rad_flux[0]),
                     BL_TO_FORTRAN_ANYD(qe[0]),
                     BL_TO_FORTRAN_ANYD(area[0][mfi]),
#if AMREX_SPACEDIM >= 2
                     BL_TO_FORTRAN_ANYD(rad_flux[1]),
                     BL_TO_FORTRAN_ANYD(qe[1]),
                     BL_TO_FORTRAN_ANYD(area[1][mfi]),
#endif
#if AMREX_SPACEDIM == 3
                     BL_TO_FORTRAN_ANYD(rad_flux[2]),
                     BL_TO_FORTRAN_ANYD(qe[2]),
                     BL_TO_FORTRAN_ANYD(area[2][mfi]),
#endif
                     &priv_nstep_fsp,
                     BL_TO_FORTRAN_ANYD(volume[mfi]),
                     AMREX_REAL_ANYD(dx), dt);

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
            if (!Geom().IsCartesian()) {
                AMREX_PARALLEL_FOR_3D(nbx, i, j, k,
                {
                    pradial_fab(i,j,k) = qex_fab(i,j,k,prescomp) * dt;
                });
            }
#endif

#if AMREX_SPACEDIM == 2
            if (!momx_flux_has_p[0]) {
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

        AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(idir), 1, i, j, k, n,
        {
            mass_fluxes_fab(i,j,k,0) = flux_fab(i,j,k,URHO);
        });

      } // idir loop

#if AMREX_SPACEDIM <= 2
      if (!Geom().IsCartesian()) {

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

#ifdef AMREX_USE_GPU
      // Check if we're going to run out of memory in the next MFIter iteration.
      // If so, do a synchronize here so that we don't oversubscribe GPU memory.
      // Note that this will capture the case where we started with more memory
      // than what the GPU has, on the logic that even in that case, it makes
      // sense to not further pile on the oversubscription demands.

      // This could (and should) be generalized in the future to operate with
      // more granularity than the MFIter loop boundary. We would have potential
      // synchronization points prior to each of the above kernel launches, and
      // we would check whether the sum of all previously allocated fabs would
      // result in oversubscription, including any contributions from a partial
      // MFIter loop. A further optimization would be to not apply a device
      // synchronize, but rather to use CUDA events to poll on a check about
      // whether enough memory has freed up to begin the next iteration, and then
      // immediately proceed to the next kernel when there's enough space for it.

      current_size += fab_size;
      if (current_size + fab_size >= Gpu::Device::totalGlobalMem()) {
          Gpu::Device::synchronize();
          current_size = starting_size;
      }
#endif

      if (oversubscribed) {
          q[mfi].prefetchToHost();
          qaux[mfi].prefetchToHost();
          src_q[mfi].prefetchToHost();
          volume[mfi].prefetchToHost();
          Sborder[mfi].prefetchToHost();
          hydro_source[mfi].prefetchToHost();
          for (int i = 0; i < AMREX_SPACEDIM; ++i) {
              area[i][mfi].prefetchToHost();
              (*fluxes[i])[mfi].prefetchToHost();
          }
#if AMREX_SPACEDIM < 3
          dLogArea[0][mfi].prefetchToHost();
          P_radial[mfi].prefetchToHost();
#endif
#ifdef RADIATION
          Erborder[mfi].prefetchToHost();
          Er_new[mfi].prefetchToHost();
#endif
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
    std::cout << "... Leaving construct_ctu_hydro_source()" << std::endl << std::endl;

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
