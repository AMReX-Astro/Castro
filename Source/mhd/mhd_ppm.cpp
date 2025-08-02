#include <Castro.H>

#include <ppm.H>

using namespace amrex;

#include <mhd_eigen.H>

void
Castro::ppm_mhd(const Box& bx,
                const int idir,
                Array4<Real const> const& q_arr,
                Array4<Real const> const& qaux,
                Array4<Real const> const& flatn,
                Array4<Real const> const& Bx,
                Array4<Real const> const& By,
                Array4<Real const> const& Bz,
                Array4<Real> const& qleft,
                Array4<Real> const& qright,
                Array4<Real const> const& srcQ,
                const Real dt) {

  // these loops are over cell-centers and for each cell-center, we find the left and
  // right interface states

  const auto dx = geom.CellSizeArray();

  Real dtdx = dt/dx[idir];

  // special care for reflecting BCs
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();

  const auto domlo = geom.Domain().loVect3d();
  const auto domhi = geom.Domain().hiVect3d();

  bool lo_symm = lo_bc[idir] == amrex::PhysBCType::symmetry;
  bool hi_symm = hi_bc[idir] == amrex::PhysBCType::symmetry;

  // these are the characteristic variables for this direction
  int cvars[NEIGN];

  if (idir == 0) {

    // component (Bx) is omitted

    cvars[IEIGN_RHO] = QRHO;
    cvars[IEIGN_U] = QU;
    cvars[IEIGN_V] = QV;
    cvars[IEIGN_W] = QW;
    cvars[IEIGN_P] = QPRES;
    cvars[IEIGN_BT] = QMAGY;
    cvars[IEIGN_BTT] = QMAGZ;

  } else if (idir == 1) {

    // component (By) is omitted

    cvars[IEIGN_RHO] = QRHO;
    cvars[IEIGN_U] = QU;
    cvars[IEIGN_V] = QV;
    cvars[IEIGN_W] = QW;
    cvars[IEIGN_P] = QPRES;
    cvars[IEIGN_BT] = QMAGX;
    cvars[IEIGN_BTT] = QMAGZ;

  } else {

    // component (Bz) is omitted

    cvars[IEIGN_RHO] = QRHO;
    cvars[IEIGN_U] = QU;
    cvars[IEIGN_V] = QV;
    cvars[IEIGN_W] = QW;
    cvars[IEIGN_P] = QPRES;
    cvars[IEIGN_BT] = QMAGX;
    cvars[IEIGN_BTT] = QMAGY;

  }


  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    // compute the eigenvectors and eigenvalues for this coordinate direction

    Array1D<Real, 0, NQ-1> q_zone;
    for (int n = 0; n < NQ; n++) {
      q_zone(n) = q_arr(i,j,k,n);
    }

    Real as = qaux(i,j,k,QC);

    Array1D<Real, 0, NEIGN-1> lam;

    evals(lam, as, q_zone, idir);

    Array2D<Real, 0, NEIGN-1, 0, NEIGN-1> leig;
    Array2D<Real, 0, NEIGN-1, 0, NEIGN-1> reig;

    if (idir == 0) {
      evecx(leig, reig, as, q_zone);

    } else if (idir == 1) {
      evecy(leig, reig, as, q_zone);

    } else {
      evecz(leig, reig, as, q_zone);
    }

    // do the parabolic reconstruction and compute the integrals under
    // the characteristic waves
    Real s[nslp];
    Real flat = flatn(i,j,k);
    Real sm;
    Real sp;

    Real Ip[NEIGN][NEIGN];  // first component is var, second is wave
    Real Im[NEIGN][NEIGN];

    Real q_ref_left[NEIGN];
    Real q_ref_right[NEIGN];


    for (int n = 0; n < NEIGN; n++) {

      int v = cvars[n];

      load_stencil(q_arr, idir, i, j, k, v, s);
      ppm_reconstruct(s, i, j, k, idir, reconstruction::Centering::zone_centered, flat, sm, sp);

      Real Ipt = 0.0;
      Real Imt = 0.0;

      // loop over waves and integrate under the parabola
      for (int ii = 0; ii < NEIGN; ii++) {
        Real speed = lam(ii);

        ppm_int_profile_single(sm, sp, s[i0], speed, dtdx, Ipt, Imt);

        // we do a special correction here for waves that are not
        // moving toward the interface.  It doesn't make sense to
        // integrate under the parabola in this case, so we simply
        // create a piecewise linear slope here.  See Stone et al. Eq
        // 44 and 45
        if (lam(ii) >= 0.0_rt) {
          Ip[n][ii] = Ipt;
        } else {
          Ip[n][ii] = 0.5_rt * dtdx * (sp - sm);
        }

        if (lam(ii) <= 0.0_rt) {
          Im[n][ii] = Imt;
        } else {
          Im[n][ii] = 0.5_rt * dtdx * (sp - sm);
        }

      }

      // and now the reference states

      // q_{i+1/2, L} uses the fastest wave moving to the right
      Real lambda_max = amrex::max(lam(NEIGN-1), 0.0_rt);
      ppm_int_profile_single(sm, sp, s[i0], lambda_max, dtdx, Ipt, Imt);
      q_ref_left[n] = Ipt;

      // q_{i+1/2, L} uses the fastest wave moving to the left
      Real lambda_min = amrex::min(lam(0), 0.0_rt);
      ppm_int_profile_single(sm, sp, s[i0], lambda_min, dtdx, Ipt, Imt);
      q_ref_right[n] = Imt;

    }

    // MHD Source Terms -- from the Miniati paper, Eq. 32 and 33
    Real smhd[NEIGN];

    smhd[IEIGN_RHO] = 0.0_rt;
    smhd[IEIGN_U] = q_zone(QMAGX) / q_zone(QRHO);
    smhd[IEIGN_V] = q_zone(QMAGY) / q_zone(QRHO);
    smhd[IEIGN_W] = q_zone(QMAGZ) / q_zone(QRHO);
    smhd[IEIGN_P] = q_zone(QMAGX) * q_zone(QU) +
                    q_zone(QMAGY) * q_zone(QV) +
                    q_zone(QMAGZ) * q_zone(QW);

    if (idir == 0) {
      smhd[IEIGN_BT] = q_zone(QV);
      smhd[IEIGN_BTT] = q_zone(QW);

      // cross-talk of normal magnetic field direction
      for (auto & source : smhd) {
          source *= (Bx(i+1,j,k) - Bx(i,j,k)) / dx[idir];
      }

    } else if (idir == 1) {
      smhd[IEIGN_BT] = q_zone(QU);
      smhd[IEIGN_BTT] = q_zone(QW);

      // cross-talk of normal magnetic field direction
      for (auto & source : smhd) {
          source *= (By(i,j+1,k) - By(i,j,k)) / dx[idir];
      }

    } else {
      smhd[IEIGN_BT] = q_zone(QU);
      smhd[IEIGN_BTT] = q_zone(QV);

      // cross-talk of normal magnetic field direction
      for (auto & source : smhd) {
          source *= (Bz(i,j,k+1) - Bz(i,j,k)) / dx[idir];
      }
    }

    // Perform the characteristic projection.  Since we are using
    // Using HLLD, we sum over all eigenvalues -- see the discussion
    // after Eq. 31

    // right state at i-1/2

    // Im is the integral from the left edge, so we take as the
    // reference state the fastest wave moving to the left

    Real summ_m[NEIGN] = {0.0_rt};

    // loop over the waves
    for (int ii = 0; ii < NEIGN; ii++) {

      Real LdQ = 0.0_rt;

      // loop over variables in Im[n][ii]
      for (int n = 0; n < NEIGN; n++) {
        if (lam(ii) <= 0.0_rt) {
          LdQ += leig(ii,n) * (q_ref_right[n] - Im[n][ii]);
        } else {
          // in this case, the integral Im is a slope
          LdQ += (lam(0) - lam(ii)) * leig(ii,n) * Im[n][ii];
        }
      }

      // add the contribution of this wave to each variable
      for (int n = 0; n < NEIGN; n++) {
        summ_m[n] += LdQ * reig(n,ii);
      }
    }

    qright(i,j,k,QRHO) = amrex::max(small_dens,
                                    q_ref_right[IEIGN_RHO] - summ_m[IEIGN_RHO] + 0.5_rt*dt*smhd[IEIGN_RHO]);
    qright(i,j,k,QU) = q_ref_right[IEIGN_U] - summ_m[IEIGN_U] + 0.5_rt*dt*smhd[IEIGN_U];
    qright(i,j,k,QV) = q_ref_right[IEIGN_V] - summ_m[IEIGN_V] + 0.5_rt*dt*smhd[IEIGN_V];
    qright(i,j,k,QW) = q_ref_right[IEIGN_W] - summ_m[IEIGN_W] + 0.5_rt*dt*smhd[IEIGN_W];
    qright(i,j,k,QPRES) = amrex::max(small_pres,
                                     q_ref_right[IEIGN_P] - summ_m[IEIGN_P] + 0.5_rt*dt*smhd[IEIGN_P]);

    if (idir == 0) {
      qright(i,j,k,QMAGX) = Bx(i,j,k); // Bx stuff
      qright(i,j,k,QMAGY) = q_ref_right[IEIGN_BT] -  summ_m[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qright(i,j,k,QMAGZ) = q_ref_right[IEIGN_BTT] - summ_m[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];

    } else if (idir == 1) {
      qright(i,j,k,QMAGX) = q_ref_right[IEIGN_BT] - summ_m[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qright(i,j,k,QMAGY) = By(i,j,k); // By stuff
      qright(i,j,k,QMAGZ) = q_ref_right[IEIGN_BTT] - summ_m[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];

    } else {
      qright(i,j,k,QMAGX) = q_ref_right[IEIGN_BT] - summ_m[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qright(i,j,k,QMAGY) = q_ref_right[IEIGN_BTT] - summ_m[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];
      qright(i,j,k,QMAGZ) = Bz(i,j,k); // Bz stuff
    }


    // left state at i+1/2

    // Ip is the integral from the right edge, so we take as the
    // reference state the fastest wave moving to the right

    Real summ_p[NEIGN] = {0.0_rt};

    // loop over the waves
    for (int ii = 0; ii < NEIGN; ii++) {

      Real LdQ = 0.0_rt;

      // loop over variables in Im[n][ii]
      for (int n = 0; n < NEIGN; n++) {
        if (lam(ii) >= 0.0_rt) {
          LdQ += leig(ii,n) * (q_ref_left[n] - Ip[n][ii]);
        } else {
          // in this case, the integral Ip is a slope
          LdQ += (lam(NEIGN-1) - lam(ii)) * leig(ii,n) * Ip[n][ii];
        }
      }

      // add the contribution of this wave to each variable
      for (int n = 0; n < NEIGN; n++) {
        summ_p[n] += LdQ * reig(n,ii);
      }
    }

    if (idir == 0) {
      qleft(i+1,j,k,QRHO) = amrex::max(small_dens,
                                       q_ref_left[IEIGN_RHO] - summ_p[IEIGN_RHO] + 0.5_rt*dt*smhd[IEIGN_RHO]);
      qleft(i+1,j,k,QU) = q_ref_left[IEIGN_U] - summ_p[IEIGN_U] + 0.5_rt*dt*smhd[IEIGN_U];
      qleft(i+1,j,k,QV) = q_ref_left[IEIGN_V] - summ_p[IEIGN_V] + 0.5_rt*dt*smhd[IEIGN_V];
      qleft(i+1,j,k,QW) = q_ref_left[IEIGN_W] - summ_p[IEIGN_W] + 0.5_rt*dt*smhd[IEIGN_W];
      qleft(i+1,j,k,QPRES) = amrex::max(small_pres,
                                        q_ref_left[IEIGN_P] - summ_p[IEIGN_P] + 0.5_rt*dt*smhd[IEIGN_P]);

      qleft(i+1,j,k,QMAGX) = Bx(i+1,j,k); // Bx stuff
      qleft(i+1,j,k,QMAGY) = q_ref_left[IEIGN_BT] - summ_p[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qleft(i+1,j,k,QMAGZ) = q_ref_left[IEIGN_BTT] - summ_p[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];

    } else if (idir == 1) {
      qleft(i,j+1,k,QRHO) = amrex::max(small_dens,
                                       q_ref_left[IEIGN_RHO] - summ_p[IEIGN_RHO] + 0.5_rt*dt*smhd[IEIGN_RHO]);
      qleft(i,j+1,k,QU) = q_ref_left[IEIGN_U] - summ_p[IEIGN_U] + 0.5_rt*dt*smhd[IEIGN_U];
      qleft(i,j+1,k,QV) = q_ref_left[IEIGN_V] - summ_p[IEIGN_V] + 0.5_rt*dt*smhd[IEIGN_V];
      qleft(i,j+1,k,QW) = q_ref_left[IEIGN_W] - summ_p[IEIGN_W] + 0.5_rt*dt*smhd[IEIGN_W];
      qleft(i,j+1,k,QPRES) = amrex::max(small_pres,
                                        q_ref_left[IEIGN_P] - summ_p[IEIGN_P] + 0.5_rt*dt*smhd[IEIGN_P]);

      qleft(i,j+1,k,QMAGX) = q_ref_left[IEIGN_BT] - summ_p[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qleft(i,j+1,k,QMAGY) = By(i,j+1,k); // By stuff
      qleft(i,j+1,k,QMAGZ) = q_ref_left[IEIGN_BTT] - summ_p[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];

    } else {
      qleft(i,j,k+1,QRHO) = amrex::max(small_dens,
                                       q_ref_left[IEIGN_RHO] - summ_p[IEIGN_RHO] + 0.5_rt*dt*smhd[IEIGN_RHO]);
      qleft(i,j,k+1,QU) = q_ref_left[IEIGN_U] - summ_p[IEIGN_U] + 0.5_rt*dt*smhd[IEIGN_U];
      qleft(i,j,k+1,QV) = q_ref_left[IEIGN_V] - summ_p[IEIGN_V] + 0.5_rt*dt*smhd[IEIGN_V];
      qleft(i,j,k+1,QW) = q_ref_left[IEIGN_W] - summ_p[IEIGN_W] + 0.5_rt*dt*smhd[IEIGN_W];
      qleft(i,j,k+1,QPRES) = amrex::max(small_pres,
                                        q_ref_left[IEIGN_P] - summ_p[IEIGN_P] + 0.5_rt*dt*smhd[IEIGN_P]);

      qleft(i,j,k+1,QMAGX) = q_ref_left[IEIGN_BT] - summ_p[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qleft(i,j,k+1,QMAGY) = q_ref_left[IEIGN_BTT] - summ_p[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];
      qleft(i,j,k+1,QMAGZ) = Bz(i,j,k+1); // Bz stuff
    }


    // species
    for (int n = 0; n < NumSpec; n++) {

      int v = QFS+n;

      Real Ips;
      Real Ims;

      Real un;

      if (idir == 0) {
        un = q_arr(i,j,k,QU);

      } else if (idir == 1) {
        un = q_arr(i,j,k,QV);

      } else {
        un = q_arr(i,j,k,QW);
      }

      load_stencil(q_arr, idir, i, j, k, v, s);
      ppm_reconstruct(s, i, j, k, idir, reconstruction::Centering::zone_centered, flat, sm, sp);
      ppm_int_profile_single(sm, sp, s[i0], un, dtdx, Ips, Ims);

      if (idir == 0) {
        qleft(i+1,j,k,QFS+n) = Ips;
      } else if (idir == 1) {
        qleft(i,j+1,k,QFS+n) = Ips;
      } else {
        qleft(i,j,k+1,QFS+n) = Ips;
      }
      qright(i,j,k,QFS+n) = Ims;
    }

    // rho e
    eos_t eos_state;

    if (idir == 0) {
      eos_state.rho = qleft(i+1,j,k,QRHO);
      eos_state.p = qleft(i+1,j,k,QPRES);
      eos_state.T = q_arr(i,j,k,QTEMP); // some initial guess?
      for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = qleft(i+1,j,k,QFS+n);
      }
#if NAUX_NET > 0
      for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = qleft(i+1,j,k,QFX+n);
      }
#endif
      eos(eos_input_rp, eos_state);
      qleft(i+1,j,k,QREINT) = eos_state.e * eos_state.rho;

    } else if (idir == 1) {
      eos_state.rho = qleft(i,j+1,k,QRHO);
      eos_state.p = qleft(i,j+1,k,QPRES);
      eos_state.T = q_arr(i,j,k,QTEMP); // some initial guess?
      for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = qleft(i,j+1,k,QFS+n);
      }
#if NAUX_NET > 0
      for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = qleft(i,j+1,k,QFX+n);
      }
#endif
      eos(eos_input_rp, eos_state);
      qleft(i,j+1,k,QREINT) = eos_state.e * eos_state.rho;

    } else {
      eos_state.rho = qleft(i,j,k+1,QRHO);
      eos_state.p = qleft(i,j,k+1,QPRES);
      eos_state.T = q_arr(i,j,k,QTEMP); // some initial guess?
      for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = qleft(i,j,k+1,QFS+n);
      }
#if NAUX_NET > 0
      for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = qleft(i,j,k+1,QFX+n);
      }
#endif
      eos(eos_input_rp, eos_state);
      qleft(i,j,k+1,QREINT) = eos_state. e * eos_state.rho;
    }

    eos_state.rho = qright(i,j,k,QRHO);
    eos_state.p = qright(i,j,k,QPRES);
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = qright(i,j,k,QFS+n);
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = qright(i,j,k,QFX+n);
    }
#endif

    eos(eos_input_rp, eos_state);
    qright(i,j,k,QREINT) = eos_state.e * eos_state.rho;

    // add source terms
    if (idir == 0) {
      qleft(i+1,j,k,QRHO) = amrex::max(small_dens,
                                       qleft(i+1,j,k,QRHO) + 0.5_rt*dt*srcQ(i,j,k,QRHO));
      qleft(i+1,j,k,QU) = qleft(i+1,j,k,QU) + 0.5_rt*dt*srcQ(i,j,k,QU);
      qleft(i+1,j,k,QV) = qleft(i+1,j,k,QV) + 0.5_rt*dt*srcQ(i,j,k,QV);
      qleft(i+1,j,k,QW) = qleft(i+1,j,k,QW) + 0.5_rt*dt*srcQ(i,j,k,QW);
      qleft(i+1,j,k,QPRES) = qleft(i+1,j,k,QPRES) + 0.5_rt*dt*srcQ(i,j,k,QPRES);
      qleft(i+1,j,k,QREINT) = qleft(i+1,j,k,QREINT) + 0.5_rt*dt*srcQ(i,j,k,QREINT);

    } else if (idir == 1) {
      qleft(i,j+1,k,QRHO) = amrex::max(small_dens,
                                       qleft(i,j+1,k,QRHO) + 0.5_rt*dt*srcQ(i,j,k,QRHO));
      qleft(i,j+1,k,QU) = qleft(i,j+1,k,QU) + 0.5_rt*dt*srcQ(i,j,k,QU);
      qleft(i,j+1,k,QV) = qleft(i,j+1,k,QV) + 0.5_rt*dt*srcQ(i,j,k,QV);
      qleft(i,j+1,k,QW) = qleft(i,j+1,k,QW) + 0.5_rt*dt*srcQ(i,j,k,QW);
      qleft(i,j+1,k,QPRES) = qleft(i,j+1,k,QPRES) + 0.5_rt*dt*srcQ(i,j,k,QPRES);
      qleft(i,j+1,k,QREINT) = qleft(i,j+1,k,QREINT) + 0.5_rt*dt*srcQ(i,j,k,QREINT);

    } else {
      qleft(i,j,k+1,QRHO) = amrex::max(small_dens, qleft(i,j,k+1,QRHO) + 0.5_rt*dt*srcQ(i,j,k,QRHO));
      qleft(i,j,k+1,QU) = qleft(i,j,k+1,QU) + 0.5_rt*dt*srcQ(i,j,k,QU);
      qleft(i,j,k+1,QV) = qleft(i,j,k+1,QV) + 0.5_rt*dt*srcQ(i,j,k,QV);
      qleft(i,j,k+1,QW) = qleft(i,j,k+1,QW) + 0.5_rt*dt*srcQ(i,j,k,QW);
      qleft(i,j,k+1,QPRES) = qleft(i,j,k+1,QPRES) + 0.5_rt*dt*srcQ(i,j,k,QPRES);
      qleft(i,j,k+1,QREINT) = qleft(i,j,k+1,QREINT) + 0.5_rt*dt*srcQ(i,j,k,QREINT);
    }

    qright(i,j,k,QRHO) = amrex::max(small_dens, qright(i,j,k,QRHO) + 0.5_rt*dt*srcQ(i,j,k,QRHO));
    qright(i,j,k,QU) = qright(i,j,k,QU) + 0.5_rt*dt*srcQ(i,j,k,QU);
    qright(i,j,k,QV) = qright(i,j,k,QV) + 0.5_rt*dt*srcQ(i,j,k,QV);
    qright(i,j,k,QW) = qright(i,j,k,QW) + 0.5_rt*dt*srcQ(i,j,k,QW);
    qright(i,j,k,QPRES) = qright(i,j,k,QPRES) + 0.5_rt*dt*srcQ(i,j,k,QPRES);
    qright(i,j,k,QREINT) = qright(i,j,k,QREINT) + 0.5_rt*dt*srcQ(i,j,k,QREINT);

  });
}
