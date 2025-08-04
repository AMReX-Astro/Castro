#include <Castro.H>

using namespace amrex;

#include <mhd_eigen.H>
#include <slope.H>
void
Castro::plm(const Box& bx,
            const int idir,
            Array4<Real const> const& s,
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

  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();

  bool lo_symm = lo_bc[idir] == amrex::PhysBCType::symmetry;
  bool hi_symm = hi_bc[idir] == amrex::PhysBCType::symmetry;

  const auto domlo = geom.Domain().loVect3d();
  const auto domhi = geom.Domain().hiVect3d();

  Real dtdx = dt/dx[idir];

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    bool lo_bc_test = lo_symm && ((idir == 0 && i == domlo[0]) ||
                                  (idir == 1 && j == domlo[1]) ||
                                  (idir == 2 && k == domlo[2]));

    bool hi_bc_test = hi_symm && ((idir == 0 && i == domhi[0]) ||
                                  (idir == 1 && j == domhi[1]) ||
                                  (idir == 2 && k == domhi[2]));

    // compute the 1-sided differences used for the slopes

    Real Q[NEIGN][nslp];

    // we use a reduced eigensystem, the normal B field component is
    // omitted

    // for reflect BCs, we need to know what the normal component of
    // the velocity is
    int QUN;

    load_stencil(s, idir, i, j, k, QRHO, Q[IEIGN_RHO]);
    load_stencil(s, idir, i, j, k, QU, Q[IEIGN_U]);
    load_stencil(s, idir, i, j, k, QV, Q[IEIGN_V]);
    load_stencil(s, idir, i, j, k, QW, Q[IEIGN_W]);
    load_stencil(s, idir, i, j, k, QPRES, Q[IEIGN_P]);

    if (idir == 0) {
        QUN = IEIGN_U;

        // component (Bx) is omitted
        load_stencil(s, idir, i, j, k, QMAGY, Q[IEIGN_BT]);
        load_stencil(s, idir, i, j, k, QMAGZ, Q[IEIGN_BTT]);

    } else if (idir == 1) {

        QUN = IEIGN_V;

        // component (By) is omitted
        load_stencil(s, idir, i, j, k, QMAGX, Q[IEIGN_BT]);
        load_stencil(s, idir, i, j, k, QMAGZ, Q[IEIGN_BTT]);

    } else {

        QUN = IEIGN_W;

        // component (Bz) is omitted
        load_stencil(s, idir, i, j, k, QMAGX, Q[IEIGN_BT]);
        load_stencil(s, idir, i, j, k, QMAGY, Q[IEIGN_BTT]);

    }

    // compute the eigenvectors and eigenvalues for this coordinate direction

    Array1D<Real, 0, NQ-1> q_zone;
    for (int n = 0; n < NQ; n++) {
      q_zone(n) = s(i,j,k,n);
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

    // compute the slopes
    Real dq[NEIGN] = {};

    if (mhd_limit_characteristic == 1) {

      // we are limiting on characteristic variables
      for (int ii = 0; ii < NEIGN; ii++) {

        // construct the ii-th primitive variables
        Real W[5] = {};
        for (int n = 0; n < NEIGN; n++) {
          W[im2] += leig(ii,n) * Q[n][im2];
          W[im1] += leig(ii,n) * Q[n][im1];
          W[i0] += leig(ii,n) * Q[n][i0];
          W[ip1] += leig(ii,n) * Q[n][ip1];
          W[ip2] += leig(ii,n) * Q[n][ip2];
        }

        // now limit
        Real dW = uslope(W, flatn(i,j,k), false, false);

        // now add this charactistic variable's contribution to the
        // primitive variable slope
        for (int n = 0; n < NEIGN; n++) {
          dq[n] += dW * reig(n, ii);
        }
      }

    } else {

      // we are limiting on primitive variables

      for (int ii = 0; ii < NEIGN; ii++) {
        bool vtest = ii == QUN;
        dq[ii] = uslope(Q[ii], flatn(i,j,k), lo_bc_test && vtest, hi_bc_test && vtest);
      }

    }


    // Perform the characteristic projection.  Since we are using
    // Using HLLD, we sum over all eigenvalues -- see the discussion after Eq. 31
    Real summ_p[NEIGN] = {0.0_rt};
    Real summ_m[NEIGN] = {0.0_rt};

    for (int ii = 0; ii < NEIGN; ii++) {
      Real Ldq = 0.0;

      for (int n = 0; n < NEIGN; n++) {
        Ldq += leig(ii,n) * dq[n];
      }

      for (int n = 0; n < NEIGN; n++) {
        summ_p[n] += (1.0_rt - dtdx * lam(ii)) * Ldq * reig(n,ii);
        summ_m[n] -= (1.0_rt + dtdx * lam(ii)) * Ldq * reig(n,ii);
      }
    }

    // left state at i+1/2

    if (idir == 0) {
      qleft(i+1,j,k,QRHO) = amrex::max(small_dens,
                                       q_zone(QRHO) + 0.5_rt*summ_p[IEIGN_RHO] + 0.5_rt*dt*smhd[IEIGN_RHO]);
      qleft(i+1,j,k,QU) = q_zone(QU) + 0.5_rt*summ_p[IEIGN_U] + 0.5_rt*dt*smhd[IEIGN_U];
      qleft(i+1,j,k,QV) = q_zone(QV) + 0.5_rt*summ_p[IEIGN_V] + 0.5_rt*dt*smhd[IEIGN_V];
      qleft(i+1,j,k,QW) = q_zone(QW) + 0.5_rt*summ_p[IEIGN_W] + 0.5_rt*dt*smhd[IEIGN_W];
      qleft(i+1,j,k,QPRES) = amrex::max(small_pres,
                                        q_zone(QPRES) + 0.5_rt*summ_p[IEIGN_P] + 0.5_rt*dt*smhd[IEIGN_P]);

      qleft(i+1,j,k,QMAGX) = Bx(i+1,j,k); // Bx stuff
      qleft(i+1,j,k,QMAGY) = q_zone(QMAGY) + 0.5_rt*summ_p[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qleft(i+1,j,k,QMAGZ) = q_zone(QMAGZ) + 0.5_rt*summ_p[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];

    } else if (idir == 1) {
      qleft(i,j+1,k,QRHO) = amrex::max(small_dens,
                                       q_zone(QRHO) + 0.5_rt*summ_p[IEIGN_RHO] + 0.5_rt*dt*smhd[IEIGN_RHO]);
      qleft(i,j+1,k,QU) = q_zone(QU) + 0.5_rt*summ_p[IEIGN_U] + 0.5_rt*dt*smhd[IEIGN_U];
      qleft(i,j+1,k,QV) = q_zone(QV) + 0.5_rt*summ_p[IEIGN_V] + 0.5_rt*dt*smhd[IEIGN_V];
      qleft(i,j+1,k,QW) = q_zone(QW) + 0.5_rt*summ_p[IEIGN_W] + 0.5_rt*dt*smhd[IEIGN_W];
      qleft(i,j+1,k,QPRES) = amrex::max(small_pres,
                                        q_zone(QPRES) + 0.5_rt*summ_p[IEIGN_P] + 0.5_rt*dt*smhd[IEIGN_P]);

      qleft(i,j+1,k,QMAGX) = q_zone(QMAGX) + 0.5_rt*summ_p[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qleft(i,j+1,k,QMAGY) = By(i,j+1,k); // By stuff
      qleft(i,j+1,k,QMAGZ) = q_zone(QMAGZ) + 0.5_rt*summ_p[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];

    } else {
      qleft(i,j,k+1,QRHO) = amrex::max(small_dens,
                                       q_zone(QRHO) + 0.5_rt*summ_p[IEIGN_RHO] + 0.5_rt*dt*smhd[IEIGN_RHO]);
      qleft(i,j,k+1,QU) = q_zone(QU) + 0.5_rt*summ_p[IEIGN_U] + 0.5_rt*dt*smhd[IEIGN_U];
      qleft(i,j,k+1,QV) = q_zone(QV) + 0.5_rt*summ_p[IEIGN_V] + 0.5_rt*dt*smhd[IEIGN_V];
      qleft(i,j,k+1,QW) = q_zone(QW) + 0.5_rt*summ_p[IEIGN_W] + 0.5_rt*dt*smhd[IEIGN_W];
      qleft(i,j,k+1,QPRES) = amrex::max(small_pres,
                                        q_zone(QPRES) + 0.5_rt*summ_p[IEIGN_P] + 0.5_rt*dt*smhd[IEIGN_P]);

      qleft(i,j,k+1,QMAGX) = q_zone(QMAGX) + 0.5_rt*summ_p[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qleft(i,j,k+1,QMAGY) = q_zone(QMAGY) + 0.5_rt*summ_p[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];
      qleft(i,j,k+1,QMAGZ) = Bz(i,j,k+1); // Bz stuff
    }

    // right state at i-1/2
    qright(i,j,k,QRHO) = amrex::max(small_dens,
                                    q_zone(QRHO) + 0.5_rt*summ_m[IEIGN_RHO] + 0.5_rt*dt*smhd[IEIGN_RHO]);
    qright(i,j,k,QU) = q_zone(QU) + 0.5_rt*summ_m[IEIGN_U] + 0.5_rt*dt*smhd[IEIGN_U];
    qright(i,j,k,QV) = q_zone(QV) + 0.5_rt*summ_m[IEIGN_V] + 0.5_rt*dt*smhd[IEIGN_V];
    qright(i,j,k,QW) = q_zone(QW) + 0.5_rt*summ_m[IEIGN_W] + 0.5_rt*dt*smhd[IEIGN_W];
    qright(i,j,k,QPRES) = amrex::max(small_pres,
                                     q_zone(QPRES) + 0.5_rt*summ_m[IEIGN_P] + 0.5_rt*dt*smhd[IEIGN_P]);

    if (idir == 0) {
      qright(i,j,k,QMAGX) = Bx(i,j,k); // Bx stuff
      qright(i,j,k,QMAGY) = q_zone(QMAGY) + 0.5_rt*summ_m[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qright(i,j,k,QMAGZ) = q_zone(QMAGZ) + 0.5_rt*summ_m[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];

    } else if (idir == 1) {
      qright(i,j,k,QMAGX) = q_zone(QMAGX) + 0.5_rt*summ_m[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qright(i,j,k,QMAGY) = By(i,j,k); // By stuff
      qright(i,j,k,QMAGZ) = q_zone(QMAGZ) + 0.5_rt*summ_m[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];

    } else {
      qright(i,j,k,QMAGX) = q_zone(QMAGX) + 0.5_rt*summ_m[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qright(i,j,k,QMAGY) = q_zone(QMAGY) + 0.5_rt*summ_m[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];
      qright(i,j,k,QMAGZ) = Bz(i,j,k); // Bz stuff
    }

    // species
    for (int n = 0; n < NumSpec; n++) {
      Real un{};
      Real X[nslp];

      load_stencil(s, idir, i, j, k, QFS+n, X);

      if (idir == 0) {
        un = s(i,j,k,QU);
      } else if (idir == 1) {
        un = s(i,j,k,QV);
      } else {
        un = s(i,j,k,QW);
      }

      Real dX = uslope(X, flatn(i,j,k), false, false);

      if (idir == 0) {
        qleft(i+1,j,k,QFS+n) = q_zone(QFS+n) + 0.5_rt*(1.0_rt - dtdx*un) * dX;
      } else if (idir == 1) {
        qleft(i,j+1,k,QFS+n) = q_zone(QFS+n) + 0.5_rt*(1.0_rt - dtdx*un) * dX;
      } else {
        qleft(i,j,k+1,QFS+n) = q_zone(QFS+n) + 0.5_rt*(1.0_rt - dtdx*un) * dX;
      }
      qright(i,j,k,QFS+n) = q_zone(QFS+n) - 0.5_rt*(1.0_rt + dtdx*un) * dX;
    }

    // rho e
    eos_t eos_state;

    if (idir == 0) {
      eos_state.rho = qleft(i+1,j,k,QRHO);
      eos_state.p = qleft(i+1,j,k,QPRES);
      eos_state.T = s(i,j,k,QTEMP); // some initial guess?
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
      eos_state.T = s(i,j,k,QTEMP); // some initial guess?
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
      eos_state.T = s(i,j,k,QTEMP); // some initial guess?
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
