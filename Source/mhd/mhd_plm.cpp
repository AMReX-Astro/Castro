#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

#include "mhd_eigen.H"
#include "mhd_slope.H"

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

  Real dtdx = dt/dx[idir];

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    // compute the 1-sided differences used for the slopes

    Real dQL[NEIGN];
    Real dQR[NEIGN];

    // we use a reduced eigensystem, the normal B field component is
    // omitted

    if (idir == 0) {

      // component (Bx) is omitted

      dQL[IEIGN_RHO] = s(i,j,k,QRHO) - s(i-1,j,k,QRHO);
      dQL[IEIGN_U] = s(i,j,k,QU) - s(i-1,j,k,QU);
      dQL[IEIGN_V] = s(i,j,k,QV) - s(i-1,j,k,QV);
      dQL[IEIGN_W] = s(i,j,k,QW) - s(i-1,j,k,QW);
      dQL[IEIGN_P] = s(i,j,k,QPRES) - s(i-1,j,k,QPRES);
      dQL[IEIGN_BT] = s(i,j,k,QMAGY) - s(i-1,j,k,QMAGY);
      dQL[IEIGN_BTT] = s(i,j,k,QMAGZ) - s(i-1,j,k,QMAGZ);

      dQR[IEIGN_RHO] = s(i+1,j,k,QRHO) - s(i,j,k,QRHO);
      dQR[IEIGN_U] = s(i+1,j,k,QU) - s(i,j,k,QU);
      dQR[IEIGN_V] = s(i+1,j,k,QV) - s(i,j,k,QV);
      dQR[IEIGN_W] = s(i+1,j,k,QW) - s(i,j,k,QW);
      dQR[IEIGN_P] = s(i+1,j,k,QPRES) - s(i,j,k,QPRES);
      dQR[IEIGN_BT] = s(i+1,j,k,QMAGY) - s(i,j,k,QMAGY);
      dQR[IEIGN_BTT] = s(i+1,j,k,QMAGZ) - s(i,j,k,QMAGZ);

    } else if (idir == 1) {

      // component (By) is omitted

      dQL[IEIGN_RHO] = s(i,j,k,QRHO) - s(i,j-1,k,QRHO);
      dQL[IEIGN_U] = s(i,j,k,QU) - s(i,j-1,k,QU);
      dQL[IEIGN_V] = s(i,j,k,QV) - s(i,j-1,k,QV);
      dQL[IEIGN_W] = s(i,j,k,QW) - s(i,j-1,k,QW);
      dQL[IEIGN_P] = s(i,j,k,QPRES) - s(i,j-1,k,QPRES);
      dQL[IEIGN_BT] = s(i,j,k,QMAGX) - s(i,j-1,k,QMAGX);
      dQL[IEIGN_BTT] = s(i,j,k,QMAGZ) - s(i,j-1,k,QMAGZ);

      dQR[IEIGN_RHO] = s(i,j+1,k,QRHO) - s(i,j,k,QRHO);
      dQR[IEIGN_U] = s(i,j+1,k,QU) - s(i,j,k,QU);
      dQR[IEIGN_V] = s(i,j+1,k,QV) - s(i,j,k,QV);
      dQR[IEIGN_W] = s(i,j+1,k,QW) - s(i,j,k,QW);
      dQR[IEIGN_P] = s(i,j+1,k,QPRES) - s(i,j,k,QPRES);
      dQR[IEIGN_BT] = s(i,j+1,k,QMAGX) - s(i,j,k,QMAGX);
      dQR[IEIGN_BTT] = s(i,j+1,k,QMAGZ) - s(i,j,k,QMAGZ);

    } else {

      // component (Bz) is omitted

      dQL[IEIGN_RHO] = s(i,j,k,QRHO) - s(i,j,k-1,QRHO);
      dQL[IEIGN_U] = s(i,j,k,QU) - s(i,j,k-1,QU);
      dQL[IEIGN_V] = s(i,j,k,QV) - s(i,j,k-1,QV);
      dQL[IEIGN_W] = s(i,j,k,QW) - s(i,j,k-1,QW);
      dQL[IEIGN_P] = s(i,j,k,QPRES) - s(i,j,k-1,QPRES);
      dQL[IEIGN_BT] = s(i,j,k,QMAGX) - s(i,j,k-1,QMAGX);
      dQL[IEIGN_BTT] = s(i,j,k,QMAGY) - s(i,j,k-1,QMAGY);

      dQR[IEIGN_RHO] = s(i,j,k+1,QRHO) - s(i,j,k,QRHO);
      dQR[IEIGN_U] = s(i,j,k+1,QU) - s(i,j,k,QU);
      dQR[IEIGN_V] = s(i,j,k+1,QV) - s(i,j,k,QV);
      dQR[IEIGN_W] = s(i,j,k+1,QW) - s(i,j,k,QW);
      dQR[IEIGN_P] = s(i,j,k+1,QPRES) - s(i,j,k,QPRES);
      dQR[IEIGN_BT] = s(i,j,k+1,QMAGX) - s(i,j,k,QMAGX);
      dQR[IEIGN_BTT] = s(i,j,k+1,QMAGY) - s(i,j,k,QMAGY);

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
      for (int n = 0; n < NEIGN; n++) {
        smhd[n] = smhd[n] * (Bx(i+1,j,k) - Bx(i,j,k)) / dx[idir];
      }

    } else if (idir == 1) {
      smhd[IEIGN_BT] = q_zone(QU);
      smhd[IEIGN_BTT] = q_zone(QW);

      // cross-talk of normal magnetic field direction
      for (int n = 0; n < NEIGN; n++) {
        smhd[n] = smhd[n] * (By(i,j+1,k) - By(i,j,k)) / dx[idir];
      }

    } else {
      smhd[IEIGN_BT] = q_zone(QU);
      smhd[IEIGN_BTT] = q_zone(QV);

      // cross-talk of normal magnetic field direction
      for (int n = 0; n < NEIGN; n++) {
        smhd[n] = smhd[n] * (Bz(i,j,k+1) - Bz(i,j,k)) / dx[idir];
      }
    }

    // Perform the characteristic projection.  Since we are using
    // Using HLLD, we sum over all eigenvalues -- see the discussion after Eq. 31

    // We'll follow the method from Colella 1990 (as also written in Stone 2008)
    // where we first compute a reference state using the fastest moving eigenvalue
    // toward the interface and the project the difference of the waves from that

    // compute the limited primitive variable slopes where we limit on
    // the characteristic variables.  Note, we need this final result
    // ahead of time to compute the reference states

    Real dq_lim[NEIGN] = {};
    for (int ii = 0; ii < NEIGN; ii++) {
      Real dL = 0.0;
      Real dR = 0.0;

      for (int n = 0; n < NEIGN; n++) {
        dL += leig(ii,n) * dQL[n];
        dR += leig(ii,n) * dQR[n];
      }

      Real dW = 0.0;
      slope(dW, dL, dR, flatn(i,j,k));

      for (int n = 0; n < NEIGN; n++) {
        dq_lim[n] += dW * reig(n,ii);
      }
    }


    Real q_ref_left[NEIGN];
    Real ref_term = 0.5_rt * (1.0_rt - amrex::max(lam(NEIGN-1), 0.0_rt) * dtdx);

    q_ref_left[IEIGN_RHO] = q_zone(QRHO) + ref_term * dq_lim[QRHO];
    q_ref_left[IEIGN_U] = q_zone(QU) + ref_term * dq_lim[QU];
    q_ref_left[IEIGN_V] = q_zone(QV) + ref_term * dq_lim[QV];
    q_ref_left[IEIGN_W] = q_zone(QW) + ref_term * dq_lim[QW];
    q_ref_left[IEIGN_P] = q_zone(QPRES) + ref_term * dq_lim[QPRES];
    if (idir == 0) {
      q_ref_left[IEIGN_BT] = q_zone(QMAGY) + ref_term * dq_lim[IEIGN_BT];
      q_ref_left[IEIGN_BTT] = q_zone(QMAGZ) + ref_term * dq_lim[IEIGN_BTT];

    } else if (idir == 1) {
      q_ref_left[IEIGN_BT] = q_zone(QMAGX) + ref_term * dq_lim[IEIGN_BT];
      q_ref_left[IEIGN_BTT] = q_zone(QMAGZ) + ref_term * dq_lim[IEIGN_BTT];

    } else if (idir == 2) {
      q_ref_left[IEIGN_BT] = q_zone(QMAGX) + ref_term * dq_lim[IEIGN_BT];
      q_ref_left[IEIGN_BTT] = q_zone(QMAGY) + ref_term * dq_lim[IEIGN_BTT];

    }

    Real q_ref_right[NEIGN];
    // note: Stone et al. 2008 Eq. 41 has a sign error, compare to
    // Colella 1990 Eq. 2.11
    ref_term = 0.5_rt * (1.0_rt + amrex::min(lam(0), 0.0_rt) * dtdx);

    q_ref_right[IEIGN_RHO] = q_zone(QRHO) - ref_term * dq_lim[QRHO];
    q_ref_right[IEIGN_U] = q_zone(QU) - ref_term * dq_lim[QU];
    q_ref_right[IEIGN_V] = q_zone(QV) - ref_term * dq_lim[QV];
    q_ref_right[IEIGN_W] = q_zone(QW) - ref_term * dq_lim[QW];
    q_ref_right[IEIGN_P] = q_zone(QPRES) - ref_term * dq_lim[QPRES];
    if (idir == 0) {
      q_ref_right[IEIGN_BT] = q_zone(QMAGY) - ref_term * dq_lim[IEIGN_BT];
      q_ref_right[IEIGN_BTT] = q_zone(QMAGZ) - ref_term * dq_lim[IEIGN_BTT];

    } else if (idir == 1) {
      q_ref_right[IEIGN_BT] = q_zone(QMAGX) - ref_term * dq_lim[IEIGN_BT];
      q_ref_right[IEIGN_BTT] = q_zone(QMAGZ) - ref_term * dq_lim[IEIGN_BTT];

    } else if (idir == 2) {
      q_ref_right[IEIGN_BT] = q_zone(QMAGX) - ref_term * dq_lim[IEIGN_BT];
      q_ref_right[IEIGN_BTT] = q_zone(QMAGY) - ref_term * dq_lim[IEIGN_BTT];

    }

    Real summ_p[NEIGN] = {0.0_rt};
    Real summ_m[NEIGN] = {0.0_rt};

    for (int ii = 0; ii < NEIGN; ii++) {
      Real Ldq = 0.0;

      for (int n = 0; n < NEIGN; n++) {
        Ldq += leig(ii,n) * dq_lim[n];
      }

      for (int n = 0; n < NEIGN; n++) {
        summ_p[n] += 0.5_rt * dtdx * (lam(NEIGN-1) - lam(ii)) * Ldq * reig(n,ii);
        summ_m[n] += 0.5_rt * dtdx * (lam(0) - lam(ii)) * Ldq * reig(n,ii);
      }
    }

    // left state at i+1/2

    if (idir == 0) {
      qleft(i+1,j,k,QRHO) = amrex::max(small_dens,
                                       q_ref_left[IEIGN_RHO] + summ_p[IEIGN_RHO] + 0.5_rt*dt*smhd[IEIGN_RHO]);
      qleft(i+1,j,k,QU) = q_ref_left[IEIGN_U] + summ_p[IEIGN_U] + 0.5_rt*dt*smhd[IEIGN_U];
      qleft(i+1,j,k,QV) = q_ref_left[IEIGN_V] + summ_p[IEIGN_V] + 0.5_rt*dt*smhd[IEIGN_V];
      qleft(i+1,j,k,QW) = q_ref_left[IEIGN_W] + summ_p[IEIGN_W] + 0.5_rt*dt*smhd[IEIGN_W];
      qleft(i+1,j,k,QPRES) = amrex::max(small_pres,
                                        q_ref_left[IEIGN_P] + summ_p[IEIGN_P] + 0.5_rt*dt*smhd[IEIGN_P]);

      qleft(i+1,j,k,QMAGX) = Bx(i+1,j,k); // Bx stuff
      qleft(i+1,j,k,QMAGY) = q_ref_left[IEIGN_BT] + summ_p[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qleft(i+1,j,k,QMAGZ) = q_ref_left[IEIGN_BTT] + summ_p[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];

    } else if (idir == 1) {
      qleft(i,j+1,k,QRHO) = amrex::max(small_dens,
                                       q_ref_left[IEIGN_RHO] + summ_p[IEIGN_RHO] + 0.5_rt*dt*smhd[IEIGN_RHO]);
      qleft(i,j+1,k,QU) = q_ref_left[IEIGN_U] + summ_p[IEIGN_U] + 0.5_rt*dt*smhd[IEIGN_U];
      qleft(i,j+1,k,QV) = q_ref_left[IEIGN_V] + summ_p[IEIGN_V] + 0.5_rt*dt*smhd[IEIGN_V];
      qleft(i,j+1,k,QW) = q_ref_left[IEIGN_W] + summ_p[IEIGN_W] + 0.5_rt*dt*smhd[IEIGN_W];
      qleft(i,j+1,k,QPRES) = amrex::max(small_pres,
                                        q_ref_left[IEIGN_P] + summ_p[IEIGN_P] + 0.5_rt*dt*smhd[IEIGN_P]);

      qleft(i,j+1,k,QMAGX) = q_ref_left[IEIGN_BT] + summ_p[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qleft(i,j+1,k,QMAGY) = By(i,j+1,k); // By stuff
      qleft(i,j+1,k,QMAGZ) = q_ref_left[IEIGN_BTT] + summ_p[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];

    } else {
      qleft(i,j,k+1,QRHO) = amrex::max(small_dens,
                                       q_ref_left[IEIGN_RHO] + summ_p[IEIGN_RHO] + 0.5_rt*dt*smhd[IEIGN_RHO]);
      qleft(i,j,k+1,QU) = q_ref_left[IEIGN_U] + summ_p[IEIGN_U] + 0.5_rt*dt*smhd[IEIGN_U];
      qleft(i,j,k+1,QV) = q_ref_left[IEIGN_V] + summ_p[IEIGN_V] + 0.5_rt*dt*smhd[IEIGN_V];
      qleft(i,j,k+1,QW) = q_ref_left[IEIGN_W] + summ_p[IEIGN_W] + 0.5_rt*dt*smhd[IEIGN_W];
      qleft(i,j,k+1,QPRES) = amrex::max(small_pres,
                                        q_ref_left[IEIGN_P] + summ_p[IEIGN_P] + 0.5_rt*dt*smhd[IEIGN_P]);

      qleft(i,j,k+1,QMAGX) = q_ref_left[IEIGN_BT] + summ_p[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qleft(i,j,k+1,QMAGY) = q_ref_left[IEIGN_BTT] + summ_p[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];
      qleft(i,j,k+1,QMAGZ) = Bz(i,j,k+1); // Bz stuff
    }

    // right state at i-1/2
    qright(i,j,k,QRHO) = amrex::max(small_dens,
                                    q_ref_right[IEIGN_RHO] + summ_m[IEIGN_RHO] + 0.5_rt*dt*smhd[IEIGN_RHO]);
    qright(i,j,k,QU) = q_ref_right[IEIGN_U] + summ_m[IEIGN_U] + 0.5_rt*dt*smhd[IEIGN_U];
    qright(i,j,k,QV) = q_ref_right[IEIGN_V] + summ_m[IEIGN_V] + 0.5_rt*dt*smhd[IEIGN_V];
    qright(i,j,k,QW) = q_ref_right[IEIGN_W] + summ_m[IEIGN_W] + 0.5_rt*dt*smhd[IEIGN_W];
    qright(i,j,k,QPRES) = amrex::max(small_pres,
                                     q_ref_right[IEIGN_P] + summ_m[IEIGN_P] + 0.5_rt*dt*smhd[IEIGN_P]);

    if (idir == 0) {
      qright(i,j,k,QMAGX) = Bx(i,j,k); // Bx stuff
      qright(i,j,k,QMAGY) = q_ref_right[IEIGN_BT] + summ_m[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qright(i,j,k,QMAGZ) = q_ref_right[IEIGN_BTT] + summ_m[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];

    } else if (idir == 1) {
      qright(i,j,k,QMAGX) = q_ref_right[IEIGN_BT] + summ_m[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qright(i,j,k,QMAGY) = By(i,j,k); // By stuff
      qright(i,j,k,QMAGZ) = q_ref_right[IEIGN_BTT] + summ_m[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];

    } else {
      qright(i,j,k,QMAGX) = q_ref_right[IEIGN_BT] + summ_m[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qright(i,j,k,QMAGY) = q_ref_right[IEIGN_BTT] + summ_m[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];
      qright(i,j,k,QMAGZ) = Bz(i,j,k); // Bz stuff
    }

    // species
    for (int n = 0; n < NumSpec; n++) {
      Real dL;
      Real dR;
      Real un;
      Real dW = 0.0;

      if (idir == 0) {
        dL = s(i,j,k,QFS+n) - s(i-1,j,k,QFS+n);
        dR = s(i+1,j,k,QFS+n) - s(i,j,k,QFS+n);
        un = s(i,j,k,QU);

      } else if (idir == 1) {
        dL = s(i,j,k,QFS+n) - s(i,j-1,k,QFS+n);
        dR = s(i,j+1,k,QFS+n) - s(i,j,k,QFS+n);
        un = s(i,j,k,QV);

      } else {
        dL = s(i,j,k,QFS+n) - s(i,j,k-1,QFS+n);
        dR = s(i,j,k+1,QFS+n) - s(i,j,k,QFS+n);
        un = s(i,j,k,QW);
      }

      slope(dW, dL, dR, flatn(i,j,k));

      if (idir == 0) {
        qleft(i+1,j,k,QFS+n) = q_zone(QFS+n) + 0.5_rt*(1.0_rt - dtdx*un) * dW;
      } else if (idir == 1) {
        qleft(i,j+1,k,QFS+n) = q_zone(QFS+n) + 0.5_rt*(1.0_rt - dtdx*un) * dW;
      } else {
        qleft(i,j,k+1,QFS+n) = q_zone(QFS+n) + 0.5_rt*(1.0_rt - dtdx*un) * dW;
      }
      qright(i,j,k,QFS+n) = q_zone(QFS+n) - 0.5_rt*(1.0_rt + dtdx*un) * dW;
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
      for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = qleft(i+1,j,k,QFX+n);
      }
      eos(eos_input_rp, eos_state);
      qleft(i+1,j,k,QREINT) = eos_state.e * eos_state.rho;

    } else if (idir == 1) {
      eos_state.rho = qleft(i,j+1,k,QRHO);
      eos_state.p = qleft(i,j+1,k,QPRES);
      eos_state.T = s(i,j,k,QTEMP); // some initial guess?
      for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = qleft(i,j+1,k,QFS+n);
      }
      for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = qleft(i,j+1,k,QFX+n);
      }
      eos(eos_input_rp, eos_state);
      qleft(i,j+1,k,QREINT) = eos_state.e * eos_state.rho;

    } else {
      eos_state.rho = qleft(i,j,k+1,QRHO);
      eos_state.p = qleft(i,j,k+1,QPRES);
      eos_state.T = s(i,j,k,QTEMP); // some initial guess?
      for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = qleft(i,j,k+1,QFS+n);
      }
      for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = qleft(i,j,k+1,QFX+n);
      }
      eos(eos_input_rp, eos_state);
      qleft(i,j,k+1,QREINT) = eos_state. e * eos_state.rho;
    }

    eos_state.rho = qright(i,j,k,QRHO);
    eos_state.p = qright(i,j,k,QPRES);
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = qright(i,j,k,QFS+n);
    }
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = qright(i,j,k,QFX+n);
    }

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

