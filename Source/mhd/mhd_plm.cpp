
void
Castro::plm(const Box& bx,
            const int idir,
            Array4<Real const> const& s,
            Array4<Real const> const& qaux,
            Array4<Real const> const& flatn,
            Array4<Real const> const& bx,
            Array4<Real const> const& by,
            Array4<Real const> const& bz,
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

    Real q_zone[NQ];
    for (int n = 0; n < NQ; n++) {
      q_zone[n] = s(i,j,k,n);
    }

    Real as = qaux(i,j,k,QC);

    Real lam[NEIGN];
    Real leig[NEIGN];
    Real reig[NEIGN];

    evals(lam, as, q_zone, idir);

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
    smhd[IEIGN_U] = q_zone[QMAGX] / q_zone[QRHO];
    smhd[IEIGN_V] = q_zone[QMAGY] / q_zone[QRHO];
    smhd[IEIGN_W] = q_zone[QMAGZ] / q_zone[QRHO];
    smhd[IEIGN_P] = q_zone[QMAGX] * q_zone[QU] +
                    q_zone[QMAGY] * q_zone[QV] +
                    q_zone[QMAGZ] * q_zone[QW];

    if (idir == 0) {
      smhd[IEIGN_BT] = q_zone[QV];
      smhd[IEIGN_BTT] = q_zone[QW];

      // cross-talk of normal magnetic field direction
      for (int n = 0; n < NEIGN; n++) {
        smhd[n] = smhd[n] * (bx(i+1,j,k) - bx(i,j,k)) / dx[idir];
      }

    } else if (idir == 1) {
      smhd[IEIGN_BT] = q_zone[QU];
      smhd[IEIGN_BTT] = q_zone[QW];

      // cross-talk of normal magnetic field direction
      for (int n = 0; n < NEIGN; n++) {
        smhd[n] = smhd[n] * (by(i,j+1,k) - by(i,j,k)) / dx[idir];
      }

    } else {
      smhd[IEIGN_BT] = q_zone[QU];
      smhd[IEIGN_BTT] = q_zone[QV];

      // cross-talk of normal magnetic field direction
      for (int n = 0; n < NEIGN; n++) {
        smhd[n] = smhd[n] * (bz(i,j,k+1) - bz(i,j,k)) / dx[idir];
      }
    }

    // Perform the characteristic projection.  Since we are using
    // Using HLLD, we sum over all eigenvalues -- see the discussion after Eq. 31
    Real summ_p[NEIGN] = {0.0_rt};
    Real summ_m[NEIGN] = {0.0_rt};

    for (int ii = 0; ii < NEIGN; ii++) {
      Real dL[NEIGN] = {0.0_rt};
      Real dR[NEIGN] = {0.0_rt};

      for (int n = 0; n < NEIGN; n++) {
        dL[n] = leig(ii,n) * dQL[n];
        dR[n] = leig(ii,n) * dQR[n];
      }

      Real dW;
      slope(dW, dL, dR, flatn(i,j,k));

      for (int n = 0; n < NEIGN; n++) {
        summ_p[n] += (1.0_rt - dtdx * lam[ii]) * dW * reig(n,ii);
        summ_m[n] -= (1.0_rt + dtdx * lam[ii]) * dW * reig(n,ii);
      }
    }

    // left state at i+1/2

    if (idir == 0) {
      qleft(i+1,j,k,QRHO) = amrex::max(small_dens,
                                       q_zone[QRHO] + 0.5_rt*summ_p[IEIGN_RHO] + 0.5_rt*dt*smhd[IEIGN_RHO]);
      qleft(i+1,j,k,QU) = q_zone[QU] + 0.5_rt*summ_p[IEIGN_U] + 0.5_rt*dt*smhd[IEIGN_U];
      qleft(i+1,j,k,QV) = q_zone[QV] + 0.5_rt*summ_p[IEIGN_V] + 0.5_rt*dt*smhd[IEIGN_V];
      qleft(i+1,j,k,QW) = q_zone[QW] + 0.5_rt*summ_p[IEIGN_W] + 0.5_rt*dt*smhd[IEIGN_W];
      qleft(i+1,j,k,QPRES) = amrex::max(small_pres,
                                        q_zone[QPRES] + 0.5_rt*summ_p[IEIGN_P] + 0.5_rt*dt*smhd[IEIGN_P]);

      qleft(i+1,j,k,QMAGX) = bx(i+1,j,k); // Bx stuff
      qleft(i+1,j,k,QMAGY) = q_zone[QMAGY] + 0.5_rt*summ_p[IEIGN_BT] + 0.5_rt*dt*smhd[IEIGN_BT];
      qleft(i+1,j,k,QMAGZ) = q_zone[QMAGZ] + 0.5_rt*summ_p[IEIGN_BTT] + 0.5_rt*dt*smhd[IEIGN_BTT];

    } else if (idir == 1) {
                qleft(i,j+1,k,QRHO) = max(small_dens, s(i,j,k,QRHO) + 0.5_rt*summ_p(IEIGN_RHO) + 0.5_rt*dt*smhd(IEIGN_RHO))
                qleft(i,j+1,k,QU) = s(i,j,k,QU) + 0.5_rt*summ_p(IEIGN_U) + 0.5_rt*dt*smhd(IEIGN_U)
                qleft(i,j+1,k,QV) = s(i,j,k,QV) + 0.5_rt*summ_p(IEIGN_V) + 0.5_rt*dt*smhd(IEIGN_V)
                qleft(i,j+1,k,QW) = s(i,j,k,QW) + 0.5_rt*summ_p(IEIGN_W) + 0.5_rt*dt*smhd(IEIGN_W)
                qleft(i,j+1,k,QPRES) = max(small_pres, s(i,j,k,QPRES) + 0.5_rt*summ_p(IEIGN_P) + 0.5_rt*dt*smhd(IEIGN_P))

                qleft(i,j+1,k,QMAGX) = s(i,j,k,QMAGX) + 0.5_rt*summ_p(IEIGN_BT) + 0.5_rt*dt*smhd(IEIGN_BT)
                qleft(i,j+1,k,QMAGY) = by(i,j+1,k) !! By stuff
                qleft(i,j+1,k,QMAGZ) = s(i,j,k,QMAGZ) + 0.5_rt*summ_p(IEIGN_BTT) + 0.5_rt*dt*smhd(IEIGN_BTT)

             else
                qleft(i,j,k+1,QRHO) = max(small_dens, s(i,j,k,QRHO) + 0.5_rt*summ_p(IEIGN_RHO) + 0.5_rt*dt*smhd(IEIGN_RHO))
                qleft(i,j,k+1,QU) = s(i,j,k,QU) + 0.5_rt*summ_p(IEIGN_U) + 0.5_rt*dt*smhd(IEIGN_U)
                qleft(i,j,k+1,QV) = s(i,j,k,QV) + 0.5_rt*summ_p(IEIGN_V) + 0.5_rt*dt*smhd(IEIGN_V)
                qleft(i,j,k+1,QW) = s(i,j,k,QW) + 0.5_rt*summ_p(IEIGN_W) + 0.5_rt*dt*smhd(IEIGN_W)
                qleft(i,j,k+1,QPRES) = max(small_pres, s(i,j,k,QPRES) + 0.5_rt*summ_p(IEIGN_P) + 0.5_rt*dt*smhd(IEIGN_P))

                qleft(i,j,k+1,QMAGX) = s(i,j,k,QMAGX) + 0.5_rt*summ_p(IEIGN_BT) + 0.5_rt*dt*smhd(IEIGN_BT)
                qleft(i,j,k+1,QMAGY) = s(i,j,k,QMAGY) + 0.5_rt*summ_p(IEIGN_BTT) + 0.5_rt*dt*smhd(IEIGN_BTT)
                qleft(i,j,k+1,QMAGZ) = bz(i,j,k+1) !! Bz stuff
             end if

             ! right state at i-1/2
             qright(i,j,k,QRHO) = max(small_dens, s(i,j,k,QRHO) + 0.5_rt*summ_m(IEIGN_RHO) + 0.5_rt*dt*smhd(IEIGN_RHO))
             qright(i,j,k,QU) = s(i,j,k,QU) + 0.5_rt*summ_m(IEIGN_U) + 0.5_rt*dt*smhd(IEIGN_U)
             qright(i,j,k,QV) = s(i,j,k,QV) + 0.5_rt*summ_m(IEIGN_V) + 0.5_rt*dt*smhd(IEIGN_V)
             qright(i,j,k,QW) = s(i,j,k,QW) + 0.5_rt*summ_m(IEIGN_W) + 0.5_rt*dt*smhd(IEIGN_W)
             qright(i,j,k,QPRES) = max(small_pres, s(i,j,k,QPRES) + 0.5_rt*summ_m(IEIGN_P) + 0.5_rt*dt*smhd(IEIGN_P))

             if (idir == 1) then
                qright(i,j,k,QMAGX) = bx(i,j,k) !! Bx stuff
                qright(i,j,k,QMAGY) = s(i,j,k,QMAGY) + 0.5_rt*summ_m(IEIGN_BT) + 0.5_rt*dt*smhd(IEIGN_BT)
                qright(i,j,k,QMAGZ) = s(i,j,k,QMAGZ) + 0.5_rt*summ_m(IEIGN_BTT) + 0.5_rt*dt*smhd(IEIGN_BTT)

             else if (idir == 2) then
                qright(i,j,k,QMAGX) = s(i,j,k,QMAGX) + 0.5_rt*summ_m(IEIGN_BT) + 0.5_rt*dt*smhd(IEIGN_BT)
                qright(i,j,k,QMAGY) = by(i,j,k) !! By stuff
                qright(i,j,k,QMAGZ) = s(i,j,k,QMAGZ) + 0.5_rt*summ_m(IEIGN_BTT) + 0.5_rt*dt*smhd(IEIGN_BTT)

             else
                qright(i,j,k,QMAGX) = s(i,j,k,QMAGX) + 0.5_rt*summ_m(IEIGN_BT) + 0.5_rt*dt*smhd(IEIGN_BT)
                qright(i,j,k,QMAGY) = s(i,j,k,QMAGY) + 0.5_rt*summ_m(IEIGN_BTT) + 0.5_rt*dt*smhd(IEIGN_BTT)
                qright(i,j,k,QMAGZ) = bz(i,j,k) !! Bz stuff
             endif

             ! species
             do ii = QFS, QFS+nspec-1
                if (idir == 1) then
                   dL = s(i,j,k,ii) - s(i-1,j,k,ii)
                   dR = s(i+1,j,k,ii) - s(i,j,k,ii)
                   un = s(i,j,k,QU)

                else if (idir == 2) then
                   dL = s(i,j,k,ii) - s(i,j-1,k,ii)
                   dR = s(i,j+1,k,ii) - s(i,j,k,ii)
                   un = s(i,j,k,QV)

                else
                   dL = s(i,j,k,ii) - s(i,j,k-1,ii)
                   dR = s(i,j,k+1,ii) - s(i,j,k,ii)
                   un = s(i,j,k,QW)
                end if

               call slope(dW, dL, dR, flatn(i,j,k))

               if (idir == 1) then
                  qleft(i+1,j,k,ii) = s(i,j,k,ii) + 0.5_rt*(1.0d0 - dtdx*un) * dW
               else if (idir == 2) then
                  qleft(i,j+1,k,ii) = s(i,j,k,ii) + 0.5_rt*(1.0d0 - dtdx*un) * dW
               else
                  qleft(i,j,k+1,ii) = s(i,j,k,ii) + 0.5_rt*(1.0d0 - dtdx*un) * dW
               endif
               qright(i,j,k,ii) = s(i,j,k,ii) - 0.5_rt*(1.0d0 + dtdx*un) * dW
             enddo

             ! rho e
             if (idir == 1) then
                eos_state % rho = qleft(i+1,j,k,QRHO)
                eos_state % p   = qleft(i+1,j,k,QPRES)
                eos_state % T   = s(i,j,k,QTEMP) !some initial guess?
                eos_state % xn  = qleft(i+1,j,k,QFS:QFS+nspec-1)

                call eos(eos_input_rp, eos_state)
                qleft(i+1,j,k,QREINT) = eos_state % e * eos_state % rho

             else if (idir == 2) then
                eos_state % rho = qleft(i,j+1,k,QRHO)
                eos_state % p   = qleft(i,j+1,k,QPRES)
                eos_state % T   = s(i,j,k,QTEMP) !some initial guess?
                eos_state % xn  = qleft(i,j+1,k,QFS:QFS+nspec-1)

                call eos(eos_input_rp, eos_state)
                qleft(i,j+1,k,QREINT) = eos_state % e * eos_state % rho

             else
                eos_state % rho = qleft(i,j,k+1,QRHO)
                eos_state % p   = qleft(i,j,k+1,QPRES)
                eos_state % T   = s(i,j,k,QTEMP) !some initial guess?
                eos_state % xn  = qleft(i,j,k+1,QFS:QFS+nspec-1)

                call eos(eos_input_rp, eos_state)
                qleft(i,j,k+1,QREINT) = eos_state % e * eos_state % rho

             end if

             eos_state % rho = qright(i,j,k,QRHO)
             eos_state % p   = qright(i,j,k,QPRES)
             eos_state % xn  = qright(i,j,k,QFS:QFS+nspec-1)

             call eos(eos_input_rp, eos_state)
             qright(i,j,k,QREINT) = eos_state % e * eos_state % rho

             ! add source terms
             if (idir == 1) then
                qleft(i+1,j,k,QRHO) = max(small_dens, qleft(i+1,j,k,QRHO) + 0.5_rt*dt*srcQ(i,j,k,QRHO))
                qleft(i+1,j,k,QU) = qleft(i+1,j,k,QU) + 0.5_rt*dt*srcQ(i,j,k,QU)
                qleft(i+1,j,k,QV) = qleft(i+1,j,k,QV) + 0.5_rt*dt*srcQ(i,j,k,QV)
                qleft(i+1,j,k,QW) = qleft(i+1,j,k,QW) + 0.5_rt*dt*srcQ(i,j,k,QW)
                qleft(i+1,j,k,QPRES) = qleft(i+1,j,k,QPRES) + 0.5_rt*dt*srcQ(i,j,k,QPRES)
                qleft(i+1,j,k,QREINT) = qleft(i+1,j,k,QREINT) + 0.5_rt*dt*srcQ(i,j,k,QREINT)

             else if (idir == 2) then
                qleft(i,j+1,k,QRHO) = max(small_dens, qleft(i,j+1,k,QRHO) + 0.5_rt*dt*srcQ(i,j,k,QRHO))
                qleft(i,j+1,k,QU) = qleft(i,j+1,k,QU) + 0.5_rt*dt*srcQ(i,j,k,QU)
                qleft(i,j+1,k,QV) = qleft(i,j+1,k,QV) + 0.5_rt*dt*srcQ(i,j,k,QV)
                qleft(i,j+1,k,QW) = qleft(i,j+1,k,QW) + 0.5_rt*dt*srcQ(i,j,k,QW)
                qleft(i,j+1,k,QPRES) = qleft(i,j+1,k,QPRES) + 0.5_rt*dt*srcQ(i,j,k,QPRES)
                qleft(i,j+1,k,QREINT) = qleft(i,j+1,k,QREINT) + 0.5_rt*dt*srcQ(i,j,k,QREINT)

             else
                qleft(i,j,k+1,QRHO) = max(small_dens, qleft(i,j,k+1,QRHO) + 0.5_rt*dt*srcQ(i,j,k,QRHO))
                qleft(i,j,k+1,QU) = qleft(i,j,k+1,QU) + 0.5_rt*dt*srcQ(i,j,k,QU)
                qleft(i,j,k+1,QV) = qleft(i,j,k+1,QV) + 0.5_rt*dt*srcQ(i,j,k,QV)
                qleft(i,j,k+1,QW) = qleft(i,j,k+1,QW) + 0.5_rt*dt*srcQ(i,j,k,QW)
                qleft(i,j,k+1,QPRES) = qleft(i,j,k+1,QPRES) + 0.5_rt*dt*srcQ(i,j,k,QPRES)
                qleft(i,j,k+1,QREINT) = qleft(i,j,k+1,QREINT) + 0.5_rt*dt*srcQ(i,j,k,QREINT)
             end if

             qright(i,j,k,QRHO) = max(small_dens, qright(i,j,k,QRHO) + 0.5_rt*dt*srcQ(i,j,k,QRHO))
             qright(i,j,k,QU) = qright(i,j,k,QU) + 0.5_rt*dt*srcQ(i,j,k,QU)
             qright(i,j,k,QV) = qright(i,j,k,QV) + 0.5_rt*dt*srcQ(i,j,k,QV)
             qright(i,j,k,QW) = qright(i,j,k,QW) + 0.5_rt*dt*srcQ(i,j,k,QW)
             qright(i,j,k,QPRES) = qright(i,j,k,QPRES) + 0.5_rt*dt*srcQ(i,j,k,QPRES)
             qright(i,j,k,QREINT) = qright(i,j,k,QREINT) + 0.5_rt*dt*srcQ(i,j,k,QREINT)

          enddo
       enddo
    enddo

  end subroutine plm

  !======================================== Minmod TVD slope limiter =========================================
  subroutine minmod(dW, WR, WL)

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    real(rt), intent(in) :: WR, WL
    real(rt), intent(out) :: dW
    dW = 0.d0

    if (abs(WR) < abs(WL) .and. WR*WL > 0.d0) then
       dW = WR
    else if (abs(WL) < abs(WR) .and. WR*WL > 0.d0) then
       dW = WL
    endif

  end subroutine minmod


  !========================================= VanLeer TVD slope limiter =======================================
  subroutine vanleer(dW, WR, WL)

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    real(rt), intent(in )       ::  WR, WL
    real(rt), intent(out)       ::  dW
    dW = 0.0d0

    if( WR*WL .gt. 0.0d0 ) then
       dW = 2.0d0*WR*WL/(WR + WL)
    endif

  end subroutine vanleer

  !========================================== centered difference ===========================================
  subroutine centerdif(dW, WR, WL)

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    real(rt), intent(in )    :: WR, WL
    real(rt), intent(out)    :: dW

    dW = (WR+WL)/2.0d0

  end subroutine centerdif

  !================================== second order MC ==============================
  subroutine secondMC(dW, WR, WL)

    use amrex_fort_module, only : rt => amrex_real
    use amrex_constants_module, only : ZERO, HALF, ONE, TWO

    implicit none

    real(rt), intent(in  )    :: WR, WL
    real(rt), intent(out )    :: dW

    real(rt)  :: dlim

    if (WR * WL .ge. ZERO) then
       dlim = TWO * min(abs(WR), abs(WL))
    else
       dlim = ZERO
    endif

    dW = min(HALF * abs(WR + WL), dlim ) * sign(ONE, WR + WL)


  end subroutine secondMC



  !================================================================
  subroutine slope(dW, WR, WL, flat)
    use amrex_fort_module, only : rt => amrex_real

    implicit none
    real(rt), intent(in )    :: flat
    real(rt), intent(in )    :: WR, WL
    real(rt), intent(out)    :: dW

    if  (mhd_plm_slope == 0)  then
       dW = 0.0
    elseif (mhd_plm_slope == 1) then
       call vanleer(dW,WR,WL)
    elseif (mhd_plm_slope == 2) then
       call centerdif(dW,WR,WL)
    elseif (mhd_plm_slope == 3) then
       call secondMC(dW,WR,WL)
    endif

    if (use_flattening == 1) then
        dW = flat * dW
    endif

  end subroutine slope

  !=========================================== Evals =========================================================

  subroutine evals(lam, qaux, Q, dir)

    use amrex_fort_module, only : rt => amrex_real
    use network, only: nspec
    implicit none

    real(rt), intent(in)  :: Q(NQ), qaux(NQAUX)
    real(rt), intent(out) :: lam(NEIGN) !7 waves
    integer, intent(in)   :: dir !Choose direction, 1 for x, 2 for y, 3 for z

    !The characteristic speeds of the system
    real(rt)  :: cfx, cfy, cfz, cax, cay, caz, csx, csy, csz, ca, as

    ! speeds
    as = qaux(QC) * qaux(QC)

    ! Alfven
    ca  = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
    cax = (Q(QMAGX)**2)/Q(QRHO)
    cay = (Q(QMAGY)**2)/Q(QRHO)
    caz = (Q(QMAGZ)**2)/Q(QRHO)

    ! Slow
    csx = 0.5_rt*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cax))
    csy = 0.5_rt*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cay))
    csz = 0.5_rt*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*caz))

    ! Fast
    cfx = 0.5_rt*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cax))
    cfy = 0.5_rt*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cay))
    cfz = 0.5_rt*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*caz))

    if(dir.eq.1) then

       !Ax eigenvalues
       lam(1) = Q(QU) - sqrt(cfx)
       lam(2) = Q(QU) - sqrt(cax)
       lam(3) = Q(QU) - sqrt(csx)
       lam(4) = Q(QU)
       lam(5) = Q(QU) + sqrt(csx)
       lam(6) = Q(QU) + sqrt(cax)
       lam(7) = Q(QU) + sqrt(cfx)

    elseif(dir.eq.2) then

       !Ay eigenvalues
       lam(1) = Q(QV) - sqrt(cfy)
       lam(2) = Q(QV) - sqrt(cay)
       lam(3) = Q(QV) - sqrt(csy)
       lam(4) = Q(QV)
       lam(5) = Q(QV) + sqrt(csy)
       lam(6) = Q(QV) + sqrt(cay)
       lam(7) = Q(QV) + sqrt(cfy)

    else

       !Az eigenvalues
       lam(1) = Q(QW) - sqrt(cfz)
       lam(2) = Q(QW) - sqrt(caz)
       lam(3) = Q(QW) - sqrt(csz)
       lam(4) = Q(QW)
       lam(5) = Q(QW) + sqrt(csz)
       lam(6) = Q(QW) + sqrt(caz)
       lam(7) = Q(QW) + sqrt(cfz)
    endif

  end subroutine evals

  !====================================== Left Eigenvectors ===============================================

  !x direction
  subroutine evecx(leig, reig, qaux, Q)
    use amrex_fort_module, only : rt => amrex_real
    use network, only : nspec

    implicit none

    !returnes Leig, where the rows are the left eigenvectors of the characteristic matrix Ax
    real(rt), intent(in)  ::Q(NQ), qaux(NQAUX)
    real(rt), intent(out) :: leig(NEIGN,NEIGN), reig(NEIGN,NEIGN)

    !The characteristic speeds of the system
    real(rt) :: cfx, cax, csx, ca, as, S, N
    real(rt) :: cff, css, Qf, Qs, AAf, AAs, alf, als, bety, betz

    ! Speeds

    as = qaux(QC) * qaux(QC)

    ! Alfven
    ca = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
    cax = (Q(QMAGX)**2)/Q(QRHO)

    ! Slow
    csx = 0.5_rt*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cax))

    ! Fast
    cfx = 0.5_rt*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cax))

    !useful constants
    alf = sqrt((as - csx)/(cfx - csx))
    als = sqrt((cfx - as)/(cfx - csx))
    if(cfx - as .lt. 0.d0) als = 0.d0
    if(as - csx .lt. 0.d0) alf = 0.d0

    if(abs(Q(QMAGY)).le. 1.d-14 .and.abs(Q(QMAGZ)).le. 1.d-14) then
       bety = 1.d0/sqrt(2.d0)
       betz = bety
    else
       bety = Q(QMAGY)/(sqrt(Q(QMAGY)**2 + Q(QMAGZ)**2))
       betz = Q(QMAGZ)/(sqrt(Q(QMAGY)**2 + Q(QMAGZ)**2))
    endif

    cff = sqrt(cfx)*alf
    css = sqrt(csx)*als
    S = sign(1.0d0, Q(QMAGX))
    Qf = sqrt(cfx)*alf*S
    Qs = sqrt(csx)*als*S
    N = 0.5_rt/as
    AAf = sqrt(as)*alf*sqrt(Q(QRHO))
    AAs = sqrt(as)*als*sqrt(Q(QRHO))


    leig(1,:) = (/0.d0, -N*Cff, N*Qs*bety, N*Qs*betz, N*alf/Q(QRHO), N*AAs*bety/Q(QRHO), N*AAs*betz/Q(QRHO)/) !u - cf
    leig(2,:) = (/0.d0,  0.d0, -0.5_rt*betz, 0.5_rt*bety, 0.d0, -0.5_rt*S*betz/(sqrt(Q(QRHO))), 0.5_rt*bety*S/(sqrt(Q(QRHO)))/) !u - cAx
    leig(3,:) = (/0.d0, -N*Css, -N*Qf*bety, -N*Qf*betz, N*als/Q(QRHO), -N*AAf*bety/Q(QRHO), -N*AAf*betz/Q(QRHO)/) !u - cs
    leig(4,:) = (/1.d0,  0.d0,  0.d0, 0.d0, -1.d0/as, 0.d0, 0.d0/) !u
    leig(5,:) = (/0.d0,  N*Css, N*Qf*bety, N*Qf*betz, N*als/Q(QRHO), -N*AAf*bety/Q(QRHO), -N*AAf*betz/Q(QRHO)/) !u + cs
    leig(6,:) = (/0.d0,  0.d0, 0.5_rt*betz, -0.5_rt*bety, 0.d0, -0.5_rt*betz*S/(sqrt(Q(QRHO))), 0.5_rt*bety*S/(sqrt(Q(QRHO)))/) !u + cAx
    leig(7,:) = (/0.d0, N*Cff, -N*Qs*bety, -N*Qs*betz, N*alf/Q(QRHO), N*AAs*bety/Q(QRHO), N*AAs*betz/Q(QRHO)/) !u + cf

    !   u - cf       u - Cax      u - cs   u    u + cs   u + Cax   u + cf
    reig(IEIGN_RHO,:) = (/Q(QRHO)*alf, 0.d0, Q(QRHO)*als, 1.d0, Q(QRHO)*als, 0.d0, Q(QRHO)*alf/)
    reig(IEIGN_U,:) = (/-cff , 0.d0, -css, 0.d0, css, 0.d0, cff/)
    reig(IEIGN_V,:) = (/Qs*bety, -betz, -Qf*bety, 0.d0, Qf*bety, betz, -Qs*bety/)
    reig(IEIGN_W,:) = (/Qs*betz, bety, -Qf*betz, 0.d0, Qf*betz, -bety, -Qs*betz/)
    reig(IEIGN_P,:) = (/Q(QRHO)*as*alf, 0.d0, Q(QRHO)*as*als, 0.d0  , Q(QRHO)*as*als, 0.d0, Q(QRHO)*as*alf/)
    reig(IEIGN_BT,:) = (/AAs*bety, -betz*S*sqrt(Q(QRHO)), -AAf*bety, 0.d0  , -AAf*bety, -betz*S*sqrt(Q(QRHO)), AAs*bety/)
    reig(IEIGN_BTT,:) = (/AAs*betz, bety*S*sqrt(Q(QRHO)), -AAf*betz, 0.d0, -AAf*betz, bety*S*sqrt(Q(QRHO)), AAs*betz/)


  end subroutine evecx

  !y direction
  subroutine evecy(leig, reig, qaux, Q)
    use amrex_fort_module, only : rt => amrex_real
    use network, only : nspec

    implicit none

    !returnes Leig, where the rows are the left eigenvectors of the characteristic matrix Ay
    real(rt), intent(in) :: Q(NQ), qaux(NQAUX)
    real(rt), intent(out) :: leig(NEIGN,NEIGN), reig(NEIGN,NEIGN)

    !The characteristic speeds of the system
    real(rt) :: cfy, cay, csy, ca, as, S, N
    real(rt) :: cff, css, Qf, Qs, AAf, AAs, alf, als, betx, betz


    ! Speeds
    as = qaux(QC) * qaux(QC)

    ! Alfven
    ca = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
    cay = (Q(QMAGY)**2)/Q(QRHO)

    ! Slow
    csy = 0.5_rt*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cay))

    ! Fast
    cfy = 0.5_rt*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cay))

    !useful constants
    alf = sqrt((as - csy)/(cfy - csy))
    als = sqrt((cfy - as)/(cfy - csy))
    if(as - csy .lt. 0.d0) alf = 0.d0
    if(cfy - as .lt. 0.d0) als = 0.d0

    if(abs(Q(QMAGX)).le. 1.d-14 .and.abs(Q(QMAGZ)).le. 1.d-14) then
       betx = 1.d0/sqrt(2.d0)
       betz = betx
    else
       betx = Q(QMAGX)/(sqrt(Q(QMAGX)**2 + Q(QMAGZ)**2))
       betz = Q(QMAGZ)/(sqrt(Q(QMAGX)**2 + Q(QMAGZ)**2))
    endif

    cff = sqrt(cfy)*alf
    css = sqrt(csy)*als
    S = sign(1.0d0, Q(QMAGY))
    Qf = sqrt(cfy)*alf*S
    Qs = sqrt(csy)*als*S
    AAf = sqrt(as)*alf*sqrt(Q(QRHO))
    AAs = sqrt(as)*als*sqrt(Q(QRHO))
    N = 0.5_rt/as

    !Need to double check the rows
    leig(1,:) = (/0.d0, N*Qs*betx, -N*Cff , N*Qs*betz, N*alf/Q(QRHO), N*AAs*betx/Q(QRHO), N*AAs*betz/Q(QRHO)/) ! v - cf
    leig(2,:) = (/0.d0, -0.5_rt*betz, 0.d0, 0.5_rt*betx, 0.d0, -0.5_rt*betz*S/(sqrt(Q(QRHO))), 0.5_rt*betx*S/(sqrt(Q(QRHO)))/) ! v - cAy
    leig(3,:) = (/0.d0, -N*Qf*betx, -N*Css, -N*Qf*betz, N*als/Q(QRHO), -N*AAf*betx/Q(QRHO), -N*AAf*betz/Q(QRHO)/) ! v - cs
    leig(4,:) = (/1.d0,  0.d0, 0.d0, 0.d0, -1.d0/as, 0.d0, 0.d0/) ! v
    leig(5,:) = (/0.d0, N*Qf*betx, N*Css, N*Qf*betz, N*als/Q(QRHO), -N*AAf*betx/Q(QRHO), -N*AAf*betz/Q(QRHO)/) ! v + cs
    leig(6,:) = (/0.d0, 0.5_rt*betz, 0.d0,  -0.5_rt*betx, 0.d0, -0.5_rt*betz*S/(sqrt(Q(QRHO))), 0.5_rt*betx*S/(sqrt(Q(QRHO)))/) ! v + cAy
    leig(7,:) = (/0.d0, -N*Qs*betx, N*Cff, -N*Qs*betz, N*alf/Q(QRHO), N*AAs*betx/Q(QRHO), N*AAs*betz/Q(QRHO)/) ! v + cf

    !   v - cf   v - Cay   v - cs   v   v + cs   v + Cay   v + cf
    reig(IEIGN_RHO,:) = (/Q(QRHO)*alf, 0.d0, Q(QRHO)*als, 1.d0, Q(QRHO)*als, 0.d0, Q(QRHO)*alf/)
    reig(IEIGN_V,:) = (/-cff, 0.d0, -css, 0.d0  , css, 0.d0, cff/)
    reig(IEIGN_U,:) = (/Qs*betx, -betz, -Qf*betx, 0.d0  , Qf*betx, betz, -Qs*betx/)
    reig(IEIGN_W,:) = (/Qs*betz, betx, -Qf*betz, 0.d0  , Qf*betz, -betx , -Qs*betz/)
    reig(IEIGN_P,:) = (/Q(QRHO)*as*alf, 0.d0, Q(QRHO)*as*als, 0.d0  , Q(QRHO)*as*als, 0.d0, Q(QRHO)*as*alf/)
    reig(IEIGN_BT,:) = (/AAs*betx, -betz*S*sqrt(Q(QRHO)), -AAf*betx, 0.d0  , -AAf*betx, -betz*S*sqrt(Q(QRHO)) , AAs*betx/)
    reig(IEIGN_BTT,:) = (/AAs*betz, betx*S*sqrt(Q(QRHO)), -AAf*betz, 0.d0, -AAf*betz, betx*S*sqrt(Q(QRHO)), AAs*betz/)

  end subroutine evecy

  !z direction
  subroutine evecz(leig, reig, qaux, Q)
    use amrex_fort_module, only : rt => amrex_real
    use network, only : nspec

    implicit none

    !returnes Leig, where the rows are the left eigenvectors of the characteristic matrix Az
    real(rt), intent(in)  :: Q(NQ), qaux(NQAUX)
    real(rt), intent(out) :: leig(NEIGN,NEIGN), reig(NEIGN,NEIGN)

    !The characteristic speeds of the system
    real(rt)  :: cfz, caz, csz, ca, as, S, N
    real(rt)  :: cff, css, Qf, Qs, AAf, AAs, alf, als, betx, bety

    ! Speeds
    as = qaux(QC) * qaux(QC)

    ! Alfven
    ca = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
    caz = (Q(QMAGZ)**2)/Q(QRHO)

    ! Slow
    csz = 0.5_rt*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*caz))

    ! Fast
    cfz = 0.5_rt*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*caz))

    !useful constants
    alf = sqrt((as - csz)/(cfz - csz))
    als = sqrt((cfz - as)/(cfz - csz))

    if(cfz - as .lt. 0.d0) als = 0.d0

    if(abs(Q(QMAGX)).le. 1.d-14 .and.abs(Q(QMAGY)).le. 1.d-14) then
       betx = 1.d0/sqrt(2.d0)
       bety = betx
    else
       betx = Q(QMAGX)/(sqrt(Q(QMAGX)**2 + Q(QMAGY)**2))
       bety = Q(QMAGY)/(sqrt(Q(QMAGX)**2 + Q(QMAGY)**2))
    endif

    cff = sqrt(cfz)*alf
    css = sqrt(csz)*als
    S = sign(1.0d0, Q(QMAGZ))
    Qf = sqrt(cfz)*alf*S
    Qs = sqrt(csz)*als*S
    AAf = sqrt(as)*alf*sqrt(Q(QRHO))
    AAs = sqrt(as)*als*sqrt(Q(QRHO))
    N = 0.5_rt/as

    !Need to double check the order
    leig(1,:) = (/0.d0, N*Qs*betx, N*Qs*bety, -N*Cff, N*alf/Q(QRHO) , N*AAs*betx/Q(QRHO), N*AAs*bety/Q(QRHO)/) !w - cf
    leig(2,:) = (/0.d0, -0.5_rt*bety, 0.5_rt*betx, 0.d0, 0.d0, -0.5_rt*S*bety/(sqrt(Q(QRHO))) , 0.5_rt*betx*S/(sqrt(Q(QRHO)))/) !w - cAz
    leig(3,:) = (/0.d0, -N*Qf*betx, -N*Qf*bety, -N*Css, N*als/Q(QRHO) , -N*AAf*betx/Q(QRHO), -N*AAf*bety/Q(QRHO)/) !w - cs
    leig(4,:) = (/1.d0,  0.d0 ,  0.d0, 0.d0, -1.d0/as, 0.d0, 0.d0/) !w
    leig(5,:) = (/0.d0, N*Qf*betx, N*Qf*bety, N*Css, N*als/Q(QRHO) , -N*AAf*betx/Q(QRHO), -N*AAf*bety/Q(QRHO)/) !w + cs
    leig(6,:) = (/0.d0, 0.5_rt*bety , -0.5_rt*betx, 0.0d0, 0.d0 , -0.5_rt*bety*S/(sqrt(Q(QRHO))) , 0.5_rt*betx*S/(sqrt(Q(QRHO)))  /) !w + cAz
    leig(7,:) = (/0.d0, -N*Qs*betx, -N*Qs*bety, N*Cff , N*alf/Q(QRHO) , N*AAs*betx/Q(QRHO) , N*AAs*bety/Q(QRHO) /) !w + cf

    !   w - cf    w - Caz     w - cs     w    w + cs    w + Caz     w + cf
    reig(IEIGN_RHO,:) = (/Q(QRHO)*alf, 0.d0, Q(QRHO)*als, 1.d0, Q(QRHO)*als, 0.d0, Q(QRHO)*alf/)
    reig(IEIGN_W,:) = (/-cff , 0.d0, -css, 0.d0, css, 0.d0 , cff/)
    reig(IEIGN_U,:) = (/Qs*betx, -bety, -Qf*betx, 0.d0, Qf*betx, bety, -Qs*betx/)
    reig(IEIGN_V,:) = (/Qs*bety, betx, -Qf*bety, 0.d0, Qf*bety, -betx , -Qs*bety/)
    reig(IEIGN_P,:) = (/Q(QRHO)*as*alf, 0.d0, Q(QRHO)*as*als, 0.d0, Q(QRHO)*as*als, 0.d0, Q(QRHO)*as*alf/)
    reig(IEIGN_BT,:) = (/AAs*betx, -bety*S*sqrt(Q(QRHO)), -AAf*betx, 0.d0, -AAf*betx, -bety*S*sqrt(Q(QRHO)), AAs*betx/)
    reig(IEIGN_BTT,:) = (/AAs*bety, betx*S*sqrt(Q(QRHO)), -AAf*bety, 0.d0, -AAf*bety, betx*S*sqrt(Q(QRHO)), AAs*bety/)

  end subroutine evecz

end module mhd_plm_module
