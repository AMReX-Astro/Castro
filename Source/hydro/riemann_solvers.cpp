void
Castro::riemanncg(const Box& bx,
                  Array4<Real> const ql,
                  Array4<Real> const qr,
                  Array4<Real const> const qaux,
                  Array4<Real> const qint,
                  const int idir,
                  const int* domlo, const int* domhi) {

  // this implements the approximate Riemann solver of Colella & Glaz
  // (1985)
  //

  constexpr Real weakwv = 1.e-3_rt;

#ifndef AMREX_USE_CUDA
  GpuArray<Real, iter_max> pstar_hist;
  GpuArray<Real, 2*iter_max> pstar_hist_extra;
#endif

#ifndef AMREX_USE_CUDA
  if (cg_blend == 2 && cg_maxiter < 5) {
    amrex::Error("Error: need cg_maxiter >= 5 to do a bisection search on secant iteration failure.");
  }
#endif

  int iu, iv1, iv2;
  int im1, im2, im3;
  int sx, sy, sz;

  if (idir == 0) {
    iu = QU;
    iv1 = QV;
    iv2 = QW;
    im1 = UMX;
    im2 = UMY;
    im3 = UMZ;
    sx = 1;
    sy = 0;
    sz = 0;

  } else if (idir == 1) {
    iu = QV;
    iv1 = QU;
    iv2 = QW;
    im1 = UMY;
    im2 = UMX;
    im3 = UMZ;
    sx = 0;
    sy = 1;
    sz = 0;

  } else {
    iu = QW;
    iv1 = QU;
    iv2 = QV;
    im1 = UMZ;
    im2 = UMX;
    im3 = UMY;
    sx = 0;
    sy = 0;
    sz = 1;
  }

  // do we want to force the flux to zero at the boundary?
  bool special_bnd_lo = (physbc_lo[idir] == Symmetry ||
                         physbc_lo[idir] == SlipWall ||
                         physbc_lo(idir) == NoSlipWall);
  bool special_bnd_hi = (physbc_hi(idir) == Symmetry ||
                         physbc_hi(idir) == SlipWall ||
                         physbc_hi(idir) == NoSlipWall);


  AMREX_PARALLEL_FOR_3D(bx, i, j, k,
  {

    // deal with hard walls
    bnd_fac = 1.0_rt;

    if (idir == 0) {
      if ((i == domlo[0] && special_bnd_lo) ||
          (i == domhi[0]+1 && special_bnd_hi)) {
        bnd_fac = 0.0_rt;
      }

    } else if (idir == 1) {
      if ((j == domlo[1] && special_bnd_lo) ||
          (j == domhi[1]+1 && special_bnd_hi)) {
        bnd_fac = 0.0_rt;
      }
    } else {
      if ((k == domlo[2] && special_bnd_lo) ||
          (k == domhi[2]+1 && special_bnd_hi)) {
        bnd_fac = 0.0_rt;
      }
    }


    // left state
    Real rl = amrex::max(ql(i,j,k,QRHO), small_dens);

    Real pl = ql(i,j,k,QPRES);
    Real rel = ql(i,j,k,QREINT);
    Real gcl = qaux(i-sx,j-sy,k-sz,QGAMC);
    if (use_reconstructed_gamma1 == 1) {
      gcl = ql(i,j,k,QGC);
    }

    // pick left velocities based on direction
    Real ul = ql(i,j,k,iu);
    Real v1l = ql(i,j,k,iv1);
    Real v2l = ql(i,j,k,iv2);

    // sometime we come in here with negative energy or pressure
    // note: reset both in either case, to remain thermo
    // consistent
    if (rel <= ZERO || pl < small_pres) {
#ifndef AMREX_USE_CUDA
      std::cout <<  "WARNING: (rho e)_l < 0 or pl < small_pres in Riemann: " << rel << " " pl << " " small_pres << std::endl;
#endif

      eos_t eos_state;
      eos_state.T = small_temp;
      eos_state.rho = rl;
      for (int n = 0; n < nspec; n++) {
        eos_state.xn[n] = ql(i,j,k,QFS+n);
      }
      for (int n = 0; n < naux; n++) {
        eos_state.aux[n] = ql(i,j,k,QFX+n);
      }

      eos(eos_input_rt, eos_state);

      rel = rl*eos_state.e;
      pl = eos_state.p;
      gcl = eos_state.gam1;
    }

    // right state
    Real rr = amrex::max(qr(i,j,k,QRHO), small_dens);

    Real pr = qr(i,j,k,QPRES);
    Real rer = qr(i,j,k,QREINT);
    Real gcr = qaux(i,j,k,QGAMC);
    if (use_reconstructed_gamma1 == 1) {
      gcr = qr(i,j,k,QGC);
    }

    // pick right velocities based on direction
    Real ur = qr(i,j,k,iu);
    Real v1r = qr(i,j,k,iv1);
    Real v2r = qr(i,j,k,iv2);

    if (rer <= ZERO || pr < small_pres) {
#ifndef AMREX_USE_CUDA
      std::cout << "WARNING: (rho e)_r < 0 or pr < small_pres in Riemann: " << rer << " " << pr << " " << small_pres << std::endl;
#endif

      eos_state.T = small_temp;
      eos_state.rho = rr;
      for (int n = 0; n < nspec; n++) {
        eos_state.xn[n] = qr(i,j,k,QFS+n);
      }
      for (int n = 0; n < naux; n++) {
        eos_state.aux[n] = qr(i,j,k,QFX+n);
      }

      eos(eos_input_rt, eos_state);

      rer = rr*eos_state.e;
      pr = eos_state.p;
      gcr = eos_state.gam1;
    }

    // common quantities
    Real taul = 1.0_rt/rl;
    Real taur = 1.0_rt/rr;

    // lagrangian sound speeds
    Real clsql = gcl*pl*rl;
    Real clsqr = gcr*pr*rr;

    Real csmall = amrex::max(small, amrex::max(small * qaux(i,j,k,QC),
                                               small * qaux(i-sx,j-sy,k-sz,QC)));

    Real cavg = 0.5_rt*(qaux(i,j,k,QC) + qaux(i-sx,j-sy,k-sz,QC));

    // Note: in the original Colella & Glaz paper, they predicted
    // gamma_e to the interfaces using a special (non-hyperbolic)
    // evolution equation.  In Castro, we instead bring (rho e)
    // to the edges, so we construct the necessary gamma_e here from
    // what we have on the interfaces.
    Real gamel = pl/rel + 1.0_rt;
    Real gamer = pr/rer + 1.0_rt;

    // these should consider a wider average of the cell-centered
    // gammas
    Real gmin = amrex::min(amrex::min(gamel, gamer), ONE);
    Real gmax = amrex::max(amrex::max(gamel, gamer), TWO);

    Real game_bar = 0.5_rt*(gamel + gamer);
    gamc_bar = 0.5_rt*(gcl + gcr);

    Real gdot = 2.0_rt*(1.0_rt - game_bar/gamc_bar)*(game_bar - 1.0_rt);

    Real wsmall = small_dens*csmall;
    Real wl = amrex::max(wsmall, std::sqrt(std::abs(clsql)));
    Real wr = amrex::max(wsmall, std::sqrt(std::abs(clsqr)));

    // make an initial guess for pstar -- this is a two-shock
    // approximation
    //pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
    Real pstar = pl + ( (pr - pl) - wr*(ur - ul) )*wl/(wl+wr);
    pstar = amrex::max(pstar, small_pres);

    // get the shock speeds -- this computes W_s from CG Eq. 34
    wsqge(pl, taul, gamel, gdot, gamstar,
          gmin, gmax, clsql, pstar, wlsq);

    wsqge(pr, taur, gamer, gdot, gamstar,
          gmin, gmax, clsqr, pstar, wrsq);

    Real pstar_old = pstar;

    Real wl = std::sqrt(wlsq);
    Real wr = std::sqrt(wrsq);

    // R-H jump conditions give ustar across each wave -- these
    // should be equal when we are done iterating.  Our notation
    // here is a little funny, comparing to CG, ustar_l = u*_L and
    // ustar_r = u*_R.
    Real ustar_l = ul - (pstar-pl)/wl;
    Real ustar_r = ur + (pstar-pr)/wr;

    // revise our pstar guess
    // pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
    pstar = pl + ( (pr - pl) - wr*(ur - ul) )*wl/(wl+wr);
    pstar = amrex::max(pstar, small_pres);

    // secant iteration
    bool converged = false;

    int iter = 0;
    while ((iter < iter_max && !converged) || iter < 2) {

      wsqge(pl, taul, gamel, gdot, gamstar,
            gmin, gmax, clsql, pstar, wlsq);

      wsqge(pr, taur, gamer, gdot, gamstar,
            gmin, gmax, clsqr, pstar, wrsq);

      // NOTE: these are really the inverses of the wave speeds!
      wl = 1.0_rt / std::sqrt(wlsq);
      wr = 1.0_rt / std::sqrt(wrsq);

      Real ustar_r_old = ustar_r;
      Real ustar_l_old = ustar_l;

      ustar_r = ur - (pr-pstar)*wr;
      ustar_l = ul + (pl-pstar)*wl;

      Real dpditer = std::abs(pstar_old-pstar);

      // Here we are going to do the Secant iteration version in
      // CG.  Note that what we call zp and zm here are not
      // actually the Z_p = |dp*/du*_p| defined in CG, by rather
      // simply |du*_p| (or something that looks like dp/Z!).
      Real zp = std::abs(ustar_l - ustar_l_old);
      if (zp - weakwv*cavg <= 0.0_rt) {
        zp = dpditer*wl;
      }

      Real zm = std::abs(ustar_r - ustar_r_old);
      if (zm - weakwv*cavg <= 0.0_rt) {
        zm = dpditer*wr;
      }

      // the new pstar is found via CG Eq. 18
      Real denom = dpditer/amrex::max(zp+zm, small*cavg);
      Real pstar_old = pstar;
      pstar = pstar - denom*(ustar_r - ustar_l);
      pstar = amrex::max(pstar, small_pres);

      Real err = std::abs(pstar - pstar_old);
      if (err < tol*pstar) {
        converged = true;
      }

#ifndef AMREX_USE_CUDA
      pstar_hist[iter] = pstar;
#endif

      iter++;
    }

    // If we failed to converge using the secant iteration, we
    // can either stop here; or, revert to the original
    // two-shock estimate for pstar; or do a bisection root
    // find using the bounds established by the most recent
    // iterations.

    if (!converged) {

      if (cg_blend == 0) {

#ifndef AMREX_USE_CUDA
        std::cout <<  "pstar history: " << std::endl;
        for (int iter=0; iter < iter_max; iter++) {
          std::cout << iter << " " << pstar_hist[iter];
        }

        std::cout << std::endl;
        std::cout << "left state  (r,u,p,re,gc): " << rl << " " << ul << " " << pl << " " << rel << " " << gcl << std::endl;
        std::cout << "right state (r,u,p,re,gc): " << rr << " " << ur << " " << pr << " " << rer << " " << gcr << std::endl;
        std::cout << "cavg, smallc: " << cavg << " " << csmall;

        amrex::Error("ERROR: non-convergence in the Riemann solver");
#endif

      } else if (cg_blend == 1) {

        pstar = pl + ( (pr - pl) - wr*(ur - ul) )*wl/(wl+wr);

      } else if (cg_blend == 2) {

        // we don't store the history if we are in CUDA, so
        // we can't do this
#ifndef AMREX_USE_CUDA
        // first try to find a reasonable bounds
        pstarl = 1.e200;
        pstaru = -1.e200;
        for (int q = iter_max-6, q < iter_max; q++) {
          pstarl = amrex::min(pstarl, pstar_hist[q]);
          pstaru = amrex::max(pstaru, pstar_hist[q]);
        }
        pstar_bisection(pstarl, pstaru,
                        ul, pl, taul, gamel, clsql,
                        ur, pr, taur, gamer, clsqr,
                        gdot, gmin, gmax,
                        pstar, gamstar, converged, pstar_hist_extra);

        if (!converged) {

          std::cout << "pstar history: ";
          for (int iter = 0; iter < iter_max; iter++) {
            std::cout << iter << " " << pstar_hist[iter];
          }
          for (int iter = 0; iter < 2*iter_max; iter++) {
            std::cout << iter << " " << pstar_hist_extra[iter];
          }

          std::cout << std::endl;
          std::cout << "left state  (r,u,p,re,gc): " << rl << " " << ul << " " << pl << " " << rel << " " << gcl << std::endl;
          std::cout << "right state (r,u,p,re,gc): " << rr << " " << ur << " " << pr << " " << rer << " " << gcr << std::endl;
          std::cout << "cavg, smallc: " << cavg << " " csmall << std::endl;

          amrex::Error("ERROR: non-convergence in the Riemann solver");
        }

#endif
      } else {

#ifndef AMREX_USE_CUDA
        amrex::Error("ERROR: unrecognized cg_blend option.");
#endif
      }

    }

    // we converged!  construct the single ustar for the region
    // between the left and right waves, using the updated wave speeds
    ustar_r = ur - (pr-pstar)*wr;  // careful -- here wl, wr are 1/W
    ustar_l = ul + (pl-pstar)*wl;

    ustar = 0.5_rt * (ustar_l + ustar_r);

    // for symmetry preservation, if ustar is really small, then we
    // set it to zero
    if (std::abs(ustar) < smallu*HALF*(std::abs(ul) + std::abs(ur))) {
      ustar = 0.0_rt;
    }

    // sample the solution -- here we look first at the direction
    // that the contact is moving.  This tells us if we need to
    // worry about the L/L* states or the R*/R states.
    Real ro;
    Real uo;
    Real po;
    Real tauo;
    Real gamco;
    Real gameo;

    if (ustar > 0.0_rt) {
      ro = rl;
      uo = ul;
      po = pl;
      tauo = taul;
      gamco = gcl;
      gameo = gamel;

    } else if (ustar < 0.0_rt) {
      ro = rr;
      uo = ur;
      po = pr;
      tauo = taur;
      gamco = gcr;
      gameo = gamer;

    } else {
      ro = 0.5_rt*(rl+rr);
      uo = 0.5_rt*(ul+ur);
      po = 0.5_rt*(pl+pr);
      tauo = 0.5_rt*(taul+taur);
      gamco = 0.5_rt*(gcl+gcr);
      gameo = 0.5_rt*(gamel + gamer);
    }

    // use tau = 1/rho as the independent variable here
    ro = amrex::max(small_dens, 1.0_rt/tauo);
    tauo = 1.0_rt/ro;

    Real co = std::sqrt(std::abs(gamco*po*tauo));
    co = amrex::max(csmall, co);
    Real clsq = std::pow(co*ro, 2);

    // now that we know which state (left or right) we need to worry
    // about, get the value of gamstar and wosq across the wave we
    // are dealing with.
    wsqge(po, tauo, gameo, gdot, gamstar,
          gmin, gmax, clsq, pstar, wosq);

    Real sgnm = std::copysign(1.0_rt, ustar);

    Real wo = std::sqrt(wosq);
    Real dpjmp = pstar - po;

    // is this max really necessary?
    //rstar=max(ONE-ro*dpjmp/wosq, (gameo-ONE)/(gameo+ONE))
    Real rstar = 1.0_rt - ro*dpjmp/wosq;
    rstar = ro/rstar;
    rstar = amrex::max(small_dens, rstar);

    Real cstar = std::sqrt(std::abs(gamco*pstar/rstar));
    cstar = amrex::max(cstar, csmall);

    Real spout = co - sgnm*uo;
    Real spin = cstar - sgnm*ustar;

    //ushock = HALF*(spin + spout)
    Real ushock = wo*tauo - sgnm*uo;

    if (pstar-po >= 0.0_rt) {
      spin = ushock;
      spout = ushock;
    }

    Real frac = 0.5_rt*(1.0_rt + (spin + spout)/amrex::max(amrex::max(spout-spin, spin+spout), small*cavg));

    // the transverse velocity states only depend on the
    // direction that the contact moves
    if (ustar > 0.0_rt) {
      qint(i,j,k,iv1) = v1l;
      qint(i,j,k,iv2) = v2l;
    } else if (ustar < 0.0_rt) {
      qint(i,j,k,iv1) = v1r;
      qint(i,j,k,iv2) = v2r;
    } else {
      qint(i,j,k,iv1) = 0.5_rt*(v1l+v1r);
      qint(i,j,k,iv2) = 0.5_rt*(v2l+v2r);
    }

    // linearly interpolate between the star and normal state -- this covers the
    // case where we are inside the rarefaction fan.
    qint(i,j,k,QRHO) = frac*rstar + (1.0_rt - frac)*ro;
    qint(i,j,k,iu) = frac*ustar + (1.0_rt - frac)*uo;
    qint(i,j,k,QPRES) = frac*pstar + (1.0_rt - frac)*po;
    game_int = frac*gamstar + (1.0_rt-frac)*gameo;

    // now handle the cases where instead we are fully in the
    // star or fully in the original (l/r) state
    if (spout < 0.0_rt) {
      qint(i,j,k,QRHO) = ro;
      qint(i,j,k,iu) = uo;
      qint(i,j,k,QPRES) = po;
      game_int = gameo;
    }

    if (spin >= 0.0_rt) {
      qint(i,j,k,QRHO) = rstar;
      qint(i,j,k,iu) = ustar;
      qint(i,j,k,QPRES) = pstar;
      game_int = gamstar;
    }

    qint(i,j,k,QPRES) = amrex::max(qint(i,j,k,QPRES), small_pres);

    qint(i,j,k,iu) = qint(i,j,k,iu) * bnd_fac;

    // Compute fluxes, order as conserved state (not q)

    // compute the total energy from the internal, p/(gamma - 1), and the kinetic
    qint(i,j,k,QREINT) = qint(i,j,k,QPRES)/(game_int - 1.0_rt);

    // advected quantities -- only the contact matters
    for (int ipassive = 0; ipassive < npassive; ipassive++) {
      int n  = upass_map[ipassive];
      int nqp = qpass_map[ipassive];

      if (ustar > 0.0_rt) {
        qint(i,j,k,nqp) = ql(i,j,k,nqp);
      } else if (ustar < 0.0_rt) {
        qint(i,j,k,nqp) = qr(i,j,k,nqp);
      } else {
        qint(i,j,k,nqp) = 0.5_rt * (ql(i,j,k,nqp) + qr(i,j,k,nqp));
      }
    }

  });
}

void
Castro::riemannus(const Box& bx,
                  Array4<Real> const ql,
                  Array4<Real> const qr,
                  Array4<Real const> const qaux,
                  Array4<Real> const qint,
#ifdef RADIATION
                  Array4<Real> const lambda_int,
#endif
                  const int idir, const int compute_gammas,
                  const int* domlo, comst int& domhi) {

  // Colella, Glaz, and Ferguson solver
  //
  // this is a 2-shock solver that uses a very simple approximation for the
  // star state, and carries an auxiliary jump condition for (rho e) to
  // deal with a real gas

  // set integer pointers for the normal and transverse velocity and
  // momentum

  int iu, iv1, iv2;
  int im1, im2, im3;

  if (idir == 0) {
    iu = QU;
    iv1 = QV;
    iv2 = QW;
    im1 = UMX;
    im2 = UMY;
    im3 = UMZ;

  } else if (idir == 1) {
    iu = QV;
    iv1 = QU;
    iv2 = QW;
    im1 = UMY;
    im2 = UMX;
    im3 = UMZ;

  } else {
    iu = QW;
    iv1 = QU;
    iv2 = QV;
    im1 = UMZ;
    im2 = UMX;
    im3 = UMY;
  }

  // do we want to force the flux to zero at the boundary?
  bool special_bnd_lo = (physbc_lo[idir] == Symmetry ||
                         physbc_lo[idir] == SlipWall ||
                         physbc_lo(idir) == NoSlipWall);
  bool special_bnd_hi = (physbc_hi(idir) == Symmetry ||
                         physbc_hi(idir) == SlipWall ||
                         physbc_hi(idir) == NoSlipWall);


  AMREX_PARALLEL_FOR_3D(bx, i, j, k,
  {

    // deal with hard walls
    bnd_fac = 1.0_rt;

    if (idir == 0) {
      if ((i == domlo[0] && special_bnd_lo) ||
          (i == domhi[0]+1 && special_bnd_hi)) {
        bnd_fac = 0.0_rt;
      }

    } else if (idir == 1) {
      if ((j == domlo[1] && special_bnd_lo) ||
          (j == domhi[1]+1 && special_bnd_hi)) {
        bnd_fac = 0.0_rt;
      }
    } else {
      if ((k == domlo[2] && special_bnd_lo) ||
          (k == domhi[2]+1 && special_bnd_hi)) {
        bnd_fac = 0.0_rt;
      }
    }


    // set the left and right states for this interface

#ifdef RADIATION
    Real laml[ngroups];
    Real lamr[ngroups];

    for (int g = 0; g < ngroups; g++) {
      if (idir == 0) {
        laml[g] = qaux(i-1,j,k,QLAMS+g);
      } else if (idir == 1) {
        laml[g] = qaux(i,j-1,k,QLAMS+g);
      } else {
        laml[g] = qaux(i,j,k-1,QLAMS+g);
      }
      lamr[g] = qaux(i,j,k,QLAMS+g);
    }
#endif

    Real rl = amrex::max(ql(i,j,k,QRHO), small_dens);

    // pick left velocities based on direction
    Real ul  = ql(i,j,k,iu);
    Real v1l = ql(i,j,k,iv1);
    Real v2l = ql(i,j,k,iv2);

#ifdef RADIATION
    Real pl = ql(i,j,k,qptot);
    Real rel = ql(i,j,k,qreitot);
    Real erl[ngroups];
    for (int g = 0; g < ngroups; g++) {
      erl[g] = ql(i,j,k,qrad+g);
    }
    Real pl_g = ql(i,j,k,QPRES);
    Real rel_g = ql(i,j,k,QREINT);
#else
    Real pl = amrex::max(ql(i,j,k,QPRES), small_pres);
    Real rel = ql(i,j,k,QREINT);
#endif

    Real rr = amrex::max(qr(i,j,k,QRHO), small_dens);

    // pick right velocities based on direction
    Real ur  = qr(i,j,k,iu);
    Real v1r = qr(i,j,k,iv1);
    Real v2r = qr(i,j,k,iv2);

#ifdef RADIATION
    Real pr = qr(i,j,k,qptot);
    Real rer = qr(i,j,k,qreitot);
    Real err[ngroups];
    for (int g = 0; g < ngroups; g++) {
      err[g] = qr(i,j,k,qrad:qrad-1+ngroups);
    }
    Real pr_g = qr(i,j,k,QPRES);
    Real rer_g = qr(i,j,k,QREINT);
#else
    Real pr = amrex::max(qr(i,j,k,QPRES), small_pres);
    Real rer = qr(i,j,k,QREINT);
#endif

    // estimate the star state: pstar, ustar

    Real csmall;
    Real cavg;
    Real gamcl;
    Real gamcr;
#ifdef RADIATION
    Real gamcgl;
    Real gamcgr;
#endif

    if (idir == 0) {
      csmall = amrex::max(small, small * amrex::max(qaux(i,j,k,QC), qaux(i-1,j,k,QC)));
      cavg = 0.5_rt*(qaux(i,j,k,QC) + qaux(i-1,j,k,QC));
      gamcl = qaux(i-1,j,k,QGAMC);
      gamcr = qaux(i,j,k,QGAMC);
#ifdef RADIATION
      gamcgl = qaux(i-1,j,k,QGAMCG);
      gamcgr = qaux(i,j,k,QGAMCG);
#endif

    } else if (idir == 1) {
      csmall = amrex::max(small, small * amrex::max(qaux(i,j,k,QC), qaux(i,j-1,k,QC)));
      cavg = 0.5_rt*(qaux(i,j,k,QC) + qaux(i,j-1,k,QC));
      gamcl = qaux(i,j-1,k,QGAMC);
      gamcr = qaux(i,j,k,QGAMC);
#ifdef RADIATION
      gamcgl = qaux(i,j-1,k,QGAMCG);
      gamcgr = qaux(i,j,k,QGAMCG);
#endif

    } else {
      csmall = amrex::max(small, small * amrex::max(qaux(i,j,k,QC), qaux(i,j,k-1,QC)));
      cavg = 0.5_rt*(qaux(i,j,k,QC) + qaux(i,j,k-1,QC));
      gamcl = qaux(i,j,k-1,QGAMC);
      gamcr = qaux(i,j,k,QGAMC);
#ifdef RADIATION
      gamcgl = qaux(i,j,k-1,QGAMCG);
      gamcgr = qaux(i,j,k,QGAMCG);
#endif
    }

#ifndef RADIATION
    if (use_reconstructed_gamma1 == 1) {
      gamcl = ql(i,j,k,QGC);
      gamcr = qr(i,j,k,QGC);

    } else if (compute_gammas == 1) {

      // we come in with a good p, rho, and X on the interfaces
      // -- use this to find the gamma used in the sound speed
      eos_t eos_state;
      eos_state.p = pl;
      eos_state.rho = rl;
      for (int n = 0; n < nspec; n++) {
        eos_state.xn[n] = ql(i,j,k,QFS+n);
      }
      eos_state.T = T_guess; // initial guess
      for (int n = 0; n < naux; n++) {
        eos_state.aux[:] = ql(i,j,k,QFX+n);
      }

      eos(eos_input_rp, eos_state);

      gamcl = eos_state.gam1;

      eos_state.p = pr;
      eos_state.rho = rr;
      for (int n = 0; n < nspec; n++) {
        eos_state.xn[n] = qr(i,j,k,QFS+n);
      }
      eos_state.T = T_guess; // initial guess
      for (int n = 0; n < naux; n++) {
        eos_state.aux[n] = qr(i,j,k,QFX+n);
      }

      eos(eos_input_rp, eos_state);

      gamcr = eos_state.gam1;
    }
#endif

    Real wsmall = small_dens*csmall;

    // this is Castro I: Eq. 33
    Real wl = amrex::max(wsmall, std::sqrt(std::abs(gamcl*pl*rl)));
    Real wr = amrex::max(wsmall, std::sqrt(std::abs(gamcr*pr*rr)));

    Real wwinv = 1.0_rt/(wl + wr);
    Real pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))*wwinv;
    Real ustar = ((wl*ul + wr*ur) + (pl - pr))*wwinv;

    pstar = amrex::max(pstar, small_pres);

    // for symmetry preservation, if ustar is really small, then we
    // set it to zero
    if (std::abs(ustar) < smallu*0.5_rt*(std::abs(ul) + std::abs(ur))) {
      ustar = 0.0_rt;
    }

    // look at the contact to determine which region we are in

    // this just determines which of the left or right states is still
    // in play.  We still need to look at the other wave to determine
    // if the star state or this state is on the interface.
    Real sgnm = std::copysign(1.0_rt, ustar);
    if (ustar == 0.0_rt) {
      sgnm = 0.0_rt;
    }

    Real fp = 0.5_rt*(1.0_rt + sgnm);
    Real fm = 0.5_rt*(1.0_rt - sgnm);

    Real ro = fp*rl + fm*rr;
    Real uo = fp*ul + fm*ur;
    Real po = fp*pl + fm*pr;
    Real reo = fp*rel + fm*rer;
    Real gamco = fp*gamcl + fm*gamcr;
#ifdef RADIATION
    Real lambda[ngroups];
    for (int g = 0; g < ngroups; g++) {
      lambda[g] = fp*laml[g] + fm*lamr[g];
    }

    if (ustar == 0) {
      // harmonic average
      for (int g = 0, g < ngroups; g++) {
        lambda[g] = 2.0_rt*(laml[g]*lamr[g])/(laml[g] + lamr[g] + 1.e-50_rt);
      }
    }

    Real po_g = fp*pl_g + fm*pr_g;
    Real reo_r[ngroups];
    Real po_r[ngroups];
    for (int g = 0; g < ngroups; g++) {
      reo_r[g] = fp*erl[g] + fm*err[g];
      po_r[g] = lambda[g]*reo_r[g];
    }
    Real reo_g = fp*rel_g + fm*rer_g;
    Real gamco_g = fp*gamcgl + fm*gamcgr;
#endif

    ro = amrex::max(small_dens, ro);

    Real roinv = 1.0_rt/ro;

    Real co = std::sqrt(std::abs(gamco*po*roinv));
    co = amrex::max(csmall, co);
    Real co2inv = 1.0_rt/(co*co);

    // we can already deal with the transverse velocities -- they
    // only jump across the contact
    qint(i,j,k,iv1) = fp*v1l + fm*v1r;
    qint(i,j,k,iv2) = fp*v2l + fm*v2r;

    // compute the rest of the star state

    Real drho = (pstar - po)*co2inv;
    Real rstar = ro + drho;
    rstar = amrex::max(small_dens, rstar);

#ifdef RADIATION
    Real estar_g = reo_g + drho*(reo_g + po_g)*roinv;

    Real co_g = std::sqrt(std::abs(gamco_g*po_g*roinv));
    co_g = amrex::max(csmall, co_g);

    Real pstar_g = po_g + drho*co_g*co_g;
    pstar_g = amrex::max(pstar_g, small_pres);

    for (int g = 0; g < ngroups; g++) {
      estar_r[g] = reo_r[g] + drho*(reo_r[g] + po_r[g])*roinv;
    }
#else
    Real entho = (reo + po)*roinv*co2inv;
    Real estar = reo + (pstar - po)*entho;
#endif

    Real cstar = std::sqrt(std::abs(gamco*pstar/rstar));
    cstar = amrex::max(cstar, csmall);

    // finish sampling the solution

    // look at the remaining wave to determine if the star state or the
    // 'o' state above is on the interface

    // the values of u +/- c on either side of the non-contact wave
    Real spout = co - sgnm*uo;
    Real spin = cstar - sgnm*ustar;

    // a simple estimate of the shock speed
    Real ushock = 0.5_rt*(spin + spout);

    if (pstar-po > 0.0_rt) {
      spin = ushock;
      spout = ushock;
    }

    if (spout-spin == 0.0_rt) {
      scr = small*cavg;
    } else {
      scr = spout - spin;
    }

    // interpolate for the case that we are in a rarefaction
    Real frac = (1.0_rt + (spout + spin)/scr)*0.5_rt;
    frac = amrex::max(0.0_rt, amrex::min(1.0_rt, frac));

    qint(i,j,k,QRHO) = frac*rstar + (1.0_rt - frac)*ro;
    qint(i,j,k,iu  ) = frac*ustar + (1.0_rt - frac)*uo;

#ifdef RADIATION
    Real pgdnv_t = frac*pstar + (1.0_rt - frac)*po;
    Real pgdnv_g = frac*pstar_g + (1.0_rt - frac)*po_g;
    Real regdnv_g = frac*estar_g + (1.0_rt - frac)*reo_g;
    Real regdnv_r[ngroups];
    for (int g = 0; g < ngroups; g++) {
      regdnv_r[g] = frac*estar_r[g] + (1.0_rt - frac)*reo_r[g];
    }
#else
    qint(i,j,k,QPRES) = frac*pstar + (1.0_rt - frac)*po;
    Real regdnv = frac*estar + (1.0_rt - frac)*reo;
#endif

    // as it stands now, we set things assuming that the rarefaction
    // spans the interface.  We overwrite that here depending on the
    // wave speeds

    // look at the speeds on either side of the remaining wave
    // to determine which region we are in
    if (spout < 0.0_rt) {
      // the l or r state is on the interface
      qint(i,j,k,QRHO) = ro;
      qint(i,j,k,iu  ) = uo;
#ifdef RADIATION
      pgdnv_t = po;
      pgdnv_g = po_g;
      regdnv_g = reo_g;
      for (int g = 0; g < ngroups; g++) {
        regdnv_r[g] = reo_r[g];
      }
#else
      qint(i,j,k,QPRES) = po;
      regdnv = reo;
#endif
    }

    if (spin >= 0.0_rt) {
      // the star state is on the interface
      qint(i,j,k,QRHO) = rstar;
      qint(i,j,k,iu  ) = ustar;
#ifdef RADIATION
      pgdnv_t = pstar;
      pgdnv_g = pstar_g;
      regdnv_g = estar_g;
      for (int g = 0; g < ngroups; g++) {
        regdnv_r[g] = estar_r[g];
      }
#else
      qint(i,j,k,QPRES) = pstar;
      regdnv = estar;
#endif
    }

#ifdef RADIATION
    for (int g = 0, g < ngroups; g++) {
      qint(i,j,k,QRAD+g) = amrex::max(regdnv_[g], 0.0_rt);
    }

    qint(i,j,k,QPRES) = pgdnv_g;
    qint(i,j,k,QPTOT) = pgdnv_t;
    qint(i,j,k,QREINT) = regdnv_g;

    qint(i,j,k,QREITOT) = regdnv_g;
    for (int g = 0; g < ngroups; g++) {
      qint(i,j,k,QREITOT) += regdnv_r[g];
    }

    for (int g = 0; g < ngroups; g++) {
      lambda_int(i,j,k,g) = lambda[g];
    }

#else
    qint(i,j,k,QPRES) = amrex::max(qint(i,j,k,QPRES), small_pres);
    qint(i,j,k,QREINT) = regdnv;
#endif

    // we are potentially thermodynamically inconsistent, fix that
    // here
    if (use_eos_in_riemann == 1) {
      // we need to know the species -- they only jump across
      // the contact
      eos_t eos_state;

      eos_state.rho = qint(i,j,k,QRHO);
      eos_state.p = qint(i,j,k,QPRES);

      for (int n = 0; n < nspec; n++) {
        eos_state.xn[n] = fp*ql(i,j,k,QFS+n) + fm*qr(i,j,k,QFS+n);
      }

      eos_state.T = T_guess;

      for (int n = 0; n < naux; n++) {
        eos_state.aux[n] = fp*ql(i,j,k,QFX+n) + fm*qr(i,j,k,QFX+n);
      }

      eos(eos_input_rp, eos_state);

      qint(i,j,k,QREINT) = eos_state.rho * eos_state.e;
    }

    // Enforce that fluxes through a symmetry plane or wall are hard zero.
    qint(i,j,k,iu) = qint(i,j,k,iu) * bnd_fac;

    // passively advected quantities
    for (int ipassive = 0; ipassive < npassive; ipassive++) {
      int nqp = qpass_map[ipassive];
      qint(i,j,k,nqp) = fp*ql(i,j,k,nqp) + fm*qr(i,j,k,nqp);
    }

  });
}


  subroutine HLLC(ql, ql_lo, ql_hi, &
                  qr, qr_lo, qr_hi, &
                  qaux, qa_lo, qa_hi, &
                  uflx, uflx_lo, uflx_hi, &
                  qint, q_lo, q_hi, &
                  idir, lo, hi, &
                  domlo, domhi)
    ! this is an implementation of the HLLC solver described in Toro's
    ! book.  it uses the simplest estimate of the wave speeds, since
    ! those should work for a general EOS.  We also initially do the
    ! CGF Riemann construction to get pstar and ustar, since we'll
    ! need to know the pressure and velocity on the interface for the
    ! pdV term in the internal energy update.

    use prob_params_module, only : physbc_lo, physbc_hi, &
         Symmetry, SlipWall, NoSlipWall

    implicit none

    integer, intent(in) :: ql_lo(3), ql_hi(3)
    integer, intent(in) :: qr_lo(3), qr_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: uflx_lo(3), uflx_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: idir, lo(3), hi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: ql(ql_lo(1):ql_hi(1),ql_lo(2):ql_hi(2),ql_lo(3):ql_hi(3),NQ)
    real(rt), intent(in) :: qr(qr_lo(1):qr_hi(1),qr_lo(2):qr_hi(2),qr_lo(3):qr_hi(3),NQ)

    ! note: qaux comes in dimensioned as the fully box, so use k3d to
    ! index in z
    real(rt), intent(in) :: qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3),NQAUX)

    real(rt), intent(inout) :: uflx(uflx_lo(1):uflx_hi(1),uflx_lo(2):uflx_hi(2),uflx_lo(3):uflx_hi(3),NVAR)
    real(rt), intent(inout) :: qint(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)

    integer :: i, j, k

    real(rt) :: rgdnv, regdnv
    real(rt) :: rl, ul, v1l, v2l, pl, rel
    real(rt) :: rr, ur, v1r, v2r, pr, rer
    real(rt) :: wl, wr, scr
    real(rt) :: rstar, cstar, estar, pstar, ustar
    real(rt) :: ro, uo, po, reo, co, gamco, entho
    real(rt) :: sgnm, spin, spout, ushock, frac
    real(rt) :: wsmall, csmall
    real(rt) :: cavg, gamcl, gamcr

    integer :: iu, iv1, iv2, sx, sy, sz
    logical :: special_bnd_lo, special_bnd_hi, special_bnd_lo_x, special_bnd_hi_x
    integer :: bnd_fac_x, bnd_fac_y, bnd_fac_z, bnd_fac
    real(rt) :: wwinv, roinv, co2inv

    real(rt) :: U_hllc_state(nvar), U_state(nvar), F_state(nvar)
    real(rt) :: S_l, S_r, S_c

    real(rt) :: q_zone(NQ)

    !$gpu

    if (idir == 1) then
       iu = QU
       iv1 = QV
       iv2 = QW
       sx = 1
       sy = 0
       sz = 0
    else if (idir == 2) then
       iu = QV
       iv1 = QU
       iv2 = QW
       sx = 0
       sy = 1
       sz = 0
    else
       iu = QW
       iv1 = QU
       iv2 = QV
       sx = 0
       sy = 0
       sz = 1
    end if

    special_bnd_lo = (physbc_lo(idir) == Symmetry &
         .or.         physbc_lo(idir) == SlipWall &
         .or.         physbc_lo(idir) == NoSlipWall)
    special_bnd_hi = (physbc_hi(idir) == Symmetry &
         .or.         physbc_hi(idir) == SlipWall &
         .or.         physbc_hi(idir) == NoSlipWall)

    if (idir == 1) then
       special_bnd_lo_x = special_bnd_lo
       special_bnd_hi_x = special_bnd_hi
    else
       special_bnd_lo_x = .false.
       special_bnd_hi_x = .false.
    end if

    do k = lo(3), hi(3)

       bnd_fac_z = 1
       if (idir == 3) then
          if ( k == domlo(3)   .and. special_bnd_lo .or. &
               k == domhi(3)+1 .and. special_bnd_hi ) then
             bnd_fac_z = 0
          end if
       end if

       do j = lo(2), hi(2)

          bnd_fac_y = 1
          if (idir == 2) then
             if ( j == domlo(2)   .and. special_bnd_lo .or. &
                  j == domhi(2)+1 .and. special_bnd_hi ) then
                bnd_fac_y = 0
             end if
          end if

          !dir$ ivdep
          do i = lo(1), hi(1)

             rl = max(ql(i,j,k,QRHO), small_dens)

             ! pick left velocities based on direction
             ul  = ql(i,j,k,iu)
             v1l = ql(i,j,k,iv1)
             v2l = ql(i,j,k,iv2)

             pl  = max(ql(i,j,k,QPRES), small_pres)
             rel = ql(i,j,k,QREINT)

             rr = max(qr(i,j,k,QRHO), small_dens)

             ! pick right velocities based on direction
             ur  = qr(i,j,k,iu)
             v1r = qr(i,j,k,iv1)
             v2r = qr(i,j,k,iv2)

             pr  = max(qr(i,j,k,QPRES), small_pres)
             rer = qr(i,j,k,QREINT)

             ! now we essentially do the CGF solver to get p and u on the
             ! interface, but we won't use these in any flux construction.
             csmall = max( small, max(small * qaux(i,j,k,QC) , small * qaux(i-sx,j-sy,k-sz,QC)) )
             cavg = 0.5_rt*(qaux(i,j,k,QC) + qaux(i-sx,j-sy,k-sz,QC))

             if (use_reconstructed_gamma1 == 1) then
                gamcl = ql(i,j,k,QGC)
                gamcr = qr(i,j,k,QGC)
             else
                gamcl = qaux(i-sx,j-sy,k-sz,QGAMC)
                gamcr = qaux(i,j,k,QGAMC)
             endif

             wsmall = small_dens*csmall
             wl = max(wsmall, sqrt(abs(gamcl*pl*rl)))
             wr = max(wsmall, sqrt(abs(gamcr*pr*rr)))

             wwinv = 1.0_rt/(wl + wr)
             pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))*wwinv
             ustar = ((wl*ul + wr*ur) + (pl - pr))*wwinv

             pstar = max(pstar, small_pres)
             ! for symmetry preservation, if ustar is really small, then we
             ! set it to zero
             if (abs(ustar) < smallu*0.5_rt*(abs(ul) + abs(ur))) then
                ustar = 0.0_rt
             endif

             if (ustar > 0.0_rt) then
                ro = rl
                uo = ul
                po = pl
                reo = rel
                gamco = gamcl

             else if (ustar < 0.0_rt) then
                ro = rr
                uo = ur
                po = pr
                reo = rer
                gamco = gamcr
             else
                ro = 0.5_rt*(rl + rr)
                uo = 0.5_rt*(ul + ur)
                po = 0.5_rt*(pl + pr)
                reo = 0.5_rt*(rel + rer)
                gamco = 0.5_rt*(gamcl + gamcr)
             endif
             ro = max(small_dens, ro)

             roinv = 1.0_rt/ro
             co = sqrt(abs(gamco*po*roinv))
             co = max(csmall, co)
             co2inv = 1.0_rt/(co*co)

             rstar = ro + (pstar - po)*co2inv
             rstar = max(small_dens, rstar)

             entho = (reo + po)*co2inv * roinv
             estar = reo + (pstar - po)*entho

             cstar = sqrt(abs(gamco*pstar/rstar))
             cstar = max(cstar, csmall)

             sgnm = sign(1.0_rt, ustar)
             spout = co - sgnm*uo
             spin = cstar - sgnm*ustar
             ushock = 0.5_rt*(spin + spout)

             if (pstar-po > 0.0_rt) then
                spin = ushock
                spout = ushock
             endif
             if (spout-spin == 0.0_rt) then
                scr = small*cavg
             else
                scr = spout-spin
             endif
             frac = (1.0_rt + (spout + spin)/scr)*0.5_rt
             frac = max(0.0_rt, min(1.0_rt, frac))

             rgdnv = frac*rstar + (1.0_rt - frac)*ro
             regdnv = frac*estar + (1.0_rt - frac)*reo

             qint(i,j,k,iu) = frac*ustar + (1.0_rt - frac)*uo
             qint(i,j,k,QPRES) = frac*pstar + (1.0_rt - frac)*po


             ! now we do the HLLC construction


             ! Enforce that the fluxes through a symmetry plane or wall are hard zero.
             if ( special_bnd_lo_x .and. i== domlo(1) .or. &
                  special_bnd_hi_x .and. i== domhi(1)+1 ) then
                bnd_fac_x = 0
             else
                bnd_fac_x = 1
             end if

             bnd_fac = bnd_fac_x*bnd_fac_y*bnd_fac_z

             ! use the simplest estimates of the wave speeds
             S_l = min(ul - sqrt(gamcl*pl/rl), ur - sqrt(gamcr*pr/rr))
             S_r = max(ul + sqrt(gamcl*pl/rl), ur + sqrt(gamcr*pr/rr))

             ! estimate of the contact speed -- this is Toro Eq. 10.8
             S_c = (pr - pl + rl*ul*(S_l - ul) - rr*ur*(S_r - ur))/ &
                  (rl*(S_l - ul) - rr*(S_r - ur))

             if (S_r <= 0.0_rt) then
                ! R region
                q_zone(:) = qr(i,j,k,:)
                call cons_state(q_zone, U_state)
                call compute_flux(idir, bnd_fac, U_state, pr, F_state)

             else if (S_r > 0.0_rt .and. S_c <= 0.0_rt) then
                ! R* region
                q_zone(:) = qr(i,j,k,:)
                call cons_state(q_zone, U_state)
                call compute_flux(idir, bnd_fac, U_state, pr, F_state)

                call HLLC_state(idir, S_r, S_c, q_zone, U_hllc_state)

                ! correct the flux
                F_state(:) = F_state(:) + S_r*(U_hllc_state(:) - U_state(:))

             else if (S_c > 0.0_rt .and. S_l < 0.0_rt) then
                ! L* region
                q_zone(:) = ql(i,j,k,:)
                call cons_state(q_zone, U_state)
                call compute_flux(idir, bnd_fac, U_state, pl, F_state)

                call HLLC_state(idir, S_l, S_c, q_zone, U_hllc_state)

                ! correct the flux
                F_state(:) = F_state(:) + S_l*(U_hllc_state(:) - U_state(:))

             else
                ! L region
                q_zone(:) = ql(i,j,k,:)
                call cons_state(q_zone, U_state)
                call compute_flux(idir, bnd_fac, U_state, pl, F_state)

             endif

             uflx(i,j,k,:) = F_state(:)
          end do
       end do
    end do

  end subroutine HLLC
