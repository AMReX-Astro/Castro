#ifndef riemann_solvers_H
#define riemann_solvers_H

#include <Castro_util.H>
#include <riemann.H>

///
/// The Colella-Glaz Riemann solver for pure hydrodynamics.  This is a
/// two shock approximate state Riemann solver.
///
/// @param bx         the box to operate over
/// @param ql         the left interface state
/// @param qr         the right interface state
/// @param qaux_arr   the auxillary state
/// @param qint       the full Godunov state on the interface
/// @param idir       coordinate direction for the solve (0 = x, 1 = y, 2 = z)
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
riemanncg(const int i, const int j, const int k,
          const RiemannState ql, const RiemannState qr, const RiemannAux raux,
          Array4<Real> const& qint,
          const int idir) {

  // this implements the approximate Riemann solver of Colella & Glaz
  // (1985)
  //

  constexpr Real weakwv = 1.e-3_rt;

  int iu, iv1, iv2;

  if (idir == 0) {
    iu = QU;
    iv1 = QV;
    iv2 = QW;
  } else if (idir == 1) {
    iu = QV;
    iv1 = QU;
    iv2 = QW;
  } else {
    iu = QW;
    iv1 = QU;
    iv2 = QV;
  }



#ifndef AMREX_USE_GPU
  GpuArray<Real, HISTORY_SIZE> pstar_hist;
#endif


  // common quantities
  Real taul = 1.0_rt / ql.rho;
  Real taur = 1.0_rt / qr.rho;

  // lagrangian sound speeds
  Real clsql = ql.gamc * ql.p * ql.rho;
  Real clsqr = qr.gamc * qr.p * qr.rho;

  // Note: in the original Colella & Glaz paper, they predicted
  // gamma_e to the interfaces using a special (non-hyperbolic)
  // evolution equation.  In Castro, we instead bring (rho e)
  // to the edges, so we construct the necessary gamma_e here from
  // what we have on the interfaces.
  Real gamel = ql.p / ql.rhoe + 1.0_rt;
  Real gamer = qr.p / qr.rhoe + 1.0_rt;

  // these should consider a wider average of the cell-centered
  // gammas
  Real gmin = amrex::min(amrex::min(gamel, gamer), 1.0_rt);
  Real gmax = amrex::max(amrex::max(gamel, gamer), 2.0_rt);

  Real game_bar = 0.5_rt*(gamel + gamer);
  Real gamc_bar = 0.5_rt*(ql.gamc + qr.gamc);

  Real gdot = 2.0_rt*(1.0_rt - game_bar/gamc_bar)*(game_bar - 1.0_rt);

  Real wsmall = small_dens * raux.csmall;
  Real wl = amrex::max(wsmall, std::sqrt(std::abs(clsql)));
  Real wr = amrex::max(wsmall, std::sqrt(std::abs(clsqr)));

  // make an initial guess for pstar -- this is a two-shock
  // approximation
  //pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
  Real pstar = ql.p + ( (qr.p - ql.p) - wr*(qr.un - ql.un) ) * wl / (wl + wr);
  pstar = amrex::max(pstar, small_pres);

  // get the shock speeds -- this computes W_s from CG Eq. 34
  Real gamstar = 0.0;

  Real wlsq = 0.0;
  wsqge(ql.p, taul, gamel, gdot, gamstar,
        gmin, gmax, clsql, pstar, wlsq);

  Real wrsq = 0.0;
  wsqge(qr.p, taur, gamer, gdot, gamstar,
        gmin, gmax, clsqr, pstar, wrsq);

  Real pstar_old = pstar;

  wl = std::sqrt(wlsq);
  wr = std::sqrt(wrsq);

  // R-H jump conditions give ustar across each wave -- these
  // should be equal when we are done iterating.  Our notation
  // here is a little funny, comparing to CG, ustar_l = u*_L and
  // ustar_r = u*_R.
  Real ustar_l = ql.un - (pstar - ql.p) / wl;
  Real ustar_r = qr.un + (pstar - qr.p) / wr;

  // revise our pstar guess
  // pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
  pstar = ql.p + ( (qr.p - ql.p) - wr * (qr.un - ql.un) ) * wl / (wl + wr);
  pstar = amrex::max(pstar, small_pres);

  // secant iteration
  bool converged = false;

  int iter = 0;
  while ((iter < cg_maxiter && !converged) || iter < 2) {

      wsqge(ql.p, taul, gamel, gdot, gamstar,
            gmin, gmax, clsql, pstar, wlsq);

      wsqge(qr.p, taur, gamer, gdot, gamstar,
            gmin, gmax, clsqr, pstar, wrsq);


      // NOTE: these are really the inverses of the wave speeds!
      wl = 1.0_rt / std::sqrt(wlsq);
      wr = 1.0_rt / std::sqrt(wrsq);

      Real ustar_r_old = ustar_r;
      Real ustar_l_old = ustar_l;

      ustar_r = qr.un - (qr.p - pstar) * wr;
      ustar_l = ql.un + (ql.p - pstar) * wl;

      Real dpditer = std::abs(pstar_old - pstar);

      // Here we are going to do the Secant iteration version in
      // CG.  Note that what we call zp and zm here are not
      // actually the Z_p = |dp*/du*_p| defined in CG, by rather
      // simply |du*_p| (or something that looks like dp/Z!).
      Real zp = std::abs(ustar_l - ustar_l_old);
      if (zp - weakwv * raux.cavg <= 0.0_rt) {
          zp = dpditer * wl;
      }

      Real zm = std::abs(ustar_r - ustar_r_old);
      if (zm - weakwv * raux.cavg <= 0.0_rt) {
          zm = dpditer * wr;
      }

      // the new pstar is found via CG Eq. 18
      Real denom = dpditer / amrex::max(zp + zm, riemann_constants::small * raux.cavg);
      pstar_old = pstar;
      pstar = pstar - denom*(ustar_r - ustar_l);
      pstar = amrex::max(pstar, small_pres);

      Real err = std::abs(pstar - pstar_old);
      if (err < cg_tol*pstar) {
          converged = true;
      }

#ifndef AMREX_USE_GPU
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

#ifndef AMREX_USE_GPU
          std::cout <<  "pstar history: " << std::endl;
          for (int iter_l=0; iter_l < cg_maxiter; iter_l++) {
              std::cout << iter_l << " " << pstar_hist[iter_l] << std::endl;
          }

          std::cout << std::endl;
          std::cout << "left state: " << std::endl << ql << std::endl;
          std::cout << "right state: " << std::endl << qr << std::endl;
          std::cout << "aux information: " << std::endl << raux << std::endl;

          amrex::Error("ERROR: non-convergence in the Riemann solver");
#endif

      } else if (cg_blend == 1) {

          pstar = ql.p + ( (qr.p - ql.p) - wr * (qr.un - ql.un) ) * wl / (wl + wr);

      } else if (cg_blend == 2) {

          // we don't store the history if we are in CUDA, so
          // we can't do this
#ifndef AMREX_USE_GPU
          // first try to find a reasonable bounds
          Real pstarl = 1.e200;
          Real pstaru = -1.e200;
          for (int n = cg_maxiter-6; n < cg_maxiter; n++) {
              pstarl = amrex::min(pstarl, pstar_hist[n]);
              pstaru = amrex::max(pstaru, pstar_hist[n]);
          }

          pstarl = amrex::max(pstarl, small_pres);
          pstaru = amrex::max(pstaru, small_pres);

          GpuArray<Real, PSTAR_BISECT_FACTOR*HISTORY_SIZE> pstar_hist_extra;

          pstar_bisection(pstarl, pstaru,
                          ql.un, ql.p, taul, gamel, clsql,
                          qr.un, qr.p, taur, gamer, clsqr,
                          gdot, gmin, gmax,
                          cg_maxiter, cg_tol,
                          pstar, gamstar, converged, pstar_hist_extra);

          if (!converged) {

              std::cout << "pstar history: " << std::endl;
              for (int iter_l = 0; iter_l < cg_maxiter; iter_l++) {
                  std::cout << iter_l << " " << pstar_hist[iter_l] << std::endl;
              }
              std::cout << "pstar extra history: " << std::endl;
              for (int iter_l = 0; iter_l < PSTAR_BISECT_FACTOR*cg_maxiter; iter_l++) {
                  std::cout << iter_l << " " << pstar_hist_extra[iter_l] << std::endl;
              }

              std::cout << std::endl;
              std::cout << "left state: " << std::endl << ql << std::endl;
              std::cout << "right state: " << std::endl << qr << std::endl;
              std::cout << "aux information: " << std::endl << raux << std::endl;

              amrex::Error("ERROR: non-convergence in the Riemann solver");
          }

#endif
      } else {

#ifndef AMREX_USE_GPU
          amrex::Error("ERROR: unrecognized cg_blend option.");
#endif
      }

  }

  // we converged!  construct the single ustar for the region
  // between the left and right waves, using the updated wave speeds
  ustar_r = qr.un - (qr.p - pstar) * wr;  // careful -- here wl, wr are 1/W
  ustar_l = ql.un + (ql.p - pstar) * wl;

  Real ustar = 0.5_rt * (ustar_l + ustar_r);

  // for symmetry preservation, if ustar is really small, then we
  // set it to zero
  if (std::abs(ustar) < riemann_constants::smallu * 0.5_rt * (std::abs(ql.un) + std::abs(qr.un))) {
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
      ro = ql.rho;
      uo = ql.un;
      po = ql.p;
      tauo = taul;
      gamco = ql.gamc;
      gameo = gamel;

  } else if (ustar < 0.0_rt) {
      ro = qr.rho;
      uo = qr.un;
      po = qr.p;
      tauo = taur;
      gamco = qr.gamc;
      gameo = gamer;

  } else {
      ro = 0.5_rt * (ql.rho + qr.rho);
      uo = 0.5_rt * (ql.un + qr.un);
      po = 0.5_rt * (ql.p + qr.p);
      tauo = 0.5_rt * (taul + taur);
      gamco = 0.5_rt * (ql.gamc + qr.gamc);
      gameo = 0.5_rt * (gamel + gamer);
  }

  // use tau = 1/rho as the independent variable here
  ro = amrex::max(small_dens, 1.0_rt/tauo);
  tauo = 1.0_rt/ro;

  Real co = std::sqrt(std::abs(gamco*po*tauo));
  co = amrex::max(raux.csmall, co);
  Real clsq = std::pow(co*ro, 2);

  // now that we know which state (left or right) we need to worry
  // about, get the value of gamstar and wosq across the wave we
  // are dealing with.
  Real wosq = 0.0;
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

  Real cstar = std::sqrt(std::abs(gamco * pstar / rstar));
  cstar = amrex::max(cstar, raux.csmall);

  Real spout = co - sgnm*uo;
  Real spin = cstar - sgnm*ustar;

  //ushock = 0.5_rt*(spin + spout)
  Real ushock = wo*tauo - sgnm*uo;

  if (pstar - po >= 0.0_rt) {
      spin = ushock;
      spout = ushock;
  }

  Real frac = 0.5_rt*(1.0_rt + (spin + spout)/amrex::max(amrex::max(spout-spin, spin+spout),
                                                         riemann_constants::small * raux.cavg));

  // the transverse velocity states only depend on the
  // direction that the contact moves
  if (ustar > 0.0_rt) {
      qint(i,j,k,iv1) = ql.ut;
      qint(i,j,k,iv2) = ql.utt;
  } else if (ustar < 0.0_rt) {
      qint(i,j,k,iv1) = qr.ut;
      qint(i,j,k,iv2) = qr.utt;
  } else {
      qint(i,j,k,iv1) = 0.5_rt * (ql.ut + qr.ut);
      qint(i,j,k,iv2) = 0.5_rt * (ql.utt + qr.utt);
  }

  // linearly interpolate between the star and normal state -- this covers the
  // case where we are inside the rarefaction fan.
  qint(i,j,k,QRHO) = frac*rstar + (1.0_rt - frac)*ro;
  qint(i,j,k,iu) = frac*ustar + (1.0_rt - frac)*uo;
  qint(i,j,k,QPRES) = frac*pstar + (1.0_rt - frac)*po;
  Real game_int = frac*gamstar + (1.0_rt-frac)*gameo;

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

  qint(i,j,k,iu) = qint(i,j,k,iu) * raux.bnd_fac;

  // compute the total energy from the internal, p/(gamma - 1), and the kinetic
  qint(i,j,k,QREINT) = qint(i,j,k,QPRES)/(game_int - 1.0_rt);

  // we'll do the passive scalars separately

}


///
/// The Colella-Glaz-Ferguson Riemann solver for hydrodynamics and
/// radiation hydrodynamics.  This is a two shock approximate state
/// Riemann solver.
///
/// @param bx          the box to operate over
/// @param ql          the left interface state
/// @param qr          the right interface state
/// @param qaux_arr    the auxillary state
/// @param qint        the full Godunov state on the interface
/// @param lambda_int  the radiation flux limiter on the interface
/// @param idir        coordinate direction for the solve (0 = x, 1 = y, 2 = z)
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
riemannus(const int i, const int j, const int k,
          const RiemannState ql, const RiemannState qr, const RiemannAux raux,
          Array4<Real> const& qint,
#ifdef RADIATION
          Array4<Real> const& lambda_int,
#endif
          const int idir) {

  // Colella, Glaz, and Ferguson solver
  //
  // this is a 2-shock solver that uses a very simple approximation for the
  // star state, and carries an auxiliary jump condition for (rho e) to
  // deal with a real gas

  // set integer pointers for the normal and transverse velocity and
  // momentum

  int iu, iv1, iv2;

  if (idir == 0) {
    iu = QU;
    iv1 = QV;
    iv2 = QW;

  } else if (idir == 1) {
    iu = QV;
    iv1 = QU;
    iv2 = QW;

  } else {
    iu = QW;
    iv1 = QU;
    iv2 = QV;

  }

  // estimate the star state: pstar, ustar

  Real wsmall = small_dens * raux.csmall;

  // this is Castro I: Eq. 33

  Real wl = amrex::max(wsmall, std::sqrt(std::abs(ql.gamc * ql.p * ql.rho)));
  Real wr = amrex::max(wsmall, std::sqrt(std::abs(qr.gamc * qr.p * qr.rho)));

  Real wwinv = 1.0_rt/(wl + wr);
  Real pstar = ((wr * ql.p + wl * qr.p) + wl * wr * (ql.un - qr.un)) * wwinv;
  Real ustar = ((wl * ql.un + wr * qr.un) + (ql.p - qr.p)) * wwinv;

  pstar = amrex::max(pstar, small_pres);

  // for symmetry preservation, if ustar is really small, then we
  // set it to zero

  if (std::abs(ustar) < riemann_constants::smallu * 0.5_rt * (std::abs(ql.un) + std::abs(qr.un))) {
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

  Real ro = fp * ql.rho + fm * qr.rho;
  Real uo = fp * ql.un + fm * qr.un;
  Real po = fp * ql.p + fm * qr.p;
  Real reo = fp * ql.rhoe + fm * qr.rhoe;
  Real gamco = fp * ql.gamc + fm * qr.gamc;
#ifdef RADIATION
  Real lambda[NGROUPS];
  for (int g = 0; g < NGROUPS; g++) {
      lambda[g] = fp * ql.lam[g] + fm * qr.lam[g];
  }

  if (ustar == 0) {
      // harmonic average
      for (int g = 0; g < NGROUPS; g++) {
          lambda[g] = 2.0_rt * (ql.lam[g] * qr.lam[g]) / (ql.lam[g] + qr.lam[g] + 1.e-50_rt);
      }
  }

  Real po_g = fp * ql.p_g + fm * qr.p_g;
  Real reo_r[NGROUPS];
  Real po_r[NGROUPS];
  for (int g = 0; g < NGROUPS; g++) {
      reo_r[g] = fp * ql.er[g] + fm * qr.er[g];
      po_r[g] = lambda[g] * reo_r[g];
  }
  Real reo_g = fp * ql.rhoe_g + fm * qr.rhoe_g;
  Real gamco_g = fp * ql.gamcg + fm * qr.gamcg;
#endif

  ro = amrex::max(small_dens, ro);

  Real roinv = 1.0_rt / ro;

  Real co = std::sqrt(std::abs(gamco * po * roinv));
  co = amrex::max(raux.csmall, co);
  Real co2inv = 1.0_rt / (co*co);

  // we can already deal with the transverse velocities -- they
  // only jump across the contact

  qint(i,j,k,iv1) = fp * ql.ut + fm * qr.ut;
  qint(i,j,k,iv2) = fp * ql.utt + fm * qr.utt;

  // compute the rest of the star state

  Real drho = (pstar - po)*co2inv;
  Real rstar = ro + drho;
  rstar = amrex::max(small_dens, rstar);

#ifdef RADIATION
  Real estar_g = reo_g + drho*(reo_g + po_g)*roinv;

  Real co_g = std::sqrt(std::abs(gamco_g*po_g*roinv));
  co_g = amrex::max(raux.csmall, co_g);

  Real pstar_g = po_g + drho*co_g*co_g;
  pstar_g = amrex::max(pstar_g, small_pres);

  Real estar_r[NGROUPS];
  for (int g = 0; g < NGROUPS; g++) {
      estar_r[g] = reo_r[g] + drho*(reo_r[g] + po_r[g])*roinv;
  }
#else
  Real entho = (reo + po)*roinv*co2inv;
  Real estar = reo + (pstar - po)*entho;
#endif

  Real cstar = std::sqrt(std::abs(gamco*pstar/rstar));
  cstar = amrex::max(cstar, raux.csmall);

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

  Real scr = spout - spin;
  if (spout-spin == 0.0_rt) {
      scr = riemann_constants::small * raux.cavg;
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
  Real regdnv_r[NGROUPS];
  for (int g = 0; g < NGROUPS; g++) {
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
      for (int g = 0; g < NGROUPS; g++) {
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
      for (int g = 0; g < NGROUPS; g++) {
          regdnv_r[g] = estar_r[g];
      }
#else
      qint(i,j,k,QPRES) = pstar;
      regdnv = estar;
#endif
    }

#ifdef RADIATION
  for (int g = 0; g < NGROUPS; g++) {
      qint(i,j,k,QRAD+g) = amrex::max(regdnv_r[g], 0.0_rt);
  }

  qint(i,j,k,QPRES) = pgdnv_g;
  qint(i,j,k,QPTOT) = pgdnv_t;
  qint(i,j,k,QREINT) = regdnv_g;

  qint(i,j,k,QREITOT) = regdnv_g;
  for (int g = 0; g < NGROUPS; g++) {
      qint(i,j,k,QREITOT) += regdnv_r[g];
  }

  for (int g = 0; g < NGROUPS; g++) {
      lambda_int(i,j,k,g) = lambda[g];
  }

#else
  qint(i,j,k,QPRES) = amrex::max(qint(i,j,k,QPRES), small_pres);
  qint(i,j,k,QREINT) = regdnv;
#endif

  // Enforce that fluxes through a symmetry plane or wall are hard zero.
  qint(i,j,k,iu) = qint(i,j,k,iu) * raux.bnd_fac;

  // we'll do the passive scalars separately

}

///
/// A simple HLL Riemann solver for pure hydrodynamics.  This takes just a
/// single interface's data and returns the HLL flux
///
/// @param ql     the left interface state
/// @param qr     the right interface state
/// @param cl     sound speed on the left interface
/// @param cr     sound speed on the right interface
/// @param idir   coordinate direction for the solve (0 = x, 1 = y, 2 = z)
/// @param coord  geometry type (0 = Cartesian, 1 = axisymmetric, 2 = spherical)
/// @param f      the HLL fluxes
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
HLL(const Real* ql, const Real* qr,
    const Real cl, const Real cr,
    const int idir, const int coord,
    Real* flux_hll) {

  // This is the HLLE solver.  We should apply it to zone averages
  // (not reconstructed states) at an interface in the presence of
  // shocks to avoid the odd-even decoupling / carbuncle phenomenon.
  //
  // See: Einfeldt, B.  et al. 1991, JCP, 92, 273
  //      Einfeldt, B. 1988, SIAM J NA, 25, 294


  constexpr Real small_hll = 1.e-10_rt;

  int ivel, ivelt, iveltt;
  int imom, imomt, imomtt;

  if (idir == 0) {
    ivel = QU;
    ivelt = QV;
    iveltt = QW;

    imom = UMX;
    imomt = UMY;
    imomtt = UMZ;

  } else if (idir == 1) {
    ivel = QV;
    ivelt = QU;
    iveltt = QW;

    imom = UMY;
    imomt = UMX;
    imomtt = UMZ;

  } else {
    ivel = QW;
    ivelt = QU;
    iveltt = QV;

    imom = UMZ;
    imomt = UMX;
    imomtt = UMY;
  }

  Real rhol_sqrt = std::sqrt(ql[QRHO]);
  Real rhor_sqrt = std::sqrt(qr[QRHO]);

  Real rhod = 1.0_rt/(rhol_sqrt + rhor_sqrt);


  // compute the average sound speed. This uses an approximation from
  // E88, eq. 5.6, 5.7 that assumes gamma falls between 1
  // and 5/3
  Real cavg = std::sqrt( (rhol_sqrt*cl*cl + rhor_sqrt*cr*cr)*rhod +
                         0.5_rt*rhol_sqrt*rhor_sqrt*rhod*rhod*std::pow(qr[ivel] - ql[ivel], 2));


  // Roe eigenvalues (E91, eq. 5.3b)
  Real uavg = (rhol_sqrt*ql[ivel] + rhor_sqrt*qr[ivel])*rhod;

  Real a1 = uavg - cavg;
  Real a4 = uavg + cavg;


  // signal speeds (E91, eq. 4.5)
  Real bl = amrex::min(a1, ql[ivel] - cl);
  Real br = amrex::max(a4, qr[ivel] + cr);

  Real bm = amrex::min(0.0_rt, bl);
  Real bp = amrex::max(0.0_rt, br);

  Real bd = bp - bm;

  if (std::abs(bd) < small_hll*amrex::max(std::abs(bm), std::abs(bp))) return;

  // we'll overwrite the passed in flux with the HLL flux

  bd = 1.0_rt/bd;

  // compute the fluxes according to E91, eq. 4.4b -- note that the
  // min/max above picks the correct flux if we are not in the star
  // region

  // density flux
  Real fl_tmp = ql[QRHO]*ql[ivel];
  Real fr_tmp = qr[QRHO]*qr[ivel];

  flux_hll[URHO] = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr[QRHO] - ql[QRHO]);

  // normal momentum flux.  Note for 1-d and 2-d non cartesian
  // r-coordinate, we leave off the pressure term and handle that
  // separately in the update, to accommodate different geometries
  fl_tmp = ql[QRHO]*ql[ivel]*ql[ivel];
  fr_tmp = qr[QRHO]*qr[ivel]*qr[ivel];
  if (mom_flux_has_p(idir, idir, coord)) {
    fl_tmp = fl_tmp + ql[QPRES];
    fr_tmp = fr_tmp + qr[QPRES];
  }

  flux_hll[imom] = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr[QRHO]*qr[ivel] - ql[QRHO]*ql[ivel]);

  // transverse momentum flux
  fl_tmp = ql[QRHO]*ql[ivel]*ql[ivelt];
  fr_tmp = qr[QRHO]*qr[ivel]*qr[ivelt];

  flux_hll[imomt] = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr[QRHO]*qr[ivelt] - ql[QRHO]*ql[ivelt]);


  fl_tmp = ql[QRHO]*ql[ivel]*ql[iveltt];
  fr_tmp = qr[QRHO]*qr[ivel]*qr[iveltt];

  flux_hll[imomtt] = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr[QRHO]*qr[iveltt] - ql[QRHO]*ql[iveltt]);

  // total energy flux
  Real rhoEl = ql[QREINT] + 0.5_rt*ql[QRHO]*(ql[ivel]*ql[ivel] + ql[ivelt]*ql[ivelt] + ql[iveltt]*ql[iveltt]);
  fl_tmp = ql[ivel]*(rhoEl + ql[QPRES]);

  Real rhoEr = qr[QREINT] + 0.5_rt*qr[QRHO]*(qr[ivel]*qr[ivel] + qr[ivelt]*qr[ivelt] + qr[iveltt]*qr[iveltt]);
  fr_tmp = qr[ivel]*(rhoEr + qr[QPRES]);

  flux_hll[UEDEN] = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(rhoEr - rhoEl);


  // eint flux
  fl_tmp = ql[QREINT]*ql[ivel];
  fr_tmp = qr[QREINT]*qr[ivel];

  flux_hll[UEINT] = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr[QREINT] - ql[QREINT]);


  // passively-advected scalar fluxes
  for (int ipassive = 0; ipassive < npassive; ipassive++) {
    int n  = upassmap(ipassive);
    int nqs = qpassmap(ipassive);

    fl_tmp = ql[QRHO]*ql[nqs]*ql[ivel];
    fr_tmp = qr[QRHO]*qr[nqs]*qr[ivel];

    flux_hll[n] = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr[QRHO]*qr[nqs] - ql[QRHO]*ql[nqs]);
  }
}

#endif
