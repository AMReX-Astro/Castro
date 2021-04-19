#ifndef riemann_solvers_H
#define riemann_solvers_H

#include <Castro_util.H>
#include <riemann.H>
#ifdef HYBRID_MOMENTUM
#include <hybrid.H>
#endif
#ifdef RADIATION
#include <Radiation.H>
#include <rad_util.H>
#endif

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
compute_flux_q(const int i, const int j, const int k, const int idir,
               const GeometryData& geomdata,
               const RiemannState& qint,
               Array4<Real> const& F,
#ifdef RADIATION
               Array4<Real> const& rF,
#endif
               Array4<Real> const& qgdnv, const bool store_full_state) {

    // given a primitive state, compute the flux in direction idir
    //

    int im1, im2, im3;

    auto coord = geomdata.Coord();
    auto mom_check = mom_flux_has_p(idir, idir, coord);

    if (idir == 0) {
        im1 = UMX;
        im2 = UMY;
        im3 = UMZ;

    } else if (idir == 1) {
        im1 = UMY;
        im2 = UMX;
        im3 = UMZ;

    } else {
        im1 = UMZ;
        im2 = UMX;
        im3 = UMY;
    }

#ifdef RADIATION
    int fspace_t = Radiation::fspace_advection_type;
    int comov = Radiation::comoving;
    int limiter = Radiation::limiter;
    int closure = Radiation::closure;
#endif

    // Compute fluxes, order as conserved state (not q)
    F(i,j,k,URHO) = qint.rho * qint.un;

    F(i,j,k,im1) = F(i,j,k,URHO) * qint.un;
    F(i,j,k,im2) = F(i,j,k,URHO) * qint.ut;
    F(i,j,k,im3) = F(i,j,k,URHO) * qint.utt;

#ifndef RADIATION
    if (mom_check) {
      F(i,j,k,im1) += qint.p;
    }

    Real rhoetot = qint.rhoe + 0.5_rt * qint.rho *
        (qint.un * qint.un +
         qint.ut * qint.ut +
         qint.utt * qint.utt);

    F(i,j,k,UEDEN) = qint.un * (rhoetot + qint.p);
    F(i,j,k,UEINT) = qint.un * qint.rhoe;
#else
    if (mom_check) {
      F(i,j,k,im1) += qint.p_g;
    }

    Real rhoetot = qint.rhoe_g + 0.5_rt * qint.rho *
        (qint.un * qint.un +
         qint.ut * qint.ut +
         qint.utt * qint.utt);

    F(i,j,k,UEDEN) = qint.un * (rhoetot + qint.p_g);
    F(i,j,k,UEINT) = qint.un * qint.rhoe_g;
#endif

    F(i,j,k,UTEMP) = 0.0;
#ifdef SHOCK_VAR
    F(i,j,k,USHK) = 0.0;
#endif

#ifdef RADIATION
    if (fspace_t == 1) {
        for (int g = 0; g < NGROUPS; g++) {
            Real eddf = Edd_factor(qint.lam[g], limiter, closure);
            Real f1 = 0.5e0_rt*(1.0_rt-eddf);
            rF(i,j,k,g) = (1.0_rt + f1) * qint.er[g] * qint.un;
        }
    } else {
        // type 2
        for (int g = 0; g < NGROUPS; g++) {
            rF(i,j,k,g) = qint.er[g] * qint.un;
        }
    }
#endif

#ifdef HYBRID_MOMENTUM
    // the hybrid routine uses the Godunov indices, not the full NQ state
    GpuArray<Real, NGDNV> qgdnv_zone;
    qgdnv_zone[GDRHO] = qint.rho;
    if (idir == 0) {
        qgdnv_zone[GDU] = qint.un;
        qgdnv_zone[GDV] = qint.ut;
        qgdnv_zone[GDW] = qint.utt;
    } else if (idir == 1) {
        qgdnv_zone[GDV] = qint.un;
        qgdnv_zone[GDU] = qint.ut;
        qgdnv_zone[GDW] = qint.utt;
    } else {
        qgdnv_zone[GDW] = qint.un;
        qgdnv_zone[GDU] = qint.ut;
        qgdnv_zone[GDV] = qint.utt;
    }
    qgdnv_zone[GDPRES] = qint.p;
#ifdef RADIATION
    for (int g = 0; g < NGROUPS; g++) {
        qgdnv_zone[GDLAMS+g] = qint.lam[g];
        qgdnv_zone[GDERADS+g] = qint.er[g];
    }
#endif
    GpuArray<Real, NUM_STATE> F_zone;
    for (int n = 0; n < NUM_STATE; n++) {
        F_zone[n] = F(i,j,k,n);
    }
    compute_hybrid_flux(qgdnv_zone, geomdata, idir, i, j, k, F_zone);
    for (int n = 0; n < NUM_STATE; n++) {
        F(i,j,k,n) = F_zone[n];
    }
#endif

    // now we store the interface state in the full 4-d array.  We
    // either store just the NGDNV subset or the whole NQ hydro state,
    // depending on the value of store_full_state

    // note: we don't handle the passives here

    if (store_full_state) {

        // this is the primitive variable state on the interface

        qgdnv(i,j,k,QRHO) = qint.rho;
        if (idir == 0) {
            qgdnv(i,j,k,QU) = qint.un;
            qgdnv(i,j,k,QV) = qint.ut;
            qgdnv(i,j,k,QW) = qint.utt;
        } else if (idir == 1) {
            qgdnv(i,j,k,QV) = qint.un;
            qgdnv(i,j,k,QU) = qint.ut;
            qgdnv(i,j,k,QW) = qint.utt;
        } else {
            qgdnv(i,j,k,QW) = qint.un;
            qgdnv(i,j,k,QU) = qint.ut;
            qgdnv(i,j,k,QV) = qint.utt;
        }

        qgdnv(i,j,k,QTEMP) = 0.0_rt;

#ifdef RADIATION
        qgdnv(i,j,k,QPRES) = qint.p_g;
        qgdnv(i,j,k,QREINT) = qint.rhoe_g;
        qgdnv(i,j,k,QPTOT) = qint.p;
        qgdnv(i,j,k,QREITOT) = qint.rhoe;
        for (int g = 0; g < NGROUPS; g++) {
            qgdnv(i,j,k,QRAD) = qint.er[g];
        }
#else
        qgdnv(i,j,k,QPRES) = qint.p;
        qgdnv(i,j,k,QREINT) = qint.rhoe;
#endif

    } else {

        // this is the reduced "Godunov" state on the interface,
        // with NGDNV components

        qgdnv(i,j,k,GDRHO) = qint.rho;
        if (idir == 0) {
            qgdnv(i,j,k,GDU) = qint.un;
            qgdnv(i,j,k,GDV) = qint.ut;
            qgdnv(i,j,k,GDW) = qint.utt;
        } else if (idir == 1) {
            qgdnv(i,j,k,GDV) = qint.un;
            qgdnv(i,j,k,GDU) = qint.ut;
            qgdnv(i,j,k,GDW) = qint.utt;
        } else {
            qgdnv(i,j,k,GDW) = qint.un;
            qgdnv(i,j,k,GDU) = qint.ut;
            qgdnv(i,j,k,GDV) = qint.utt;
        }
#ifndef RADIATION
        qgdnv(i,j,k,GDPRES) = qint.p;
#else
        qgdnv(i,j,k,GDPRES) = qint.p_g;
        for (int g = 0; g < NGROUPS; g++) {
            qgdnv(i,j,k,GDLAMS+g) = qint.lam[g];
            qgdnv(i,j,k,GDERADS+g) = qint.er[g];
        }
#endif

    }

}


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
riemanncg(const RiemannState& ql, const RiemannState& qr, const RiemannAux& raux,
          RiemannState& qint,
          const int idir) {

  // this implements the approximate Riemann solver of Colella & Glaz
  // (1985)
  //

  constexpr Real weakwv = 1.e-3_rt;

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
      qint.ut = ql.ut;
      qint.utt = ql.utt;
  } else if (ustar < 0.0_rt) {
      qint.ut = qr.ut;
      qint.utt = qr.utt;
  } else {
      qint.ut = 0.5_rt * (ql.ut + qr.ut);
      qint.utt = 0.5_rt * (ql.utt + qr.utt);
  }

  // linearly interpolate between the star and normal state -- this covers the
  // case where we are inside the rarefaction fan.
  qint.rho = frac*rstar + (1.0_rt - frac)*ro;
  qint.un = frac*ustar + (1.0_rt - frac)*uo;
  qint.p = frac*pstar + (1.0_rt - frac)*po;
  Real game_int = frac*gamstar + (1.0_rt-frac)*gameo;

  // now handle the cases where instead we are fully in the
  // star or fully in the original (l/r) state
  if (spout < 0.0_rt) {
      qint.rho = ro;
      qint.un = uo;
      qint.p = po;
      game_int = gameo;
  }

  if (spin >= 0.0_rt) {
      qint.rho = rstar;
      qint.un = ustar;
      qint.p = pstar;
      game_int = gamstar;
  }

  qint.p = amrex::max(qint.p, small_pres);

  qint.un = qint.un * raux.bnd_fac;

  // compute the total energy from the internal, p/(gamma - 1), and the kinetic
  qint.rhoe = qint.p / (game_int - 1.0_rt);

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
riemannus(const RiemannState& ql, const RiemannState& qr, const RiemannAux& raux,
          RiemannState& qint,
          const int idir) {

  // Colella, Glaz, and Ferguson solver
  //
  // this is a 2-shock solver that uses a very simple approximation for the
  // star state, and carries an auxiliary jump condition for (rho e) to
  // deal with a real gas


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
  for (int g = 0; g < NGROUPS; g++) {
      qint.lam[g] = fp * ql.lam[g] + fm * qr.lam[g];
  }

  if (ustar == 0) {
      // harmonic average
      for (int g = 0; g < NGROUPS; g++) {
          qint.lam[g] = 2.0_rt * (ql.lam[g] * qr.lam[g]) / (ql.lam[g] + qr.lam[g] + 1.e-50_rt);
      }
  }

  Real po_g = fp * ql.p_g + fm * qr.p_g;
  Real reo_r[NGROUPS];
  Real po_r[NGROUPS];
  for (int g = 0; g < NGROUPS; g++) {
      reo_r[g] = fp * ql.er[g] + fm * qr.er[g];
      po_r[g] = qint.lam[g] * reo_r[g];
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

  qint.ut = fp * ql.ut + fm * qr.ut;
  qint.utt = fp * ql.utt + fm * qr.utt;

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

  qint.rho = frac*rstar + (1.0_rt - frac)*ro;
  qint.un = frac*ustar + (1.0_rt - frac)*uo;

#ifdef RADIATION
  Real pgdnv_t = frac*pstar + (1.0_rt - frac)*po;
  Real pgdnv_g = frac*pstar_g + (1.0_rt - frac)*po_g;
  Real regdnv_g = frac*estar_g + (1.0_rt - frac)*reo_g;
  Real regdnv_r[NGROUPS];
  for (int g = 0; g < NGROUPS; g++) {
      regdnv_r[g] = frac*estar_r[g] + (1.0_rt - frac)*reo_r[g];
  }
#else
  qint.p = frac*pstar + (1.0_rt - frac)*po;
  Real regdnv = frac*estar + (1.0_rt - frac)*reo;
#endif

  // as it stands now, we set things assuming that the rarefaction
  // spans the interface.  We overwrite that here depending on the
  // wave speeds

  // look at the speeds on either side of the remaining wave
  // to determine which region we are in
  if (spout < 0.0_rt) {
      // the l or r state is on the interface
      qint.rho = ro;
      qint.un = uo;
#ifdef RADIATION
      pgdnv_t = po;
      pgdnv_g = po_g;
      regdnv_g = reo_g;
      for (int g = 0; g < NGROUPS; g++) {
          regdnv_r[g] = reo_r[g];
      }
#else
      qint.p = po;
      regdnv = reo;
#endif
  }

  if (spin >= 0.0_rt) {
      // the star state is on the interface
      qint.rho = rstar;
      qint.un = ustar;
#ifdef RADIATION
      pgdnv_t = pstar;
      pgdnv_g = pstar_g;
      regdnv_g = estar_g;
      for (int g = 0; g < NGROUPS; g++) {
          regdnv_r[g] = estar_r[g];
      }
#else
      qint.p = pstar;
      regdnv = estar;
#endif
    }

#ifdef RADIATION
  for (int g = 0; g < NGROUPS; g++) {
      qint.er[g] = amrex::max(regdnv_r[g], 0.0_rt);
  }

  qint.p_g = pgdnv_g;
  qint.p = pgdnv_t;
  qint.rhoe_g = regdnv_g;

  qint.rhoe = regdnv_g;
  for (int g = 0; g < NGROUPS; g++) {
      qint.rhoe += regdnv_r[g];
  }

#else
  qint.p = amrex::max(qint.p, small_pres);
  qint.rhoe = regdnv;
#endif

  // Enforce that fluxes through a symmetry plane or wall are hard zero.
  qint.un = qint.un * raux.bnd_fac;

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


///
/// A HLLC Riemann solver for pure hydrodynamics
///
/// @param ql          the left interface state
/// @param qr          the right interface state
/// @param qaux_arr    the auxillary state
/// @param uflx        the flux through the interface
/// @param qint        an approximate Godunov state on the interface
/// @param idir        coordinate direction for the solve (0 = x, 1 = y, 2 = z)
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
HLLC(const int i, const int j, const int k, const int idir,
     Array4<Real const> const& ql,
     Array4<Real const> const& qr,
     Array4<Real const> const& qaux_arr,
     Array4<Real> const& uflx,
     Array4<Real> const& qgdnv, const bool store_full_state,
     const GeometryData& geom,
     const bool special_bnd_lo, const bool special_bnd_hi,
     GpuArray<int, 3> const& domlo, GpuArray<int, 3> const& domhi) {

    // this is an implementation of the HLLC solver described in Toro's
    // book.  it uses the simplest estimate of the wave speeds, since
    // those should work for a general EOS.  We also initially do the
    // CGF Riemann construction to get pstar and ustar, since we'll
    // need to know the pressure and velocity on the interface for the
    // pdV term in the internal energy update.

    const Real small = 1.e-8_rt;

    int iu;
    int sx, sy, sz;

    if (idir == 0) {
        iu = QU;
        sx = 1;
        sy = 0;
        sz = 0;

    } else if (idir == 1) {
        iu = QV;
        sx = 0;
        sy = 1;
        sz = 0;

    } else {
        iu = QW;
        sx = 0;
        sy = 0;
        sz = 1;

    }


    int coord = geom.Coord();

    // deal with hard walls
    Real bnd_fac = 1.0_rt;

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


    Real rl = amrex::max(ql(i,j,k,QRHO), small_dens);

    // pick left velocities based on direction
    Real ul  = ql(i,j,k,iu);

    Real pl = amrex::max(ql(i,j,k,QPRES), small_pres);

    Real rr = amrex::max(qr(i,j,k,QRHO), small_dens);

    // pick right velocities based on direction
    Real ur  = qr(i,j,k,iu);

    Real pr = amrex::max(qr(i,j,k,QPRES), small_pres);

    // now we essentially do the CGF solver to get p and u on the
    // interface, but we won't use these in any flux construction.
    Real csmall = amrex::max(small, amrex::max(small * qaux_arr(i,j,k,QC),
                                                small * qaux_arr(i-sx,j-sy,k-sz,QC)));
    Real cavg = 0.5_rt*(qaux_arr(i,j,k,QC) + qaux_arr(i-sx,j-sy,k-sz,QC));

    Real gamcl = qaux_arr(i-sx,j-sy,k-sz,QGAMC);
    Real gamcr = qaux_arr(i,j,k,QGAMC);

#ifdef TRUE_SDC
    if (use_reconstructed_gamma1 == 1) {
        gamcl = ql(i,j,k,QGC);
        gamcr = qr(i,j,k,QGC);
    }
#endif

    Real wsmall = small_dens*csmall;
    Real wl = amrex::max(wsmall, std::sqrt(std::abs(gamcl*pl*rl)));
    Real wr = amrex::max(wsmall, std::sqrt(std::abs(gamcr*pr*rr)));

    Real wwinv = 1.0_rt/(wl + wr);
    Real pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))*wwinv;
    Real ustar = ((wl*ul + wr*ur) + (pl - pr))*wwinv;

    pstar = amrex::max(pstar, small_pres);
    // for symmetry preservation, if ustar is really small, then we
    // set it to zero
    if (std::abs(ustar) < riemann_constants::smallu*0.5_rt*(std::abs(ul) + std::abs(ur))){
        ustar = 0.0_rt;
    }

    Real ro;
    Real uo;
    Real po;
    Real gamco;

    if (ustar > 0.0_rt) {
        ro = rl;
        uo = ul;
        po = pl;
        gamco = gamcl;

    } else if (ustar < 0.0_rt) {
        ro = rr;
        uo = ur;
        po = pr;
        gamco = gamcr;

    } else {
        ro = 0.5_rt*(rl + rr);
        uo = 0.5_rt*(ul + ur);
        po = 0.5_rt*(pl + pr);
        gamco = 0.5_rt*(gamcl + gamcr);
    }

    ro = amrex::max(small_dens, ro);

    Real roinv = 1.0_rt/ro;
    Real co = std::sqrt(std::abs(gamco*po*roinv));
    co = amrex::max(csmall, co);
    Real co2inv = 1.0_rt/(co*co);

    Real rstar = ro + (pstar - po)*co2inv;
    rstar = amrex::max(small_dens, rstar);

    Real cstar = std::sqrt(std::abs(gamco*pstar/rstar));
    cstar = max(cstar, csmall);

    Real sgnm = std::copysign(1.0_rt, ustar);
    Real spout = co - sgnm*uo;
    Real spin = cstar - sgnm*ustar;
    Real ushock = 0.5_rt*(spin + spout);

    if (pstar-po > 0.0_rt) {
        spin = ushock;
        spout = ushock;
    }

    Real scr = spout-spin;
    if (spout-spin == 0.0_rt) {
        scr = small*cavg;
    }

    Real frac = (1.0_rt + (spout + spin)/scr)*0.5_rt;
    frac = amrex::max(0.0_rt, amrex::min(1.0_rt, frac));

    Real qint[NQ] = {0.0};

    qint[QRHO] = frac*rstar + (1.0_rt - frac)*ro;
    qint[iu] = frac*ustar + (1.0_rt - frac)*uo;
    qint[QPRES] = frac*pstar + (1.0_rt - frac)*po;


    // now we do the HLLC construction

    // use the simplest estimates of the wave speeds
    Real S_l = amrex::min(ul - std::sqrt(gamcl*pl/rl), ur - std::sqrt(gamcr*pr/rr));
    Real S_r = amrex::max(ul + std::sqrt(gamcl*pl/rl), ur + std::sqrt(gamcr*pr/rr));

    // estimate of the contact speed -- this is Toro Eq. 10.8
    Real S_c = (pr - pl + rl*ul*(S_l - ul) - rr*ur*(S_r - ur))/
        (rl*(S_l - ul) - rr*(S_r - ur));

    Real q_zone[NQ];
    Real U_state[NUM_STATE];
    Real U_hllc_state[NUM_STATE];
    Real F_state[NUM_STATE];

    if (S_r <= 0.0_rt) {
        // R region
        for (int n = 0; n < NQ; n++) {
            q_zone[n] = qr(i,j,k,n);
        }
        cons_state(q_zone, U_state);
        compute_flux(idir, bnd_fac, coord,
                     U_state, pr, F_state);

    } else if (S_r > 0.0_rt && S_c <= 0.0_rt) {
        // R* region
        for (int n = 0; n < NQ; n++) {
            q_zone[n] = qr(i,j,k,n);
        }
        cons_state(q_zone, U_state);
        compute_flux(idir, bnd_fac, coord,
                     U_state, pr, F_state);
        HLLC_state(idir, S_r, S_c, q_zone, U_hllc_state);

        // correct the flux
        for (int n = 0; n < NUM_STATE; n++) {
            F_state[n] = F_state[n] + S_r*(U_hllc_state[n] - U_state[n]);
        }

    } else if (S_c > 0.0_rt && S_l < 0.0_rt) {
        // L* region
        for (int n = 0; n < NQ; n++) {
            q_zone[n] = ql(i,j,k,n);
        }
        cons_state(q_zone, U_state);
        compute_flux(idir, bnd_fac, coord,
                     U_state, pl, F_state);
        HLLC_state(idir, S_l, S_c, q_zone, U_hllc_state);

        // correct the flux
        for (int n = 0; n < NUM_STATE; n++) {
            F_state[n] = F_state[n] + S_l*(U_hllc_state[n] - U_state[n]);
        }

    } else {
        // L region
        for (int n = 0; n < NQ; n++) {
            q_zone[n] = ql(i,j,k,n);
        }
        cons_state(q_zone, U_state);
        compute_flux(idir, bnd_fac, coord,
                     U_state, pl, F_state);
    }

    for (int n = 0; n < NUM_STATE; n++) {
        uflx(i,j,k,n) = F_state[n];
    }


    // now store the state

    if (store_full_state) {

        for (int n = 0; n < NQ; n++) {
            qgdnv(i,j,k,n) = qint[n];
        }

    } else {

        qgdnv(i,j,k,GDRHO) = qint[QRHO];
        qgdnv(i,j,k,GDU) = qint[QU];
        qgdnv(i,j,k,GDV) = qint[QV];
        qgdnv(i,j,k,GDW) = qint[QW];
        qgdnv(i,j,k,GDPRES) = qint[QPRES];

    }


}



AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
riemann_state(const int i, const int j, const int k, const int idir,
              Array4<Real> const& qm,
              Array4<Real> const& qp,
              Array4<Real const> const& qaux_arr,
              RiemannState& qint,
              const GeometryData& geom,
              const bool special_bnd_lo, const bool special_bnd_hi,
              GpuArray<int, 3> const& domlo, GpuArray<int, 3> const& domhi) {

  // just compute the hydrodynamic state on the interfaces
  // don't compute the fluxes

  // note: bx is not necessarily the limits of the valid (no ghost
  // cells) domain, but could be hi+1 in some dimensions.  We rely on
  // the caller to specify the interfaces over which to solve the
  // Riemann problems


  if (ppm_temp_fix == 2) {
      // recompute the thermodynamics on the interface to make it
      // all consistent

      // we want to take the edge states of rho, e, and X, and get
      // new values for p on the edges that are
      // thermodynamically consistent.

      eos_rep_t eos_state;

      // this is an initial guess for iterations, since we
      // can't be certain what temp is on interfaces
      eos_state.T = T_guess;

      // minus state
      eos_state.rho = qm(i,j,k,QRHO);
      eos_state.e = qm(i,j,k,QREINT)/qm(i,j,k,QRHO);
      for (int n = 0; n < NumSpec; n++) {
          eos_state.xn[n] = qm(i,j,k,QFS+n);
      }
#if NAUX_NET > 0
      for (int n = 0; n < NumAux; n++) {
          eos_state.aux[n] = qm(i,j,k,QFX+n);
      }
#endif

      eos(eos_input_re, eos_state);

      qm(i,j,k,QREINT) = eos_state.e * eos_state.rho;
      qm(i,j,k,QPRES) = eos_state.p;

      // plus state
      eos_state.rho = qp(i,j,k,QRHO);
      eos_state.e = qp(i,j,k,QREINT)/qp(i,j,k,QRHO);
      for (int n = 0; n < NumSpec; n++) {
          eos_state.xn[n] = qp(i,j,k,QFS+n);
      }
#if NAUX_NET > 0
      for (int n = 0; n < NumAux; n++) {
          eos_state.aux[n] = qp(i,j,k,QFX+n);
      }
#endif

      eos(eos_input_re, eos_state);

      qp(i,j,k,QREINT) = eos_state.e * eos_state.rho;
      qp(i,j,k,QPRES) = eos_state.p;
  }


  RiemannState ql;
  RiemannState qr;
  RiemannAux raux;


  load_input_states(i, j, k, idir,
                    qm, qp, qaux_arr,
                    ql, qr, raux);

  // deal with hard walls
  raux.bnd_fac = 1.0_rt;

  if (idir == 0) {
      if ((i == domlo[0] && special_bnd_lo) ||
          (i == domhi[0]+1 && special_bnd_hi)) {
          raux.bnd_fac = 0.0_rt;
      }

  } else if (idir == 1) {
      if ((j == domlo[1] && special_bnd_lo) ||
          (j == domhi[1]+1 && special_bnd_hi)) {
          raux.bnd_fac = 0.0_rt;
      }

  } else {
      if ((k == domlo[2] && special_bnd_lo) ||
          (k == domhi[2]+1 && special_bnd_hi)) {
          raux.bnd_fac = 0.0_rt;
      }
  }


  // Solve Riemann problem
  if (riemann_solver == 0) {
      // Colella, Glaz, & Ferguson solver

      riemannus(ql, qr, raux,
                qint,
                idir);

  } else if (riemann_solver == 1) {
      // Colella & Glaz solver

#ifndef RADIATION
      riemanncg(ql, qr, raux,
                qint,
                idir);
#endif

#ifndef AMREX_USE_GPU
  } else {
      amrex::Error("ERROR: invalid value of riemann_solver");
#endif
  }


}


#endif
