#ifndef RIEMANN_2SHOCK_SOLVERS_H
#define RIEMANN_2SHOCK_SOLVERS_H

#include <AMReX_REAL.H>

using namespace amrex::literals;

#include <riemann_type.H>
#ifdef RADIATION
#include <Radiation.H>
#include <rad_util.H>
#endif

namespace TwoShock {

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void
    wsqge(const amrex::Real p, const amrex::Real v,
          const amrex::Real gam, const amrex::Real gdot, amrex::Real& gstar,
          const amrex::Real gmin, const amrex::Real gmax, const amrex::Real csq,
          const amrex::Real pstar, amrex::Real& wsq) {

      // compute the lagrangian wave speeds -- this is the approximate
      // version for the Colella & Glaz algorithm


      // First predict a value of game across the shock

      // CG Eq. 31
      gstar = (pstar-p)*gdot/(pstar+p) + gam;
      gstar = amrex::max(gmin, amrex::min(gmax, gstar));

      // Now use that predicted value of game with the R-H jump conditions
      // to compute the wave speed.

      // this is CG Eq. 34
      amrex::Real alpha = pstar - (gstar - 1.0_rt)*p/(gam - 1.0_rt);
      if (alpha == 0.0_rt) {
        alpha = riemann_constants::smlp1*(pstar + p);
      }

      amrex::Real beta = pstar + 0.5_rt*(gstar - 1.0_rt)*(pstar+p);

      wsq = (pstar-p)*beta/(v*alpha);

      if (std::abs(pstar - p) < riemann_constants::smlp1*(pstar + p)) {
        wsq = csq;
      }
      wsq = amrex::max(wsq, (0.5_rt * (gam - 1.0_rt)/gam)*csq);
    }


    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void
    pstar_bisection(amrex::Real& pstar_lo, amrex::Real& pstar_hi,
                    const amrex::Real ul, const amrex::Real pl, const amrex::Real taul,
                    const amrex::Real gamel, const amrex::Real clsql,
                    const amrex::Real ur, const amrex::Real pr, const amrex::Real taur,
                    const amrex::Real gamer, const amrex::Real clsqr,
                    const amrex::Real gdot, const amrex::Real gmin, const amrex::Real gmax,
                    const int lriemann_shock_maxiter, const amrex::Real lriemann_pstar_tol,
                    amrex::Real& pstar, amrex::Real& gamstar,
                    bool& converged, amrex::GpuArray<amrex::Real, riemann_constants::PSTAR_BISECT_FACTOR * riemann_constants::HISTORY_SIZE>& pstar_hist_extra) {

      // we want to zero
      // f(p*) = u*_l(p*) - u*_r(p*)
      // we'll do bisection
      //
      // this version is for the approximate Colella & Glaz
      // version


      // lo bounds
      amrex::Real wlsq = 0.0;
      wsqge(pl, taul, gamel, gdot,
             gamstar, gmin, gmax, clsql, pstar_lo, wlsq);

      amrex::Real wrsq = 0.0;
      wsqge(pr, taur, gamer, gdot,
             gamstar, gmin, gmax, clsqr, pstar_lo, wrsq);

      amrex::Real wl = 1.0_rt / std::sqrt(wlsq);
      amrex::Real wr = 1.0_rt / std::sqrt(wrsq);

      amrex::Real ustar_l = ul - (pstar_lo - pstar)*wl;
      amrex::Real ustar_r = ur + (pstar_lo - pstar)*wr;

      amrex::Real f_lo = ustar_l - ustar_r;

      // hi bounds
      wsqge(pl, taul, gamel, gdot,
            gamstar, gmin, gmax, clsql, pstar_hi, wlsq);

      wsqge(pr, taur, gamer, gdot,
            gamstar, gmin, gmax, clsqr, pstar_hi, wrsq);

      // wl = 1.0_rt / std::sqrt(wlsq);
      // wr = 1.0_rt / std::sqrt(wrsq);

      //ustar_l = ul - (pstar_hi - pstar)*wl;
      //ustar_r = ur + (pstar_hi - pstar)*wr;

      //amrex::Real f_hi = ustar_l - ustar_r;

      // bisection
      converged = false;
      amrex::Real pstar_c = 0.0;

      for (int iter = 0; iter < riemann_constants::PSTAR_BISECT_FACTOR * lriemann_shock_maxiter; iter++) {

        pstar_c = 0.5_rt * (pstar_lo + pstar_hi);
        pstar_hist_extra[iter] = pstar_c;

        wsqge(pl, taul, gamel, gdot,
              gamstar, gmin, gmax, clsql, pstar_c, wlsq);

        wsqge(pr, taur, gamer, gdot,
              gamstar, gmin, gmax, clsqr, pstar_c, wrsq);

        wl = 1.0_rt / std::sqrt(wlsq);
        wr = 1.0_rt / std::sqrt(wrsq);

        ustar_l = ul - (pstar_c - pl)*wl;
        ustar_r = ur - (pstar_c - pr)*wr;

        amrex::Real f_c = ustar_l - ustar_r;

        if ( 0.5_rt * std::abs(pstar_lo - pstar_hi) < lriemann_pstar_tol * pstar_c ) {
          converged = true;
          break;
        }

        if (f_lo * f_c < 0.0_rt) {
          // root is in the left half
          pstar_hi = pstar_c;
          //f_hi = f_c;
        } else {
          pstar_lo = pstar_c;
          f_lo = f_c;
        }
      }

      pstar = pstar_c;
    }


    ///
    /// The Colella-Glaz Riemann solver for pure hydrodynamics.  This is a
    /// two shock approximate state Riemann solver.
    ///
    /// @param bx         the box to operate over
    /// @param ql         the left interface state
    /// @param qr         the right interface state
    /// @param qaux_arr   the auxiliary state
    /// @param qint       the full Godunov state on the interface
    /// @param idir       coordinate direction for the solve (0 = x, 1 = y, 2 = z)
    ///
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void
    riemanncg(const RiemannState& ql, const RiemannState& qr, const RiemannAux& raux,
              RiemannState& qint) {

      // this implements the approximate Riemann solver of Colella & Glaz
      // (1985)
      //

      constexpr amrex::Real weakwv = 1.e-3_rt;

    #ifndef AMREX_USE_GPU
      amrex::GpuArray<amrex::Real, riemann_constants::HISTORY_SIZE> pstar_hist;
    #endif


      // common quantities
      amrex::Real taul = 1.0_rt / ql.rho;
      amrex::Real taur = 1.0_rt / qr.rho;

      // lagrangian sound speeds
      amrex::Real clsql = ql.gamc * ql.p * ql.rho;
      amrex::Real clsqr = qr.gamc * qr.p * qr.rho;

      // Note: in the original Colella & Glaz paper, they predicted
      // gamma_e to the interfaces using a special (non-hyperbolic)
      // evolution equation.  In Castro, we instead bring (rho e)
      // to the edges, so we construct the necessary gamma_e here from
      // what we have on the interfaces.
      amrex::Real gamel = ql.p / ql.rhoe + 1.0_rt;
      amrex::Real gamer = qr.p / qr.rhoe + 1.0_rt;

      // these should consider a wider average of the cell-centered
      // gammas
      amrex::Real gmin = amrex::min(amrex::min(gamel, gamer), 1.0_rt);
      amrex::Real gmax = amrex::max(amrex::max(gamel, gamer), 2.0_rt);

      amrex::Real game_bar = 0.5_rt*(gamel + gamer);
      amrex::Real gamc_bar = 0.5_rt*(ql.gamc + qr.gamc);

      amrex::Real gdot = 2.0_rt*(1.0_rt - game_bar/gamc_bar)*(game_bar - 1.0_rt);

      amrex::Real wsmall = small_dens * raux.csmall;
      amrex::Real wl = amrex::max(wsmall, std::sqrt(std::abs(clsql)));
      amrex::Real wr = amrex::max(wsmall, std::sqrt(std::abs(clsqr)));

      // make an initial guess for pstar -- this is a two-shock
      // approximation
      //pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
      amrex::Real pstar = ql.p + ( (qr.p - ql.p) - wr*(qr.un - ql.un) ) * wl / (wl + wr);
      pstar = amrex::max(pstar, small_pres);

      // get the shock speeds -- this computes W_s from CG Eq. 34
      amrex::Real gamstar = 0.0;

      amrex::Real wlsq = 0.0;
      wsqge(ql.p, taul, gamel, gdot, gamstar,
            gmin, gmax, clsql, pstar, wlsq);

      amrex::Real wrsq = 0.0;
      wsqge(qr.p, taur, gamer, gdot, gamstar,
            gmin, gmax, clsqr, pstar, wrsq);

      amrex::Real pstar_old = pstar;

      wl = std::sqrt(wlsq);
      wr = std::sqrt(wrsq);

      // R-H jump conditions give ustar across each wave -- these
      // should be equal when we are done iterating.  Our notation
      // here is a little funny, comparing to CG, ustar_l = u*_L and
      // ustar_r = u*_R.
      amrex::Real ustar_l = ql.un - (pstar - ql.p) / wl;
      amrex::Real ustar_r = qr.un + (pstar - qr.p) / wr;

      // revise our pstar guess
      // pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
      pstar = ql.p + ( (qr.p - ql.p) - wr * (qr.un - ql.un) ) * wl / (wl + wr);
      pstar = amrex::max(pstar, small_pres);

      // secant iteration
      bool converged = false;

      int iter = 0;
      while ((iter < riemann_shock_maxiter && !converged) || iter < 2) {

          wsqge(ql.p, taul, gamel, gdot, gamstar,
                gmin, gmax, clsql, pstar, wlsq);

          wsqge(qr.p, taur, gamer, gdot, gamstar,
                gmin, gmax, clsqr, pstar, wrsq);


          // NOTE: these are really the inverses of the wave speeds!
          wl = 1.0_rt / std::sqrt(wlsq);
          wr = 1.0_rt / std::sqrt(wrsq);

          amrex::Real ustar_r_old = ustar_r;
          amrex::Real ustar_l_old = ustar_l;

          ustar_r = qr.un - (qr.p - pstar) * wr;
          ustar_l = ql.un + (ql.p - pstar) * wl;

          amrex::Real dpditer = std::abs(pstar_old - pstar);

          // Here we are going to do the Secant iteration version in
          // CG.  Note that what we call zp and zm here are not
          // actually the Z_p = |dp*/du*_p| defined in CG, by rather
          // simply |du*_p| (or something that looks like dp/Z!).
          amrex::Real zp = std::abs(ustar_l - ustar_l_old);
          if (zp - weakwv * raux.cavg <= 0.0_rt) {
              zp = dpditer * wl;
          }

          amrex::Real zm = std::abs(ustar_r - ustar_r_old);
          if (zm - weakwv * raux.cavg <= 0.0_rt) {
              zm = dpditer * wr;
          }

          // the new pstar is found via CG Eq. 18
          amrex::Real denom = dpditer / amrex::max(zp + zm, riemann_constants::small * raux.cavg);
          pstar_old = pstar;
          pstar = pstar - denom*(ustar_r - ustar_l);
          pstar = amrex::max(pstar, small_pres);

          amrex::Real err = std::abs(pstar - pstar_old);
          if (err < riemann_pstar_tol * pstar) {
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

          if (riemann_cg_blend == 0) {

    #ifndef AMREX_USE_GPU
              std::cout <<  "pstar history: " << std::endl;
              for (int iter_l=0; iter_l < riemann_shock_maxiter; iter_l++) {
                  std::cout << iter_l << " " << pstar_hist[iter_l] << std::endl;
              }

              std::cout << std::endl;
              std::cout << "left state: " << std::endl << ql << std::endl;
              std::cout << "right state: " << std::endl << qr << std::endl;
              std::cout << "aux information: " << std::endl << raux << std::endl;

              amrex::Error("ERROR: non-convergence in the Riemann solver");
    #endif

          } else if (riemann_cg_blend == 1) {

              pstar = ql.p + ( (qr.p - ql.p) - wr * (qr.un - ql.un) ) * wl / (wl + wr);

          } else if (riemann_cg_blend == 2) {

              // we don't store the history if we are in CUDA, so
              // we can't do this
    #ifndef AMREX_USE_GPU
              // first try to find a reasonable bounds
              amrex::Real pstarl = 1.e200;
              amrex::Real pstaru = -1.e200;
              for (int n = riemann_shock_maxiter-6; n < riemann_shock_maxiter; n++) {
                  pstarl = amrex::min(pstarl, pstar_hist[n]);
                  pstaru = amrex::max(pstaru, pstar_hist[n]);
              }

              pstarl = amrex::max(pstarl, small_pres);
              pstaru = amrex::max(pstaru, small_pres);

              amrex::GpuArray<amrex::Real, riemann_constants::PSTAR_BISECT_FACTOR * riemann_constants::HISTORY_SIZE> pstar_hist_extra;

              pstar_bisection(pstarl, pstaru,
                              ql.un, ql.p, taul, gamel, clsql,
                              qr.un, qr.p, taur, gamer, clsqr,
                              gdot, gmin, gmax,
                              riemann_shock_maxiter, riemann_pstar_tol,
                              pstar, gamstar, converged, pstar_hist_extra);

              if (!converged) {

                  std::cout << "pstar history: " << std::endl;
                  for (int iter_l = 0; iter_l < riemann_shock_maxiter; iter_l++) {
                      std::cout << iter_l << " " << pstar_hist[iter_l] << std::endl;
                  }
                  std::cout << "pstar extra history: " << std::endl;
                  for (int iter_l = 0; iter_l < riemann_constants::PSTAR_BISECT_FACTOR * riemann_shock_maxiter; iter_l++) {
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
              amrex::Error("ERROR: unrecognized riemann_cg_blend option.");
    #endif
          }

      }

      // we converged!  construct the single ustar for the region
      // between the left and right waves, using the updated wave speeds
      ustar_r = qr.un - (qr.p - pstar) * wr;  // careful -- here wl, wr are 1/W
      ustar_l = ql.un + (ql.p - pstar) * wl;

      amrex::Real ustar = 0.5_rt * (ustar_l + ustar_r);

      // for symmetry preservation, if ustar is really small, then we
      // set it to zero
      if (std::abs(ustar) < riemann_constants::smallu * 0.5_rt * (std::abs(ql.un) + std::abs(qr.un))) {
          ustar = 0.0_rt;
      }

      // sample the solution -- here we look first at the direction
      // that the contact is moving.  This tells us if we need to
      // worry about the L/L* states or the R*/R states.
      amrex::Real ro;
      amrex::Real uo;
      amrex::Real po;
      amrex::Real tauo;
      amrex::Real gamco;
      amrex::Real gameo;

      if (ustar > 0.0_rt) {
          //ro = ql.rho;
          uo = ql.un;
          po = ql.p;
          tauo = taul;
          gamco = ql.gamc;
          gameo = gamel;

      } else if (ustar < 0.0_rt) {
          //ro = qr.rho;
          uo = qr.un;
          po = qr.p;
          tauo = taur;
          gamco = qr.gamc;
          gameo = gamer;

      } else {
          //ro = 0.5_rt * (ql.rho + qr.rho);
          uo = 0.5_rt * (ql.un + qr.un);
          po = 0.5_rt * (ql.p + qr.p);
          tauo = 0.5_rt * (taul + taur);
          gamco = 0.5_rt * (ql.gamc + qr.gamc);
          gameo = 0.5_rt * (gamel + gamer);
      }

      // use tau = 1/rho as the independent variable here
      ro = amrex::max(small_dens, 1.0_rt/tauo);
      tauo = 1.0_rt/ro;

      amrex::Real co = std::sqrt(std::abs(gamco*po*tauo));
      co = amrex::max(raux.csmall, co);
      amrex::Real clsq = std::pow(co*ro, 2);

      // now that we know which state (left or right) we need to worry
      // about, get the value of gamstar and wosq across the wave we
      // are dealing with.
      amrex::Real wosq = 0.0;
      wsqge(po, tauo, gameo, gdot, gamstar,
            gmin, gmax, clsq, pstar, wosq);

      amrex::Real sgnm = std::copysign(1.0_rt, ustar);

      amrex::Real wo = std::sqrt(wosq);
      amrex::Real dpjmp = pstar - po;

      // is this max really necessary?
      //rstar=max(ONE-ro*dpjmp/wosq, (gameo-ONE)/(gameo+ONE))
      amrex::Real rstar = 1.0_rt - ro*dpjmp/wosq;
      rstar = ro/rstar;
      rstar = amrex::max(small_dens, rstar);

      amrex::Real cstar = std::sqrt(std::abs(gamco * pstar / rstar));
      cstar = amrex::max(cstar, raux.csmall);

      amrex::Real spout = co - sgnm*uo;
      amrex::Real spin = cstar - sgnm*ustar;

      //ushock = 0.5_rt*(spin + spout)
      amrex::Real ushock = wo*tauo - sgnm*uo;

      if (pstar - po >= 0.0_rt) {
          spin = ushock;
          spout = ushock;
      }

      amrex::Real frac = 0.5_rt*(1.0_rt + (spin + spout)/amrex::max(amrex::max(spout-spin, spin+spout),
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
      amrex::Real game_int = frac*gamstar + (1.0_rt-frac)*gameo;

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
    /// @param qaux_arr    the auxiliary state
    /// @param qint        the full Godunov state on the interface
    /// @param lambda_int  the radiation flux limiter on the interface
    ///
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void
    riemannus(const RiemannState& ql, const RiemannState& qr, const RiemannAux& raux,
              RiemannState& qint) {

      // Colella, Glaz, and Ferguson solver
      //
      // this is a 2-shock solver that uses a very simple approximation for the
      // star state, and carries an auxiliary jump condition for (rho e) to
      // deal with a real gas


      // estimate the star state: pstar, ustar

      amrex::Real wsmall = small_dens * raux.csmall;

      // this is Castro I: Eq. 33

      amrex::Real wl = amrex::max(wsmall, std::sqrt(std::abs(ql.gamc * ql.p * ql.rho)));
      amrex::Real wr = amrex::max(wsmall, std::sqrt(std::abs(qr.gamc * qr.p * qr.rho)));

      amrex::Real wwinv = 1.0_rt/(wl + wr);
      amrex::Real pstar = ((wr * ql.p + wl * qr.p) + wl * wr * (ql.un - qr.un)) * wwinv;
      amrex::Real ustar = ((wl * ql.un + wr * qr.un) + (ql.p - qr.p)) * wwinv;

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
      amrex::Real sgnm = std::copysign(1.0_rt, ustar);
      if (ustar == 0.0_rt) {
          sgnm = 0.0_rt;
      }

      amrex::Real fp = 0.5_rt*(1.0_rt + sgnm);
      amrex::Real fm = 0.5_rt*(1.0_rt - sgnm);

      amrex::Real ro = fp * ql.rho + fm * qr.rho;
      amrex::Real uo = fp * ql.un + fm * qr.un;
      amrex::Real po = fp * ql.p + fm * qr.p;
      amrex::Real reo = fp * ql.rhoe + fm * qr.rhoe;
      amrex::Real gamco = fp * ql.gamc + fm * qr.gamc;
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

      amrex::Real po_g = fp * ql.p_g + fm * qr.p_g;
      amrex::Real reo_r[NGROUPS];
      amrex::Real po_r[NGROUPS];
      for (int g = 0; g < NGROUPS; g++) {
          reo_r[g] = fp * ql.er[g] + fm * qr.er[g];
          po_r[g] = qint.lam[g] * reo_r[g];
      }
      amrex::Real reo_g = fp * ql.rhoe_g + fm * qr.rhoe_g;
      amrex::Real gamco_g = fp * ql.gamcg + fm * qr.gamcg;
    #endif

      ro = amrex::max(small_dens, ro);

      amrex::Real roinv = 1.0_rt / ro;

      amrex::Real co = std::sqrt(std::abs(gamco * po * roinv));
      co = amrex::max(raux.csmall, co);
      amrex::Real co2inv = 1.0_rt / (co*co);

      // we can already deal with the transverse velocities -- they
      // only jump across the contact

      qint.ut = fp * ql.ut + fm * qr.ut;
      qint.utt = fp * ql.utt + fm * qr.utt;

      // compute the rest of the star state

      amrex::Real drho = (pstar - po)*co2inv;
      amrex::Real rstar = ro + drho;
      rstar = amrex::max(small_dens, rstar);

    #ifdef RADIATION
      amrex::Real estar_g = reo_g + drho*(reo_g + po_g)*roinv;

      amrex::Real co_g = std::sqrt(std::abs(gamco_g*po_g*roinv));
      co_g = amrex::max(raux.csmall, co_g);

      amrex::Real pstar_g = po_g + drho*co_g*co_g;
      pstar_g = amrex::max(pstar_g, small_pres);

      amrex::Real estar_r[NGROUPS];
      for (int g = 0; g < NGROUPS; g++) {
          estar_r[g] = reo_r[g] + drho*(reo_r[g] + po_r[g])*roinv;
      }
    #else
      amrex::Real entho = (reo + po)*roinv*co2inv;
      amrex::Real estar = reo + (pstar - po)*entho;
    #endif

      amrex::Real cstar = std::sqrt(std::abs(gamco*pstar/rstar));
      cstar = amrex::max(cstar, raux.csmall);

      // finish sampling the solution

      // look at the remaining wave to determine if the star state or the
      // 'o' state above is on the interface

      // the values of u +/- c on either side of the non-contact wave
      amrex::Real spout = co - sgnm*uo;
      amrex::Real spin = cstar - sgnm*ustar;

      // a simple estimate of the shock speed
      amrex::Real ushock = 0.5_rt*(spin + spout);

      if (pstar-po > 0.0_rt) {
          spin = ushock;
          spout = ushock;
      }

      amrex::Real scr = spout - spin;
      if (spout-spin == 0.0_rt) {
          scr = riemann_constants::small * raux.cavg;
      }

      // interpolate for the case that we are in a rarefaction
      amrex::Real frac = (1.0_rt + (spout + spin)/scr)*0.5_rt;
      frac = amrex::max(0.0_rt, amrex::min(1.0_rt, frac));

      qint.rho = frac*rstar + (1.0_rt - frac)*ro;
      qint.un = frac*ustar + (1.0_rt - frac)*uo;

    #ifdef RADIATION
      amrex::Real pgdnv_t = frac*pstar + (1.0_rt - frac)*po;
      amrex::Real pgdnv_g = frac*pstar_g + (1.0_rt - frac)*po_g;
      amrex::Real regdnv_g = frac*estar_g + (1.0_rt - frac)*reo_g;
      amrex::Real regdnv_r[NGROUPS];
      for (int g = 0; g < NGROUPS; g++) {
          regdnv_r[g] = frac*estar_r[g] + (1.0_rt - frac)*reo_r[g];
      }
    #else
      qint.p = frac*pstar + (1.0_rt - frac)*po;
      amrex::Real regdnv = frac*estar + (1.0_rt - frac)*reo;
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
}
#endif
