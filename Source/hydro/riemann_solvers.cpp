#include <Castro.H>
#include <Castro_F.H>
#include <Castro_util.H>

#ifdef RADIATION
#include <Radiation.H>
#endif

#include <cmath>

#include <eos.H>
#include <riemann.H>

using namespace amrex;

void
Castro::riemanncg(const Box& bx,
                  Array4<Real> const& qleft_arr,
                  Array4<Real> const& qright_arr,
                  Array4<Real const> const& qaux_arr,
                  Array4<Real> const& qint,
                  const int idir) {

  // this implements the approximate Riemann solver of Colella & Glaz
  // (1985)
  //

  constexpr Real weakwv = 1.e-3_rt;

#ifndef AMREX_USE_GPU
  if (cg_maxiter > HISTORY_SIZE) {
    amrex::Error("error in riemanncg: cg_maxiter > HISTORY_SIZE");
  }
#endif

#ifndef AMREX_USE_GPU
  if (cg_blend == 2 && cg_maxiter < 5) {
    amrex::Error("Error: need cg_maxiter >= 5 to do a bisection search on secant iteration failure.");
  }
#endif

  const auto domlo = geom.Domain().loVect3d();
  const auto domhi = geom.Domain().hiVect3d();

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


  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();

  // do we want to force the flux to zero at the boundary?
  const bool special_bnd_lo = (lo_bc[idir] == Symmetry ||
                               lo_bc[idir] == SlipWall ||
                               lo_bc[idir] == NoSlipWall);
  const bool special_bnd_hi = (hi_bc[idir] == Symmetry ||
                               hi_bc[idir] == SlipWall ||
                               hi_bc[idir] == NoSlipWall);

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
  {

#ifndef AMREX_USE_GPU
    GpuArray<Real, HISTORY_SIZE> pstar_hist;
#endif


    RiemannState ql;
    RiemannState qr;
    RiemannAux raux;

    int compute_gammas = 0;

    load_input_states(i, j, k, idir,
                      qleft_arr, qright_arr, qaux_arr,
                      compute_gammas,
                      ql, qr, raux);

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

    qint(i,j,k,iu) = qint(i,j,k,iu) * bnd_fac;

    // Compute fluxes, order as conserved state (not q)

    // compute the total energy from the internal, p/(gamma - 1), and the kinetic
    qint(i,j,k,QREINT) = qint(i,j,k,QPRES)/(game_int - 1.0_rt);

    // we'll do the passive scalars separately

  });

}

void
Castro::riemannus(const Box& bx,
                  Array4<Real> const& qleft_arr,
                  Array4<Real> const& qright_arr,
                  Array4<Real const> const& qaux_arr,
                  Array4<Real> const& qint,
#ifdef RADIATION
                  Array4<Real> const& lambda_int,
#endif
                  const int idir, const int compute_gammas) {

  // Colella, Glaz, and Ferguson solver
  //
  // this is a 2-shock solver that uses a very simple approximation for the
  // star state, and carries an auxiliary jump condition for (rho e) to
  // deal with a real gas

  // set integer pointers for the normal and transverse velocity and
  // momentum

  const auto domlo = geom.Domain().loVect3d();
  const auto domhi = geom.Domain().hiVect3d();

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

  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();

  // do we want to force the flux to zero at the boundary?
  const bool special_bnd_lo = (lo_bc[idir] == Symmetry ||
                               lo_bc[idir] == SlipWall ||
                               lo_bc[idir] == NoSlipWall);
  const bool special_bnd_hi = (hi_bc[idir] == Symmetry ||
                               hi_bc[idir] == SlipWall ||
                               hi_bc[idir] == NoSlipWall);

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
  {

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


    // set the left and right states for this interface

#ifdef RADIATION
    Real laml[NGROUPS];
    Real lamr[NGROUPS];
#endif

    RiemannState ql;
    RiemannState qr;
    RiemannAux raux;

    load_input_states(i, j, k, idir,
                      qleft_arr, qright_arr, qaux_arr,
                      compute_gammas,
                      ql, qr, raux);



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
    qint(i,j,k,iu) = qint(i,j,k,iu) * bnd_fac;


    // we'll do the passive scalars separately


  });
}


void
Castro::HLLC(const Box& bx,
             Array4<Real const> const& ql,
             Array4<Real const> const& qr,
             Array4<Real const> const& qaux_arr,
             Array4<Real> const& uflx,
             Array4<Real> const& qint,
             const int idir) {

  // this is an implementation of the HLLC solver described in Toro's
  // book.  it uses the simplest estimate of the wave speeds, since
  // those should work for a general EOS.  We also initially do the
  // CGF Riemann construction to get pstar and ustar, since we'll
  // need to know the pressure and velocity on the interface for the
  // pdV term in the internal energy update.

  const auto domlo = geom.Domain().loVect3d();
  const auto domhi = geom.Domain().hiVect3d();

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

  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();

  // do we want to force the flux to zero at the boundary?
  const bool special_bnd_lo = (lo_bc[idir] == Symmetry ||
                               lo_bc[idir] == SlipWall ||
                               lo_bc[idir] == NoSlipWall);
  const bool special_bnd_hi = (hi_bc[idir] == Symmetry ||
                               hi_bc[idir] == SlipWall ||
                               hi_bc[idir] == NoSlipWall);

  const Real lsmall_dens = small_dens;
  const Real lsmall_pres = small_pres;
  const Real lsmall = riemann_constants::small;

  int coord = geom.Coord();

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
  {

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


    Real rl = amrex::max(ql(i,j,k,QRHO), lsmall_dens);

    // pick left velocities based on direction
    Real ul  = ql(i,j,k,iu);

    Real pl = amrex::max(ql(i,j,k,QPRES), lsmall_pres);

    Real rr = amrex::max(qr(i,j,k,QRHO), lsmall_dens);

    // pick right velocities based on direction
    Real ur  = qr(i,j,k,iu);

    Real pr = amrex::max(qr(i,j,k,QPRES), lsmall_pres);

    // now we essentially do the CGF solver to get p and u on the
    // interface, but we won't use these in any flux construction.
    Real csmall = amrex::max(lsmall, amrex::max(lsmall * qaux_arr(i,j,k,QC), lsmall * qaux_arr(i-sx,j-sy,k-sz,QC)));
    Real cavg = 0.5_rt*(qaux_arr(i,j,k,QC) + qaux_arr(i-sx,j-sy,k-sz,QC));

    Real gamcl = qaux_arr(i-sx,j-sy,k-sz,QGAMC);
    Real gamcr = qaux_arr(i,j,k,QGAMC);

#ifdef TRUE_SDC
    if (use_reconstructed_gamma1 == 1) {
      gamcl = ql(i,j,k,QGC);
      gamcr = qr(i,j,k,QGC);
    }
#endif

    Real wsmall = lsmall_dens*csmall;
    Real wl = amrex::max(wsmall, std::sqrt(std::abs(gamcl*pl*rl)));
    Real wr = amrex::max(wsmall, std::sqrt(std::abs(gamcr*pr*rr)));

    Real wwinv = 1.0_rt/(wl + wr);
    Real pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))*wwinv;
    Real ustar = ((wl*ul + wr*ur) + (pl - pr))*wwinv;

    pstar = amrex::max(pstar, lsmall_pres);
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

    ro = amrex::max(lsmall_dens, ro);

    Real roinv = 1.0_rt/ro;
    Real co = std::sqrt(std::abs(gamco*po*roinv));
    co = amrex::max(csmall, co);
    Real co2inv = 1.0_rt/(co*co);

    Real rstar = ro + (pstar - po)*co2inv;
    rstar = amrex::max(lsmall_dens, rstar);

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
      scr = lsmall*cavg;
    }

    Real frac = (1.0_rt + (spout + spin)/scr)*0.5_rt;
    frac = amrex::max(0.0_rt, amrex::min(1.0_rt, frac));

    qint(i,j,k,iu) = frac*ustar + (1.0_rt - frac)*uo;
    qint(i,j,k,QPRES) = frac*pstar + (1.0_rt - frac)*po;


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
  });
}
