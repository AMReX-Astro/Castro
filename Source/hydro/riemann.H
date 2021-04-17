#ifndef CASTRO_RIEMANN_H
#define CASTRO_RIEMANN_H

using namespace amrex;

namespace riemann_constants {
    const Real smlp1 = 1.e-10_rt;
    const Real small = 1.e-8_rt;
    const Real smallu = 1.e-12_rt;
}


struct RiemannState
{
    Real rho;
    Real p;
    Real rhoe;
    Real gamc;
    Real un;
    Real ut;
    Real utt;

#ifdef RADIATION
    Real lam[NGROUPS];
    Real er[NGROUPS];
    Real p_g;
    Real rhoe_g;
    Real gamcg;
#endif

};


struct RiemannAux
{
    Real csmall;
    Real cavg;
    Real bnd_fac;
};


AMREX_INLINE
std::ostream& operator<< (std::ostream& o, RiemannState const& rs)
{
    o << "rho  = " << rs.rho << std::endl;
    o << "p    = " << rs.p << std::endl;
    o << "rhoe = " << rs.rhoe << std::endl;
    o << "gamc = " << rs.gamc << std::endl;
    o << "un   = " << rs.un << std::endl;
    o << "ut   = " << rs.ut << std::endl;
    o << "utt  = " << rs.utt << std::endl;
    return o;
}


AMREX_INLINE
std::ostream& operator<< (std::ostream& o, RiemannAux const& ra)
{
    o << "csmall  = " << ra.csmall << std::endl;
    o << "cavg    = " << ra.cavg << std::endl;
    o << "bnd_fac = " << ra.bnd_fac << std::endl;
    return o;
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
load_input_states(const int i, const int j, const int k, const int idir,
                  Array4<Real const> const& qleft_arr,
                  Array4<Real const> const& qright_arr,
                  Array4<Real const> const& qaux_arr,
                  RiemannState& ql, RiemannState& qr, RiemannAux& raux) {

    const Real small = 1.e-8_rt;

    // fill the states directly from inputs

    ql.rho = amrex::max(qleft_arr(i,j,k,QRHO), small_dens);
    qr.rho = amrex::max(qright_arr(i,j,k,QRHO), small_dens);

    if (idir == 0) {
        ql.un = qleft_arr(i,j,k,QU);
        ql.ut = qleft_arr(i,j,k,QV);
        ql.utt = qleft_arr(i,j,k,QW);

        qr.un = qright_arr(i,j,k,QU);
        qr.ut = qright_arr(i,j,k,QV);
        qr.utt = qright_arr(i,j,k,QW);

        raux.csmall = amrex::max(small,
                                 small * amrex::max(qaux_arr(i,j,k,QC), qaux_arr(i-1,j,k,QC)));
        raux.cavg = 0.5_rt * (qaux_arr(i,j,k,QC) + qaux_arr(i-1,j,k,QC));

        ql.gamc = qaux_arr(i-1,j,k,QGAMC);
        qr.gamc = qaux_arr(i,j,k,QGAMC);

#ifdef RADIATION
        ql.gamcg = qaux_arr(i-1,j,k,QGAMCG);
        qr.gamcg = qaux_arr(i,j,k,QGAMCG);

        for (int g = 0; g < NGROUPS; g++) {
            ql.lam[g] = qaux_arr(i-1,j,k,QLAMS+g);
            qr.lam[g] = qaux_arr(i,j,k,QLAMS+g);
        }
#endif

    } else if (idir == 1) {
        ql.un = qleft_arr(i,j,k,QV);
        ql.ut = qleft_arr(i,j,k,QU);
        ql.utt = qleft_arr(i,j,k,QW);

        qr.un = qright_arr(i,j,k,QV);
        qr.ut = qright_arr(i,j,k,QU);
        qr.utt = qright_arr(i,j,k,QW);

        raux.csmall = amrex::max(small,
                                 small * amrex::max(qaux_arr(i,j,k,QC), qaux_arr(i,j-1,k,QC)));
        raux.cavg = 0.5_rt * (qaux_arr(i,j,k,QC) + qaux_arr(i,j-1,k,QC));

        ql.gamc = qaux_arr(i,j-1,k,QGAMC);
        qr.gamc = qaux_arr(i,j,k,QGAMC);

#ifdef RADIATION
        ql.gamcg = qaux_arr(i,j-1,k,QGAMCG);
        qr.gamcg = qaux_arr(i,j,k,QGAMCG);

        for (int g = 0; g < NGROUPS; g++) {
            ql.lam[g] = qaux_arr(i,j-1,k,QLAMS+g);
            qr.lam[g] = qaux_arr(i,j,k,QLAMS+g);
        }
#endif

    } else {
        ql.un = qleft_arr(i,j,k,QW);
        ql.ut = qleft_arr(i,j,k,QU);
        ql.utt = qleft_arr(i,j,k,QV);

        qr.un = qright_arr(i,j,k,QW);
        qr.ut = qright_arr(i,j,k,QU);
        qr.utt = qright_arr(i,j,k,QV);

        raux.csmall = amrex::max(small,
                                 small * amrex::max(qaux_arr(i,j,k,QC), qaux_arr(i,j,k-1,QC)));
        raux.cavg = 0.5_rt * (qaux_arr(i,j,k,QC) + qaux_arr(i,j,k-1,QC));

        ql.gamc = qaux_arr(i,j,k-1,QGAMC);
        qr.gamc = qaux_arr(i,j,k,QGAMC);

#ifdef RADIATION
        ql.gamcg = qaux_arr(i,j,k-1,QGAMCG);
        qr.gamcg = qaux_arr(i,j,k,QGAMCG);

        for (int g = 0; g < NGROUPS; g++) {
            ql.lam[g] = qaux_arr(i,j,k-1,QLAMS+g);
            qr.lam[g] = qaux_arr(i,j,k,QLAMS+g);
        }
#endif

    }

#ifndef RADIATION
    ql.p = qleft_arr(i,j,k,QPRES);
    qr.p = qright_arr(i,j,k,QPRES);

    ql.rhoe = qleft_arr(i,j,k,QREINT);
    qr.rhoe = qright_arr(i,j,k,QREINT);

#else
    ql.p = qleft_arr(i,j,k,QPTOT);
    ql.rhoe = qleft_arr(i,j,k,QREITOT);
    for (int g = 0; g < NGROUPS; g++) {
      ql.er[g] = qleft_arr(i,j,k,QRAD+g);
    }
    ql.p_g = qleft_arr(i,j,k,QPRES);
    ql.rhoe_g = qleft_arr(i,j,k,QREINT);

    qr.p = qright_arr(i,j,k,QPTOT);
    qr.rhoe = qright_arr(i,j,k,QREITOT);
    for (int g = 0; g < NGROUPS; g++) {
      qr.er[g] = qright_arr(i,j,k,QRAD+g);
    }
    qr.p_g = qright_arr(i,j,k,QPRES);
    qr.rhoe_g = qright_arr(i,j,k,QREINT);
#endif

    // cleaning

#ifdef TRUE_SDC
    if (use_reconstructed_gamma1 == 1) {
        ql.gamc = qleft_arr(i,j,k,QGC);
        qr.gamc = qright_arr(i,j,k,QGC);
    }
#endif

    // sometime we come in here with negative energy or pressure
    // note: reset both in either case, to remain thermo
    // consistent
    if (ql.rhoe <= 0.0_rt || ql.p < small_pres) {
#ifndef AMREX_USE_GPU
        std::cout <<  "WARNING: (rho e)_l < 0 or pl < small_pres in Riemann: " << ql.rhoe << " " << ql.p << " " << small_pres << std::endl;
#endif

        eos_rep_t eos_state;
        eos_state.T = small_temp;
        eos_state.rho = ql.rho;
        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = qleft_arr(i,j,k,QFS+n);
        }
#if NAUX_NET > 0
        for (int n = 0; n < NumAux; n++) {
            eos_state.aux[n] = qleft_arr(i,j,k,QFX+n);
        }
#endif

        eos(eos_input_rt, eos_state);

        ql.rhoe = ql.rho * eos_state.e;
        ql.p = eos_state.p;
        ql.gamc = eos_state.gam1;
    }

    if (qr.rhoe <= 0.0_rt || qr.p < small_pres) {
#ifndef AMREX_USE_GPU
        std::cout << "WARNING: (rho e)_r < 0 or pr < small_pres in Riemann: " << qr.rhoe << " " << qr.p << " " << small_pres << std::endl;
#endif

        eos_rep_t eos_state;
        eos_state.T = small_temp;
        eos_state.rho = qr.rho;
        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = qright_arr(i,j,k,QFS+n);
        }
#if NAUX_NET > 0
        for (int n = 0; n < NumAux; n++) {
            eos_state.aux[n] = qright_arr(i,j,k,QFX+n);
        }
#endif

        eos(eos_input_rt, eos_state);

        qr.rhoe = qr.rho * eos_state.e;
        qr.p = eos_state.p;
        qr.gamc = eos_state.gam1;
    }

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
wsqge(const Real p, const Real v,
      const Real gam, const Real gdot, Real& gstar,
      const Real gmin, const Real gmax, const Real csq,
      const Real pstar, Real& wsq) {

  // compute the lagrangian wave speeds -- this is the approximate
  // version for the Colella & Glaz algorithm


  // First predict a value of game across the shock

  // CG Eq. 31
  gstar = (pstar-p)*gdot/(pstar+p) + gam;
  gstar = amrex::max(gmin, amrex::min(gmax, gstar));

  // Now use that predicted value of game with the R-H jump conditions
  // to compute the wave speed.

  // this is CG Eq. 34
  Real alpha = pstar - (gstar - 1.0_rt)*p/(gam - 1.0_rt);
  if (alpha == 0.0_rt) {
    alpha = riemann_constants::smlp1*(pstar + p);
  }

  Real beta = pstar + 0.5_rt*(gstar - 1.0_rt)*(pstar+p);

  wsq = (pstar-p)*beta/(v*alpha);

  if (std::abs(pstar - p) < riemann_constants::smlp1*(pstar + p)) {
    wsq = csq;
  }
  wsq = amrex::max(wsq, (0.5_rt * (gam - 1.0_rt)/gam)*csq);
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
pstar_bisection(Real& pstar_lo, Real& pstar_hi,
                const Real ul, const Real pl, const Real taul,
                const Real gamel, const Real clsql,
                const Real ur, const Real pr, const Real taur,
                const Real gamer, const Real clsqr,
                const Real gdot, const Real gmin, const Real gmax,
                const int lcg_maxiter, const Real lcg_tol,
                Real& pstar, Real& gamstar,
                bool& converged, GpuArray<Real, PSTAR_BISECT_FACTOR*HISTORY_SIZE>& pstar_hist_extra) {

  // we want to zero
  // f(p*) = u*_l(p*) - u*_r(p*)
  // we'll do bisection
  //
  // this version is for the approximate Colella & Glaz
  // version


  // lo bounds
  Real wlsq = 0.0;
  wsqge(pl, taul, gamel, gdot,
         gamstar, gmin, gmax, clsql, pstar_lo, wlsq);

  Real wrsq = 0.0;
  wsqge(pr, taur, gamer, gdot,
         gamstar, gmin, gmax, clsqr, pstar_lo, wrsq);

  Real wl = 1.0_rt / std::sqrt(wlsq);
  Real wr = 1.0_rt / std::sqrt(wrsq);

  Real ustar_l = ul - (pstar_lo - pstar)*wl;
  Real ustar_r = ur + (pstar_lo - pstar)*wr;

  Real f_lo = ustar_l - ustar_r;

  // hi bounds
  wsqge(pl, taul, gamel, gdot,
        gamstar, gmin, gmax, clsql, pstar_hi, wlsq);

  wsqge(pr, taur, gamer, gdot,
        gamstar, gmin, gmax, clsqr, pstar_hi, wrsq);

  wl = 1.0_rt / std::sqrt(wlsq);
  wr = 1.0_rt / std::sqrt(wrsq);

  ustar_l = ul - (pstar_hi - pstar)*wl;
  ustar_r = ur + (pstar_hi - pstar)*wr;

  //Real f_hi = ustar_l - ustar_r;

  // bisection
  converged = false;
  Real pstar_c = 0.0;

  for (int iter = 0; iter < PSTAR_BISECT_FACTOR*lcg_maxiter; iter++) {

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

    Real f_c = ustar_l - ustar_r;

    if ( 0.5_rt * std::abs(pstar_lo - pstar_hi) < lcg_tol * pstar_c ) {
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


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
cons_state(const Real* qstate, Real* U) {

  U[URHO] = qstate[QRHO];

  // since we advect all 3 velocity components regardless of dimension, this
  // will be general
  U[UMX]  = qstate[QRHO]*qstate[QU];
  U[UMY]  = qstate[QRHO]*qstate[QV];
  U[UMZ]  = qstate[QRHO]*qstate[QW];

  U[UEDEN] = qstate[QREINT] + 0.5_rt*qstate[QRHO]*(qstate[QU]*qstate[QU] + qstate[QV]*qstate[QV] + qstate[QW]*qstate[QW]);
  U[UEINT] = qstate[QREINT];

  // we don't care about T here, but initialize it to make NaN
  // checking happy
  U[UTEMP] = 0.0;

#ifdef SHOCK_VAR
  U[USHK] = 0.0;
#endif

  for (int ipassive = 0; ipassive < npassive; ipassive++) {
    int n  = upassmap(ipassive);
    int nqs = qpassmap(ipassive);
    U[n] = qstate[QRHO]*qstate[nqs];
  }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
HLLC_state(const int idir, const Real S_k, const Real S_c,
           const Real* qstate, Real* U) {

  Real u_k = 0.0;
  if (idir == 0) {
    u_k = qstate[QU];
  } else if (idir == 1) {
    u_k = qstate[QV];
  } else if (idir == 2) {
    u_k = qstate[QW];
  }

  Real hllc_factor = qstate[QRHO]*(S_k - u_k)/(S_k - S_c);
  U[URHO] = hllc_factor;

  if (idir == 0) {
    U[UMX]  = hllc_factor*S_c;
    U[UMY]  = hllc_factor*qstate[QV];
    U[UMZ]  = hllc_factor*qstate[QW];

  } else if (idir == 1) {
    U[UMX]  = hllc_factor*qstate[QU];
    U[UMY]  = hllc_factor*S_c;
    U[UMZ]  = hllc_factor*qstate[QW];

  } else {
    U[UMX]  = hllc_factor*qstate[QU];
    U[UMY]  = hllc_factor*qstate[QV];
    U[UMZ]  = hllc_factor*S_c;
  }

  U[UEDEN] = hllc_factor*(qstate[QREINT]/qstate[QRHO] +
                          0.5_rt*(qstate[QU]*qstate[QU] + qstate[QV]*qstate[QV] + qstate[QW]*qstate[QW]) +
                          (S_c - u_k)*(S_c + qstate[QPRES]/(qstate[QRHO]*(S_k - u_k))));
  U[UEINT] = hllc_factor*qstate[QREINT]/qstate[QRHO];

  U[UTEMP] = 0.0; // we don't evolve T

#ifdef SHOCK_VAR
  U[USHK] = 0.0;
#endif

  for (int ipassive = 0; ipassive < npassive; ipassive++) {
    int n  = upassmap(ipassive);
    int nqs = qpassmap(ipassive);
    U[n] = hllc_factor*qstate[nqs];
  }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
compute_flux(const int idir, const Real bnd_fac, const int coord,
             const Real* U, const Real p,
             Real* F) {

  // given a conserved state, compute the flux in direction idir

  Real u_flx = U[UMX+idir]/U[URHO];

  if (bnd_fac == 0) {
    u_flx = 0.0;
  }

  F[URHO] = U[URHO]*u_flx;

  F[UMX] = U[UMX]*u_flx;
  F[UMY] = U[UMY]*u_flx;
  F[UMZ] = U[UMZ]*u_flx;

  auto mom_check = mom_flux_has_p(idir, idir, coord);

  if (mom_check) {
    // we do not include the pressure term in any non-Cartesian
    // coordinate directions
    F[UMX+idir] = F[UMX+idir] + p;
  }

  F[UEINT] = U[UEINT]*u_flx;
  F[UEDEN] = (U[UEDEN] + p)*u_flx;

  F[UTEMP] = 0.0;

#ifdef SHOCK_VAR
  F[USHK] = 0.0;
#endif

  for (int ipassive=0; ipassive < npassive; ipassive++) {
    int n = upassmap(ipassive);
    F[n] = U[n]*u_flx;
  }
}

#endif
