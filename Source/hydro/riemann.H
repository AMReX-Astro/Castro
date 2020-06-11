#ifndef CASTRO_RIEMANN_H
#define CASTRO_RIEMANN_H

using namespace amrex;

const Real smlp1 = 1.e-10_rt;
const Real small = 1.e-8_rt;
const Real smallu = 1.e-12_rt;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
Castro::wsqge(const Real p, const Real v,
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
    alpha = smlp1*(pstar + p);
  }

  Real beta = pstar + 0.5_rt*(gstar - 1.0_rt)*(pstar+p);

  wsq = (pstar-p)*beta/(v*alpha);

  if (std::abs(pstar - p) < smlp1*(pstar + p)) {
    wsq = csq;
  }
  wsq = amrex::max(wsq, (0.5_rt * (gam - 1.0_rt)/gam)*csq);
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
Castro::pstar_bisection(Real& pstar_lo, Real& pstar_hi,
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
Castro::cons_state(const Real* qstate, Real* U,
                   const GpuArray<int, npassive>& qpass_map_p,
                   const GpuArray<int, npassive>& upass_map_p) {

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
    int n  = upass_map_p[ipassive];
    int nqs = qpass_map_p[ipassive];
    U[n] = qstate[QRHO]*qstate[nqs];
  }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
Castro::HLLC_state(const int idir, const Real S_k, const Real S_c,
                   const Real* qstate, Real* U,
                   const GpuArray<int, npassive>& qpass_map_p,
                   const GpuArray<int, npassive>& upass_map_p) {

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
    int n  = upass_map_p[ipassive];
    int nqs = qpass_map_p[ipassive];
    U[n] = hllc_factor*qstate[nqs];
  }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
Castro::compute_flux(const int idir, const Real bnd_fac, const int coord,
                     const Real* U, const Real p,
                     const GpuArray<int, npassive>& upass_map_p,
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
    int n = upass_map_p[ipassive];
    F[n] = U[n]*u_flx;
  }
}

#endif
