#ifndef _riemann_H_
#define _riemann_H_

#include <castro_params.H>
#include <Castro_util.H>

using namespace amrex;

const Real smlp1 = 1.e-10_rt;
#define SMALL 1.e-8_rt
const Real smallu = 1.e-12_rt;




AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
Castro::get_riemann_input_states(Array4<Real const> const ql_arr, Array4<Real const> const qr_arr,
                                 Array4<Real const> const qaux_arr,
                                 const int i, const int j, const int k, const int idir,
                                 const int compute_gammas,
                                 Real& gcl, Real& gcr,
                                 Real& cl, Real& cr,
#ifdef RADIATION
                                 Real* laml_int, Real* lamr_int,
                                 Real& gamcgl, Real& gamcgr,
#endif
                                 Real* ql_int, Real* qr_int) {


  // copy the initial interface states

  for (int n = 0; n < NQ; n++) {
    ql_int[n] = ql_arr(i,j,k,n);
  }

  for (int n = 0; n < NQ; n++) {
    qr_int[n] = qr_arr(i,j,k,n);
  }

  if (idir == 0) {
    gcl = qaux_arr(i-1,j,k,QGAMC);
    cl = qaux_arr(i-1,j,k,QC);

  } else if (idir == 1) {
    gcl = qaux_arr(i,j-1,k,QGAMC);
    cl = qaux_arr(i,j-1,k,QC);

  } else {
    gcl = qaux_arr(i,j,k-1,QGAMC);
    cl = qaux_arr(i,j,k-1,QC);

  }

  gcr = qaux_arr(i,j,k,QGAMC);
  cr = qaux_arr(i,j,k,QC);

#ifdef TRUE_SDC
  if (castro::use_reconstructed_gamma1 == 1) {
    gcl = ql_arr(i,j,k,QGC);
    gcr = qr_arr(i,j,k,QGC);
  }
#endif

#ifndef RADIATION
  if (compute_gammas == 1) {

    // we come in with a good p, rho, and X on the interfaces
    // -- use this to find the gamma used in the sound speed
    eos_t eos_state;
    eos_state.p = ql_int[QPRES];
    eos_state.rho = ql_int[QRHO];
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = ql_int[QFS+n];
      }
      eos_state.T = castro::T_guess; // initial guess
      for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = ql_int[QFX+n];
      }

      eos(eos_input_rp, eos_state);

      gcl = eos_state.gam1;

      eos_state.p = qr_int[QPRES];
      eos_state.rho = qr_int[QRHO];
      for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = qr_int[QFS+n];
      }
      eos_state.T = castro::T_guess; // initial guess
      for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = qr_int[QFX+n];
      }

      eos(eos_input_rp, eos_state);

      gcr = eos_state.gam1;

    }
#endif

#ifdef RADIATION
  if (idir == 0) {
    for (int g = 0; g < NGROUPS; g++) {
      laml_int[g] = qaux_arr(i-1,j,k,QLAMS+g);
    }
    lamr_int[g] = qaux_arr(i,j,k,QLAMS+g);

    gamcgl = qaux_arr(i-1,j,k,QGAMCG);
    gamcgr = qaux_arr(i,j,k,QGAMCG);

  } else if (idir == 1) {
    for (int g = 0; g < NGROUPS; g++) {
      laml_int[g] = qaux_arr(i,j-1,k,QLAMS+g);
    }
    lamr_int[g] = qaux_arr(i,j,k,QLAMS+g);

    gamcgl = qaux_arr(i,j-1,k,QGAMCG);
    gamcgr = qaux_arr(i,j,k,QGAMCG);

  } else {
    for (int g = 0; g < NGROUPS; g++) {
      laml_int[g] = qaux_arr(i,j,k-1,QLAMS+g);
    }
    lamr_int[g] = qaux_arr(i,j,k,QLAMS+g);

    gamcgl = qaux_arr(i,j,k-1,QGAMCG);
    gamcgr = qaux_arr(i,j,k,QGAMCG);

  }
#endif

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
Castro::riemann_state_interface(Real* qm, Real* qp,
                                Real& gcl, Real& gcr, const Real cl, const Real cr,
#ifdef RADIATION
                                Real* lambda_int,
#endif
                                Real* qint,
                                const Real bnd_fac, const int idir) {

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

    eos_t eos_state;

    // this is an initial guess for iterations, since we
    // can't be certain what temp is on interfaces
    eos_state.T = castro::T_guess;

    // minus state
    eos_state.rho = qm[QRHO];
    eos_state.p = qm[QPRES];
    eos_state.e = qm[QREINT]/qm[QRHO];
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = qm[QFS+n];
    }
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = qm[QFX+n];
    }

    eos(eos_input_re, eos_state);

    qm[QREINT] = eos_state.e * eos_state.rho;
    qm[QPRES] = eos_state.p;

    // plus state
    eos_state.rho = qp[QRHO];
    eos_state.p = qp[QPRES];
    eos_state.e = qp[QREINT]/qp[QRHO];
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = qp[QFS+n];
    }
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = qp[QFX+n];
    }

    eos(eos_input_re, eos_state);

    qp[QREINT] = eos_state.e * eos_state.rho;
    qp[QPRES] = eos_state.p;
  }

  // Solve Riemann problem
  if (riemann_solver == 0) {
    // Colella, Glaz, & Ferguson solver

    riemannus(qm, qp,
              gcl, gcr, cl, cr,
#ifdef RADIATION
              lambda_int,
#endif
              qint,
              bnd_fac, idir);

  } else if (riemann_solver == 1) {
    // Colella & Glaz solver

#ifndef RADIATION
    riemanncg(qm, qp,
              gcl, gcr, cl, cr,
              qint,
              bnd_fac, idir);
#endif

#ifndef AMREX_USE_CUDA
  } else {
    amrex::Error("ERROR: invalid value of riemann_solver");
#endif
  }
}

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
Castro::cons_state(const Real* qstate, Real* U) {

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
Castro::HLLC_state(const int idir, const Real S_k, const Real S_c,
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
Castro::compute_flux(const int idir, const Real bnd_fac, const int coord,
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
