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
Castro::compute_flux_q(Real* q_riemann,
                       Real* F,
#ifdef RADIATION
                       Array4<Real const> const lambda,
                       Array4<Real> const rF,
#endif
#ifdef HYBRID_MOMEMTUM
                       GeometryData const& geomdata, GpuArray<Real, 3>& center,
#endif
                       const int coord_type,
                       const int idir, const int enforce_eos) {

  // given a primitive state, compute the flux in direction idir
  //

  int iu, iv1, iv2;
  int im1, im2, im3;

  auto mom_check = mom_flux_has_p(idir, idir, coord_type);

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

#ifdef RADIATION
  int fspace_t = Radiation::fspace_advection_type;
  int comov = Radiation::comoving;
  int limiter = Radiation::limiter;
  int closure = Radiation::closure;
#endif

  Real& u_adv = q_riemann[iu];
  Real& rhoeint = q_riemann[QREINT];

  // if we are enforcing the EOS, then take rho, p, and X, and
  // compute rhoe
  if (enforce_eos == 1) {
    eos_t eos_state;
    eos_state.rho = q_riemann[QRHO];
    eos_state.p = q_riemann[QPRES];
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = q_riemann[QFS+n];
    }
    eos_state.T = castro::T_guess;  // initial guess
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = q_riemann[QFX+n];
    }

    eos(eos_input_rp, eos_state);

    rhoeint = q_riemann[QRHO] * eos_state.e;
  }

  // Compute fluxes, order as conserved state (not q)
  F[URHO] = q_riemann[QRHO]*u_adv;

  F[im1] = F[URHO] * u_adv;
  if (mom_check) {
    F[im1] += q_riemann[QPRES];
  }
  F[im2] = F[URHO] * q_riemann[iv1];
  F[im3] = F[URHO] * q_riemann[iv2];

  Real rhoetot = rhoeint + 0.5_rt * q_riemann[QRHO] *
      (q_riemann[iu] * q_riemann[iu] +
       q_riemann[iv1] * q_riemann[iv1] +
       q_riemann[iv2] * q_riemann[iv2]);

  F[UEDEN] = u_adv * (rhoetot + q_riemann[QPRES]);
  F[UEINT] = u_adv * rhoeint;

  F[UTEMP] = 0.0;
#ifdef SHOCK_VAR
  F[USHK] = 0.0;
#endif

#ifdef RADIATION
  if (fspace_t == 1) {
    for (int g = 0; g < NGROUPS; g++) {
      Real eddf = Edd_factor(lambda(i,j,k,g), limiter, closure);
      Real f1 = 0.5e0_rt*(1.0_rt-eddf);
      rF(i,j,k,g) = (1.0_rt + f1) * q_riemann[QRAD+g] * u_adv;
    }
  } else {
    // type 2
    for (int g = 0; g < NGROUPS; g++) {
      rF(i,j,k,g) = q_riemann[QRAD+g] * u_adv;
    }
  }
#endif

  // passively advected quantities
  for (int ipassive = 0; ipassive < npassive; ipassive++) {
    int n  = upassmap(ipassive);
    int nqp = qpassmap(ipassive);

    F[n] = F[URHO] * q_riemann[nqp];
  }

#ifdef HYBRID_MOMENTUM
  // the hybrid routine uses the Godunov indices, not the full NQ state
  GpuArray<Real, NGDNV> qgdnv_zone;
  qgdnv_zone[GDRHO] = q_riemann[QRHO];
  qgdnv_zone[GDU] = q_riemann[QU];
  qgdnv_zone[GDV] = q_riemann[QV];
  qgdnv_zone[GDW] = q_riemann[QW];
  qgdnv_zone[GDPRES] = q_riemann[QPRES];
#ifdef RADIATION
  for (int g = 0; g < NGROUPS; g++) {
    qgdnv_zone[GDLAMS+g] = lambda(i,j,k,g);
    qgdnv_zone[GDERADS+g] = qint(i,j,k,QRAD+g);
  }
#endif
  compute_hybrid_flux(qgdnv_zone, geomdata, center, idir, i, j, k, F);
#endif

}


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
    qr_int[n] = qr_arr(i,j,k,n);
  }

  // some limits
  ql_int[QRHO] = amrex::max(ql_int[QRHO], castro::small_dens);
  qr_int[QRHO] = amrex::max(qr_int[QRHO], castro::small_dens);

  ql_int[QPRES] = amrex::max(ql_int[QPRES], castro::small_pres);
  qr_int[QPRES] = amrex::max(qr_int[QPRES], castro::small_pres);

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


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
Castro::riemanncg(Real* ql, Real* qr,
                  Real& gcl, Real& gcr,
                  const Real cl, const Real cr,
                  Real* qint,
                  const Real bnd_fac, const int idir) {

  // this implements the approximate Riemann solver of Colella & Glaz
  // (1985)
  //

  constexpr Real weakwv = 1.e-3_rt;

#ifndef AMREX_USE_CUDA
  GpuArray<Real, HISTORY_SIZE> pstar_hist;
  GpuArray<Real, PSTAR_BISECT_FACTOR*HISTORY_SIZE> pstar_hist_extra;

  if (castro::cg_maxiter > HISTORY_SIZE) {
    amrex::Error("error in riemanncg: cg_maxiter > HISTORY_SIZE");
  }
#endif

#ifndef AMREX_USE_CUDA
  if (castro::cg_blend == 2 && castro::cg_maxiter < 5) {
    amrex::Error("Error: need cg_maxiter >= 5 to do a bisection search on secant iteration failure.");
  }
#endif

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

  // left state
  Real rl = ql[QRHO];

  Real pl = ql[QPRES];
  Real rel = ql[QREINT];

  // pick left velocities based on direction
  Real ul = ql[iu];
  Real v1l = ql[iv1];
  Real v2l = ql[iv2];


  // sometime we come in here with negative energy or pressure
  // note: reset both in either case, to remain thermo
  // consistent
  if (rel <= 0.0_rt || pl < castro::small_pres) {
#ifndef AMREX_USE_CUDA
    std::cout <<  "WARNING: (rho e)_l < 0 or pl < small_pres in Riemann: " << rel << " " << pl << " " << castro::small_pres << std::endl;
#endif

    eos_t eos_state;
    eos_state.T = castro::small_temp;
    eos_state.rho = rl;
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = ql[QFS+n];
    }
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = ql[QFX+n];
    }

    eos(eos_input_rt, eos_state);

    rel = rl*eos_state.e;
    pl = eos_state.p;
    gcl = eos_state.gam1;
  }

  // right state
  Real rr = qr[QRHO];

  Real pr = qr[QPRES];
  Real rer = qr[QREINT];

  // pick right velocities based on direction
  Real ur = qr[iu];
  Real v1r = qr[iv1];
  Real v2r = qr[iv2];

  if (rer <= 0.0_rt || pr < castro::small_pres) {
#ifndef AMREX_USE_CUDA
    std::cout << "WARNING: (rho e)_r < 0 or pr < small_pres in Riemann: " << rer << " " << pr << " " << castro::small_pres << std::endl;
#endif
    eos_t eos_state;

    eos_state.T = castro::small_temp;
    eos_state.rho = rr;
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = qr[QFS+n];
    }
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = qr[QFX+n];
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

  Real csmall = amrex::max(SMALL, amrex::max(SMALL * cr, SMALL * cl));

  Real cavg = 0.5_rt*(cl + cr);

  // Note: in the original Colella & Glaz paper, they predicted
  // gamma_e to the interfaces using a special (non-hyperbolic)
  // evolution equation.  In Castro, we instead bring (rho e)
  // to the edges, so we construct the necessary gamma_e here from
  // what we have on the interfaces.
  Real gamel = pl/rel + 1.0_rt;
  Real gamer = pr/rer + 1.0_rt;

  // these should consider a wider average of the cell-centered
  // gammas
  Real gmin = amrex::min(amrex::min(gamel, gamer), 1.0_rt);
  Real gmax = amrex::max(amrex::max(gamel, gamer), 2.0_rt);

  Real game_bar = 0.5_rt*(gamel + gamer);
  Real gamc_bar = 0.5_rt*(gcl + gcr);

  Real gdot = 2.0_rt*(1.0_rt - game_bar/gamc_bar)*(game_bar - 1.0_rt);

  Real wsmall = castro::small_dens * csmall;
  Real wl = amrex::max(wsmall, std::sqrt(std::abs(clsql)));
  Real wr = amrex::max(wsmall, std::sqrt(std::abs(clsqr)));

  // make an initial guess for pstar -- this is a two-shock
  // approximation
  //pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
  Real pstar = pl + ( (pr - pl) - wr*(ur - ul) )*wl/(wl+wr);
  pstar = amrex::max(pstar, castro::small_pres);

  // get the shock speeds -- this computes W_s from CG Eq. 34
  Real gamstar = 0.0;
  Real wlsq = 0.0;

  wsqge(pl, taul, gamel, gdot, gamstar,
        gmin, gmax, clsql, pstar, wlsq);

  Real wrsq = 0.0;
  wsqge(pr, taur, gamer, gdot, gamstar,
        gmin, gmax, clsqr, pstar, wrsq);

  Real pstar_old = pstar;

  wl = std::sqrt(wlsq);
  wr = std::sqrt(wrsq);

  // R-H jump conditions give ustar across each wave -- these
  // should be equal when we are done iterating.  Our notation
  // here is a little funny, comparing to CG, ustar_l = u*_L and
  // ustar_r = u*_R.
  Real ustar_l = ul - (pstar-pl)/wl;
  Real ustar_r = ur + (pstar-pr)/wr;

  // revise our pstar guess
  // pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
  pstar = pl + ( (pr - pl) - wr*(ur - ul) )*wl/(wl+wr);
  pstar = amrex::max(pstar, castro::small_pres);

  // secant iteration
  bool converged = false;

  int iter = 0;
  while ((iter < castro::cg_maxiter && !converged) || iter < 2) {

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
    Real denom = dpditer/amrex::max(zp+zm, SMALL*cavg);
    pstar_old = pstar;
    pstar = pstar - denom*(ustar_r - ustar_l);
    pstar = amrex::max(pstar, castro::small_pres);

    Real err = std::abs(pstar - pstar_old);
    if (err < castro::cg_tol*pstar) {
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

    if (castro::cg_blend == 0) {

#ifndef AMREX_USE_CUDA
      std::cout <<  "pstar history: " << std::endl;
      for (int iter_l=0; iter_l < castro::cg_maxiter; iter_l++) {
        std::cout << iter_l << " " << pstar_hist[iter_l] << std::endl;
      }

      std::cout << std::endl;
      std::cout << "left state  (r,u,p,re,gc): " << rl << " " << ul << " " << pl << " " << rel << " " << gcl << std::endl;
      std::cout << "right state (r,u,p,re,gc): " << rr << " " << ur << " " << pr << " " << rer << " " << gcr << std::endl;
      std::cout << "cavg, smallc: " << cavg << " " << csmall;

      amrex::Error("ERROR: non-convergence in the Riemann solver");
#endif

    } else if (castro::cg_blend == 1) {

      pstar = pl + ( (pr - pl) - wr*(ur - ul) )*wl/(wl+wr);

    } else if (castro::cg_blend == 2) {

      // we don't store the history if we are in CUDA, so
      // we can't do this
#ifndef AMREX_USE_CUDA
      // first try to find a reasonable bounds
      Real pstarl = 1.e200;
      Real pstaru = -1.e200;
      for (int n = castro::cg_maxiter-6; n < castro::cg_maxiter; n++) {
        pstarl = amrex::min(pstarl, pstar_hist[n]);
        pstaru = amrex::max(pstaru, pstar_hist[n]);
      }

      pstarl = amrex::max(pstarl, castro::small_pres);
      pstaru = amrex::max(pstaru, castro::small_pres);

      pstar_bisection(pstarl, pstaru,
                      ul, pl, taul, gamel, clsql,
                      ur, pr, taur, gamer, clsqr,
                      gdot, gmin, gmax,
                      castro::cg_maxiter, castro::cg_tol, 
                      pstar, gamstar, converged, pstar_hist_extra);

      if (!converged) {

        std::cout << "pstar history: " << std::endl;
        for (int iter_l = 0; iter_l < castro::cg_maxiter; iter_l++) {
          std::cout << iter_l << " " << pstar_hist[iter_l] << std::endl;
        }
        std::cout << "pstar extra history: " << std::endl;
        for (int iter_l = 0; iter_l < PSTAR_BISECT_FACTOR*castro::cg_maxiter; iter_l++) {
          std::cout << iter_l << " " << pstar_hist_extra[iter_l] << std::endl;
        }

        std::cout << std::endl;
        std::cout << "left state  (r,u,p,re,gc): " << rl << " " << ul << " " << pl << " " << rel << " " << gcl << std::endl;
        std::cout << "right state (r,u,p,re,gc): " << rr << " " << ur << " " << pr << " " << rer << " " << gcr << std::endl;
        std::cout << "cavg, smallc: " << cavg << " " << csmall << std::endl;

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

  Real ustar = 0.5_rt * (ustar_l + ustar_r);

  // for symmetry preservation, if ustar is really small, then we
  // set it to zero
  if (std::abs(ustar) < smallu*0.5_rt*(std::abs(ul) + std::abs(ur))) {
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
  ro = amrex::max(castro::small_dens, 1.0_rt/tauo);
  tauo = 1.0_rt/ro;

  Real co = std::sqrt(std::abs(gamco*po*tauo));
  co = amrex::max(csmall, co);
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
  rstar = amrex::max(castro::small_dens, rstar);

  Real cstar = std::sqrt(std::abs(gamco*pstar/rstar));
  cstar = amrex::max(cstar, csmall);

  Real spout = co - sgnm*uo;
  Real spin = cstar - sgnm*ustar;

  //ushock = 0.5_rt*(spin + spout)
  Real ushock = wo*tauo - sgnm*uo;

  if (pstar-po >= 0.0_rt) {
    spin = ushock;
    spout = ushock;
  }

  Real frac = 0.5_rt*(1.0_rt + (spin + spout)/amrex::max(amrex::max(spout-spin, spin+spout),
                                                         SMALL*cavg));

  // the transverse velocity states only depend on the
  // direction that the contact moves
  if (ustar > 0.0_rt) {
    qint[iv1] = v1l;
    qint[iv2] = v2l;
  } else if (ustar < 0.0_rt) {
    qint[iv1] = v1r;
    qint[iv2] = v2r;
  } else {
    qint[iv1] = 0.5_rt*(v1l+v1r);
    qint[iv2] = 0.5_rt*(v2l+v2r);
  }

  // linearly interpolate between the star and normal state -- this covers the
  // case where we are inside the rarefaction fan.
  qint[QRHO] = frac*rstar + (1.0_rt - frac)*ro;
  qint[iu] = frac*ustar + (1.0_rt - frac)*uo;
  qint[QPRES] = frac*pstar + (1.0_rt - frac)*po;
  Real game_int = frac*gamstar + (1.0_rt-frac)*gameo;

  // now handle the cases where instead we are fully in the
  // star or fully in the original (l/r) state
  if (spout < 0.0_rt) {
    qint[QRHO] = ro;
    qint[iu] = uo;
    qint[QPRES] = po;
    game_int = gameo;
  }

  if (spin >= 0.0_rt) {
    qint[QRHO] = rstar;
    qint[iu] = ustar;
    qint[QPRES] = pstar;
    game_int = gamstar;
  }

  qint[QPRES] = amrex::max(qint[QPRES], castro::small_pres);

  qint[iu] = qint[iu] * bnd_fac;

  // Compute fluxes, order as conserved state (not q)

  // compute the total energy from the internal, p/(gamma - 1), and the kinetic
  qint[QREINT] = qint[QPRES]/(game_int - 1.0_rt);

  // advected quantities -- only the contact matters
  for (int ipassive = 0; ipassive < npassive; ipassive++) {
    int nqp = qpassmap(ipassive);

    if (ustar > 0.0_rt) {
      qint[nqp] = ql[nqp];
    } else if (ustar < 0.0_rt) {
      qint[nqp] = qr[nqp];
    } else {
      qint[nqp] = 0.5_rt * (ql[nqp] + qr[nqp]);
    }
  }

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
Castro::riemannus(Real* ql, Real* qr,
                  const Real gamcl, const Real gamcr,
                  const Real cl, const Real cr,
#ifdef RADIATION
                  const Real gamcgl, const Real gamcgr,
                  Real* laml, Real* lamr,
#endif
                  Real* qint,
                  const Real bnd_fac, const int idir) {

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

  // set the left and right states for this interface

  Real& rl = ql[QRHO];

  // pick left velocities based on direction
  Real& ul  = ql[iu];
  Real& v1l = ql[iv1];
  Real& v2l = ql[iv2];

#ifdef RADIATION
  Real pl = ql[QPTOT];
  Real rel = ql[QREITOT];
  Real erl[NGROUPS];
  for (int g = 0; g < NGROUPS; g++) {
    erl[g] = ql[QRAD+g];
  }
  Real pl_g = ql[QPRES];
  Real rel_g = ql[QREINT];
#else
  Real& pl = ql[QPRES];
  Real& rel = ql[QREINT];
#endif

  Real& rr = qr[QRHO];

  // pick right velocities based on direction
  Real& ur  = qr[iu];
  Real& v1r = qr[iv1];
  Real& v2r = qr[iv2];

#ifdef RADIATION
  Real pr = qr[QPTOT];
  Real rer = qr[QREITOT];
  Real err[NGROUPS];
  for (int g = 0; g < NGROUPS; g++) {
    err[g] = qr[QRAD+g];
  }
  Real pr_g = qr[QPRES];
  Real rer_g = qr[QREINT];
#else
  Real& pr = qr[QPRES];
  Real& rer = qr[QREINT];
#endif

  // estimate the star state: pstar, ustar

  Real csmall;
  Real cavg;
#ifdef RADIATION
  Real gamcgl;
  Real gamcgr;
#endif

  csmall = amrex::max(SMALL, SMALL * amrex::max(cl, cr));
  cavg = 0.5_rt*(cl + cr);

  Real wsmall = castro::small_dens*csmall;

  // this is Castro I: Eq. 33
  Real wl = amrex::max(wsmall, std::sqrt(std::abs(gamcl*pl*rl)));
  Real wr = amrex::max(wsmall, std::sqrt(std::abs(gamcr*pr*rr)));

  Real wwinv = 1.0_rt/(wl + wr);
  Real pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))*wwinv;
  Real ustar = ((wl*ul + wr*ur) + (pl - pr))*wwinv;

  pstar = amrex::max(pstar, castro::small_pres);

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
  Real lambda[NGROUPS];
  for (int g = 0; g < NGROUPS; g++) {
    lambda[g] = fp*laml[g] + fm*lamr[g];
  }

  if (ustar == 0) {
    // harmonic average
    for (int g = 0; g < NGROUPS; g++) {
      lambda[g] = 2.0_rt*(laml[g]*lamr[g])/(laml[g] + lamr[g] + 1.e-50_rt);
    }
  }

  Real po_g = fp*pl_g + fm*pr_g;
  Real reo_r[NGROUPS];
  Real po_r[NGROUPS];
  for (int g = 0; g < NGROUPS; g++) {
    reo_r[g] = fp*erl[g] + fm*err[g];
    po_r[g] = lambda[g]*reo_r[g];
  }
  Real reo_g = fp*rel_g + fm*rer_g;
  Real gamco_g = fp*gamcgl + fm*gamcgr;
#endif

  ro = amrex::max(castro::small_dens, ro);

  Real roinv = 1.0_rt/ro;

  Real co = std::sqrt(std::abs(gamco*po*roinv));
  co = amrex::max(csmall, co);
  Real co2inv = 1.0_rt/(co*co);

  // we can already deal with the transverse velocities -- they
  // only jump across the contact
  qint[iv1] = fp*v1l + fm*v1r;
  qint[iv2] = fp*v2l + fm*v2r;

  // compute the rest of the star state

  Real drho = (pstar - po)*co2inv;
  Real rstar = ro + drho;
  rstar = amrex::max(castro::small_dens, rstar);

#ifdef RADIATION
  Real estar_g = reo_g + drho*(reo_g + po_g)*roinv;

  Real co_g = std::sqrt(std::abs(gamco_g*po_g*roinv));
  co_g = amrex::max(csmall, co_g);

  Real pstar_g = po_g + drho*co_g*co_g;
  pstar_g = amrex::max(pstar_g, castro::small_pres);

  Real estar_r[NGROUPS];
  for (int g = 0; g < NGROUPS; g++) {
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

  Real scr = spout - spin;
  if (spout-spin == 0.0_rt) {
    scr = SMALL*cavg;
  }

  // interpolate for the case that we are in a rarefaction
  Real frac = (1.0_rt + (spout + spin)/scr)*0.5_rt;
  frac = amrex::max(0.0_rt, amrex::min(1.0_rt, frac));

  qint[QRHO] = frac*rstar + (1.0_rt - frac)*ro;
  qint[iu] = frac*ustar + (1.0_rt - frac)*uo;

#ifdef RADIATION
  Real pgdnv_t = frac*pstar + (1.0_rt - frac)*po;
  Real pgdnv_g = frac*pstar_g + (1.0_rt - frac)*po_g;
  Real regdnv_g = frac*estar_g + (1.0_rt - frac)*reo_g;
  Real regdnv_r[NGROUPS];
  for (int g = 0; g < NGROUPS; g++) {
    regdnv_r[g] = frac*estar_r[g] + (1.0_rt - frac)*reo_r[g];
  }
#else
  qint[QPRES] = frac*pstar + (1.0_rt - frac)*po;
  Real regdnv = frac*estar + (1.0_rt - frac)*reo;
#endif

  // as it stands now, we set things assuming that the rarefaction
  // spans the interface.  We overwrite that here depending on the
  // wave speeds

  // look at the speeds on either side of the remaining wave
  // to determine which region we are in
  if (spout < 0.0_rt) {
    // the l or r state is on the interface
    qint[QRHO] = ro;
    qint[iu] = uo;
#ifdef RADIATION
    pgdnv_t = po;
    pgdnv_g = po_g;
    regdnv_g = reo_g;
    for (int g = 0; g < NGROUPS; g++) {
      regdnv_r[g] = reo_r[g];
    }
#else
    qint[QPRES] = po;
    regdnv = reo;
#endif
  }

  if (spin >= 0.0_rt) {
    // the star state is on the interface
    qint[QRHO] = rstar;
    qint[iu] = ustar;
#ifdef RADIATION
    pgdnv_t = pstar;
    pgdnv_g = pstar_g;
    regdnv_g = estar_g;
    for (int g = 0; g < NGROUPS; g++) {
      regdnv_r[g] = estar_r[g];
    }
#else
    qint[QPRES] = pstar;
    regdnv = estar;
#endif
  }

#ifdef RADIATION
  for (int g = 0; g < NGROUPS; g++) {
    qint[QRAD+g] = amrex::max(regdnv_r[g], 0.0_rt);
  }

  qint[QPRES] = pgdnv_g;
  qint[QPTOT] = pgdnv_t;
  qint[QREINT] = regdnv_g;

  qint[QREITOT] = regdnv_g;
  for (int g = 0; g < NGROUPS; g++) {
    qint[QREITOT] += regdnv_r[g];
  }

  for (int g = 0; g < NGROUPS; g++) {
    lambda_int(g) = lambda[g];
  }

#else
  qint[QPRES] = amrex::max(qint[QPRES], castro::small_pres);
  qint[QREINT] = regdnv;
#endif

  // we are potentially thermodynamically inconsistent, fix that
  // here
  if (castro::use_eos_in_riemann == 1) {
    // we need to know the species -- they only jump across
    // the contact
    eos_t eos_state;

    eos_state.rho = qint[QRHO];
    eos_state.p = qint[QPRES];

    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = fp*ql[QFS+n] + fm*qr[QFS+n];
    }

    eos_state.T = castro::T_guess;

    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = fp*ql[QFX+n] + fm*qr[QFX+n];
    }

    eos(eos_input_rp, eos_state);

    qint[QREINT] = eos_state.rho * eos_state.e;
  }

  // Enforce that fluxes through a symmetry plane or wall are hard zero.
  qint[iu] = qint[iu] * bnd_fac;

  // passively advected quantities
  for (int ipassive = 0; ipassive < npassive; ipassive++) {
    int nqp = qpassmap(ipassive);
    qint[nqp] = fp*ql[nqp] + fm*qr[nqp];
  }
}



#endif
