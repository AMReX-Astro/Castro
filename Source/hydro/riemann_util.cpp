#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#include "fluxlimiter.H"
#endif

#include <eos.H>

#include <cmath>

using namespace amrex;

#include <riemann.H>

AMREX_GPU_HOST_DEVICE
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


AMREX_GPU_HOST_DEVICE
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

  Real f_hi = ustar_l - ustar_r;

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
      f_hi = f_c;
    } else {
      pstar_lo = pstar_c;
      f_lo = f_c;
    }
  }

  pstar = pstar_c;
}


AMREX_GPU_HOST_DEVICE
void
Castro::cons_state(const Real* q, Real* U,
                   const GpuArray<int, npassive>& qpass_map_p,
                   const GpuArray<int, npassive>& upass_map_p) {


  U[URHO] = q[QRHO];

  // since we advect all 3 velocity components regardless of dimension, this
  // will be general
  U[UMX]  = q[QRHO]*q[QU];
  U[UMY]  = q[QRHO]*q[QV];
  U[UMZ]  = q[QRHO]*q[QW];

  U[UEDEN] = q[QREINT] + 0.5_rt*q[QRHO]*(q[QU]*q[QU] + q[QV]*q[QV] + q[QW]*q[QW]);
  U[UEINT] = q[QREINT];

  // we don't care about T here, but initialize it to make NaN
  // checking happy
  U[UTEMP] = 0.0;

#ifdef SHOCK_VAR
  U[USHK] = 0.0;
#endif

  for (int ipassive = 0; ipassive < npassive; ipassive++) {
    int n  = upass_map_p[ipassive];
    int nqs = qpass_map_p[ipassive];
    U[n] = q[QRHO]*q[nqs];
  }
}

AMREX_GPU_HOST_DEVICE
void
Castro::HLLC_state(const int idir, const Real S_k, const Real S_c,
                   const Real* q, Real* U) {

  Real u_k;
  if (idir == 0) {
    u_k = q[QU];
  } else if (idir == 1) {
    u_k = q[QV];
  } else if (idir == 2) {
    u_k = q[QW];
  }

  Real hllc_factor = q[QRHO]*(S_k - u_k)/(S_k - S_c);
  U[URHO] = hllc_factor;

  if (idir == 0) {
    U[UMX]  = hllc_factor*S_c;
    U[UMY]  = hllc_factor*q[QV];
    U[UMZ]  = hllc_factor*q[QW];

  } else if (idir == 1) {
    U[UMX]  = hllc_factor*q[QU];
    U[UMY]  = hllc_factor*S_c;
    U[UMZ]  = hllc_factor*q[QW];

  } else {
    U[UMX]  = hllc_factor*q[QU];
    U[UMY]  = hllc_factor*q[QV];
    U[UMZ]  = hllc_factor*S_c;
  }

  U[UEDEN] = hllc_factor*(q[QREINT]/q[QRHO] +
                          0.5_rt*(q[QU]*q[QU] + q[QV]*q[QV] + q[QW]*q[QW]) +
                          (S_c - u_k)*(S_c + q[QPRES]/(q[QRHO]*(S_k - u_k))));
  U[UEINT] = hllc_factor*q[QREINT]/q[QRHO];

  U[UTEMP] = 0.0; // we don't evolve T

#ifdef SHOCK_VAR
  U[USHK] = 0.0;
#endif

  for (int ipassive = 0; ipassive < npassive; ipassive++) {
    int n  = upass_map[ipassive];
    int nqs = qpass_map[ipassive];
    U[n] = hllc_factor*q[nqs];
  }
}

AMREX_GPU_HOST_DEVICE
void
Castro::compute_flux_q(const Box& bx,
                       Array4<Real const> const qint,
                       Array4<Real> const F,
#ifdef RADIATION
                       Array4<Real const> const lambda,
                       Array4<Real> const rF,
#endif
                       const int idir, const int enforce_eos) {

  // given a primitive state, compute the flux in direction idir
  //

  int iu, iv1, iv2;
  int im1, im2, im3;
  int mom_check = 0;

  if (idir == 0) {
    iu = QU;
    iv1 = QV;
    iv2 = QW;
    im1 = UMX;
    im2 = UMY;
    im3 = UMZ;
    mom_check = momx_flux_has_p[idir];

  } else if (idir == 1) {
    iu = QV;
    iv1 = QU;
    iv2 = QW;
    im1 = UMY;
    im2 = UMX;
    im3 = UMZ;
    mom_check = momy_flux_has_p[idir];

  } else {
    iu = QW;
    iv1 = QU;
    iv2 = QV;
    im1 = UMZ;
    im2 = UMX;
    im3 = UMY;
    mom_check = momz_flux_has_p[idir];
  }

#ifdef RADIATION
    int fspace_t = Radiation::fspace_advection_type;
    int comov = Radiation::comoving;
    int limiter = Radiation::limiter;
    int closure = Radiation::closure;
#endif

  AMREX_PARALLEL_FOR_3D(bx, i, j, k,
  {

    Real u_adv = qint(i,j,k,iu);
    Real rhoeint = qint(i,j,k,QREINT);

    // if we are enforcing the EOS, then take rho, p, and X, and
    // compute rhoe
    if (enforce_eos == 1) {
      eos_t eos_state;
      eos_state.rho = qint(i,j,k,QRHO);
      eos_state.p = qint(i,j,k,QPRES);
      for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = qint(i,j,k,QFS+n);
      }
      eos_state.T = T_guess;  // initial guess
      for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = qint(i,j,k,QFX+n);
      }

      eos(eos_input_rp, eos_state);

      rhoeint = qint(i,j,k,QRHO) * eos_state.e;
    }

    // Compute fluxes, order as conserved state (not q)
    F(i,j,k,URHO) = qint(i,j,k,QRHO)*u_adv;

    F(i,j,k,im1) = F(i,j,k,URHO)*qint(i,j,k,iu);
    if (mom_check) {
      F(i,j,k,im1) += qint(i,j,k,QPRES);
    }
    F(i,j,k,im2) = F(i,j,k,URHO)*qint(i,j,k,iv1);
    F(i,j,k,im3) = F(i,j,k,URHO)*qint(i,j,k,iv2);

    Real rhoetot = rhoeint + 0.5_rt * qint(i,j,k,QRHO)*
      (qint(i,j,k,iu)*qint(i,j,k,iu) +
       qint(i,j,k,iv1)*qint(i,j,k,iv1) +
       qint(i,j,k,iv2)*qint(i,j,k,iv2));

    F(i,j,k,UEDEN) = u_adv*(rhoetot + qint(i,j,k,QPRES));
    F(i,j,k,UEINT) = u_adv*rhoeint;

    F(i,j,k,UTEMP) = 0.0;
#ifdef SHOCK_VAR
    F(i,j,k,USHK) = 0.0;
#endif

#ifdef RADIATION
    if (fspace_t == 1) {
      for (int g = 0; g < NGROUPS; g++) {
        Real eddf = Edd_factor(lambda(i,j,k,g), limiter, closure);
        Real f1 = 0.5e0_rt*(1.0_rt-eddf);
        rF(i,j,k,g) = (1.0_rt + f1) * qint(i,j,k,QRAD+g) * u_adv;
      }
    } else {
      // type 2
      for (int g = 0; g < NGROUPS; g++) {
        rF(i,j,k,g) = qint(i,j,k,QRAD+g) * u_adv;
      }
    }
#endif

    // passively advected quantities
    for (int ipassive = 0; ipassive < npassive; ipassive++) {
      int n  = upass_map[ipassive];
      int nqp = qpass_map[ipassive];

      F(i,j,k,n) = F(i,j,k,URHO)*qint(i,j,k,nqp);
    }

#ifdef HYBRID_MOMENTUM
    // the hybrid routine uses the Godunov indices, not the full NQ state
    Real qgdnv_zone[NGDNV];
    qgdnv_zone[GDRHO] = qint(i,j,k,QRHO);
    qgdnv_zone[GDU] = qint(i,j,k,QU);
    qgdnv_zone[GDV] = qint(i,j,k,QV);
    qgdnv_zone[GDW] = qint(i,j,k,QW);
    qgdnv_zone[GDPRES] = qint(i,j,k,QPRES);
#ifdef RADIATION
    for (int g = 0; g < NGROUPS; g++) {
      qgdnv_zone[GDLAMS+g] = lambda(i,j,k,g);
      qgdnv_zone(GDERADS+g] = qint(i,j,k,QRAD+g);
    }
#endif
    Real F_zone[NUM_STATE];
    for (int n = 0; n < NUM_STATE; n++) {
      F_zone[n] = F(i,j,k,n);
    }
    compute_hybrid_flux(qgdnv_zone, F_zone, idir, i, j, k);
    for (int n = 0; n < NUM_STATE; n++) {
      F(i,j,k,n) = F_zone[n];
    }
#endif
  });
}


AMREX_GPU_HOST_DEVICE
void
Castro::store_godunov_state(const Box& bx,
                            Array4<Real const> const qint,
#ifdef RADIATION
                            Array4<Real const> const lambda,
#endif
                            Array4<Real> const qgdnv) {

  // this copies the full interface state (NQ -- one for each primitive
  // variable) over to a smaller subset of size NGDNV for use later in the
  // hydro advancement.

  AMREX_PARALLEL_FOR_3D(bx, i, j, k,
  {

    // the hybrid routine uses the Godunov indices, not the full NQ state
    qgdnv(i,j,k,GDRHO) = qint(i,j,k,QRHO);
    qgdnv(i,j,k,GDU) = qint(i,j,k,QU);
    qgdnv(i,j,k,GDV) = qint(i,j,k,QV);
    qgdnv(i,j,k,GDW) = qint(i,j,k,QW);
    qgdnv(i,j,k,GDPRES) = qint(i,j,k,QPRES);
#ifdef RADIATION
    for (int g = 0; g < NGROUPS; g++) {
      qgdnv(i,j,k,GDLAMS+g) = lambda(i,j,k,g);
      qgdnv(i,j,k,GDERADS+g) = qint(i,j,k,QRAD+g);
    }
#endif
  });
}

AMREX_GPU_HOST_DEVICE
void
Castro::compute_flux(const int idir, const Real bnd_fac,
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

  int mom_check = 0;
  if (idir == 0) {
    mom_check = momx_flux_has_p[idir];
  } else if (idir == 1) {
    mom_check = momy_flux_has_p[idir];
  } else {
    mom_check = momz_flux_has_p[idir];
  }

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
    int n = upass_map[ipassive];
    F[n] = U[n]*u_flx;
  }
}


