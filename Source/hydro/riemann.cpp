#include <Castro.H>
#include <Castro_F.H>

#include <riemann_solvers.H>

#ifdef RADIATION
#include <Radiation.H>
#endif

#include <cmath>

#include <eos.H>
using namespace amrex;

void
Castro::cmpflx_plus_godunov(const Box& bx,
                            Array4<Real> const& qm,
                            Array4<Real> const& qp,
                            Array4<Real> const& flx,
                            Array4<Real> const& qint,
#ifdef RADIATION
                            Array4<Real> const& rflx,
                            Array4<Real> const& lambda_int,
#endif
                            Array4<Real> const& qgdnv,
                            Array4<Real const> const& qaux_arr,
                            Array4<Real const> const& shk,
                            const int idir) {

  // note: bx is not necessarily the limits of the valid (no ghost
  // cells) domain, but could be hi+1 in some dimensions.  We rely on
  // the caller to specify the interfaces over which to solve the
  // Riemann problems

  // Solve Riemann problem to get the fluxes

  if (riemann_solver == 0 || riemann_solver == 1) {
    // approximate state Riemann solvers

    riemann_state(bx,
                  qm, qp,
                  qint,
#ifdef RADIATION
                  lambda_int,
#endif
                  qaux_arr,
                  idir, 0);

    compute_flux_q(bx,
                   qint, flx,
#ifdef RADIATION
                   lambda_int, rflx,
#endif
                   idir);

  } else if (riemann_solver == 2) {
    // HLLC
    HLLC(bx,
         qm, qp,
         qaux_arr,
         flx, qint,
         idir);

#ifndef AMREX_USE_GPU
  } else {
    amrex::Error("ERROR: invalid value of riemann_solver");
#endif
  }

  if (hybrid_riemann == 1) {
    // correct the fluxes using an HLL scheme if we are in a shock
    // and doing the hybrid approach

    auto coord = geom.Coord();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

      int is_shock = 0;

      if (idir == 0) {
        is_shock = static_cast<int>(shk(i-1,j,k) + shk(i,j,k));
      } else if (idir == 1) { 
        is_shock = static_cast<int>(shk(i,j-1,k) + shk(i,j,k));
      } else {
        is_shock = static_cast<int>(shk(i,j,k-1) + shk(i,j,k));
      }

      if (is_shock >= 1) {

        Real cl;
        Real cr;
        if (idir == 0) {
          cl = qaux_arr(i-1,j,k,QC);
          cr = qaux_arr(i,j,k,QC);
        } else if (idir == 1) { 
          cl = qaux_arr(i,j-1,k,QC);
          cr = qaux_arr(i,j,k,QC);
        } else {
          cl = qaux_arr(i,j,k-1,QC);
          cr = qaux_arr(i,j,k,QC);
        }

        Real ql_zone[NQ];
        Real qr_zone[NQ];
        Real flx_zone[NUM_STATE];

        for (int n = 0; n < NQ; n++) {
          ql_zone[n] = qm(i,j,k,n);
          qr_zone[n] = qp(i,j,k,n);
        }

        // pass in the current flux -- the
        // HLL solver will overwrite this
        // if necessary
        for (int n = 0; n < NUM_STATE; n++) {
          flx_zone[n] = flx(i,j,k,n);
        }

        HLL(ql_zone, qr_zone, cl, cr,
            idir, coord,
            flx_zone);

        for (int n = 0; n < NUM_STATE; n++) {
          flx(i,j,k,n) = flx_zone[n];
        }
      }
    });
  }


  store_godunov_state(bx,
                      qint,
#ifdef RADIATION
                      lambda_int,
#endif
                      qgdnv);
}


void
Castro::riemann_state(const Box& bx,
                      Array4<Real> const& qm,
                      Array4<Real> const& qp,
                      Array4<Real> const& qint,
#ifdef RADIATION
                      Array4<Real> const& lambda_int,
#endif
                      Array4<Real const> const& qaux_arr,
                      const int idir, const int compute_gammas) {

  // just compute the hydrodynamic state on the interfaces
  // don't compute the fluxes

  // note: bx is not necessarily the limits of the valid (no ghost
  // cells) domain, but could be hi+1 in some dimensions.  We rely on
  // the caller to specify the interfaces over which to solve the
  // Riemann problems

#ifdef RADIATION
#ifndef AMREX_USE_GPU
  if (hybrid_riemann == 1) {
    amrex::Error("ERROR: hybrid Riemann not supported for radiation");
  }

  if (riemann_solver > 0) {
    amrex::Error("ERROR: only the CGF Riemann solver is supported for radiation");
  }
#endif
#endif

#if AMREX_SPACEDIM == 1
#ifndef AMREX_USE_GPU
  if (riemann_solver > 1) {
    amrex::Error("ERROR: HLLC not implemented for 1-d");
  }
#endif
#endif

  if (ppm_temp_fix == 2) {
    // recompute the thermodynamics on the interface to make it
    // all consistent

    // we want to take the edge states of rho, e, and X, and get
    // new values for p on the edges that are
    // thermodynamically consistent.

    const Real lT_guess = T_guess;

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

     eos_t eos_state;

     // this is an initial guess for iterations, since we
     // can't be certain what temp is on interfaces
     eos_state.T = lT_guess;

     // minus state
     eos_state.rho = qm(i,j,k,QRHO);
     eos_state.p = qm(i,j,k,QPRES);
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
     eos_state.p = qp(i,j,k,QPRES);
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
    });
  }

  // Solve Riemann problem
  if (riemann_solver == 0) {
    // Colella, Glaz, & Ferguson solver

    riemannus(bx,
              qm, qp,
              qaux_arr, qint,
#ifdef RADIATION
              lambda_int,
#endif
              idir, compute_gammas);

  } else if (riemann_solver == 1) {
    // Colella & Glaz solver

#ifndef RADIATION
    riemanncg(bx,
              qm, qp,
              qaux_arr, qint,
              idir);
#endif

#ifndef AMREX_USE_GPU
  } else {
    amrex::Error("ERROR: invalid value of riemann_solver");
#endif
  }
}
