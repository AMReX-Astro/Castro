#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

#include <cmath>

#include <eos.H>
using namespace amrex;

void
Castro::cmpflx_plus_godunov(const Box& bx,
                            Array4<Real> const qm,
                            Array4<Real> const qp,
                            Array4<Real> const flx,
                            Array4<Real> const qint,
#ifdef RADIATION
                            Array4<Real> const rflx,
                            Array4<Real> const lambda_int,
#endif
                            Array4<Real> const qgdnv,
                            Array4<Real const> const qaux_arr,
                            Array4<Real const> const shk,
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
                   idir, 0);

  } else if (riemann_solver == 2) {
    // HLLC
    HLLC(bx,
         qm, qp,
         qaux_arr,
         flx, qint,
         idir);

#ifndef AMREX_USE_CUDA
  } else {
    amrex::Error("ERROR: invalid value of riemann_solver");
#endif
  }

  if (hybrid_riemann == 1) {
    // correct the fluxes using an HLL scheme if we are in a shock
    // and doing the hybrid approach

    GpuArray<int, npassive> upass_map_p;
    GpuArray<int, npassive> qpass_map_p;
    for (int n = 0; n < npassive; ++n) {
      upass_map_p[n] = upass_map[n];
      qpass_map_p[n] = qpass_map[n];
    }

    auto coord = geom.Coord();

    AMREX_PARALLEL_FOR_3D(bx, i, j, k,
    {

      int is_shock = 0;

      if (idir == 0) {
        is_shock = shk(i-1,j,k) + shk(i,j,k);
      } else if (idir == 1) { 
        is_shock = shk(i,j-1,k) + shk(i,j,k);
      } else {
        is_shock = shk(i,j,k-1) + shk(i,j,k);
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
            upass_map_p, qpass_map_p,
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
                      Array4<Real> const qm,
                      Array4<Real> const qp,
                      Array4<Real> const qint,
#ifdef RADIATION
                      Array4<Real> lambda_int,
#endif
                      Array4<Real const> const qaux_arr,
                      const int idir, const int compute_gammas) {


  AMREX_PARALLEL_FOR_3D(bx, i, j, k,
  {

    GpuArray<Real, NQ> qm_int;
    for (int n = 0; n < NQ; n++) {
      qm_int[n] = qm(i,j,k,n);
    }

    GpuArray<Real, NQ> qp_int;
    for (int n = 0; n < NQ; n++) {
      qp_int[n] = qp(i,j,k,n);
    }

    GpuArray<Real, NQ> q_riemann;

    Real gcl = qaux_arr(i-sx,j-sy,k-sz,QGAMC);
#ifdef TRUE_SDC
    if (luse_reconstructed_gamma1 == 1) {
      gcl = ql(i,j,k,QGC);
    }
#endif

    Real gcr = qaux_arr(i,j,k,QGAMC);
#ifdef TRUE_SDC
    if (luse_reconstructed_gamma1 == 1) {
      gcr = qr(i,j,k,QGC);
    }
#endif

    Real cl = qaux_arr(i-sx,j-sy,k-sz,QC);
    Real cr = qaux_arr(i,j,k,QC);

#ifdef RADIATION
    Real laml[NGROUPS];
    Real lamr[NGROUPS];

    if (idir == 0) {
      for (int g = 0; g < NGROUPS; g++) {
        laml[g] = qaux_arr(i-1,j,k,QLAMS+g);
      }
      lamr[g] = qaux_arr(i,j,k,QLAMS+g);

      gamcgl = qaux_arr(i-1,j,k,QGAMCG);
      gamcgr = qaux_arr(i,j,k,QGAMCG);

    } else if (idir == 1) {
      for (int g = 0; g < NGROUPS; g++) {
        laml[g] = qaux_arr(i,j-1,k,QLAMS+g);
      }
      lamr[g] = qaux_arr(i,j,k,QLAMS+g);

      gamcgl = qaux_arr(i,j-1,k,QGAMCG);
      gamcgr = qaux_arr(i,j,k,QGAMCG);

    } else {
      for (int g = 0; g < NGROUPS; g++) {
        laml[g] = qaux_arr(i,j,k-1,QLAMS+g);
      }
      lamr[g] = qaux_arr(i,j,k,QLAMS+g);

      gamcgl = qaux_arr(i,j,k-1,QGAMCG);
      gamcgr = qaux_arr(i,j,k,QGAMCG);

    }
#endif

    riemann_state_interface(qm_int, qp_int,
                            qcl, qcr, cl, cr,
                            q_riemann, idir);

    for (int n = 0; n < NQ; n++) {
      qint(i,j,k,n) = q_riemann[n];
    }

  });

}


void
Castro::riemann_state_interface(GpuArray<Real, NQ>& qm, GpuArray<Real, NQ>& qp,
                                const Real qcl, const Real qcr, const Real cl, const Real cr,
#ifdef RADIATION
                                GpuArray<Real, Radiation::ngroups> lambda_int,
#endif
                                GpuArray<Real, NQ>& qint,
                                const int idir) {

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

    const Real lT_guess = T_guess;

    eos_t eos_state;

    // this is an initial guess for iterations, since we
    // can't be certain what temp is on interfaces
    eos_state.T = lT_guess;

    // minus state
    eos_state.rho = qm(QRHO);
    eos_state.p = qm(QPRES);
    eos_state.e = qm(QREINT)/qm(QRHO);
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = qm(QFS+n);
    }
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = qm(QFX+n);
    }

    eos(eos_input_re, eos_state);

    qm(QREINT) = eos_state.e * eos_state.rho;
    qm(QPRES) = eos_state.p;

    // plus state
    eos_state.rho = qp(QRHO);
    eos_state.p = qp(QPRES);
    eos_state.e = qp(QREINT)/qp(QRHO);
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = qp(QFS+n);
    }
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = qp(QFX+n);
    }

    eos(eos_input_re, eos_state);

     qp(QREINT) = eos_state.e * eos_state.rho;
     qp(QPRES) = eos_state.p;
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
    riemanncg(qm, qp,
              qcl, qcr, cl, cr,
              qint,
              idir);
#endif

#ifndef AMREX_USE_CUDA
  } else {
    amrex::Error("ERROR: invalid value of riemann_solver");
#endif
  }
}
