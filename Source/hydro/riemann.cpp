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
                  idir);

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
                      const int idir) {

  // just compute the hydrodynamic state on the interfaces
  // don't compute the fluxes

  // note: bx is not necessarily the limits of the valid (no ghost
  // cells) domain, but could be hi+1 in some dimensions.  We rely on
  // the caller to specify the interfaces over which to solve the
  // Riemann problems

#ifndef AMREX_USE_GPU

#ifdef RADIATION
  if (hybrid_riemann == 1) {
    amrex::Error("ERROR: hybrid Riemann not supported for radiation");
  }

  if (riemann_solver > 0) {
    amrex::Error("ERROR: only the CGF Riemann solver is supported for radiation");
  }
#endif

#if AMREX_SPACEDIM == 1
  if (riemann_solver > 1) {
    amrex::Error("ERROR: HLLC not implemented for 1-d");
  }
#endif

  if (riemann_solver == 1) {
      if (cg_maxiter > HISTORY_SIZE) {
          amrex::Error("error in riemanncg: cg_maxiter > HISTORY_SIZE");
      }

      if (cg_blend == 2 && cg_maxiter < 5) {
          amrex::Error("Error: need cg_maxiter >= 5 to do a bisection search on secant iteration failure.");
      }
  }
#endif


  const auto domlo = geom.Domain().loVect3d();
  const auto domhi = geom.Domain().hiVect3d();

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


      if (ppm_temp_fix == 2) {
          // recompute the thermodynamics on the interface to make it
          // all consistent

          // we want to take the edge states of rho, e, and X, and get
          // new values for p on the edges that are
          // thermodynamically consistent.

          eos_t eos_state;

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

          riemannus(i, j, k,
                    ql, qr, raux,
                    qint,
#ifdef RADIATION
                    lambda_int,
#endif
                    idir);

      } else if (riemann_solver == 1) {
          // Colella & Glaz solver

#ifndef RADIATION
          riemanncg(i, j, k,
                    ql, qr, raux,
                    qint,
                    idir);
#endif

#ifndef AMREX_USE_GPU
      } else {
          amrex::Error("ERROR: invalid value of riemann_solver");
#endif
      }

      // the passives are always just upwinded, so we do that here
      // regardless of the solver

      Real sgnm = std::copysign(1.0_rt, qint(i,j,k,QU+idir));
      if (qint(i,j,k,QU+idir) == 0.0_rt) {
          sgnm = 0.0_rt;
      }

      Real fp = 0.5_rt*(1.0_rt + sgnm);
      Real fm = 0.5_rt*(1.0_rt - sgnm);

      for (int ipassive = 0; ipassive < npassive; ipassive++) {
          int nqp = qpassmap(ipassive);
          qint(i,j,k,nqp) = fp * qm(i,j,k,nqp) + fm * qp(i,j,k,nqp);
      }

  });

}
