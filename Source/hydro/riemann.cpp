#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

#include <cmath>

#include <eos.H>
#include <riemann.H>

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
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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


    GpuArray<Real, NQ> qm_int;
    for (int n = 0; n < NQ; n++) {
      qm_int[n] = qm(i,j,k,n);
    }

    GpuArray<Real, NQ> qp_int;
    for (int n = 0; n < NQ; n++) {
      qp_int[n] = qp(i,j,k,n);
    }

    GpuArray<Real, NQ> q_riemann;

    Real gcl;
    Real cl;

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

    Real gcr = qaux_arr(i,j,k,QGAMC);
    Real cr = qaux_arr(i,j,k,QC);

#ifdef TRUE_SDC
    if (luse_reconstructed_gamma1 == 1) {
      gcl = qm(i,j,k,QGC);
      gcr = qp(i,j,k,QGC);
    }
#endif

#ifndef RADIATION
    if (compute_gammas == 1) {

      // we come in with a good p, rho, and X on the interfaces
      // -- use this to find the gamma used in the sound speed
      eos_t eos_state;
      eos_state.p = qm_int[QPRES];
      eos_state.rho = qm_int[QRHO];
      for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = qm_int[QFS+n];
      }
      eos_state.T = castro::T_guess; // initial guess
      for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = qm_int[QFX+n];
      }

      eos(eos_input_rp, eos_state);

      gcl = eos_state.gam1;

      eos_state.p = qp_int[QPRES];
      eos_state.rho = qp_int[QRHO];
      for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = qp_int[QFS+n];
      }
      eos_state.T = castro::T_guess; // initial guess
      for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = qp_int[QFX+n];
      }

      eos(eos_input_rp, eos_state);

      gcr = eos_state.gam1;

    }
#endif

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
                            gcl, gcr, cl, cr,
                            q_riemann, bnd_fac, idir);

    for (int n = 0; n < NQ; n++) {
      qint(i,j,k,n) = q_riemann[n];
    }

  });

}
