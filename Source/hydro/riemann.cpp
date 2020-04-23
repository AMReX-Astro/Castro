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
                            Array4<Real> const ql_arr,
                            Array4<Real> const qr_arr,
                            Array4<Real> const flx,
#ifdef RADIATION
                            Array4<Real> const rflx,
#endif
                            Array4<Real> const qgdnv,
                            Array4<Real const> const qaux_arr,
                            Array4<Real const> const shk,
                            const int idir) {

  // note: bx is not necessarily the limits of the valid (no ghost
  // cells) domain, but could be hi+1 in some dimensions.  We rely on
  // the caller to specify the interfaces over which to solve the
  // Riemann problems

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

  auto coord = geom.Coord();

#ifdef HYBRID_MOMENTUM
  GeometryData geomdata = geom.data();

  GpuArray<Real, 3> center;
  ca_get_center(center.begin());
#endif

  if (castro::riemann_solver == 0 || castro::riemann_solver == 1) {
    // approximate state Riemann solvers

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {

      // Solve Riemann problem to get the fluxes

      GpuArray<Real, NQ> ql_int;
      GpuArray<Real, NQ> qr_int;

      Real gcl, gcr;
      Real cl, cr;

#ifdef RADIATION
      GpuArray<Real, Radiation::ngroups> laml_int;
      GpuArray<Real, Radiation::ngroups> lamr_int;
      Real gamcgl, gamcgr;
#endif

      get_riemann_input_states(ql_arr, qr_arr,
                               qaux_arr,
                               i, j, k, idir, 0,
                               gcl, gcr, cl, cr,
#ifdef RADIATION
                               laml_int, lamr_int,
                               gamcgl, gamcgr,
#endif
                               ql_int, qr_int);


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


      // solve the Riemann problem

      GpuArray<Real, NQ> q_riemann;

      riemann_state_interface(ql_int, qr_int,
                              gcl, gcr, cl, cr,
                              q_riemann, bnd_fac, idir);


      // compute the flux

      GpuArray<Real, NUM_STATE> F;

      compute_flux_q(q_riemann, F,
#ifdef RADIATION
                     lambda_int, rflx,
#endif
#ifdef HYBRID_MOMENTUM
                     geomdata, center,
#endif
                     coord,
                     idir, 0);

      // store the godunov state
      qgdnv(i,j,k,GDRHO) = q_riemann[QRHO];
      qgdnv(i,j,k,GDU) = q_riemann[QU];
      qgdnv(i,j,k,GDV) = q_riemann[QV];
      qgdnv(i,j,k,GDW) = q_riemann[QW];
      qgdnv(i,j,k,GDPRES) = q_riemann[QPRES];
#ifdef RADIATION
      for (int g = 0; g < NGROUPS; g++) {
        qgdnv(i,j,k,GDLAMS+g) = lambda(i,j,k,g);
        qgdnv(i,j,k,GDERADS+g) = qint(i,j,k,QRAD+g);
      }
#endif

      // store the flux
      for (int n = 0; n < NUM_STATE; n++) {
        flx(i,j,k,n) = F[n];
      }

    });

  } else if (riemann_solver == 2) {
    // HLLC
    HLLC(bx,
         ql_arr, qr_arr,
         qaux_arr,
         flx, qgdnv,  // wrong
         idir);

    // need to do this zone-by-zone and store the godunov state

#ifndef AMREX_USE_CUDA
  } else {
    amrex::Error("ERROR: invalid value of riemann_solver");
#endif
  }

  if (hybrid_riemann == 1) {
    // correct the fluxes using an HLL scheme if we are in a shock
    // and doing the hybrid approach


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
          ql_zone[n] = ql_arr(i,j,k,n);
          qr_zone[n] = qr_arr(i,j,k,n);
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

}


void
Castro::riemann_state(const Box& bx,
                      Array4<Real> const ql_arr,
                      Array4<Real> const qr_arr,
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

    GpuArray<Real, NQ> ql_int;
    GpuArray<Real, NQ> qr_int;

    Real gcl, gcr;
    Real cl, cr;

#ifdef RADIATION
    GpuArray<Real, Radiation::ngroups> laml_int;
    GpuArray<Real, Radiation::ngroups> lamr_int;
    Real gamcgl, gamcgr;
#endif

    get_riemann_input_states(ql_arr, qr_arr,
                             qaux_arr,
                             i, j, k, idir, 0,
                             gcl, gcr, cl, cr,
#ifdef RADIATION
                             laml_int, lamr_int,
                             gamcgl, gamcgr,
#endif
                             ql_int, qr_int);

    GpuArray<Real, NQ> q_riemann;

    riemann_state_interface(ql_int, qr_int,
                            gcl, gcr, cl, cr,
                            q_riemann, bnd_fac, idir);

    for (int n = 0; n < NQ; n++) {
      qint(i,j,k,n) = q_riemann[n];
    }

  });

}
