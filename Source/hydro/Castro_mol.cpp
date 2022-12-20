#include <Castro.H>
#include <Castro_F.H>

#ifdef RADIATION
#include <Radiation.H>
#endif

#ifdef GRAVITY
#include <Gravity.H>
#endif

#include <ppm.H>
#include <slope.H>

using namespace amrex;

void
Castro::mol_plm_reconstruct(const Box& bx,
                            const int idir,
                            Array4<Real const> const& q_arr,
                            Array4<Real const> const& flatn_arr,
                            Array4<Real const> const& src_q_arr,
                            Array4<Real> const& dq,
                            Array4<Real> const& qm,
                            Array4<Real> const& qp) {

  const auto dx = geom.CellSizeArray();

  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();

  bool lo_symm = lo_bc[idir] == Symmetry;
  bool hi_symm = hi_bc[idir] == Symmetry;

  const auto domlo = geom.Domain().loVect3d();
  const auto domhi = geom.Domain().hiVect3d();

  // piecewise linear slopes
  amrex::ParallelFor(bx, NQ,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int n)
  {

    bool lo_bc_test = lo_symm && ((idir == 0 && i == domlo[0]) ||
                                  (idir == 1 && j == domlo[1]) ||
                                  (idir == 2 && k == domlo[2]));

    bool hi_bc_test = hi_symm && ((idir == 0 && i == domhi[0]) ||
                                  (idir == 1 && j == domhi[1]) ||
                                  (idir == 2 && k == domhi[2]));

    Real s[5];
    Real flat = flatn_arr(i,j,k);

    load_stencil(q_arr, idir, i, j, k, n, s);

    // normal velocity?
    bool vtest = n == QU+idir;

    dq(i,j,k,n) = uslope(s, flat, lo_bc_test && vtest, hi_bc_test && vtest);
  });

  if (use_pslope == 1) {

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

      Real s[5];
      Real flat = flatn_arr(i,j,k);

      Real trho[5];
      Real src[5];

      bool lo_bc_test = lo_symm && ((idir == 0 && i == domlo[0]) ||
                                    (idir == 1 && j == domlo[1]) ||
                                    (idir == 2 && k == domlo[2]));

      bool hi_bc_test = hi_symm && ((idir == 0 && i == domhi[0]) ||
                                    (idir == 1 && j == domhi[1]) ||
                                    (idir == 2 && k == domhi[2]));

      load_stencil(q_arr, idir, i, j, k, QPRES, s);
      load_stencil(q_arr, idir, i, j, k, QRHO, trho);
      load_stencil(src_q_arr, idir, i, j, k, QU+idir, src);

      Real dp = dq(i,j,k,QPRES);
      pslope(trho, s, src, flat, lo_bc_test, hi_bc_test, dx[idir], dp);
      dq(i,j,k,QPRES) = dp;

    });
  }

  amrex::ParallelFor(bx, NQ,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int n)
  {


   // this is a loop over zones.  For each slope in the zone, fill the
   // two adjacent edge states (e.g., the right state at i-1/2 and the
   // left state at i+1/2

   if (idir == 0) {

     // left state at i+1/2 interface
     qm(i+1,j,k,n) = q_arr(i,j,k,n) + 0.5_rt*dq(i,j,k,n);

     // right state at i-1/2 interface
     qp(i,j,k,n) = q_arr(i,j,k,n) - 0.5_rt*dq(i,j,k,n);


#if AMREX_SPACEDIM >= 2
   } else if (idir == 1) {

     // left state at j+1/2 interface
     qm(i,j+1,k,n) = q_arr(i,j,k,n) + 0.5_rt*dq(i,j,k,n);

     // right state at j-1/2 interface
     qp(i,j,k,n) = q_arr(i,j,k,n) - 0.5_rt*dq(i,j,k,n);
#endif

#if AMREX_SPACEDIM == 3
   } else {

     // left state at k+1/2 interface
     qm(i,j,k+1,n) = q_arr(i,j,k,n) + 0.5_rt*dq(i,j,k,n);

     // right state at k-1/2 interface
     qp(i,j,k,n) = q_arr(i,j,k,n) - 0.5_rt*dq(i,j,k,n);

#endif
   }

  });


  // special care for reflecting BCs

  // we have to do this after the loops above, since here we will
  // consider interfaces, not zones

  enforce_reflect_states(bx, idir, qm, qp);

}


void
Castro::mol_ppm_reconstruct(const Box& bx,
                            const int idir,
                            Array4<Real const> const& q_arr,
                            Array4<Real const> const& flatn_arr,
                            Array4<Real> const& qm,
                            Array4<Real> const& qp) {

  amrex::ParallelFor(bx, NQ,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int n)
  {

    Real s[5];
    Real flat = flatn_arr(i,j,k);
    Real sm;
    Real sp;

    load_stencil(q_arr, idir, i, j, k, n, s);
    ppm_reconstruct(s, flat, sm, sp);

    if (idir == 0) {
      // right state at i-1/2
      qp(i,j,k,n) = sm;

      // left state at i+1/2
      qm(i+1,j,k,n) = sp;

    } else if (idir == 1) {
      // right state at j-1/2
      qp(i,j,k,n) = sm;

      // left state at j+1/2
      qm(i,j+1,k,n) = sp;

    } else {
      // right state at k-1/2
      qp(i,j,k,n) = sm;

      // left state at k+1/2
      qm(i,j,k+1,n) = sp;

    }

  });

  // special care for reflecting BCs

  // we have to do this after the loops above, since here we will
  // consider interfaces, not zones

  enforce_reflect_states(bx, idir, qm, qp);

}


void
Castro::mol_consup(const Box& bx,
#ifdef SHOCK_VAR
                   Array4<Real const> const& shk,
#endif
                   Array4<Real const> const& srcU,
                   Array4<Real> const& update,
                   const Real dt,
                   Array4<Real const> const& flux0,
#if AMREX_SPACEDIM >= 2
                   Array4<Real const> const& flux1,
#endif
#if AMREX_SPACEDIM == 3
                   Array4<Real const> const& flux2,
#endif
                   Array4<Real const> const& area0,
#if AMREX_SPACEDIM >= 2
                   Array4<Real const> const& area1,
#endif
#if AMREX_SPACEDIM == 3
                   Array4<Real const> const& area2,
#endif
#if AMREX_SPACEDIM <= 2
                   Array4<Real const> const& q0,
#endif
                   Array4<Real const> const& vol) {


  // For hydro, we will create an update source term that is
  // essentially the flux divergence.  This can be added with dt to
  // get the update

#if AMREX_SPACEDIM <= 2
  const auto dx = geom.CellSizeArray();
#endif

#if AMREX_SPACEDIM == 2
  auto coord = geom.Coord();
#endif

  amrex::ParallelFor(bx, NUM_STATE,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int n)
  {

#if AMREX_SPACEDIM == 1
    update(i,j,k,n) += (flux0(i,j,k,n) * area0(i,j,k) - flux0(i+1,j,k,n) * area0(i+1,j,k) ) / vol(i,j,k);

#elif AMREX_SPACEDIM == 2
    update(i,j,k,n) += (flux0(i,j,k,n) * area0(i,j,k) - flux0(i+1,j,k,n) * area0(i+1,j,k) +
                        flux1(i,j,k,n) * area1(i,j,k) - flux1(i,j+1,k,n) * area1(i,j+1,k) ) / vol(i,j,k);

#else
    update(i,j,k,n) += (flux0(i,j,k,n) * area0(i,j,k) - flux0(i+1,j,k,n) * area0(i+1,j,k) +
                        flux1(i,j,k,n) * area1(i,j,k) - flux1(i,j+1,k,n) * area1(i,j+1,k) +
                        flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j,k+1,n) * area2(i,j,k+1) ) / vol(i,j,k);
#endif

#if AMREX_SPACEDIM == 1
    if (do_hydro == 1) {
      if (n == UMX) {
        update(i,j,k,UMX) -= (q0(i+1,j,k,GDPRES) - q0(i,j,k,GDPRES)) / dx[0];
      }
    }
#endif

#if AMREX_SPACEDIM == 2
    if (do_hydro == 1) {
      if (n == UMX) {
        // add the pressure source term for axisymmetry
        if (coord > 0) {
          update(i,j,k,n) -= (q0(i+1,j,k,GDPRES) - q0(i,j,k,GDPRES)) / dx[0];
        }
      }
    }
#endif

    // this assumes that the species are at the end of the conserved state
    if (n < NSRC) {
      update(i,j,k,n) += srcU(i,j,k,n);
    }

  });

#ifdef SHOCK_VAR
  // We'll update the shock data for future use in the burning step.
  // For the update, we are starting from USHK == 0 (set at the
  // beginning of the timestep) and we need to divide by dt since
  // we'll be multiplying that for the update calculation.

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
  {
    update(i,j,k,USHK) = shk(i,j,k) / dt;
  });
#endif

}


void
Castro::mol_diffusive_flux(const Box& bx,
                           const int idir,
                           Array4<Real const> const& U,
                           Array4<Real const> const& cond,
                           Array4<Real> const& flux) {

  const auto dx = geom.CellSizeArray();

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
  {

    Real cond_int;
    Real diff_term;

    if (idir == 0) {
      cond_int = 0.5_rt * (cond(i,j,k) + cond(i-1,j,k));
      diff_term = -cond_int * (U(i,j,k,UTEMP) - U(i-1,j,k,UTEMP)) / dx[0];

    } else if (idir == 1) {
      cond_int = 0.5_rt * (cond(i,j,k) + cond(i,j-1,k));
      diff_term = -cond_int * (U(i,j,k,UTEMP) - U(i,j-1,k,UTEMP)) / dx[1];

    } else {
      cond_int = 0.5_rt * (cond(i,j,k) + cond(i,j,k-1));
      diff_term = -cond_int * (U(i,j,k,UTEMP) - U(i,j,k-1,UTEMP)) / dx[2];
    }

    flux(i,j,k,UEINT) += diff_term;
    flux(i,j,k,UEDEN) += diff_term;

  });
}
