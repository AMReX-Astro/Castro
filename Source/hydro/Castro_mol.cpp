#include <Castro.H>

#ifdef RADIATION
#include <Radiation.H>
#endif

#ifdef GRAVITY
#include <Gravity.H>
#endif

#ifdef HYBRID_MOMENTUM
#include <hybrid.H>
#endif

#include <ppm.H>
#include <slope.H>

using namespace amrex;
using namespace reconstruction;

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

  bool lo_symm = lo_bc[idir] == amrex::PhysBCType::symmetry;
  bool hi_symm = hi_bc[idir] == amrex::PhysBCType::symmetry;

  const auto domlo = geom.Domain().loVect3d();
  const auto domhi = geom.Domain().hiVect3d();

  // piecewise linear slopes
  amrex::ParallelFor(bx, NQ,
  [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {

    bool lo_bc_test = lo_symm && ((idir == 0 && i == domlo[0]) ||
                                  (idir == 1 && j == domlo[1]) ||
                                  (idir == 2 && k == domlo[2]));

    bool hi_bc_test = hi_symm && ((idir == 0 && i == domhi[0]) ||
                                  (idir == 1 && j == domhi[1]) ||
                                  (idir == 2 && k == domhi[2]));

    Real s[nslp];
    Real flat = flatn_arr(i,j,k);

    load_stencil(q_arr, idir, i, j, k, n, s);

    // normal velocity?
    bool vtest = n == QU+idir;

    dq(i,j,k,n) = uslope(s, flat, lo_bc_test && vtest, hi_bc_test && vtest);
  });

  if (use_pslope == 1) {

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {

      Real s[nslp];
      Real flat = flatn_arr(i,j,k);

      Real trho[nslp];
      Real src[nslp];

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
  [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
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
  [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {

    Real s[nslp];
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
Castro::mol_consup(const Box& bx,  // NOLINT(readability-convert-member-functions-to-static)
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
#if AMREX_SPACEDIM == 2
                   Array4<Real const> const& q1,
#endif
#endif
                   Array4<Real const> const& vol) {

  amrex::ignore_unused(dt);

  // For hydro, we will create an update source term that is
  // essentially the flux divergence.  This can be added with dt to
  // get the update

#if AMREX_SPACEDIM <= 2
  const auto dx = geom.CellSizeArray();
  auto coord = geom.Coord();
  auto prob_lo = geom.ProbLoArray();
#endif

  amrex::ParallelFor(bx, NUM_STATE,
  [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
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

#if AMREX_SPACEDIM <= 2
    if (do_hydro == 1) {
        if (n == UMX && !mom_flux_has_p(0, 0, coord)) {
            // Add gradp term to radial momentum equation -- only for axisymmetric
            // coords.

            update(i,j,k,UMX) -= (q0(i+1,j,k,GDPRES) - q0(i,j,k,GDPRES)) / dx[0];

#if AMREX_SPACEDIM == 2
        } else if (n == UMY && !mom_flux_has_p(1, 1, coord)) {
            // Add gradp term to polar(theta) momentum equation for Spherical 2D geometry

            Real r = prob_lo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0];
            update(i,j,k,UMY) -= (q1(i,j+1,k,GDPRES) - q1(i,j,k,GDPRES)) / (r * dx[1]);
#endif
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
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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


void
Castro::compute_flux_from_q(const Box& bx,
                            Array4<Real const> const& qint,
                            Array4<Real> const& F,
                            const int idir) {

  // given a primitive state, compute the flux in direction idir
  //

  int iu, iv1, iv2;
  int im1, im2, im3;

  auto coord = geom.Coord();
  auto mom_check = mom_flux_has_p(idir, idir, coord);

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

#ifdef HYBRID_MOMENTUM
  GeometryData geomdata = geom.data();
#endif

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    Real u_adv = qint(i,j,k,iu);
    Real rhoeint = qint(i,j,k,QREINT);

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

#ifdef NSE_NET
    F(i,j,k,UMUP) = 0.0;
    F(i,j,k,UMUN) = 0.0;
#endif
    // passively advected quantities
    for (int ipassive = 0; ipassive < npassive; ipassive++) {
      int n  = upassmap(ipassive);
      int nqp = qpassmap(ipassive);

      F(i,j,k,n) = F(i,j,k,URHO)*qint(i,j,k,nqp);
    }

#ifdef HYBRID_MOMENTUM
    // the hybrid routine uses the Godunov indices, not the full NQ state
    GpuArray<Real, NGDNV> qgdnv_zone;
    qgdnv_zone[GDRHO] = qint(i,j,k,QRHO);
    qgdnv_zone[GDU] = qint(i,j,k,QU);
    qgdnv_zone[GDV] = qint(i,j,k,QV);
    qgdnv_zone[GDW] = qint(i,j,k,QW);
    qgdnv_zone[GDPRES] = qint(i,j,k,QPRES);
    GpuArray<Real, NUM_STATE> F_zone;
    for (int n = 0; n < NUM_STATE; n++) {
        F_zone[n] = F(i,j,k,n);
    }
    compute_hybrid_flux(qgdnv_zone, geomdata, idir, i, j, k, F_zone);
    for (int n = 0; n < NUM_STATE; n++) {
        F(i,j,k,n) = F_zone[n];
    }
#endif
  });
}

void
Castro::store_godunov_state(const Box& bx,
                            Array4<Real const> const& qint,
#ifdef RADIATION
                            Array4<Real const> const& lambda,
#endif
                            Array4<Real> const& qgdnv) {

  // this copies the full interface state (NQ -- one for each primitive
  // variable) over to a smaller subset of size NGDNV for use later in the
  // hydro advancement.

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {


#ifdef HYBRID_MOMENTUM
    qgdnv(i,j,k,GDRHO) = qint(i,j,k,QRHO);
#endif
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
