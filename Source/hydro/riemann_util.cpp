#include <Castro.H>
#include <Castro_F.H>
#include <Castro_util.H>

#ifdef RADIATION
#include <Radiation.H>
#include <fluxlimiter.H>
#endif

#ifdef HYBRID_MOMENTUM
#include <hybrid.H>
#endif

#include <cmath>

using namespace amrex;

#include <riemann.H>

void
Castro::compute_flux_q(const Box& bx,
                       Array4<Real const> const& qint,
                       Array4<Real> const& F,
#ifdef RADIATION
                       Array4<Real const> const& lambda,
                       Array4<Real> const& rF,
#endif
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

#ifdef RADIATION
  int fspace_t = Radiation::fspace_advection_type;
  int comov = Radiation::comoving;
  int limiter = Radiation::limiter;
  int closure = Radiation::closure;
#endif

  const Real lT_guess = T_guess;

#ifdef HYBRID_MOMENTUM
  GeometryData geomdata = geom.data();
#endif

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
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
#ifdef RADIATION
    for (int g = 0; g < NGROUPS; g++) {
        qgdnv_zone[GDLAMS+g] = lambda(i,j,k,g);
        qgdnv_zone[GDERADS+g] = qint(i,j,k,QRAD+g);
    }
#endif
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
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
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
