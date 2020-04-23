#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"
#include "Castro_util.H"

#ifdef RADIATION
#include "Radiation.H"
#include "fluxlimiter.H"
#endif

#ifdef HYBRID_MOMENTUM
#include "hybrid.H"
#endif

#include <eos.H>

#include <cmath>

using namespace amrex;

#include <riemann.H>

AMREX_GPU_HOST_DEVICE
void
Castro::compute_flux_q(GpuArray<Real, NQ>& q_riemann,
                       GpuArray<Real, NUM_STATE>& F,
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

  Real u_adv = q_riemann[iu];
  Real rhoeint = q_riemann[QREINT];

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
