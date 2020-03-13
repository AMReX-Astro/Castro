#include "Castro.H"
#include "Castro_F.H"
#include "Castro_util.H"
#include "Castro_hydro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

using namespace amrex;

void
Castro::consup_hydro(const Box& bx,
                     Array4<Real const> const shk,
                     Array4<Real> const update,
                     Array4<Real> const flux0,
                     Array4<Real const> const qx,
                     Array4<Real const> const area0,
#if AMREX_SPACEDIM >= 2
                     Array4<Real> const flux1,
                     Array4<Real const> const qy,
                     Array4<Real const> const area1,
#endif
#if AMREX_SPACEDIM == 3
                     Array4<Real> const flux2,
                     Array4<Real const> const qz,
                     Array4<Real const> const area2,
#endif
                     Array4<Real const> const vol,
                     const Real dt)
{


  const auto dx = geom.CellSizeArray();

  // For hydro, we will create an update source term that is
  // essentially the flux divergence.  This can be added with dt to
  // get the update

  int coord = geom.Coord();

  GpuArray<Real, 3> center;
  ca_get_center(center.begin());

  AMREX_PARALLEL_FOR_4D(bx, NUM_STATE, i, j, k, n,
  {

    Real volinv = 1.0 / vol(i,j,k);

    update(i,j,k,n) = update(i,j,k,n) +
      ( flux0(i,j,k,n) * area0(i,j,k) - flux0(i+1,j,k,n) * area0(i+1,j,k)
#if AMREX_SPACEDIM >= 2
      + flux1(i,j,k,n) * area1(i,j,k) - flux1(i,j+1,k,n) * area1(i,j+1,k)
#endif
#if AMREX_SPACEDIM == 3
      + flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j,k+1,n) * area2(i,j,k+1)
#endif
        ) * volinv;

   // Add the p div(u) source term to (rho e).
    if (n == UEINT) {

      Real pdu = (qx(i+1,j,k,GDPRES) + qx(i,j,k,GDPRES)) *
                 (qx(i+1,j,k,GDU) * area0(i+1,j,k) - qx(i,j,k,GDU) * area0(i,j,k));

#if AMREX_SPACEDIM >= 2
      pdu += (qy(i,j+1,k,GDPRES) + qy(i,j,k,GDPRES)) *
             (qy(i,j+1,k,GDV) * area1(i,j+1,k) - qy(i,j,k,GDV) * area1(i,j,k));
#endif

#if AMREX_SPACEDIM == 3
      pdu += (qz(i,j,k+1,GDPRES) + qz(i,j,k,GDPRES)) *
             (qz(i,j,k+1,GDW) * area2(i,j,k+1) - qz(i,j,k,GDW) * area2(i,j,k));
#endif

      pdu = 0.5 * pdu * volinv;

      update(i,j,k,n) = update(i,j,k,n) - pdu;

#ifdef SHOCK_VAR
    } else if (n == USHK) {
      update(i,j,k,USHK) = shk(i,j,k) / dt;
#endif

#ifndef RADIATION
    } else if (n == UMX) {
      // Add gradp term to momentum equation -- only for axisymmetric
      // coords (and only for the radial flux).

      if (!mom_flux_has_p(0, 0, coord)) {
        update(i,j,k,UMX) += - (qx(i+1,j,k,GDPRES) - qx(i,j,k,GDPRES)) / dx[0];
      }
#endif
    }
  });
}


void
Castro::ctu_ppm_states(const Box& bx, const Box& vbx,
                       Array4<Real const> const q_arr,
                       Array4<Real const> const flatn,
                       Array4<Real const> const qaux_arr,
                       Array4<Real const> const srcQ,
                       Array4<Real> const qxm,
                       Array4<Real> const qxp,
#if AMREX_SPACEDIM >= 2
                       Array4<Real> const qym,
                       Array4<Real> const qyp,
#endif
#if AMREX_SPACEDIM == 3
                       Array4<Real> const qzm,
                       Array4<Real> const qzp,
#endif
#if AMREX_SPACEDIM < 3
                       Array4<Real const> const dloga,
#endif
                       const Real dt) {

  // Compute the normal interface states by reconstructing
  // the primitive variables using the piecewise parabolic method
  // and doing characteristic tracing.  We do not apply the
  // transverse terms here.

  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {

    if (idir == 0) {
      trace_ppm(bx,
                idir,
                q_arr, qaux_arr, srcQ, flatn,
                qxm, qxp,
#if AMREX_SPACEDIM <= 2
                dloga,
#endif
                vbx, dt);

#if AMREX_SPACEDIM >= 2
    } else if (idir == 1) {
      trace_ppm(bx,
                idir,
                q_arr, qaux_arr, srcQ, flatn,
                qym, qyp,
#if AMREX_SPACEDIM <= 2
                dloga,
#endif
                vbx, dt);
#endif

#if AMREX_SPACEDIM == 3
    } else {
      trace_ppm(bx,
                idir,
                q_arr, qaux_arr, srcQ, flatn,
                qzm, qzp,
                vbx, dt);

#endif
    }
  }
}


#ifdef RADIATION
void
Castro::ctu_ppm_rad_states(const Box& bx, const Box& vbx,
                           Array4<Real const> const q_arr,
                           Array4<Real const> const flatn,
                           Array4<Real const> const qaux_arr,
                           Array4<Real const> const srcQ,
                           Array4<Real> const qxm,
                           Array4<Real> const qxp,
#if AMREX_SPACEDIM >= 2
                           Array4<Real> const qym,
                           Array4<Real> const qyp,
#endif
#if AMREX_SPACEDIM == 3
                           Array4<Real> const qzm,
                           Array4<Real> const qzp,
#endif
#if AMREX_SPACEDIM < 3
                           Array4<Real const> const dloga,
#endif
                           const Real dt) {


  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {

    if (idir == 0) {

      trace_ppm_rad(bx,
                    idir,
                    q_arr, qaux_arr, srcQ, flatn,
                    qxm, qxp,
#if AMREX_SPACEDIM <= 2
                    dloga,
#endif
                    vbx, dt);

#if AMREX_SPACEDIM >= 2
    } else if (idir == 1) {
      trace_ppm_rad(bx,
                    idir,
                    q_arr, qaux_arr, srcQ, flatn,
                    qym, qyp,
#if AMREX_SPACEDIM <= 2
                    dloga,
#endif
                    vbx, dt);
#endif

#if AMREX_SPACEDIM == 3
    } else {
      trace_ppm_rad(bx,
                    idir,
                    q_arr, qaux_arr, srcQ, flatn,
                    qzm, qzp,
                    vbx, dt);

#endif
    }
  }
}
#endif


void
Castro::ctu_plm_states(const Box& bx, const Box& vbx,
                       Array4<Real const> const q_arr,
                       Array4<Real const> const flatn_arr,
                       Array4<Real const> const qaux_arr,
                       Array4<Real const> const srcQ,
                       Array4<Real> const dq,
                       Array4<Real> const qxm,
                       Array4<Real> const qxp,
#if AMREX_SPACEDIM >= 2
                       Array4<Real> const qym,
                       Array4<Real> const qyp,
#endif
#if AMREX_SPACEDIM == 3
                       Array4<Real> const qzm,
                       Array4<Real> const qzp,
#endif
#if AMREX_SPACEDIM < 3
                       Array4<Real const> const dloga,
#endif
                       const Real dt) {

  // Compute the normal interface states by reconstructing
  // the primitive variables using piecewise linear slopes and doing
  // characteristic tracing.  We do not apply the transverse terms here.
  //
  // .. todo::
  //    we can get rid of the the different temporary q Godunov
  //    state arrays
  //


  Real hdt = 0.5_rt * dt;

#ifdef RADIATION
#ifndef AMREX_USE_CUDA
  amrex::Error("ppm_type <=0 is not supported in with radiation");
#endif
#endif

  // Compute all slopes
  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {

    for (int n = 0; n < NQ; n++) {
      if (n == QTEMP)
        continue;

      uslope(bx, idir,
             q_arr, n,
             flatn_arr, dq);
    }

    if (use_pslope == 1) {
      pslope(bx, idir,
             q_arr,
             flatn_arr, dq, srcQ);
    }

    // compute the interface states

    if (idir == 0) {
      trace_plm(bx, 0,
                q_arr, qaux_arr, dq,
                qxm, qxp,
#if AMREX_SPACEDIM < 3
                dloga,
#endif
                srcQ, vbx, dt);

#if AMREX_SPACEDIM >= 2
    } else if (idir == 1) {
      trace_plm(bx, 1,
                q_arr, qaux_arr, dq,
                qym, qyp,
#if AMREX_SPACEDIM < 3
                dloga,
#endif
                srcQ, vbx, dt);
#endif

#if AMREX_SPACEDIM == 3
    } else {
      trace_plm(bx, 2,
                q_arr, qaux_arr, dq,
                qzm, qzp,
                srcQ, vbx, dt);
#endif
    }
  }
}
