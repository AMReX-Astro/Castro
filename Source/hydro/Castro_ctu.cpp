#include <Castro.H>
#include <Castro_F.H>
#include <Castro_util.H>

#ifdef RADIATION
#include <Radiation.H>
#endif

using namespace amrex;

void
Castro::consup_hydro(const Box& bx,
#ifdef SHOCK_VAR
                     Array4<Real const> const& shk,
#endif
                     Array4<Real> const& update,
                     Array4<Real> const& flux0,
                     Array4<Real const> const& qx,
                     Array4<Real const> const& area0,
#if AMREX_SPACEDIM >= 2
                     Array4<Real> const& flux1,
                     Array4<Real const> const& qy,
                     Array4<Real const> const& area1,
#endif
#if AMREX_SPACEDIM == 3
                     Array4<Real> const& flux2,
                     Array4<Real const> const& qz,
                     Array4<Real const> const& area2,
#endif
                     Array4<Real const> const& vol,
                     const Real dt)
{


  const auto dx = geom.CellSizeArray();

  // For hydro, we will create an update source term that is
  // essentially the flux divergence.  This can be added with dt to
  // get the update

  int coord = geom.Coord();

  amrex::ParallelFor(bx, NUM_STATE,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int n)
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

    } else if (n == UMX) {
      // Add gradp term to momentum equation -- only for axisymmetric
      // coords (and only for the radial flux).

      if (!mom_flux_has_p(0, 0, coord)) {
        update(i,j,k,UMX) += - (qx(i+1,j,k,GDPRES) - qx(i,j,k,GDPRES)) / dx[0];
      }
    }
  });
}


void
Castro::ctu_ppm_states(const Box& bx, const Box& vbx,
                       Array4<Real const> const& q_arr,
                       Array4<Real const> const& flatn,
                       Array4<Real const> const& qaux_arr,
                       Array4<Real const> const& srcQ,
                       Array4<Real> const& qxm,
                       Array4<Real> const& qxp,
#if AMREX_SPACEDIM >= 2
                       Array4<Real> const& qym,
                       Array4<Real> const& qyp,
#endif
#if AMREX_SPACEDIM == 3
                       Array4<Real> const& qzm,
                       Array4<Real> const& qzp,
#endif
#if AMREX_SPACEDIM < 3
                       Array4<Real const> const& dloga,
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
                           Array4<Real const> const& q_arr,
                           Array4<Real const> const& flatn,
                           Array4<Real const> const& qaux_arr,
                           Array4<Real const> const& srcQ,
                           Array4<Real> const& qxm,
                           Array4<Real> const& qxp,
#if AMREX_SPACEDIM >= 2
                           Array4<Real> const& qym,
                           Array4<Real> const& qyp,
#endif
#if AMREX_SPACEDIM == 3
                           Array4<Real> const& qzm,
                           Array4<Real> const& qzp,
#endif
#if AMREX_SPACEDIM < 3
                           Array4<Real const> const& dloga,
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
                       Array4<Real const> const& q_arr,
                       Array4<Real const> const& flatn_arr,
                       Array4<Real const> const& qaux_arr,
                       Array4<Real const> const& srcQ,
                       Array4<Real> const& qxm,
                       Array4<Real> const& qxp,
#if AMREX_SPACEDIM >= 2
                       Array4<Real> const& qym,
                       Array4<Real> const& qyp,
#endif
#if AMREX_SPACEDIM == 3
                       Array4<Real> const& qzm,
                       Array4<Real> const& qzp,
#endif
#if AMREX_SPACEDIM < 3
                       Array4<Real const> const& dloga,
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


#ifdef RADIATION
#ifndef AMREX_USE_GPU
  amrex::Error("ppm_type <=0 is not supported in with radiation");
#endif
#endif

  // Compute all slopes
  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {

    // compute the interface states

    if (idir == 0) {
      trace_plm(bx, 0,
                q_arr, qaux_arr, flatn_arr,
                qxm, qxp,
#if AMREX_SPACEDIM < 3
                dloga,
#endif
                srcQ, vbx, dt);

#if AMREX_SPACEDIM >= 2
    } else if (idir == 1) {
      trace_plm(bx, 1,
                q_arr, qaux_arr, flatn_arr,
                qym, qyp,
#if AMREX_SPACEDIM < 3
                dloga,
#endif
                srcQ, vbx, dt);
#endif

#if AMREX_SPACEDIM == 3
    } else {
      trace_plm(bx, 2,
                q_arr, qaux_arr, flatn_arr,
                qzm, qzp,
                srcQ, vbx, dt);
#endif
    }

    // special care for reflecting BCs
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();

    const auto domlo = geom.Domain().loVect3d();
    const auto domhi = geom.Domain().hiVect3d();

    bool lo_bc_test = lo_bc[idir] == Symmetry;
    bool hi_bc_test = hi_bc[idir] == Symmetry;

    // we have to do this after the loops above, since here we will
    // consider interfaces, not zones

    if (idir == 0) {
      if (lo_bc_test) {

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {

          // reset the left state at domlo(0) if needed -- it is outside the domain
          if (i == domlo[0]) {
            for (int n = 0; n < NQ; n++) {
              if (n == QU) {
                qxm(i,j,k,QU) = -qxp(i,j,k,QU);
              } else {
                qxm(i,j,k,n) = qxp(i,j,k,n);
              }
            }
          }
        });

      }


      if (hi_bc_test) {

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {

          // reset the right state at domhi(0)+1 if needed -- it is outside the domain
          if (i == domhi[0]+1) {
            for (int n = 0; n < NQ; n++) {
              if (n == QU) {
                qxp(i,j,k,QU) = -qxm(i,j,k,QU);
              } else {
                qxp(i,j,k,n) = qxm(i,j,k,n);
              }
            }
          }
        });

      }

#if AMREX_SPACEDIM >= 2
    } else if (idir == 1) {

      if (lo_bc_test) {

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {

          // reset the left state at domlo(0) if needed -- it is outside the domain
          if (j == domlo[1]) {
            for (int n = 0; n < NQ; n++) {
              if (n == QV) {
                qym(i,j,k,QV) = -qyp(i,j,k,QV);
              } else {
                qym(i,j,k,n) = qyp(i,j,k,n);
              }
            }
          }
        });

      }


      if (hi_bc_test) {

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {

          // reset the right state at domhi(0)+1 if needed -- it is outside the domain
          if (j == domhi[1]+1) {
            for (int n = 0; n < NQ; n++) {
              if (n == QV) {
                qyp(i,j,k,QV) = -qym(i,j,k,QV);
              } else {
                qyp(i,j,k,n) = qym(i,j,k,n);
              }
            }
          }
        });

      }

#endif
#if AMREX_SPACEDIM == 3
    } else {

      if (lo_bc_test) {

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {

          // reset the left state at domlo(0) if needed -- it is outside the domain
          if (k == domlo[2]) {
            for (int n = 0; n < NQ; n++) {
              if (n == QW) {
                qzm(i,j,k,QW) = -qzp(i,j,k,QW);
              } else {
                qzm(i,j,k,n) = qzp(i,j,k,n);
              }
            }
          }
        });

      }


      if (hi_bc_test) {

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {

          // reset the right state at domhi(0)+1 if needed -- it is outside the domain
          if (k == domhi[2]+1) {
            for (int n = 0; n < NQ; n++) {
              if (n == QW) {
                qzp(i,j,k,QW) = -qzm(i,j,k,QW);
              } else {
                qzp(i,j,k,n) = qzm(i,j,k,n);
              }
            }
          }
        });

      }
#endif

    }
  }
}


void
Castro::src_to_prim(const Box& bx,
                    Array4<Real const> const& q_arr,
                    Array4<Real const> const& src,
                    Array4<Real> const& srcQ)
{

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
  {


      for (int n = 0; n < NQSRC; ++n) {
        srcQ(i,j,k,n) = 0.0_rt;
      }

      Real rhoinv = 1.0_rt / q_arr(i,j,k,QRHO);

      // get the needed derivatives
      eos_rep_t eos_state;
      eos_state.T = q_arr(i,j,k,QTEMP);
      eos_state.rho = q_arr(i,j,k,QRHO);
      eos_state.e = q_arr(i,j,k,QREINT) * rhoinv;
      for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n]  = q_arr(i,j,k,QFS+n);
      }
#if NAUX_NET > 0
      for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = q_arr(i,j,k,QFX+n);
      }
#endif

      eos(eos_input_re, eos_state);


      srcQ(i,j,k,QRHO) = src(i,j,k,URHO);
      srcQ(i,j,k,QU) = (src(i,j,k,UMX) - q_arr(i,j,k,QU) * srcQ(i,j,k,QRHO)) * rhoinv;
      srcQ(i,j,k,QV) = (src(i,j,k,UMY) - q_arr(i,j,k,QV) * srcQ(i,j,k,QRHO)) * rhoinv;
      srcQ(i,j,k,QW) = (src(i,j,k,UMZ) - q_arr(i,j,k,QW) * srcQ(i,j,k,QRHO)) * rhoinv;
      srcQ(i,j,k,QREINT) = src(i,j,k,UEINT);
      srcQ(i,j,k,QPRES ) = eos_state.dpde *
        (srcQ(i,j,k,QREINT) - q_arr(i,j,k,QREINT) * srcQ(i,j,k,QRHO)*rhoinv) *
        rhoinv + eos_state.dpdr_e * srcQ(i,j,k,QRHO);

#ifdef PRIM_SPECIES_HAVE_SOURCES
      for (int ipassive = 0; ipassive < npassive; ++ipassive) {
        int n = upassmap(ipassive);
        int iq = qpassmap(ipassive);

       // we may not be including the ability to have species sources,
       //  so check to make sure that we are < NQSRC
        srcQ(i,j,k,iq) = (src(i,j,k,n) - q_arr(i,j,k,iq) * srcQ(i,j,k,QRHO) ) /
          q_arr(i,j,k,QRHO);
      }
#endif
  });

}
