#include <Castro.H>
#include <Castro_util.H>

#ifdef RADIATION
#include <Radiation.H>
#endif

using namespace amrex;

void
Castro::consup_hydro(const Box& bx,
                     Array4<Real> const& U_new,
                     Array4<Real> const& flux0,
                     Array4<Real const> const& qx,
#if AMREX_SPACEDIM >= 2
                     Array4<Real> const& flux1,
                     Array4<Real const> const& qy,
#endif
#if AMREX_SPACEDIM == 3
                     Array4<Real> const& flux2,
                     Array4<Real const> const& qz,
#endif
                     const Real dt)
{
  auto geomdata = geom.data();

  amrex::ParallelFor(bx, NUM_STATE,
  [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
  {
    Real volinv = 1.0 / geometry_util::volume(i, j, k, geomdata);

    if (n == UMY && i == 4 && j == 0) {
        std::cout << "original velocities = " << U_new(i,j,k,UMX) << " " << U_new(i,j,k,UMY) << std::endl;
    }

    U_new(i,j,k,n) = U_new(i,j,k,n) + dt *
        ( flux0(i,  j,k,n) * geometry_util::area(i,   j, k, 0, geomdata)
        - flux0(i+1,j,k,n) * geometry_util::area(i+1, j, k, 0, geomdata)
#if AMREX_SPACEDIM >= 2
        + flux1(i,j,  k,n) * geometry_util::area(i, j,   k, 1, geomdata)
        - flux1(i,j+1,k,n) * geometry_util::area(i, j+1, k, 1, geomdata)
#endif
#if AMREX_SPACEDIM == 3
        + flux2(i,j,k,  n) * geometry_util::area(i, j, k,   2, geomdata)
        - flux2(i,j,k+1,n) * geometry_util::area(i, j, k+1, 2, geomdata)
#endif
        ) * volinv;

    // Add the p div(u) source term to (rho e).
    if (n == UEINT) {

      Real pdu = (qx(i+1,j,k,GDPRES) + qx(i,j,k,GDPRES)) *
          (qx(i+1,j,k,GDU) * geometry_util::area(i+1, j, k, 0, geomdata) -
           qx(i,  j,k,GDU) * geometry_util::area(i,   j, k, 0, geomdata));

#if AMREX_SPACEDIM >= 2
      pdu += (qy(i,j+1,k,GDPRES) + qy(i,j,k,GDPRES)) *
          (qy(i,j+1,k,GDV) * geometry_util::area(i, j+1, k, 1, geomdata) -
           qy(i,j,  k,GDV) * geometry_util::area(i, j,   k, 1, geomdata));
#endif

#if AMREX_SPACEDIM == 3
      pdu += (qz(i,j,k+1,GDPRES) + qz(i,j,k,GDPRES)) *
          (qz(i,j,k+1,GDW) * geometry_util::area(i, j, k+1, 2, geomdata) -
           qz(i,j,k  ,GDW) * geometry_util::area(i, j, k,   2, geomdata));
#endif

      pdu = 0.5 * pdu * volinv;

      U_new(i,j,k,n) = U_new(i,j,k,n) - dt * pdu;

    } else if (n == UMX && !mom_flux_has_p(0, 0, geomdata.Coord())) {
        // Add gradp term to momentum equation -- only for axisymmetric
        // coords (and only for the radial flux).

        U_new(i,j,k,UMX) += - dt * (qx(i+1,j,k,GDPRES) - qx(i,j,k,GDPRES)) / geomdata.CellSize(0);

#if AMREX_SPACEDIM >= 2
    } else if (n == UMY && !mom_flux_has_p(1, 1, geomdata.Coord())) {
        // Add gradp term to polar(theta) momentum equation for Spherical 2D geometry

        Real r = geomdata.ProbLo(0) + (static_cast<Real>(i) + 0.5_rt) * geomdata.CellSize(0);
        //U_new(i,j,k,UMY) += - dt * (qy(i,j+1,k,GDPRES) - qy(i,j,k,GDPRES)) / (r * geomdata.CellSize(1));
#endif
    }

    if (n == UMY && i == 4 && j == 0) {
        std::cout << "new velocities = " << U_new(i,j,k,UMX) << " " << U_new(i,j,k,UMY) << std::endl;
        std::cout << "x fluxes: " << flux0(i, j, k, UMY) << " " << flux1(i+1, j, k, UMY) << std::endl;
        std::cout << "y fluxes: " << flux1(i, j, k, UMY) << " " << flux1(i, j+1, k, UMY) << std::endl;
        std::cout << "y pressures: " << qy(i,j,k,GDPRES) << " " << qy(i,j+1,k,GDPRES) << std::endl;
        //amrex::Error("stop");
    }

  });
}


void
Castro::ctu_ppm_states(const Box& bx, const Box& vbx,
                       Array4<Real const> const& U_arr,
                       Array4<Real const> const& rho_inv_arr,
                       Array4<Real const> const& q_arr,
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
                       Array4<Real const> const& dlogaX,
#endif
#if AMREX_SPACEDIM == 2
                       Array4<Real const> const& dlogaY,
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
                U_arr, rho_inv_arr, q_arr, qaux_arr, srcQ,
                qxm, qxp,
#if AMREX_SPACEDIM <= 2
                dlogaX,
#endif
                vbx, dt);

      enforce_reflect_states(bx, 0, qxm, qxp);

#if AMREX_SPACEDIM >= 2
    } else if (idir == 1) {
      trace_ppm(bx,
                idir,
                U_arr, rho_inv_arr, q_arr, qaux_arr, srcQ,
                qym, qyp,
#if AMREX_SPACEDIM <= 2
                dlogaY,
#endif
                vbx, dt);

      enforce_reflect_states(bx, 1, qym, qyp);
#endif


#if AMREX_SPACEDIM == 3
    } else {
      trace_ppm(bx,
                idir,
                U_arr, rho_inv_arr, q_arr, qaux_arr, srcQ,
                qzm, qzp,
                vbx, dt);

      enforce_reflect_states(bx, 2, qzm, qzp);
#endif
    }
  }
}


#ifdef RADIATION
void
Castro::ctu_ppm_rad_states(const Box& bx, const Box& vbx,
                           Array4<Real const> const& U_arr,
                           Array4<Real const> const& rho_inv_arr,
                           Array4<Real const> const& q_arr,
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
                           Array4<Real const> const& dlogaX,
#endif
#if AMREX_SPACEDIM == 2
                           Array4<Real const> const& dlogaY,
#endif
                           const Real dt) {


  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {

    if (idir == 0) {

      trace_ppm_rad(bx,
                    idir,
                    U_arr, rho_inv_arr, q_arr, qaux_arr, srcQ,
                    qxm, qxp,
#if AMREX_SPACEDIM <= 2
                    dlogaX,
#endif
                    vbx, dt);

      enforce_reflect_states(bx, 0, qxm, qxp);

#if AMREX_SPACEDIM >= 2
    } else if (idir == 1) {
      trace_ppm_rad(bx,
                    idir,
                    U_arr, rho_inv_arr, q_arr, qaux_arr, srcQ,
                    qym, qyp,
#if AMREX_SPACEDIM <= 2
                    dlogaY,
#endif
                    vbx, dt);

      enforce_reflect_states(bx, 1, qym, qyp);
#endif


#if AMREX_SPACEDIM == 3
    } else {
      trace_ppm_rad(bx,
                    idir,
                    U_arr, rho_inv_arr, q_arr, qaux_arr, srcQ,
                    qzm, qzp,
                    vbx, dt);

      enforce_reflect_states(bx, 2, qzm, qzp);

#endif
    }
  }
}
#endif


void
Castro::ctu_plm_states(const Box& bx, const Box& vbx,
                       Array4<Real const> const& U_arr,
                       Array4<Real const> const& rho_inv_arr,
                       Array4<Real const> const& q_arr,
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
                       Array4<Real const> const& dlogaX,
#endif
#if AMREX_SPACEDIM == 2
                       Array4<Real const> const& dlogaY,
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
                U_arr, rho_inv_arr, q_arr, qaux_arr,
                qxm, qxp,
#if AMREX_SPACEDIM < 3
                dlogaX,
#endif
                srcQ, vbx, dt);

      enforce_reflect_states(bx, 0, qxm, qxp);

#if AMREX_SPACEDIM >= 2
    } else if (idir == 1) {
      trace_plm(bx, 1,
                U_arr, rho_inv_arr, q_arr, qaux_arr,
                qym, qyp,
#if AMREX_SPACEDIM < 3
                dlogaY,
#endif
                srcQ, vbx, dt);

      enforce_reflect_states(bx, 1, qym, qyp);
#endif

#if AMREX_SPACEDIM == 3
    } else {
      trace_plm(bx, 2,
                U_arr, rho_inv_arr, q_arr, qaux_arr,
                qzm, qzp,
                srcQ, vbx, dt);

      enforce_reflect_states(bx, 2, qzm, qzp);
#endif
    }
  }
}


void
Castro::add_sdc_source_to_states(const Box& bx, const int idir, const Real dt,
                                 Array4<Real> const& qleft,
                                 Array4<Real> const& qright,
                                 Array4<Real const> const& sdc_src)
{

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {

        // of the state variables, only pressure, rhoe, and
        // composition have source terms

        Real p_old = qleft(i,j,k,QPRES);
        Real rhoe_old = qleft(i,j,k,QREINT);

        if (idir == 0) {
            qleft(i,j,k,QPRES) += 0.5 * dt * sdc_src(i-1,j,k,QPRES);
            qleft(i,j,k,QREINT) += 0.5 * dt * sdc_src(i-1,j,k,QREINT);
        } else if (idir == 1) {
            qleft(i,j,k,QPRES) += 0.5 * dt * sdc_src(i,j-1,k,QPRES);
            qleft(i,j,k,QREINT) += 0.5 * dt * sdc_src(i,j-1,k,QREINT);
        } else {
            qleft(i,j,k,QPRES) += 0.5 * dt * sdc_src(i,j,k-1,QPRES);
            qleft(i,j,k,QREINT) += 0.5 * dt * sdc_src(i,j,k-1,QREINT);
        }

        if (qleft(i,j,k,QPRES) < small_pres || qleft(i,j,k,QREINT) < small_dens * small_ener) {
            qleft(i,j,k,QPRES) = p_old;
            qleft(i,j,k,QREINT) = rhoe_old;
        }

        p_old = qright(i,j,k,QPRES);
        rhoe_old = qright(i,j,k,QREINT);

        qright(i,j,k,QPRES) += 0.5 * dt * sdc_src(i,j,k,QPRES);
        qright(i,j,k,QREINT) += 0.5 * dt * sdc_src(i,j,k,QREINT);

        if (qright(i,j,k,QPRES) < small_pres || qright(i,j,k,QREINT) < small_dens * small_ener) {
            qright(i,j,k,QPRES) = p_old;
            qright(i,j,k,QREINT) = rhoe_old;
        }


        for (int ipassive = 0; ipassive < npassive; ipassive++) {

            int n = qpassmap(ipassive);

            if (idir == 0) {
                qleft(i,j,k,n) += 0.5 * dt * sdc_src(i-1,j,k,n);
            } else if (idir == 1) {
                qleft(i,j,k,n) += 0.5 * dt * sdc_src(i,j-1,k,n);
            } else {
                qleft(i,j,k,n) += 0.5 * dt * sdc_src(i,j,k-1,n);
            }
            qright(i,j,k,n) += 0.5 * dt * sdc_src(i,j,k,n);

            if (n >= QFS && n <= QFS-1+NumSpec) {
                // mass fractions should be in [0, 1]
                qleft(i,j,k,n) = std::clamp(qleft(i,j,k,n), 0.0_rt, 1.0_rt);
                qright(i,j,k,n) = std::clamp(qright(i,j,k,n), 0.0_rt, 1.0_rt);
            }
        }

    });

}
