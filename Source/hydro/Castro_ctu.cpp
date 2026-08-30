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
        U_new(i,j,k,UMY) += - dt * (qy(i,j+1,k,GDPRES) - qy(i,j,k,GDPRES)) / (r * geomdata.CellSize(1));
#endif
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

        int isrc = i;
        int jsrc = j;
        int ksrc = k;

        // for left face, move one cell left in that dir
        if (idir == 0) {
            isrc -= 1;
        } else if (idir == 1) {
            jsrc -= 1;
        } else {
            ksrc -= 1;
        }

        const Real half_dt = 0.5_rt * dt;
        const Real rhoe_floor = small_dens * small_ener; // floor for rhoe in castro

        // no damping is default unless needed
        Real damping_factor_left = 1.0_rt;
        Real damping_factor_right = 1.0_rt;

        int limiting_comp_left = -1;
        int limiting_comp_right = -1;

        // helper function to update the current damping factor to keep the
        // smallest allowed value only
        auto update_damping = [] AMREX_GPU_DEVICE (Real& damping_factor,
                                                   int& limiting_comp,
                                                   Real trial,
                                                   int comp) noexcept
        {
            if (trial < damping_factor) {
                damping_factor = trial;
                limiting_comp = comp;
            }
        };

        // current state pressure and rhoe
        Real p_old = qleft(i,j,k,QPRES);
        Real rhoe_old = qleft(i,j,k,QREINT);

        // source terms 
        Real p_src = half_dt * sdc_src(isrc,jsrc,ksrc,QPRES);
        Real rhoe_src = half_dt * sdc_src(isrc,jsrc,ksrc,QREINT);

        // damping factor for pressure source term if needed
        if (p_src < 0.0_rt) {
            Real trial = amrex::max(0.0_rt, (p_old - small_pres) / (-p_src));
            update_damping(damping_factor_left, limiting_comp_left, trial, QPRES);
        }

        // damping factor for rho source term if needed
        if (rhoe_src < 0.0_rt) {
            Real trial = amrex::max(0.0_rt, (rhoe_old - rhoe_floor) / (-rhoe_src));
            update_damping(damping_factor_left, limiting_comp_left, trial, QREINT);
        }

        for (int n = 0; n < NumSpec; ++n) {
            Real X_old = qleft(i,j,k,QFS+n);
            Real X_src = half_dt * sdc_src(isrc,jsrc,ksrc,QFS+n);

            // daping factor for species source term if needed 
            if (X_src < 0.0_rt) {
                Real trial = amrex::max(0.0_rt, X_old / (-X_src));
                update_damping(damping_factor_left, limiting_comp_left, trial, QFS+n);
            } else if (X_src > 0.0_rt) {
                Real trial = amrex::max(0.0_rt, (1.0_rt - X_old) / X_src);
                update_damping(damping_factor_left, limiting_comp_left, trial, QFS+n);
            }
        }

        damping_factor_left = amrex::Clamp(damping_factor_left, 0.0_rt, 1.0_rt);

        // new state with the damping factor for left face
        qleft(i,j,k,QPRES) = p_old + damping_factor_left * p_src;
        qleft(i,j,k,QREINT) = rhoe_old + damping_factor_left * rhoe_src;

        // we apply the same mechanisms to the right face
        p_old = qright(i,j,k,QPRES);
        rhoe_old = qright(i,j,k,QREINT);

        p_src = half_dt * sdc_src(i,j,k,QPRES);
        rhoe_src = half_dt * sdc_src(i,j,k,QREINT);

        if (p_src < 0.0_rt) {
            Real trial = amrex::max(0.0_rt, (p_old - small_pres) / (-p_src));
            update_damping(damping_factor_right, limiting_comp_right, trial, QPRES);
        }

        if (rhoe_src < 0.0_rt) {
            Real trial = amrex::max(0.0_rt, (rhoe_old - rhoe_floor) / (-rhoe_src));
            update_damping(damping_factor_right, limiting_comp_right, trial, QREINT);
        }

        for (int n = 0; n < NumSpec; ++n) {
            Real X_old = qright(i,j,k,QFS+n);
            Real X_src = half_dt * sdc_src(i,j,k,QFS+n);

            if (X_src < 0.0_rt) {
                Real trial = amrex::max(0.0_rt, X_old / (-X_src));
                update_damping(damping_factor_right, limiting_comp_right, trial, QFS+n);
            } else if (X_src > 0.0_rt) {
                Real trial = amrex::max(0.0_rt, (1.0_rt - X_old) / X_src);
                update_damping(damping_factor_right, limiting_comp_right, trial, QFS+n);
            }
        }

        damping_factor_right = amrex::Clamp(damping_factor_right, 0.0_rt, 1.0_rt);

        qright(i,j,k,QPRES) = p_old + damping_factor_right * p_src;
        qright(i,j,k,QREINT) = rhoe_old + damping_factor_right * rhoe_src;


        for (int ipassive = 0; ipassive < npassive; ipassive++) {
            int n = qpassmap(ipassive);

            qleft(i,j,k,n) += damping_factor_left * half_dt * sdc_src(isrc,jsrc,ksrc,n);
            qright(i,j,k,n) += damping_factor_right * half_dt * sdc_src(i,j,k,n);

            if (n >= QFS && n < QFS + NumSpec) {
                qleft(i,j,k,n) = amrex::Clamp(qleft(i,j,k,n), 0.0_rt, 1.0_rt);
                qright(i,j,k,n) = amrex::Clamp(qright(i,j,k,n), 0.0_rt, 1.0_rt);
            }
        }

        Real X_sum_left = 0.0_rt;
        Real X_sum_right = 0.0_rt;

        for (int n = 0; n < NumSpec; ++n) {
            X_sum_left += qleft(i,j,k,QFS+n);
            X_sum_right += qright(i,j,k,QFS+n);
        }

        if (X_sum_left > 0.0_rt) {
            Real X_sum_inv = 1.0_rt / X_sum_left;
            for (int n = 0; n < NumSpec; ++n) {
                qleft(i,j,k,QFS+n) *= X_sum_inv;
            }
        }

        if (X_sum_right > 0.0_rt) {
            Real X_sum_inv = 1.0_rt / X_sum_right;
            for (int n = 0; n < NumSpec; ++n) {
                qright(i,j,k,QFS+n) *= X_sum_inv;
            }
        }

        // printing only for debug
#ifndef AMREX_USE_GPU
        if (damping_factor_left < 1.0_rt) {
            std::cout << "SDC damping applied on LEFT face: "
                      << "idir=" << idir
                      << " face=(" << i << "," << j << "," << k << ")"
                      << " src=(" << isrc << "," << jsrc << "," << ksrc << ")"
                      << " damping_factor=" << damping_factor_left
                      << " limiting_comp=" << limiting_comp_left
                      << " rho=" << qleft(i,j,k,QRHO)
                      << " p=" << qleft(i,j,k,QPRES)
                      << " rhoe=" << qleft(i,j,k,QREINT)
                      << std::endl;
        }

        if (damping_factor_right < 1.0_rt) {
            std::cout << "SDC damping applied on RIGHT face: "
                      << "idir=" << idir
                      << " face=(" << i << "," << j << "," << k << ")"
                      << " src=(" << i << "," << j << "," << k << ")"
                      << " damping_factor=" << damping_factor_right
                      << " limiting_comp=" << limiting_comp_right
                      << " rho=" << qright(i,j,k,QRHO)
                      << " p=" << qright(i,j,k,QPRES)
                      << " rhoe=" << qright(i,j,k,QREINT)
                      << std::endl;
        }
#elif defined(ALLOW_GPU_PRINTF)
        if (damping_factor_left < 1.0_rt) {
            AMREX_DEVICE_PRINTF(
                "SDC damping LEFT idir=%d face=(%d,%d,%d) src=(%d,%d,%d) damping_factor=%g limiting_comp=%d rho=%g p=%g rhoe=%g\n",
                idir, i, j, k, isrc, jsrc, ksrc, damping_factor_left, limiting_comp_left,
                qleft(i,j,k,QRHO), qleft(i,j,k,QPRES), qleft(i,j,k,QREINT));
        }

        if (damping_factor_right < 1.0_rt) {
            AMREX_DEVICE_PRINTF(
                "SDC damping RIGHT idir=%d face=(%d,%d,%d) src=(%d,%d,%d) damping_factor=%g limiting_comp=%d rho=%g p=%g rhoe=%g\n",
                idir, i, j, k, i, j, k, damping_factor_right, limiting_comp_right,
                qright(i,j,k,QRHO), qright(i,j,k,QPRES), qright(i,j,k,QREINT));
        }
#endif

    });

}

