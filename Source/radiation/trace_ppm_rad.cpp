#include <Castro.H>
#include <Castro_util.H>

#ifdef RADIATION
#include <Radiation.H>
#endif

#include <cmath>

#include <ppm.H>
#include <flatten.H>

using namespace amrex;

void
Castro::trace_ppm_rad(const Box& bx,
                      const int idir,
                      Array4<Real const> const& U_arr,
                      Array4<Real const> const& rho_inv_arr,
                      Array4<Real const> const& q_arr,
                      Array4<Real const> const& qaux_arr,
                      Array4<Real const> const& srcQ,
                      Array4<Real> const& qm,
                      Array4<Real> const& qp,
#if (AMREX_SPACEDIM < 3)
                      Array4<Real const> const& dloga,
#endif
                      const Box& vbx,
                      const Real dt) {

  // These routines do the characteristic tracing under the parabolic
  // profiles in each zone to the edge / half-time.

  // To allow for easy integration of radiation, we adopt the
  // following conventions:
  //
  // rho : mass density
  // u, v, w : velocities
  // p : gas (hydro) pressure
  // ptot : total pressure (note for pure hydro, this is
  //        just the gas pressure)
  // rhoe_g : gas specific internal energy
  // cgas : sound speed for just the gas contribution
  // cc : total sound speed (including radiation)
  // h_g : gas specific enthalpy / cc**2
  // gam_g : the gas Gamma_1
  //
  // for pure hydro, we will only consider:
  //   rho, u, v, w, ptot, rhoe_g, cc, h_g

  const auto dx = geom.CellSizeArray();

  Real hdt = 0.5_rt * dt;
  Real dtdx = dt / dx[idir];

  auto lo = bx.loVect3d();
  auto hi = bx.hiVect3d();

  auto vlo = vbx.loVect3d();
  auto vhi = vbx.hiVect3d();

  // special care for reflecting BCs
  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();

  const auto domlo = geom.Domain().loVect3d();
  const auto domhi = geom.Domain().hiVect3d();

  bool lo_symm = lo_bc[idir] == amrex::PhysBCType::symmetry;
  bool hi_symm = hi_bc[idir] == amrex::PhysBCType::symmetry;

#ifndef AMREX_USE_GPU

  // if we're on the CPU, we preprocess the sources over the whole
  // tile up front -- we don't want to trace under a source that is
  // empty. This check only needs to be done over the tile we're
  // working on, since the PPM reconstruction and integration done
  // here is only local to this tile.

  GpuArray<int, NQSRC> do_source_trace;


  for (int n = 0; n < NQSRC; n++) {
    do_source_trace[n] = 0;

    for (int k = lo[2]-2*dg2; k <= hi[2]+2*dg2; k++) {
      for (int j = lo[1]-2*dg1; j <= hi[1]+2*dg1; j++) {
        for (int i = lo[0]-2; i <= hi[0]+2; i++) {
          if (std::abs(srcQ(i,j,k,n)) > 0.0_rt) {
            do_source_trace[n] = 1;
            break;
          }
        }
        if (do_source_trace[n] == 1) break;
      }
      if (do_source_trace[n] == 1) break;
    }
  }
#endif

  // This does the characteristic tracing to build the interface
  // states using the normal predictor only (no transverse terms).
  //
  // For each zone, we construct Im and Ip arrays -- these are the averages
  // of the various primitive state variables under the parabolic
  // interpolant over the region swept out by one of the 3 different
  // characteristic waves.
  //
  // Im is integrating to the left interface of the current zone
  // (which will be used to build the right ("p") state at that interface)
  // and Ip is integrating to the right interface of the current zone
  // (which will be used to build the left ("m") state at that interface).
  //
  //
  // The choice of reference state is designed to minimize the
  // effects of the characteristic projection.  We subtract the I's
  // off of the reference state, project the quantity such that it is
  // in terms of the characteristic variables, and then add all the
  // jumps that are moving toward the interface to the reference
  // state to get the full state on that interface.

  int QUN, QUT, QUTT;

  if (idir == 0) {
    QUN = QU;
    QUT = QV;
    QUTT = QW;
  } else if (idir == 1) {
    QUN = QV;
    QUT = QW;
    QUTT = QU;
  } else if (idir == 2) {
    QUN = QW;
    QUT = QU;
    QUTT = QV;
  }

  Real lsmall_dens = small_dens;
  Real lsmall_pres = small_pres;

  // Trace to left and right edges using upwind PPM
  AMREX_PARALLEL_FOR_3D(bx, i, j, k,
  {

   Real lam0[NGROUPS];
   Real lamp[NGROUPS];
   Real lamm[NGROUPS];

    for (int g = 0; g < NGROUPS; g++) {
      lam0[g] = qaux_arr(i,j,k,QLAMS+g);
      lamp[g] = qaux_arr(i,j,k,QLAMS+g);
      lamm[g] = qaux_arr(i,j,k,QLAMS+g);
    }

    Real rho = q_arr(i,j,k,QRHO);

    // cgassq is the gas soundspeed **2
    // cc is the total soundspeed **2 (gas + radiation)
    Real cgassq = qaux_arr(i,j,k,QCG)*qaux_arr(i,j,k,QCG);
    Real cc = qaux_arr(i,j,k,QC);
    Real csq = cc*cc;

    Real un = q_arr(i,j,k,QUN);

    Real p = q_arr(i,j,k,QPRES);
    Real rhoe_g = q_arr(i,j,k,QREINT);
    Real h_g = ( (p+rhoe_g)/rho)/csq;

    Real gam_g = qaux_arr(i,j,k,QGAMCG);

    Real ptot = q_arr(i,j,k,QPTOT);

    Real er[NGROUPS];
    Real hr[NGROUPS];
    for (int g = 0; g < NGROUPS; g++) {
      er[g] = q_arr(i,j,k,QRAD+g);
      hr[g] = (lam0[g] + 1.0_rt)*er[g]/rho;
    }

    // do the parabolic reconstruction and compute the
    // integrals under the characteristic waves
    Real s[5];

    Real flat = 1.0;

    if (castro::first_order_hydro) {
        flat = 0.0;
    }
    else if (castro::use_flattening) {
        flat = hydro::flatten(i, j, k, q_arr, QPRES);

#ifdef RADIATION
        flat *= hydro::flatten(i, j, k, q_arr, QPTOT);

        if (radiation::flatten_pp_threshold > 0.0) {
            if ( q_arr(i-1,j,k,QU) + q_arr(i,j-1,k,QV) + q_arr(i,j,k-1,QW) >
                 q_arr(i+1,j,k,QU) + q_arr(i,j+1,k,QV) + q_arr(i,j,k+1,QW) ) {

                if (q_arr(i,j,k,QPRES) < radiation::flatten_pp_threshold * q_arr(i,j,k,QPTOT)) {
                    flat = 0.0;
                }
            }
        }
#endif
    }

    Real sm;
    Real sp;

    Real Ip[NQ][3];
    Real Im[NQ][3];


    // do the non-passives first

    for (int n = 0; n < NQTHERM; n++) {
      if (n == QTEMP) continue;

      load_stencil(q_arr, idir, i, j, k, n, s);
      ppm_reconstruct(s, i, j, k, idir, reconstruction::Centering::zone_centered, flat, sm, sp);
      ppm_int_profile(sm, sp, s[i0], un, cc, dtdx, Ip[n], Im[n]);

    }

    // now the passives

    for (int ipassive = 0; ipassive < npassive; ++ipassive) {

      const int n = qpassmap(ipassive);
      const int nc = upassmap(ipassive);

      load_passive_stencil(U_arr, rho_inv_arr, idir, i, j, k, nc, s);
      ppm_reconstruct(s, i, j, k, idir, reconstruction::Centering::zone_centered, flat, sm, sp);
      ppm_int_profile(sm, sp, s[i0], un, cc, dtdx, Ip[n], Im[n]);

    }


    // source terms
    Real Ip_src[NQSRC][3];
    Real Im_src[NQSRC][3];

    for (int n = 0; n < NQSRC; n++) {

      // do we even need to trace (non-zero source?)
#ifndef AMREX_USE_GPU
      int do_trace = do_source_trace[n];
#else
      int do_trace = 0;
      if (idir == 0) {
        for (int b = i-2; b <= i+2; b++) {
          if (std::abs(srcQ(b,j,k,n)) > 0.0_rt) {
            do_trace = 1;
            break;
          }
        }
      } else if (idir == 1) {
        for (int b = j-2; b <= j+2; b++) {
          if (std::abs(srcQ(i,b,k,n)) > 0.0_rt) {
            do_trace = 1;
            break;
          }
        }
      } else {
        for (int b = k-2; b <= k+2; b++) {
          if (std::abs(srcQ(i,j,b,n)) > 0.0_rt) {
            do_trace = 1;
            break;
          }
        }
      }
#endif

      if (do_trace) {

          load_stencil(srcQ, idir, i, j, k, n, s);
          ppm_reconstruct(s, i, j, k, idir, lo_symm, hi_symm, domlo, domhi, flat, sm, sp);
          ppm_int_profile(sm, sp, s[i0], un, cc, dtdx, Ip_src[n], Im_src[n]);

      } else {
        Ip_src[n][0] = 0.0_rt;
        Ip_src[n][1] = 0.0_rt;
        Ip_src[n][2] = 0.0_rt;

        Im_src[n][0] = 0.0_rt;
        Im_src[n][1] = 0.0_rt;
        Im_src[n][2] = 0.0_rt;
      }

    }


    // do the passives separately

    // the passive stuff is the same regardless of the tracing

    for (int ipassive = 0; ipassive < npassive; ipassive++) {

      int n = qpassmap(ipassive);

      // Plus state on face i
      if ((idir == 0 && i >= vlo[0]) ||
          (idir == 1 && j >= vlo[1]) ||
          (idir == 2 && k >= vlo[2])) {

        // We have
        //
        // q_l = q_ref - Proj{(q_ref - I)}
        //
        // and Proj{} represents the characteristic projection.
        // But for these, there is only 1-wave that matters, the u
        // wave, so no projection is needed.  Since we are not
        // projecting, the reference state doesn't matter

        qp(i,j,k,n) = Im[n][1];
      }

      // Minus state on face i+1
      if (idir == 0 && i <= vhi[0]) {
        qm(i+1,j,k,n) = Ip[n][1];
      } else if (idir == 1 && j <= vhi[1]) {
        qm(i,j+1,k,n) = Ip[n][1];
      } else if (idir == 2 && k <= vhi[2]) {
        qm(i,j,k+1,n) = Ip[n][1];
      }
    }

    // plus state on face i

    if ((idir == 0 && i >= vlo[0]) ||
        (idir == 1 && j >= vlo[1]) ||
        (idir == 2 && k >= vlo[2])) {

      // Set the reference state
      // This will be the fastest moving state to the left --
      // this is the method that Miller & Colella and Colella &
      // Woodward use
      Real rho_ref = Im[QRHO][0];
      Real un_ref = Im[QUN][0];

      Real p_ref = Im[QPRES][0];
      Real rhoe_g_ref = Im[QREINT][0];

      Real ptot_ref = Im[QPTOT][0];
      Real er_ref[NGROUPS];
      for (int g = 0; g < NGROUPS; g++) {
        er_ref[g] = Im[QRAD+g][0];
      }

      rho_ref = amrex::max(rho_ref, lsmall_dens);
      p_ref = amrex::max(p_ref, lsmall_pres);

      // *m are the jumps carried by u-c
      // *p are the jumps carried by u+c

      // Note: for the transverse velocities, the jump is carried
      //       only by the u wave (the contact)

      // we also add the sources here so they participate in the tracing
      Real dum = un_ref - Im[QUN][0] - hdt*Im_src[QUN][0];
      Real dptotm = ptot_ref - Im[QPTOT][0] - hdt*Im_src[QPRES][0];

      Real drho = rho_ref - Im[QRHO][1] - hdt*Im_src[QRHO][1];
      Real dptot = ptot_ref - Im[QPTOT][1] - hdt*Im_src[QPRES][1];
      Real drhoe_g = rhoe_g_ref - Im[QREINT][1] - hdt*Im_src[QREINT][1];

      Real der[NGROUPS];
      for (int g = 0; g < NGROUPS; g++) {
        der[g] = er_ref[g] - Im[QRAD+g][1];
      }

      Real dup = un_ref - Im[QUN][2] - hdt*Im_src[QUN][2];
      Real dptotp = ptot_ref - Im[QPTOT][2] - hdt*Im_src[QPRES][2];

      // {rho, u, p, (rho e)} eigensystem

      // These are analogous to the beta's from the original PPM
      // paper (except we work with rho instead of tau).  This is
      // simply (l . dq), where dq = qref - I(q)

      Real alpham = 0.5_rt*(dptotm/(rho*cc) - dum)*rho/cc;
      Real alphap = 0.5_rt*(dptotp/(rho*cc) + dup)*rho/cc;
      Real alpha0r = drho - dptot/csq;
      Real alpha0e_g = drhoe_g - dptot*h_g;

      Real alphar[NGROUPS];
      for (int g = 0; g < NGROUPS; g++) {
        alphar[g] = der[g] - dptot/csq*hr[g];
      }

      alpham = un-cc > 0.0_rt ? 0.0_rt : -alpham;
      alphap = un+cc > 0.0_rt ? 0.0_rt : -alphap;
      alpha0r = un > 0.0_rt ? 0.0_rt : -alpha0r;
      alpha0e_g = un > 0.0_rt ? 0.0_rt : -alpha0e_g;

      for (int g = 0; g < NGROUPS; g++) {
        alphar[g] = un > 0.0_rt ? 0.0_rt : -alphar[g];
      }

      // The final interface states are just
      // q_s = q_ref - sum(l . dq) r
      // note that the a{mpz}right as defined above have the minus already
      qp(i,j,k,QRHO) = amrex::max(lsmall_dens, rho_ref + alphap + alpham + alpha0r);
      qp(i,j,k,QUN) = un_ref + (alphap - alpham)*cc/rho;
      qp(i,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g;

      qp(i,j,k,QPRES) = p_ref + (alphap + alpham)*cgassq;
      for (int g = 0; g < NGROUPS; g++) {
        qp(i,j,k,QPRES) += -lamp[g]*alphar[g];
      }
      qp(i,j,k,QPRES) = amrex::max(lsmall_pres, qp(i,j,k,QPRES));

      qp(i,j,k,QPTOT) = ptot_ref + (alphap + alpham)*csq;
      qp(i,j,k,QREITOT) = qp(i,j,k,QREINT);

      Real qrtmp;
      for (int g = 0; g < NGROUPS; g++) {
        qrtmp = er_ref[g] + (alphap + alpham)*hr[g] + alphar[g];
        qp(i,j,k,QRAD+g) = qrtmp;
        qp(i,j,k,QREITOT) += qrtmp;
      }

      for (int g = 0; g < NGROUPS; g++) {
        if (qp(i,j,k,QRAD+g) < 0.0_rt) {
          Real er_foo = -qp(i,j,k,QRAD+g);
          qp(i,j,k,QRAD+g) = 0.0_rt;
          qp(i,j,k,QPTOT) += lamp[g] * er_foo;
          qp(i,j,k,QREITOT) += er_foo;
        }
      }

      // Transverse velocities -- there's no projection here, so
      // we don't need a reference state.  We only care about
      // the state traced under the middle wave

      // Recall that I already takes the limit of the parabola
      // in the event that the wave is not moving toward the
      // interface
      qp(i,j,k,QUT) = Im[QUT][1] + hdt*Im_src[QUT][1];
      qp(i,j,k,QUTT) = Im[QUTT][1] + hdt*Im_src[QUTT][1];

    }

    // minus state on face i + 1

    if ((idir == 0 && i <= vhi[0]) ||
        (idir == 1 && j <= vhi[1]) ||
        (idir == 2 && k <= vhi[2])) {

      // Set the reference state
      // This will be the fastest moving state to the right
      Real rho_ref = Ip[QRHO][2];
      Real un_ref = Ip[QUN][2];

      Real p_ref = Ip[QPRES][2];
      Real rhoe_g_ref = Ip[QREINT][2];

      Real ptot_ref = Ip[QPTOT][2];
      Real er_ref[NGROUPS];
      for (int g = 0; g < NGROUPS; g++) {
        er_ref[g] = Ip[QRAD+g][2];
      }

      rho_ref = amrex::max(rho_ref, lsmall_dens);
      p_ref = amrex::max(p_ref, lsmall_pres);

      // *m are the jumps carried by u-c
      // *p are the jumps carried by u+c

      Real dum = un_ref - Ip[QUN][0] - hdt*Ip_src[QUN][0];
      Real dptotm = ptot_ref - Ip[QPTOT][0] - hdt*Ip_src[QPRES][0];

      Real drho = rho_ref - Ip[QRHO][1] - hdt*Ip_src[QRHO][1];
      Real dptot = ptot_ref - Ip[QPTOT][1] - hdt*Ip_src[QPRES][1];
      Real drhoe_g = rhoe_g_ref - Ip[QREINT][1] - hdt*Ip_src[QREINT][1];

      Real der[NGROUPS];
      for (int g = 0; g < NGROUPS; g++) {
        der[g]  = er_ref[g]  - Ip[QRAD+g][1];
      }

      Real dup = un_ref - Ip[QUN][2] - hdt*Ip_src[QUN][2];
      Real dptotp = ptot_ref - Ip[QPTOT][2] - hdt*Ip_src[QPRES][2];

      // {rho, u, p, (rho e)} eigensystem

      // These are analogous to the beta's from the original PPM
      // paper (except we work with rho instead of tau).  This is
      // simply (l . dq), where dq = qref - I(q)

      Real alpham = 0.5_rt*(dptotm/(rho*cc) - dum)*rho/cc;
      Real alphap = 0.5_rt*(dptotp/(rho*cc) + dup)*rho/cc;
      Real alpha0r = drho - dptot/csq;
      Real alpha0e_g = drhoe_g - dptot*h_g;

      Real alphar[NGROUPS];
      for (int g = 0; g < NGROUPS; g++) {
        alphar[g] = der[g] - dptot/csq*hr[g];
      }

      alpham = un-cc > 0.0_rt ? -alpham : 0.0_rt;
      alphap = un+cc > 0.0_rt ? -alphap : 0.0_rt;
      alpha0r = un > 0.0_rt ? -alpha0r : 0.0_rt;
      alpha0e_g = un > 0.0_rt ? -alpha0e_g : 0.0_rt;

      for (int g = 0; g < NGROUPS; g++) {
        alphar[g] = un > 0.0_rt ? -alphar[g] : 0.0_rt;
      }

      // The final interface states are just
      // q_s = q_ref - sum (l . dq) r
      // note that the a{mpz}left as defined above have the minus already

      if (idir == 0) {
        qm(i+1,j,k,QRHO) = amrex::max(lsmall_dens, rho_ref + alphap + alpham + alpha0r);
        qm(i+1,j,k,QUN) = un_ref + (alphap - alpham)*cc/rho;
        qm(i+1,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g;

        qm(i+1,j,k,QPRES) = p_ref + (alphap + alpham)*cgassq;
        for (int g = 0; g < NGROUPS; g++) {
          qm(i+1,j,k,QPRES) += -lamm[g]*alphar[g];
        }

        qm(i+1,j,k,QPTOT) = ptot_ref + (alphap + alpham)*csq;
        qm(i+1,j,k,QREITOT) = qm(i+1,j,k,QREINT);

        Real qrtmp;
        for (int g = 0; g < NGROUPS; g++) {
          qrtmp = er_ref[g] + (alphap + alpham)*hr[g] + alphar[g];
          qm(i+1,j,k,QRAD+g) = qrtmp;
          qm(i+1,j,k,QREITOT) += qrtmp;
        }

        for (int g = 0; g < NGROUPS; g++) {
          if (qm(i+1,j,k,QRAD+g) < 0.0_rt) {
            Real er_foo = -qm(i+1,j,k,QRAD+g);
            qm(i+1,j,k,QRAD+g) = 0.0_rt;
            qm(i+1,j,k,QPTOT) += lamm[g] * er_foo;
            qm(i+1,j,k,QREITOT) += er_foo;
          }
        }

        // transverse velocities
        qm(i+1,j,k,QUT) = Ip[QUT][1] + hdt*Ip_src[QUT][1];
        qm(i+1,j,k,QUTT) = Ip[QUTT][1] + hdt*Ip_src[QUTT][1];

      } else if (idir == 1) {
        qm(i,j+1,k,QRHO) = amrex::max(lsmall_dens, rho_ref + alphap + alpham + alpha0r);
        qm(i,j+1,k,QUN) = un_ref + (alphap - alpham)*cc/rho;
        qm(i,j+1,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g;

        qm(i,j+1,k,QPRES) = p_ref + (alphap + alpham)*cgassq;
        for (int g = 0; g < NGROUPS; g++) {
          qm(i,j+1,k,QPRES) += -lamm[g]*alphar[g];
        }

        qm(i,j+1,k,QPTOT) = ptot_ref + (alphap + alpham)*csq;
        qm(i,j+1,k,QREITOT) = qm(i,j+1,k,QREINT);

        Real qrtmp;
        for (int g = 0; g < NGROUPS; g++) {
          qrtmp = er_ref[g] + (alphap + alpham)*hr[g] + alphar[g];
          qm(i,j+1,k,QRAD+g) = qrtmp;
          qm(i,j+1,k,QREITOT) += qrtmp;
        }

        for (int g = 0; g < NGROUPS; g++) {
          qrtmp = er_ref[g] + (alphap + alpham)*hr[g] + alphar[g];
          qm(i,j+1,k,QRAD+g) = qrtmp;
          qm(i,j+1,k,QREITOT) += qrtmp;
        }

        for (int g = 0; g < NGROUPS; g++) {
          if (qm(i,j+1,k,QRAD+g) < 0.0_rt) {
            Real er_foo = -qm(i,j+1,k,QRAD+g);
            qm(i,j+1,k,QRAD+g) = 0.0_rt;
            qm(i,j+1,k,QPTOT) += lamm[g] * er_foo;
            qm(i,j+1,k,QREITOT) += er_foo;
          }
        }

        // transverse velocities
        qm(i,j+1,k,QUT) = Ip[QUT][1] + hdt*Ip_src[QUT][1];
        qm(i,j+1,k,QUTT) = Ip[QUTT][1] + hdt*Ip_src[QUTT][1];


      } else {

        qm(i,j,k+1,QRHO) = amrex::max(lsmall_dens, rho_ref + alphap + alpham + alpha0r);
        qm(i,j,k+1,QUN) = un_ref + (alphap - alpham)*cc/rho;
        qm(i,j,k+1,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g;

        qm(i,j,k+1,QPRES) = p_ref + (alphap + alpham)*cgassq;
        for (int g = 0; g < NGROUPS; g++) {
          qm(i,j,k+1,QPRES) += -lamm[g]*alphar[g];
        }

        qm(i,j,k+1,QPTOT) = ptot_ref + (alphap + alpham)*csq;
        qm(i,j,k+1,QREITOT) = qm(i,j,k+1,QREINT);

        Real qrtmp;
        for (int g = 0; g < NGROUPS; g++) {
          qrtmp = er_ref[g] + (alphap + alpham)*hr[g] + alphar[g];
          qm(i,j,k+1,QRAD+g) = qrtmp;
          qm(i,j,k+1,QREITOT) += qrtmp;
        }

        for (int g = 0; g < NGROUPS; g++) {
          qrtmp = er_ref[g] + (alphap + alpham)*hr[g] + alphar[g];
          qm(i,j,k+1,QRAD+g) = qrtmp;
          qm(i,j,k+1,QREITOT) += qrtmp;
        }

        for (int g = 0; g < NGROUPS; g++) {
          if (qm(i,j,k+1,QRAD+g) < 0.0_rt) {
            Real er_foo = -qm(i,j,k+1,QRAD+g);
            qm(i,j,k+1,QRAD+g) = 0.0_rt;
            qm(i,j,k+1,QPTOT) += lamm[g] * er_foo;
            qm(i,j,k+1,QREITOT) += er_foo;
          }
        }

        // transverse velocities
        qm(i,j,k+1,QUT) = Ip[QUT][1] + hdt*Ip_src[QUT][1];
        qm(i,j,k+1,QUTT) = Ip[QUTT][1] + hdt*Ip_src[QUTT][1];

      }

    }

    // geometry source terms

#if (AMREX_SPACEDIM < 3)
    // these only apply for x states (idir = 0)
    if (idir == 0 && dloga(i,j,k) != 0.0_rt) {
      Real courn = dt/dx[0]*(cc+std::abs(un));
      Real eta = (1.0_rt - courn)/(cc*dt*std::abs(dloga(i,j,k)));
      Real dlogatmp = amrex::min(eta, 1.0_rt)*dloga(i,j,k);
      Real sourcr = -0.5_rt*dt*rho*dlogatmp*un;
      Real sourcp = sourcr*cgassq;
      Real source = sourcp*h_g;
      Real sourcer[NGROUPS];
      for (int g = 0; g < NGROUPS; g++) {
        sourcer[g] = -0.5_rt*dt*dlogatmp*un*(lam0[g] + 1.0_rt)*er[g];
      }

      if (i <= vhi[0]) {
        qm(i+1,j,k,QRHO) = qm(i+1,j,k,QRHO) + sourcr;
        qm(i+1,j,k,QRHO) = amrex::max(qm(i+1,j,k,QRHO), lsmall_dens);
        qm(i+1,j,k,QPRES ) = qm(i+1,j,k,QPRES ) + sourcp;
        qm(i+1,j,k,QREINT) = qm(i+1,j,k,QREINT) + source;
        for (int g = 0; g < NGROUPS; g++) {
          qm(i+1,j,k,QRAD+g) += sourcer[g];
        }
        qm(i+1,j,k,QPTOT) = qm(i+1,j,k,QPTOT) + sourcp;
        qm(i+1,j,k,QREITOT) = qm(i+1,j,k,QREINT);
        for (int g = 0; g < NGROUPS; g++) {
          qm(i+1,j,k,QPTOT) += lamm[g]*sourcer[g];
          qm(i+1,j,k,QREITOT) += qm(i+1,j,k,QRAD+g);
          }
      }

      if (i >= vlo[0]) {
        qp(i,j,k,QRHO) = qp(i,j,k,QRHO) + sourcr;
        qp(i,j,k,QRHO) = amrex::max(qp(i,j,k,QRHO), lsmall_dens);
        qp(i,j,k,QPRES) = qp(i,j,k,QPRES) + sourcp;
        qp(i,j,k,QREINT) = qp(i,j,k,QREINT) + source;
        for (int g = 0; g < NGROUPS; g++) {
          qp(i,j,k,QRAD+g) += sourcer[g];
        }
        qp(i,j,k,QPTOT) = qp(i,j,k,QPTOT) + sourcp;
        qp(i,j,k,QREITOT) = qp(i,j,k,QREINT);
        for (int g = 0; g < NGROUPS; g++) {
          qp(i,j,k,QPTOT) += lamp[g]*sourcer[g];
          qp(i,j,k,QREITOT) += qp(i,j,k,QRAD+g);
        }
      }
    }
#endif

  });
}
