#include <Castro.H>
#include <Castro_util.H>

#ifdef RADIATION
#include <Radiation.H>
#endif

#include <cmath>

#include <ppm.H>
#include <reconstruction.H>
#include <flatten.H>

using namespace amrex;
using namespace reconstruction;

void
Castro::trace_ppm(const Box& bx,
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

  // here, lo and hi are the range we loop over -- this can include ghost cells
  // vlo and vhi are the bounds of the valid box (no ghost cells)



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
  // game : gas gamma_e
  //
  // for pure hydro, we will only consider:
  //    rho, u, v, w, ptot, rhoe_g, cc, h_g


  const auto dx = geom.CellSizeArray();
  const auto problo = geom.ProbLoArray();
  const int coord = geom.Coord();

  Real hdt = 0.5_rt * dt;

#ifndef AMREX_USE_GPU
  auto lo = bx.loVect3d();
  auto hi = bx.hiVect3d();
#endif

  auto vlo = vbx.loVect3d();
  auto vhi = vbx.hiVect3d();


#ifndef AMREX_USE_GPU

  // if we're on the CPU, we preprocess the sources over the whole
  // tile up front -- we don't want to trace under a source that is
  // empty. This check only needs to be done over the tile we're
  // working on, since the PPM reconstruction and integration done
  // here is only local to this tile.

  GpuArray<int, NQSRC> do_source_trace;


  for (int n = 0; n < NQSRC; n++) {
    do_source_trace[n] = 0;

    // geometric source terms in r-direction for nonCartesian coordinate
    // or theta-direction in spherical coordinate need tracing

    if ((coord == 2 || (idir == 0 && coord == 1))
        && (n == QRHO || n == QPRES || n == QREINT)) {
        do_source_trace[n] = 1;
        continue;
    }

    for (int k = lo[2]-2*dg2; k <= hi[2]+2*dg2; k++) {
      for (int j = lo[1]-2*dg1; j <= hi[1]+2*dg1; j++) {
        for (int i = lo[0]-2; i <= hi[0]+2; i++) {
          if (std::abs(srcQ(i,j,k,n)) > 0.0_rt) {
            do_source_trace[n] = 1;
            break;
          }
        }
        if (do_source_trace[n] == 1) {
            break;
        }
      }
      if (do_source_trace[n] == 1) {
          break;
      }
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

  int QUN = -1;
  int QUT = -1;
  int QUTT = -1;

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
  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {
    Real dtdL = dt / dx[idir];
    Real dL = dx[idir];

    // Want dt/(rdtheta) instead of dt/dtheta for 2d Spherical
    if (coord == 2 && idir == 1) {
        Real r = problo[0] + static_cast<Real>(i + 0.5_rt) * dx[0];
        dL = r * dx[1];
        dtdL = dt / dL;
    }

    Real cc = qaux_arr(i,j,k,QC);
    Real un = q_arr(i,j,k,QUN);

    // do the parabolic reconstruction and compute the
    // integrals under the characteristic waves
    Real s[nslp];

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


    // reconstruct density

    Real Ip_rho[3];
    Real Im_rho[3];

    load_stencil(q_arr, idir, i, j, k, QRHO, s);
    ppm_reconstruct(s, flat, sm, sp);
    ppm_int_profile(sm, sp, s[i0], un, cc, dtdL, Ip_rho, Im_rho);

    // reconstruct normal velocity

    Real Ip_un_0;
    Real Im_un_0;
    Real Ip_un_2;
    Real Im_un_2;

    load_stencil(q_arr, idir, i, j, k, QUN, s);
    ppm_reconstruct(s, flat, sm, sp);
    ppm_int_profile_single(sm, sp, s[i0], un-cc, dtdL, Ip_un_0, Im_un_0);
    ppm_int_profile_single(sm, sp, s[i0], un+cc, dtdL, Ip_un_2, Im_un_2);

    // reconstruct pressure

    Real Ip_p[3];
    Real Im_p[3];

    load_stencil(q_arr, idir, i, j, k, QPRES, s);

    bool in_hse{};

    // HSE pressure on interfaces -- needed if we are dealing with
    // perturbation pressure as the parabolic reconstruction

    amrex::Real p_m_hse{};
    amrex::Real p_p_hse{};

    if (castro::use_pslope || castro::ppm_well_balanced) {
        Real trho[nslp];
        Real src[nslp];

        load_stencil(q_arr, idir, i, j, k, QRHO, trho);
        load_stencil(srcQ, idir, i, j, k, QUN, src);

        in_hse = ppm_reconstruct_pslope(trho, s, src, flat, dx[idir], sm, sp);

        if (in_hse && castro::ppm_well_balanced) {
            // we are working with the perturbational pressure
            ppm_int_profile(sm, sp, 0.0_rt, un, cc, dtdL, Ip_p, Im_p);
            p_m_hse = s[i0] - 0.5_rt * dx[idir] * trho[i0] * src[i0];
            p_p_hse = s[i0] + 0.5_rt * dx[idir] * trho[i0] * src[i0];

        } else {
            ppm_int_profile(sm, sp, s[i0], un, cc, dtdL, Ip_p, Im_p);
        }

    } else {
        ppm_reconstruct(s, flat, sm, sp);
        ppm_int_profile(sm, sp, s[i0], un, cc, dtdL, Ip_p, Im_p);
    }

    // reconstruct rho e

    Real Ip_rhoe[3];
    Real Im_rhoe[3];

    load_stencil(q_arr, idir, i, j, k, QREINT, s);
    ppm_reconstruct(s, flat, sm, sp);
    ppm_int_profile(sm, sp, s[i0], un, cc, dtdL, Ip_rhoe, Im_rhoe);

    // reconstruct transverse velocities

    Real Ip_ut_1;
    Real Im_ut_1;
    Real Ip_utt_1;
    Real Im_utt_1;

    load_stencil(q_arr, idir, i, j, k, QUT, s);
    ppm_reconstruct(s, flat, sm, sp);
    ppm_int_profile_single(sm, sp, s[i0], un, dtdL, Ip_ut_1, Im_ut_1);

    load_stencil(q_arr, idir, i, j, k, QUTT, s);
    ppm_reconstruct(s, flat, sm, sp);
    ppm_int_profile_single(sm, sp, s[i0], un, dtdL, Ip_utt_1, Im_utt_1);

    // gamma_c

    Real Ip_gc_0;
    Real Im_gc_0;
    Real Ip_gc_2;
    Real Im_gc_2;

    load_stencil(qaux_arr, idir, i, j, k, QGAMC, s);
    ppm_reconstruct(s, flat, sm, sp);
    ppm_int_profile_single(sm, sp, s[i0], un-cc, dtdL, Ip_gc_0, Im_gc_0);
    ppm_int_profile_single(sm, sp, s[i0], un+cc, dtdL, Ip_gc_2, Im_gc_2);


    // source terms -- we only trace if the terms in the stencil are non-zero
    int do_trace;

    // density

    Real Ip_src_rho[3] = {0.0_rt};
    Real Im_src_rho[3] = {0.0_rt};

#ifndef AMREX_USE_GPU
    do_trace = do_source_trace[QRHO];
#else
    if (coord == 2 || (idir == 0 && coord == 1)) {
        do_trace = 1;
    } else {
        do_trace = check_trace_source(srcQ, idir, i, j, k, QRHO);
    }
#endif

    if (do_trace) {
        load_stencil(srcQ, idir, i, j, k, QRHO, s);
#if AMREX_SPACEDIM <= 2
        if (coord == 2 || (idir == 0 && coord == 1)) {
            add_geometric_rho_source(q_arr, dloga, i, j, k, QUN, s);
        }
#endif
        ppm_reconstruct(s, flat, sm, sp);
        ppm_int_profile(sm, sp, s[i0], un, cc, dtdL, Ip_src_rho, Im_src_rho);
    }

    // normal velocity

    Real Ip_src_un_0 = 0.0_rt;
    Real Im_src_un_0 = 0.0_rt;
    Real Ip_src_un_2 = 0.0_rt;
    Real Im_src_un_2 = 0.0_rt;

#ifndef AMREX_USE_GPU
    do_trace = do_source_trace[QUN];
#else
    do_trace = check_trace_source(srcQ, idir, i, j, k, QUN);
#endif

    if (do_trace) {
        load_stencil(srcQ, idir, i, j, k, QUN, s);
        ppm_reconstruct(s, flat, sm, sp);
        ppm_int_profile_single(sm, sp, s[i0], un-cc, dtdL, Ip_src_un_0, Im_src_un_0);
        ppm_int_profile_single(sm, sp, s[i0], un+cc, dtdL, Ip_src_un_2, Im_src_un_2);
    }

    // pressure

    Real Ip_src_p[3] = {0.0_rt};
    Real Im_src_p[3] = {0.0_rt};

#ifndef AMREX_USE_GPU
    do_trace = do_source_trace[QPRES];
#else
    if (coord == 2 || (idir == 0 && coord == 1)) {
        do_trace = 1;
    } else {
        do_trace = check_trace_source(srcQ, idir, i, j, k, QPRES);
    }
#endif

    if (do_trace) {
        load_stencil(srcQ, idir, i, j, k, QPRES, s);
#if AMREX_SPACEDIM <= 2
        if (coord == 2 || (idir == 0 && coord == 1)) {
            add_geometric_p_source(q_arr, qaux_arr, dloga, i, j, k, QUN, s);
        }
#endif
        ppm_reconstruct(s, flat, sm, sp);
        ppm_int_profile(sm, sp, s[i0], un, cc, dtdL, Ip_src_p, Im_src_p);
    }

    // rho e

    Real Ip_src_rhoe[3] = {0.0_rt};
    Real Im_src_rhoe[3] = {0.0_rt};

#ifndef AMREX_USE_GPU
    do_trace = do_source_trace[QREINT];
#else
    if (coord == 2 || (idir == 0 && coord == 1)) {
        do_trace = 1;
    } else {
        do_trace = check_trace_source(srcQ, idir, i, j, k, QREINT);
    }
#endif

    if (do_trace) {
        load_stencil(srcQ, idir, i, j, k, QREINT, s);
#if AMREX_SPACEDIM <= 2
        if (coord == 2 || (idir == 0 && coord == 1)) {
            add_geometric_rhoe_source(q_arr, dloga, i, j, k, QUN, s);
        }
#endif
        ppm_reconstruct(s, flat, sm, sp);
        ppm_int_profile(sm, sp, s[i0], un, cc, dtdL, Ip_src_rhoe, Im_src_rhoe);
    }

    // transverse velocities

    Real Ip_src_ut_1 = 0.0_rt;
    Real Im_src_ut_1 = 0.0_rt;

#ifndef AMREX_USE_GPU
    do_trace = do_source_trace[QUT];
#else
    do_trace = check_trace_source(srcQ, idir, i, j, k, QUT);
#endif

    if (do_trace) {
        load_stencil(srcQ, idir, i, j, k, QUT, s);
        ppm_reconstruct(s, flat, sm, sp);
        ppm_int_profile_single(sm, sp, s[i0], un, dtdL, Ip_src_ut_1, Im_src_ut_1);
    }

    Real Ip_src_utt_1 = 0.0_rt;
    Real Im_src_utt_1 = 0.0_rt;

#ifndef AMREX_USE_GPU
    do_trace = do_source_trace[QUTT];
#else
    do_trace = check_trace_source(srcQ, idir, i, j, k, QUTT);
#endif

    if (do_trace) {
        load_stencil(srcQ, idir, i, j, k, QUTT, s);
        ppm_reconstruct(s, flat, sm, sp);
        ppm_int_profile_single(sm, sp, s[i0], un, dtdL, Ip_src_utt_1, Im_src_utt_1);
    }


    // do the passives separately

    // the passive stuff is the same regardless of the tracing

    Real Ip_passive;
    Real Im_passive;

    for (int ipassive = 0; ipassive < npassive; ipassive++) {

        const int nc = upassmap(ipassive);
        const int n = qpassmap(ipassive);

        load_passive_stencil(U_arr, rho_inv_arr, idir, i, j, k, nc, s);
        ppm_reconstruct(s, flat, sm, sp);
        ppm_int_profile_single(sm, sp, s[i0], un, dtdL, Ip_passive, Im_passive);

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

            qp(i,j,k,n) = Im_passive;
        }

        // Minus state on face i+1
        if (idir == 0 && i <= vhi[0]) {
            qm(i+1,j,k,n) = Ip_passive;
        } else if (idir == 1 && j <= vhi[1]) {
            qm(i,j+1,k,n) = Ip_passive;
        } else if (idir == 2 && k <= vhi[2]) {
            qm(i,j,k+1,n) = Ip_passive;
        }
    }


    // for well-balanced, the velocity sources should not be added

    amrex::Real fac = (castro::ppm_well_balanced && in_hse) ? 0.0_rt : 1.0_rt;

    // plus state on face i

    if ((idir == 0 && i >= vlo[0]) ||
        (idir == 1 && j >= vlo[1]) ||
        (idir == 2 && k >= vlo[2])) {

      // Set the reference state
      // This will be the fastest moving state to the left --
      // this is the method that Miller & Colella and Colella &
      // Woodward use. These papers don't include the effect of
      // the source term in constructing the reference state, but
      // we do because the source term could be large relative to
      // the current state (for example, consider the gravity term
      // acting on a fluid that is initially uniformly at rest --
      // the dt * g term will be the only non-zero contribution).
      // We ignore the effect of the source term for gamma.
      Real rho_ref = Im_rho[0] + hdt * Im_src_rho[0];
      Real un_ref = Im_un_0 + fac * hdt * Im_src_un_0;

      Real p_ref = Im_p[0] + hdt * Im_src_p[0];
      Real rhoe_g_ref = Im_rhoe[0] + hdt * Im_src_rhoe[0];

      Real gam_g_ref = Im_gc_0;

      rho_ref = std::max(rho_ref, lsmall_dens);

      Real rho_ref_inv = 1.0_rt/rho_ref;
      p_ref = std::max(p_ref, lsmall_pres);

      // For tracing
      Real csq_ref = gam_g_ref * (p_m_hse + p_ref) * rho_ref_inv;
      Real cc_ref = std::sqrt(csq_ref);
      Real cc_ref_inv = 1.0_rt/cc_ref;
      Real h_g_ref = ((p_m_hse + p_ref) + rhoe_g_ref) * rho_ref_inv;

      // *m are the jumps carried by un-c
      // *p are the jumps carried by un+c

      // Note: for the transverse velocities, the jump is carried
      //       only by the u wave (the contact)


      // we also add the sources here so they participate in the tracing

      Real dum = un_ref - Im_un_0 - fac * hdt * Im_src_un_0;
      Real dptotm = p_ref - Im_p[0] - hdt * Im_src_p[0];

      Real drho = rho_ref - Im_rho[1] - hdt * Im_src_rho[1];
      Real dptot = p_ref - Im_p[1] - hdt * Im_src_p[1];
      Real drhoe_g = rhoe_g_ref - Im_rhoe[1] - hdt * Im_src_rhoe[1];

      Real dup = un_ref - Im_un_2 - fac * hdt * Im_src_un_2;
      Real dptotp = p_ref - Im_p[2] - hdt * Im_src_p[2];

      // {rho, u, p, (rho e)} eigensystem

      // These are analogous to the beta's from the original PPM
      // paper (except we work with rho instead of tau).  This is
      // simply (l . dq), where dq = qref - I(q)

      Real alpham = 0.5_rt*(dptotm*rho_ref_inv*cc_ref_inv - dum)*rho_ref*cc_ref_inv;
      Real alphap = 0.5_rt*(dptotp*rho_ref_inv*cc_ref_inv + dup)*rho_ref*cc_ref_inv;
      Real alpha0r = drho - dptot/csq_ref;
      Real alpha0e_g = drhoe_g - dptot*h_g_ref/csq_ref;

      alpham = un-cc > 0.0_rt ? 0.0_rt : -alpham;
      alphap = un+cc > 0.0_rt ? 0.0_rt : -alphap;
      alpha0r = un > 0.0_rt ? 0.0_rt : -alpha0r;
      alpha0e_g = un > 0.0_rt ? 0.0_rt : -alpha0e_g;

      // The final interface states are just
      // q_s = q_ref - sum(l . dq) r
      // note that the a{mpz}right as defined above have the minus already
      qp(i,j,k,QRHO) = std::max(lsmall_dens, rho_ref +  alphap + alpham + alpha0r);
      qp(i,j,k,QUN) = un_ref + (alphap - alpham)*cc_ref*rho_ref_inv;
      qp(i,j,k,QREINT) = std::max(castro::small_dens * castro::small_ener,
                                    rhoe_g_ref + (alphap + alpham)*h_g_ref + alpha0e_g);
      qp(i,j,k,QPRES) = std::max(lsmall_pres, p_m_hse + p_ref + (alphap + alpham)*csq_ref);

      // Transverse velocities -- there's no projection here, so we
      // don't need a reference state.  We only care about the state
      // traced under the middle wave

      // Recall that I already takes the limit of the parabola
      // in the event that the wave is not moving toward the
      // interface
      qp(i,j,k,QUT) = Im_ut_1 + hdt * Im_src_ut_1;
      qp(i,j,k,QUTT) = Im_utt_1 + hdt * Im_src_utt_1;

    }

    // minus state on face i + 1

    if ((idir == 0 && i <= vhi[0]) ||
        (idir == 1 && j <= vhi[1]) ||
        (idir == 2 && k <= vhi[2])) {

      // Set the reference state
      // This will be the fastest moving state to the right
      Real rho_ref = Ip_rho[2] + hdt * Ip_src_rho[2];
      Real un_ref = Ip_un_2 + fac * hdt * Ip_src_un_2;

      Real p_ref = Ip_p[2] + hdt * Ip_src_p[2];
      Real rhoe_g_ref = Ip_rhoe[2] + hdt * Ip_src_rhoe[2];

      Real gam_g_ref = Ip_gc_2;

      rho_ref = std::max(rho_ref, lsmall_dens);
      Real rho_ref_inv = 1.0_rt/rho_ref;
      p_ref = std::max(p_ref, lsmall_pres);

      // For tracing
      Real csq_ref = gam_g_ref * (p_p_hse + p_ref) * rho_ref_inv;
      Real cc_ref = std::sqrt(csq_ref);
      Real cc_ref_inv = 1.0_rt/cc_ref;
      Real h_g_ref = ((p_p_hse + p_ref) + rhoe_g_ref) * rho_ref_inv;

      // *m are the jumps carried by u-c
      // *p are the jumps carried by u+c

      Real dum = un_ref - Ip_un_0 - fac * hdt * Ip_src_un_0;
      Real dptotm  = p_ref - Ip_p[0] - hdt * Ip_src_p[0];

      Real drho = rho_ref - Ip_rho[1] - hdt * Ip_src_rho[1];
      Real dptot = p_ref - Ip_p[1] - hdt * Ip_src_p[1];
      Real drhoe_g = rhoe_g_ref - Ip_rhoe[1] - hdt * Ip_src_rhoe[1];

      Real dup = un_ref - Ip_un_2 - fac * hdt * Ip_src_un_2;
      Real dptotp = p_ref - Ip_p[2] - hdt * Ip_src_p[2];

      // {rho, u, p, (rho e)} eigensystem

      // These are analogous to the beta's from the original PPM
      // paper (except we work with rho instead of tau).  This is
      // simply (l . dq), where dq = qref - I(q)

      Real alpham = 0.5_rt*(dptotm*rho_ref_inv*cc_ref_inv - dum)*rho_ref*cc_ref_inv;
      Real alphap = 0.5_rt*(dptotp*rho_ref_inv*cc_ref_inv + dup)*rho_ref*cc_ref_inv;
      Real alpha0r = drho - dptot/csq_ref;
      Real alpha0e_g = drhoe_g - dptot*h_g_ref/csq_ref;

      alpham = un-cc > 0.0_rt ? -alpham : 0.0_rt;
      alphap = un+cc > 0.0_rt ? -alphap : 0.0_rt;
      alpha0r = un > 0.0_rt ? -alpha0r : 0.0_rt;
      alpha0e_g = un > 0.0_rt ? -alpha0e_g : 0.0_rt;

      // The final interface states are just
      // q_s = q_ref - sum (l . dq) r
      // note that the a{mpz}left as defined above have the minus already
      if (idir == 0) {
        qm(i+1,j,k,QRHO) = std::max(lsmall_dens, rho_ref +  alphap + alpham + alpha0r);
        qm(i+1,j,k,QUN) = un_ref + (alphap - alpham)*cc_ref*rho_ref_inv;
        qm(i+1,j,k,QREINT) = std::max(castro::small_dens * castro::small_ener,
                                        rhoe_g_ref + (alphap + alpham)*h_g_ref + alpha0e_g);
        qm(i+1,j,k,QPRES) = std::max(lsmall_pres, p_p_hse + p_ref + (alphap + alpham)*csq_ref);

        // transverse velocities
        qm(i+1,j,k,QUT) = Ip_ut_1 + hdt*Ip_src_ut_1;
        qm(i+1,j,k,QUTT) = Ip_utt_1 + hdt*Ip_src_utt_1;

      } else if (idir == 1) {
        qm(i,j+1,k,QRHO) = std::max(lsmall_dens, rho_ref +  alphap + alpham + alpha0r);
        qm(i,j+1,k,QUN) = un_ref + (alphap - alpham)*cc_ref*rho_ref_inv;
        qm(i,j+1,k,QREINT) = std::max(castro::small_dens * castro::small_ener,
                                      rhoe_g_ref + (alphap + alpham)*h_g_ref + alpha0e_g);
        qm(i,j+1,k,QPRES) = std::max(lsmall_pres, p_p_hse + p_ref + (alphap + alpham)*csq_ref);

        // transverse velocities
        qm(i,j+1,k,QUT) = Ip_ut_1 + hdt*Ip_src_ut_1;
        qm(i,j+1,k,QUTT) = Ip_utt_1 + hdt*Ip_src_utt_1;

      } else if (idir == 2) {
        qm(i,j,k+1,QRHO) = std::max(lsmall_dens, rho_ref +  alphap + alpham + alpha0r);
        qm(i,j,k+1,QUN) = un_ref + (alphap - alpham)*cc_ref*rho_ref_inv;
        qm(i,j,k+1,QREINT) = std::max(castro::small_dens * castro::small_ener,
                                        rhoe_g_ref + (alphap + alpham)*h_g_ref + alpha0e_g);
        qm(i,j,k+1,QPRES) = std::max(lsmall_pres, p_p_hse + p_ref + (alphap + alpham)*csq_ref);

        // transverse velocities
        qm(i,j,k+1,QUT) = Ip_ut_1 + hdt*Ip_src_ut_1;
        qm(i,j,k+1,QUTT) = Ip_utt_1 + hdt*Ip_src_utt_1;
      }

    }

  });
}
