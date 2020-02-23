#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

#include <cmath>

#include <ppm.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE
void
Castro::trace_ppm(const Box& bx,
                  const int idir,
                  Array4<Real const> const q,
                  Array4<Real const> const qaux,
                  Array4<Real const> const srcQ,
                  Array4<Real const> const flatn,
                  Array4<Real> const qm,
                  Array4<Real> const qp,
#if (AMREX_SPACEDIM < 3)
                  Array4<Real const> const dloga,
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

  Real hdt = 0.5_rt * dt;
  Real dtdx = dt / dx[idir];

#ifndef AMREX_USE_CUDA

  // if we're on the CPU, we preprocess the sources over the whole
  // tile up front -- we don't want to trace under a source that is
  // empty. This check only needs to be done over the tile we're
  // working on, since the PPM reconstruction and integration done
  // here is only local to this tile.

  GpuArray<int, NQSRC> do_source_trace;

  const int* lo = bx.loVect();
  const int* hi = bx.hiVect();


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
  // in terms of the characteristic varaibles, and then add all the
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

  // Trace to left and right edges using upwind PPM
  AMREX_PARALLEL_FOR_3D(bx, i, j, k,
  {

    Real rho = q(i,j,k,QRHO);


    Real cc = qaux(i,j,k,QC);
    Real csq = cc*cc;

    Real un = q(i,j,k,QUN);


    // do the parabolic reconstruction and compute the
    // integrals under the characteristic waves
    Real s[5];
    Real flat = flatn(i,j,k);
    Real sm;
    Real sp;

    Real Ip[NQ][3];
    Real Im[NQ][3];


    for (int n = 0; n < NQ; n++) {
      if (n == UTEMP) continue;

      if (idir == 0) {
        s[im2] = q(i-2,j,k,n);
        s[im1] = q(i-1,j,k,n);
        s[i0]  = q(i,j,k,n);
        s[ip1] = q(i+1,j,k,n);
        s[ip2] = q(i+2,j,k,n);

      } else if (idir == 1) {
        s[im2] = q(i,j-2,k,n);
        s[im1] = q(i,j-1,k,n);
        s[i0]  = q(i,j,k,n);
        s[ip1] = q(i,j+1,k,n);
        s[ip2] = q(i,j+2,k,n);

      } else {
        s[im2] = q(i,j,k-2,n);
        s[im1] = q(i,j,k-1,n);
        s[i0]  = q(i,j,k,n);
        s[ip1] = q(i,j,k+1,n);
        s[ip2] = q(i,j,k+2,n);

      }

      ppm_reconstruct(s, flat, sm, sp);
      ppm_int_profile(sm, sp, s[i0], un, cc, dtdx, Ip[n], Im[n]);
    }

    // gamma_c

    Real Ip_gc[3];
    Real Im_gc[3];

    if (idir == 0) {
        s[im2] = qaux(i-2,j,k,QGAMC);
        s[im1] = qaux(i-1,j,k,QGAMC);
        s[i0]  = qaux(i,j,k,QGAMC);
        s[ip1] = qaux(i+1,j,k,QGAMC);
        s[ip2] = qaux(i+2,j,k,QGAMC);

    } else if (idir == 1) {
        s[im2] = qaux(i,j-2,k,QGAMC);
        s[im1] = qaux(i,j-1,k,QGAMC);
        s[i0]  = qaux(i,j,k,QGAMC);
        s[ip1] = qaux(i,j+1,k,QGAMC);
        s[ip2] = qaux(i,j+2,k,QGAMC);

    } else {
        s[im2] = qaux(i,j,k-2,QGAMC);
        s[im1] = qaux(i,j,k-1,QGAMC);
        s[i0]  = qaux(i,j,k,QGAMC);
        s[ip1] = qaux(i,j,k+1,QGAMC);
        s[ip2] = qaux(i,j,k+2,QGAMC);

    }


    ppm_reconstruct(s, flat, sm, sp);
    ppm_int_profile(sm, sp, s[i0], un, cc, dtdx, Ip_gc, Im_gc);


    // source terms
    Real Ip_src[NQSRC][3];
    Real Im_src[NQSRC][3];

    for (int n = 0; n < NQSRC; n++) {

      // do we even need to trace (non-zero source?)
#ifndef AMREX_USE_CUDA
      int do_trace = do_source_trace[n];
#else
      int do_trace = 0;
      if (idir == 0) {
        for (int q = lo[0]-2; q <= hi[0]+2) {
          if (std::abs(srcQ(q,j,k,n) > 0.0_rt)) {
            do_trace = 1;
            break;
          }
        }
      } else if (idir == 1) {
        for (int q = lo[1]-2; q <= hi[1]+2) {
          if (std::abs(srcQ(i,q,k,n) > 0.0_rt)) {
            do_trace = 1;
            break;
          }
        }
      } else {
        for (int q = lo[2]-2; q <= hi[2]+2) {
          if (std::abs(srcQ(i,j,q,n) > 0.0_rt)) {
            do_trace = 1;
            break;
          }
        }
      }
#endif

      if (do_trace) {

        if (idir == 0) {
          s[im2] = srcQ(i-2,j,k,n);
          s[im1] = srcQ(i-1,j,k,n);
          s[i0]  = srcQ(i,j,k,n);
          s[ip1] = srcQ(i+1,j,k,n);
          s[ip2] = srcQ(i+2,j,k,n);

        } else if (idir == 1) {
          s[im2] = srcQ(i,j-2,k,n);
          s[im1] = srcQ(i,j-1,k,n);
          s[i0]  = srcQ(i,j,k,n);
          s[ip1] = srcQ(i,j+1,k,n);
          s[ip2] = srcQ(i,j+2,k,n);

        } else {
          s[im2] = srcQ(i,j,k-2,n);
          s[im1] = srcQ(i,j,k-1,n);
          s[i0]  = srcQ(i,j,k,n);
          s[ip1] = srcQ(i,j,k+1,n);
          s[ip2] = srcQ(i,j,k+2,n);

        }

        ppm_reconstruct(s, flat, sm, sp);
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

    const int* vlo = vbx.loVect();
    const int* vhi = vbx.hiVect();

    // do the passives separately
    trace_ppm_species(i, j, k,
                      idir,
                      Ip, Im, Ip_src, Im_src,
                      qm, qp,
                      vlo, vhi,
                      dt);


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

      Real gam_g_ref = Im_gc[0];

      rho_ref = std::max(rho_ref, small_dens);

      Real rho_ref_inv = 1.0_rt/rho_ref;
      p_ref = std::max(p_ref, small_pres);

      // For tracing
      Real csq_ref = gam_g_ref*p_ref*rho_ref_inv;
      Real cc_ref = std::sqrt(csq_ref);
      Real cc_ref_inv = 1.0_rt/cc_ref;
      Real h_g_ref = (p_ref + rhoe_g_ref)*rho_ref_inv;

      // *m are the jumps carried by un-c
      // *p are the jumps carried by un+c

      // Note: for the transverse velocities, the jump is carried
      //       only by the u wave (the contact)


      // we also add the sources here so they participate in the tracing
      Real dum = un_ref - Im[QUN][0] - hdt*Im_src[QUN][0];
      Real dptotm = p_ref - Im[QPRES][0] - hdt*Im_src[QPRES][0];

      Real drho = rho_ref - Im[QRHO][1] - hdt*Im_src[QRHO][1];
      Real dptot = p_ref - Im[QPRES][1] - hdt*Im_src[QPRES][1];
      Real drhoe_g = rhoe_g_ref - Im[QREINT][1] - hdt*Im_src[QREINT][1];

      Real dup = un_ref - Im[QUN][2] - hdt*Im_src[QUN][2];
      Real dptotp = p_ref - Im[QPRES][2] - hdt*Im_src[QPRES][2];

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
      qp(i,j,k,QRHO) = std::max(small_dens, rho_ref +  alphap + alpham + alpha0r);
      qp(i,j,k,QUN) = un_ref + (alphap - alpham)*cc_ref*rho_ref_inv;
      qp(i,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ref + alpha0e_g;
      qp(i,j,k,QPRES) = std::max(small_pres, p_ref + (alphap + alpham)*csq_ref);

      // Transverse velocities -- there's no projection here, so we
      // don't need a reference state.  We only care about the state
      // traced under the middle wave

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

      Real gam_g_ref = Ip_gc[2];

      rho_ref = std::max(rho_ref, small_dens);
      Real rho_ref_inv = 1.0_rt/rho_ref;
      p_ref = std::max(p_ref, small_pres);

      // For tracing
      Real csq_ref = gam_g_ref*p_ref*rho_ref_inv;
      Real cc_ref = std::sqrt(csq_ref);
      Real cc_ref_inv = 1.0_rt/cc_ref;
      Real h_g_ref = (p_ref + rhoe_g_ref)*rho_ref_inv;

      // *m are the jumps carried by u-c
      // *p are the jumps carried by u+c

      Real dum = un_ref - Ip[QUN][0] - hdt*Ip_src[QUN][0];
      Real dptotm  = p_ref - Ip[QPRES][0] - hdt*Ip_src[QPRES][0];

      Real drho = rho_ref - Ip[QRHO][1] - hdt*Ip_src[QRHO][1];
      Real dptot = p_ref - Ip[QPRES][1] - hdt*Ip_src[QPRES][1];
      Real drhoe_g = rhoe_g_ref - Ip[QREINT][1] - hdt*Ip_src[QREINT][1];

      Real dup = un_ref - Ip[QUN][2] - hdt*Ip_src[QUN][2];
      Real dptotp = p_ref - Ip[QPRES][2] - hdt*Ip_src[QPRES][2];

      // (rho, u, p, (rho e)) eigensystem

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
        qm(i+1,j,k,QRHO) = std::max(small_dens, rho_ref +  alphap + alpham + alpha0r);
        qm(i+1,j,k,QUN) = un_ref + (alphap - alpham)*cc_ref*rho_ref_inv;
        qm(i+1,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ref + alpha0e_g;
        qm(i+1,j,k,QPRES) = max(small_pres, p_ref + (alphap + alpham)*csq_ref);

        // transverse velocities
        qm(i+1,j,k,QUT) = Ip[QUT][1] + hdt*Ip_src[QUT][1];
        qm(i+1,j,k,QUTT) = Ip[QUTT][1] + hdt*Ip_src[QUTT][1];

      } else if (idir == 1) {
        qm(i,j+1,k,QRHO) = std::max(small_dens, rho_ref +  alphap + alpham + alpha0r);
        qm(i,j+1,k,QUN) = un_ref + (alphap - alpham)*cc_ref*rho_ref_inv;
        qm(i,j+1,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ref + alpha0e_g;
        qm(i,j+1,k,QPRES) = max(small_pres, p_ref + (alphap + alpham)*csq_ref);

        // transverse velocities
        qm(i,j+1,k,QUT) = Ip[QUT][1] + hdt*Ip_src[QUT][1];
        qm(i,j+1,k,QUTT) = Ip[QUTT][1] + hdt*Ip_src[QUTT][1];

      } else if (idir == 2) {
        qm(i,j,k+1,QRHO) = std::max(small_dens, rho_ref +  alphap + alpham + alpha0r);
        qm(i,j,k+1,QUN) = un_ref + (alphap - alpham)*cc_ref*rho_ref_inv;
        qm(i,j,k+1,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g_ref + alpha0e_g;
        qm(i,j,k+1,QPRES) = max(small_pres, p_ref + (alphap + alpham)*csq_ref);

        // transverse velocities
        qm(i,j,k+1,QUT) = Ip[QUT][1] + hdt*Ip_src[QUT][1];
        qm(i,j,k+1,QUTT) = Ip[QUTT][1] + hdt*Ip_src[QUTT][1];
      }

      // geometry source terms

#if (AMREX_SPACEDIM < 3)
      // these only apply for x states (idir = 0)
      if (idir == 0 && dloga(i,j,k) != 0.0_rt) {
        Real courn = dt/dx[0]*(cc+std::abs(un));
        Real eta = (1.0_rt - courn)/(cc*dt*std::abs(dloga(i,j,k)));
        Real dlogatmp = std::min(eta, 1.0_rt)*dloga(i,j,k);
        Real sourcr = -0.5_rt*dt*rho*dlogatmp*un;
        Real sourcp = sourcr*csq;
        Real source = sourcp*((q(i,j,k,QPRES) + q(i,j,k,QREINT))/rho)/csq;

        if (i <= vhi[0]) {
          qm(i+1,j,k,QRHO) = qm(i+1,j,k,QRHO) + sourcr;
          qm(i+1,j,k,QRHO) = std::max(qm(i+1,j,k,QRHO), small_dens);
          qm(i+1,j,k,QPRES) = qm(i+1,j,k,QPRES) + sourcp;
          qm(i+1,j,k,QREINT) = qm(i+1,j,k,QREINT) + source;
        }

        if (i >= vlo[1]) {
          qp(i,j,k,QRHO) = qp(i,j,k,QRHO) + sourcr;
          qp(i,j,k,QRHO) = std::max(qp(i,j,k,QRHO), small_dens);
          qp(i,j,k,QPRES) = qp(i,j,k,QPRES) + sourcp;
          qp(i,j,k,QREINT) = qp(i,j,k,QREINT) + source;
        }
      }
#endif
    }

  });

}


AMREX_GPU_HOST_DEVICE
void
Castro::trace_ppm_species(const int i, const int j, const int k,
                          const int idir,
                          Real Ip[][3], Real Im[][3],
                          Real Ip_src[][3], Real Im_src[][3],
                          Array4<Real> const qm, Array4<Real> const qp,
                          const int* vlo, const int* vhi,
                          const Real dt) {

  // here, lo and hi are the range we loop over -- this can include ghost cells
  // vlo and vhi are the bounds of the valid box (no ghost cells)


  // the passive stuff is the same regardless of the tracing
  
  for (int ipassive = 0; ipassive < npassive; ipassive++) {

    int n = qpass_map[ipassive];

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
#ifdef PRIM_SPECIES_HAVE_SOURCES
      qp(i,j,k,n) += 0.5_rt * dt * Im_src[n][1];
#endif
    }

    // Minus state on face i+1
    if (idir == 0 && i <= vhi[0]) {
      qm(i+1,j,k,n) = Ip[n][1];
#ifdef PRIM_SPECIES_HAVE_SOURCES
      qm(i+1,j,k,n) += 0.5_rt * dt * Ip_src[n][1];
#endif

    } else if (idir == 1 && j <= vhi[1]) {
      qm(i,j+1,k,n) = Ip[n][1];
#ifdef PRIM_SPECIES_HAVE_SOURCES
      qm(i,j+1,k,n) += 0.5_rt * dt * Ip_src[n][1];
#endif

    } else if (idir == 2 && k <= vhi[2]) {
      qm(i,j,k+1,n) = Ip[n][1];
#ifdef PRIM_SPECIES_HAVE_SOURCES
      qm(i,j,k+1,n) += 0.5_rt * dt * Ip_src[n][1];
#endif
    }
  }
}

