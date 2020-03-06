#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

#include <cmath>

#include <ppm.H>

using namespace amrex;

void
Castro::trace_ppm_rad(const Box& bx,
                      const int idir,
                      Array4<Real const> const q_arr,
                      Array4<Real const> const qaux_arr,
                      Array4<Real const> const srcQ,
                      Array4<Real const> const flatn,
                      Array4<Real> const qm,
                      Array4<Real> const qp,
#if (AMREX_SPACEDIM < 3)
                      Array4<Real const> const dloga,
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

#ifndef AMREX_USE_CUDA

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

  Real lsmall_dens = small_dens;
  Real lsmall_pres = small_pres;

  GpuArray<int, npassive> qpass_map_p;
  for (int n = 0; n < npassive; n++){
    qpass_map_p[n] = qpass_map[n];
  }

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


    // do the parabolic reconstruction and compute the
    // integrals under the characteristic waves
    Real s[5];
    Real flat = flatn(i,j,k);
    Real sm;
    Real sp;

    Real Ip[NQ][3];
    Real Im[NQ][3];


    for (int n = 0; n < NQ; n++) {
      if (n == QTEMP) continue;

      if (idir == 0) {
        s[im2] = q_arr(i-2,j,k,n);
        s[im1] = q_arr(i-1,j,k,n);
        s[i0]  = q_arr(i,j,k,n);
        s[ip1] = q_arr(i+1,j,k,n);
        s[ip2] = q_arr(i+2,j,k,n);

      } else if (idir == 1) {
        s[im2] = q_arr(i,j-2,k,n);
        s[im1] = q_arr(i,j-1,k,n);
        s[i0]  = q_arr(i,j,k,n);
        s[ip1] = q_arr(i,j+1,k,n);
        s[ip2] = q_arr(i,j+2,k,n);

      } else {
        s[im2] = q_arr(i,j,k-2,n);
        s[im1] = q_arr(i,j,k-1,n);
        s[i0]  = q_arr(i,j,k,n);
        s[ip1] = q_arr(i,j,k+1,n);
        s[ip2] = q_arr(i,j,k+2,n);

      }

      ppm_reconstruct(s, flat, sm, sp);
      ppm_int_profile(sm, sp, s[i0], un, cc, dtdx, Ip[n], Im[n]);

    }


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


    // do the passives separately

    // the passive stuff is the same regardless of the tracing
  
    for (int ipassive = 0; ipassive < npassive; ipassive++) {

      int n = qpass_map_p[ipassive];

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
      for (int g=0; g < NGROUPS; g++) {
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
      Real dptotm = ptot_ref - Im[qptot][0] - hdt*Im_src[QPRES][0];

      Real drho = rho_ref - Im[QRHO][1] - hdt*Im_src[QRHO][1];
      Real dptot = ptot_ref - Im[qptot][1] - hdt*Im_src[QPRES][1];
      Real drhoe_g = rhoe_g_ref - Im[QREINT][1] - hdt*Im_src[QREINT][1];

      Real der[NGROUPS];
      for (int g-0; g < NGROUPS; g++) {
        der[g] = er_ref[g] - Im[QRAD+g][1];
      }

      Real dup = un_ref - Im[QUN][2] - hdt*Im_src[QUN][2];
      Real dptotp = ptot_ref - Im[QPTOT][2] - hdt*Im_src[QPRES][2];

      // {rho, u, p, (rho e)} eigensystem

      // These are analogous to the beta's from the original PPM
      // paper (except we work with rho instead of tau).  This is
      // simply (l . dq), where dq = qref - I(q)

                alpham = HALF*(dptotm/(rho*cc) - dum)*rho/cc
                alphap = HALF*(dptotp/(rho*cc) + dup)*rho/cc
                alpha0r = drho - dptot/csq
                alpha0e_g = drhoe_g - dptot*h_g

                alphar(:) = der(:) - dptot/csq*hr

                if (un-cc > ZERO) then
                   alpham = ZERO
                else
                   alpham = -alpham
                end if

                if (un+cc > ZERO) then
                   alphap = ZERO
                else
                   alphap = -alphap
                end if

                if (un > ZERO) then
                   alpha0r = ZERO
                else
                   alpha0r = -alpha0r
                end if

                if (un > ZERO) then
                   alphar(:) = ZERO
                else
                   alphar(:) = -alphar(:)
                end if

                if (un > ZERO) then
                   alpha0e_g = ZERO
                else
                   alpha0e_g = -alpha0e_g
                end if


                ! The final interface states are just
                ! q_s = q_ref - sum(l . dq) r
                ! note that the a{mpz}right as defined above have the minus already
                qp(i,j,k,QRHO) = rho_ref + alphap + alpham + alpha0r
                qp(i,j,k,QUN) = un_ref + (alphap - alpham)*cc/rho
                qp(i,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
                qp(i,j,k,QPRES) = p_ref + (alphap + alpham)*cgassq - sum(lamp(:)*alphar(:))

                qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
                qp(i,j,k,qrad:qrad-1+ngroups) = qrtmp

                qp(i,j,k,qptot) = ptot_ref + (alphap + alpham)*csq
                qp(i,j,k,qreitot) = qp(i,j,k,QREINT) + sum(qrtmp)


                ! Enforce small_*
                qp(i,j,k,QRHO) = max(qp(i,j,k,QRHO), small_dens)
                qp(i,j,k,QPRES) = max(qp(i,j,k,QPRES),small_pres)

                do g = 0, ngroups-1
                   if (qp(i,j,k,qrad+g) < ZERO) then
                      er_foo = - qp(i,j,k,qrad+g)
                      qp(i,j,k,qrad+g) = ZERO
                      qp(i,j,k,qptot) = qp(i,j,k,qptot) + lamp(g) * er_foo
                      qp(i,j,k,qreitot) = qp(i,j,k,qreitot) + er_foo
                   end if
                end do


                ! Transverse velocities -- there's no projection here, so
                ! we don't need a reference state.  We only care about
                ! the state traced under the middle wave

                ! Recall that I already takes the limit of the parabola
                ! in the event that the wave is not moving toward the
                ! interface
                qp(i,j,k,QUT) = Im(2,QUT) + hdt*Im_src(2,QUT)
                qp(i,j,k,QUTT) = Im(2,QUTT) + hdt*Im_src(2,QUTT)

             end if


             !-------------------------------------------------------------------
             ! minus state on face i + 1
             !-------------------------------------------------------------------
             if ((idir == 1 .and. i <= vhi(1)) .or. &
                 (idir == 2 .and. j <= vhi(2)) .or. &
                 (idir == 3 .and. k <= vhi(3))) then

                ! Set the reference state
                ! This will be the fastest moving state to the right
                rho_ref  = Ip(3,QRHO)
                un_ref    = Ip(3,QUN)

                p_ref    = Ip(3,QPRES)
                rhoe_g_ref = Ip(3,QREINT)

                tau_ref  = ONE/Ip(3,QRHO)

                !gam_g_ref  = Ip_gc(i,j,k,3,1)

                ptot_ref = Ip(3,QPTOT)

                er_ref(:) = Ip(3,QRAD:QRAD-1+ngroups)

                rho_ref = max(rho_ref,small_dens)
                p_ref = max(p_ref,small_pres)

                ! *m are the jumps carried by u-c
                ! *p are the jumps carried by u+c

                dum    = un_ref    - Ip(1,QUN) - hdt*Ip_src(1,QUN)
                dptotm = ptot_ref - Ip(1,qptot) - hdt*Ip_src(1,QPRES)

                drho    = rho_ref    - Ip(2,QRHO) - hdt*Ip_src(2,QRHO)
                dptot   = ptot_ref   - Ip(2,qptot) - hdt*Ip_src(2,QPRES)
                drhoe_g = rhoe_g_ref - Ip(2,QREINT) - hdt*Ip_src(2,QREINT)
                dtau  = tau_ref  - ONE/Ip(2,QRHO) + hdt*Ip_src(2,QRHO)/Ip(2,QRHO)**2
                der(:)  = er_ref(:)  - Ip(2,qrad:qrad-1+ngroups)

                dup    = un_ref    - Ip(3,QUN) - hdt*Ip_src(3,QUN)
                dptotp = ptot_ref - Ip(3,qptot) - hdt*Ip_src(3,QPRES)


                ! Optionally use the reference state in evaluating the
                ! eigenvectors -- NOT YET IMPLEMENTED

                ! (rho, u, p, (rho e)) eigensystem

                ! These are analogous to the beta's from the original PPM
                ! paper (except we work with rho instead of tau).  This is
                ! simply (l . dq), where dq = qref - I(q)

                alpham = HALF*(dptotm/(rho*cc) - dum)*rho/cc
                alphap = HALF*(dptotp/(rho*cc) + dup)*rho/cc
                alpha0r = drho - dptot/csq
                alpha0e_g = drhoe_g - dptot*h_g

                alphar(:) = der(:) - dptot/csq*hr

                if (un-cc > ZERO) then
                   alpham = -alpham
                else
                   alpham = ZERO
                end if

                if (un+cc > ZERO) then
                   alphap = -alphap
                else
                   alphap = ZERO
                end if

                if (un > ZERO) then
                   alpha0r = -alpha0r
                else
                   alpha0r = ZERO
                end if

                if (un > ZERO) then
                   alphar(:) = -alphar(:)
                else
                   alphar(:) = ZERO
                end if

                if (un > ZERO) then
                   alpha0e_g = -alpha0e_g
                else
                   alpha0e_g = ZERO
                end if

                ! The final interface states are just
                ! q_s = q_ref - sum (l . dq) r
                ! note that the a{mpz}left as defined above have the minus already

                if (idir == 1) then
                   qm(i+1,j,k,QRHO) = max(small_dens, rho_ref + alphap + alpham + alpha0r)
                   qm(i+1,j,k,QUN) = un_ref + (alphap - alpham)*cc/rho
                   qm(i+1,j,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
                   qm(i+1,j,k,QPRES) = max(small_pres, p_ref + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:)))

                   qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
                   qm(i+1,j,k,qrad:qrad-1+ngroups) = qrtmp

                   qm(i+1,j,k,qptot) = ptot_ref + (alphap + alpham)*csq
                   qm(i+1,j,k,qreitot) = qm(i+1,j,k,QREINT) + sum(qrtmp)

                else if (idir == 2) then
                   qm(i,j+1,k,QRHO) = max(small_dens, rho_ref + alphap + alpham + alpha0r)
                   qm(i,j+1,k,QUN) = un_ref + (alphap - alpham)*cc/rho
                   qm(i,j+1,k,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
                   qm(i,j+1,k,QPRES) = max(small_pres, p_ref + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:)))

                   qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
                   qm(i,j+1,k,qrad:qrad-1+ngroups) = qrtmp

                   qm(i,j+1,k,qptot) = ptot_ref + (alphap + alpham)*csq
                   qm(i,j+1,k,qreitot) = qm(i,j+1,k,QREINT) + sum(qrtmp)

                else if (idir == 3) then
                   qm(i,j,k+1,QRHO) = max(small_dens, rho_ref + alphap + alpham + alpha0r)
                   qm(i,j,k+1,QUN) = un_ref + (alphap - alpham)*cc/rho
                   qm(i,j,k+1,QREINT) = rhoe_g_ref + (alphap + alpham)*h_g*csq + alpha0e_g
                   qm(i,j,k+1,QPRES) = max(small_pres, p_ref + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:)))

                   qrtmp = er_ref(:) + (alphap + alpham)*hr + alphar(:)
                   qm(i,j,k+1,qrad:qrad-1+ngroups) = qrtmp

                   qm(i,j,k+1,qptot) = ptot_ref + (alphap + alpham)*csq
                   qm(i,j,k+1,qreitot) = qm(i,j,k+1,QREINT) + sum(qrtmp)

                end if

                if (idir == 1) then
                   do g=0,ngroups-1
                      if (qm(i+1,j,k,qrad+g) < ZERO) then
                         er_foo = - qm(i+1,j,k,qrad+g)
                         qm(i+1,j,k,qrad+g) = ZERO
                         qm(i+1,j,k,qptot) = qm(i+1,j,k,qptot) + lamm(g) * er_foo
                         qm(i+1,j,k,qreitot) = qm(i+1,j,k,qreitot) + er_foo
                      end if
                   end do

                   ! transverse velocities
                   qm(i+1,j,k,QUT) = Ip(2,QUT) + hdt*Ip_src(2,QUT)
                   qm(i+1,j,k,QUTT) = Ip(2,QUTT) + hdt*Ip_src(2,QUTT)

                else if (idir == 2) then
                   do g=0,ngroups-1
                      if (qm(i,j+1,k,qrad+g) < ZERO) then
                         er_foo = - qm(i,j+1,k,qrad+g)
                         qm(i,j+1,k,qrad+g) = ZERO
                         qm(i,j+1,k,qptot) = qm(i,j+1,k,qptot) + lamm(g) * er_foo
                         qm(i,j+1,k,qreitot) = qm(i,j+1,k,qreitot) + er_foo
                      end if
                   end do

                   ! transverse velocities
                   qm(i,j+1,k,QUT) = Ip(2,QUT) + hdt*Ip_src(2,QUT)
                   qm(i,j+1,k,QUTT) = Ip(2,QUTT) + hdt*Ip_src(2,QUTT)

                else if (idir == 3) then
                   do g=0,ngroups-1
                      if (qm(i,j,k+1,qrad+g) < ZERO) then
                         er_foo = - qm(i,j,k+1,qrad+g)
                         qm(i,j,k+1,qrad+g) = ZERO
                         qm(i,j,k+1,qptot) = qm(i,j,k+1,qptot) + lamm(g) * er_foo
                         qm(i,j,k+1,qreitot) = qm(i,j,k+1,qreitot) + er_foo
                      end if
                   end do

                   ! transverse velocities
                   qm(i,j,k+1,QUT) = Ip(2,QUT) + hdt*Ip_src(2,QUT)
                   qm(i,j,k+1,QUTT) = Ip(2,QUTT) + hdt*Ip_src(2,QUTT)

                end if

             end if


             !-------------------------------------------------------------------
             ! geometry source terms
             !-------------------------------------------------------------------

#if (AMREX_SPACEDIM < 3)
             if (idir == 1 .and. dloga(i,j,k) /= 0) then
                courn = dt/dx(1)*(cc+abs(un))
                eta = (ONE-courn)/(cc*dt*abs(dloga(i,j,k)))
                dlogatmp = min(eta, ONE)*dloga(i,j,k)
                sourcr = -HALF*dt*rho*dlogatmp*un
                sourcp = sourcr*cgassq
                source = sourcp*h_g
                sourcer(:) = -HALF*dt*dlogatmp*un*(lam0(:)+ONE)*er(:)

                if (i <= vhi(1)) then
                   qm(i+1,j,k,QRHO  ) = qm(i+1,j,k,QRHO  ) + sourcr
                   qm(i+1,j,k,QRHO  ) = max(qm(i+1,j,k,QRHO), small_dens)
                   qm(i+1,j,k,QPRES ) = qm(i+1,j,k,QPRES ) + sourcp
                   qm(i+1,j,k,QREINT) = qm(i+1,j,k,QREINT) + source
                   qm(i+1,j,k,qrad:qrad-1+ngroups) = qm(i+1,j,k,qrad:qrad-1+ngroups) + sourcer(:)
                   ! qm(i+1,j,k,qptot ) = sum(lamm(:)*qm(i+1,j,k,qrad:qradhi)) + qm(i+1,j,k,QPRES)
                   qm(i+1,j,k,qptot) = qm(i+1,j,k,qptot) + sum(lamm(:)*sourcer(:)) + sourcp
                   qm(i+1,j,k,qreitot) = sum(qm(i+1,j,k,qrad:qrad-1+ngroups))  + qm(i+1,j,k,QREINT)
                end if

                if (i >= vlo(1)) then
                   qp(i,j,k,QRHO  ) = qp(i,j,k,QRHO  ) + sourcr
                   qp(i,j,k,QRHO  ) = max(qp(i,j,k,QRHO), small_dens)
                   qp(i,j,k,QPRES ) = qp(i,j,k,QPRES ) + sourcp
                   qp(i,j,k,QREINT) = qp(i,j,k,QREINT) + source
                   qp(i,j,k,qrad:qrad-1+ngroups) = qp(i,j,k,qrad:qrad-1+ngroups) + sourcer(:)
                   ! qp(i  ,qptot ) = sum(lamp(:)*qp(i,qrad:qradhi)) + qp(i,QPRES)
                   qp(i,j,k,qptot) = qp(i,j,k,qptot) + sum(lamp(:)*sourcer(:)) + sourcp
                   qp(i,j,k,qreitot) = sum(qp(i,j,k,qrad:qrad-1+ngroups))  + qp(i,j,k,QREINT)
                end if
             endif
#endif

          end do
       end do
    end do

  end subroutine trace_ppm_rad

  subroutine trace_ppm_species(i, j, k, &
                               idir, &
                               q, qd_lo, qd_hi, &
                               Ip, Im, &
                               Ip_src, Im_src, &
                               qm, qm_lo, qm_hi, &
                               qp, qp_lo, qp_hi, &
                               vlo, vhi, domlo, domhi, &
                               dx, dt)
    ! here, lo and hi are the range we loop over -- this can include ghost cells
    ! vlo and vhi are the bounds of the valid box (no ghost cells)

    use network, only : nspec, naux
    use meth_params_module, only : NQ, NQSRC, npassive, qpass_map

    implicit none

    integer, intent(in) :: idir
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qp_lo(3), qp_hi(3)
    integer, intent(in) :: qm_lo(3), qm_hi(3)

    integer, intent(in) :: vlo(3), vhi(3)
    integer, intent(in) :: domlo(3), domhi(3)

    real(rt), intent(in) :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),NQ)

    real(rt), intent(in) :: Ip(1:3,NQ)
    real(rt), intent(in) :: Im(1:3,NQ)

    real(rt), intent(in) :: Ip_src(1:3,NQSRC)
    real(rt), intent(in) :: Im_src(1:3,NQSRC)

    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),NQ)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQ)
    real(rt), intent(in) :: dt, dx(3)

    integer, intent(in) :: i, j, k

    integer :: ipassive, n

    !$gpu

    ! the passive stuff is the same regardless of the tracing
    do ipassive = 1, npassive
       n = qpass_map(ipassive)

       ! Plus state on face i
       if ((idir == 1 .and. i >= vlo(1)) .or. &
           (idir == 2 .and. j >= vlo(2)) .or. &
           (idir == 3 .and. k >= vlo(3))) then

          ! We have
          !
          ! q_l = q_ref - Proj{(q_ref - I)}
          !
          ! and Proj{} represents the characteristic projection.
          ! But for these, there is only 1-wave that matters, the u
          ! wave, so no projection is needed.  Since we are not
          ! projecting, the reference state doesn't matter

          qp(i,j,k,n) = Im(2,n)
          if (n <= NQSRC) qp(i,j,k,n) = qp(i,j,k,n) + HALF*dt*Im_src(2,n)

       end if

       ! Minus state on face i+1
       if (idir == 1 .and. i <= vhi(1)) then
          qm(i+1,j,k,n) = Ip(2,n)
          if (n <= NQSRC) qm(i+1,j,k,n) = qm(i+1,j,k,n) + HALF*dt*Ip_src(2,n)

       else if (idir == 2 .and. j <= vhi(2)) then
          qm(i,j+1,k,n) = Ip(2,n)
          if (n <= NQSRC) qm(i,j+1,k,n) = qm(i,j+1,k,n) + HALF*dt*Ip_src(2,n)

       else if (idir == 3 .and. k <= vhi(3)) then
          qm(i,j,k+1,n) = Ip(2,n)
          if (n <= NQSRC) qm(i,j,k+1,n) = qm(i,j,k+1,n) + HALF*dt*Ip_src(2,n)
       end if

    end do
  end subroutine trace_ppm_species


end module trace_ppm_rad_module
