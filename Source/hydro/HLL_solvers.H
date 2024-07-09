#ifndef HLL_SOLVERS_H
#define HLL_SOLVERS_H

#include <AMReX_REAL.H>

namespace HLL {

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void
    HLLC_state(const int idir, const Real S_k, const Real S_c,
               const Real* qstate, Real* U) {

      Real u_k = 0.0;
      if (idir == 0) {
        u_k = qstate[QU];
      } else if (idir == 1) {
        u_k = qstate[QV];
      } else if (idir == 2) {
        u_k = qstate[QW];
      }

      Real hllc_factor = qstate[QRHO]*(S_k - u_k)/(S_k - S_c);
      U[URHO] = hllc_factor;

      if (idir == 0) {
        U[UMX]  = hllc_factor*S_c;
        U[UMY]  = hllc_factor*qstate[QV];
        U[UMZ]  = hllc_factor*qstate[QW];

      } else if (idir == 1) {
        U[UMX]  = hllc_factor*qstate[QU];
        U[UMY]  = hllc_factor*S_c;
        U[UMZ]  = hllc_factor*qstate[QW];

      } else {
        U[UMX]  = hllc_factor*qstate[QU];
        U[UMY]  = hllc_factor*qstate[QV];
        U[UMZ]  = hllc_factor*S_c;
      }

      U[UEDEN] = hllc_factor*(qstate[QREINT]/qstate[QRHO] +
                              0.5_rt*(qstate[QU]*qstate[QU] + qstate[QV]*qstate[QV] + qstate[QW]*qstate[QW]) +
                              (S_c - u_k)*(S_c + qstate[QPRES]/(qstate[QRHO]*(S_k - u_k))));
      U[UEINT] = hllc_factor*qstate[QREINT]/qstate[QRHO];

      U[UTEMP] = 0.0; // we don't evolve T

    #ifdef SHOCK_VAR
      U[USHK] = 0.0;
    #endif

    #ifdef NSE_NET
      U[UMUP] = 0.0;
      U[UMUN] = 0.0;
    #endif

      for (int ipassive = 0; ipassive < npassive; ipassive++) {
          const int n  = upassmap(ipassive);
          const int nqs = qpassmap(ipassive);
          U[n] = hllc_factor*qstate[nqs];
      }
    }


    ///
    /// given a conserved state, compute the flux in direction idir
    ///
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void
    compute_flux(const int idir, const Real bnd_fac, const int coord,
                 const Real* U, const Real p,
                 Real* F) {

      Real u_flx = U[UMX+idir]/U[URHO];

      if (bnd_fac == 0) {
        u_flx = 0.0;
      }

      F[URHO] = U[URHO]*u_flx;

      F[UMX] = U[UMX]*u_flx;
      F[UMY] = U[UMY]*u_flx;
      F[UMZ] = U[UMZ]*u_flx;

      auto mom_check = mom_flux_has_p(idir, idir, coord);

      if (mom_check) {
        // we do not include the pressure term in any non-Cartesian
        // coordinate directions
        F[UMX+idir] = F[UMX+idir] + p;
      }

      F[UEINT] = U[UEINT]*u_flx;
      F[UEDEN] = (U[UEDEN] + p)*u_flx;

      F[UTEMP] = 0.0;

    #ifdef SHOCK_VAR
      F[USHK] = 0.0;
    #endif

    #ifdef NSE_NET
      F[UMUP] = 0.0;
      F[UMUN] = 0.0;
    #endif

      for (int ipassive=0; ipassive < npassive; ipassive++) {
          const int n = upassmap(ipassive);
          F[n] = U[n]*u_flx;
      }
    }

    ///
    /// convert the primitive variable state to the conserved state
    ///
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void
    cons_state(const Real* qstate, Real* U) {

      U[URHO] = qstate[QRHO];

      // since we advect all 3 velocity components regardless of dimension, this
      // will be general
      U[UMX]  = qstate[QRHO]*qstate[QU];
      U[UMY]  = qstate[QRHO]*qstate[QV];
      U[UMZ]  = qstate[QRHO]*qstate[QW];

      U[UEDEN] = qstate[QREINT] + 0.5_rt*qstate[QRHO]*(qstate[QU]*qstate[QU] + qstate[QV]*qstate[QV] + qstate[QW]*qstate[QW]);
      U[UEINT] = qstate[QREINT];

      // we don't care about T here, but initialize it to make NaN
      // checking happy
      U[UTEMP] = 0.0;

    #ifdef SHOCK_VAR
      U[USHK] = 0.0;
    #endif

    #ifdef NSE_NET
      U[UMUP] = 0.0;
      U[UMUN] = 0.0;
    #endif

      for (int ipassive = 0; ipassive < npassive; ipassive++) {
          const int n  = upassmap(ipassive);
          const int nqs = qpassmap(ipassive);
          U[n] = qstate[QRHO]*qstate[nqs];
      }
    }


    ///
    /// A simple HLL Riemann solver for pure hydrodynamics.  This takes just a
    /// single interface's data and returns the HLL flux
    ///
    /// @param ql     the left interface state
    /// @param qr     the right interface state
    /// @param cl     sound speed on the left interface
    /// @param cr     sound speed on the right interface
    /// @param idir   coordinate direction for the solve (0 = x, 1 = y, 2 = z)
    /// @param coord  geometry type (0 = Cartesian, 1 = axisymmetric, 2 = spherical)
    /// @param f      the HLL fluxes
    ///
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void
    HLL(const Real* ql, const Real* qr,
        const Real cl, const Real cr,
        const int idir, const int coord,
        Real* flux_hll) {

      // This is the HLLE solver.  We should apply it to zone averages
      // (not reconstructed states) at an interface in the presence of
      // shocks to avoid the odd-even decoupling / carbuncle phenomenon.
      //
      // See: Einfeldt, B.  et al. 1991, JCP, 92, 273
      //      Einfeldt, B. 1988, SIAM J NA, 25, 294


      constexpr Real small_hll = 1.e-10_rt;

      int ivel, ivelt, iveltt;
      int imom, imomt, imomtt;

      if (idir == 0) {
        ivel = QU;
        ivelt = QV;
        iveltt = QW;

        imom = UMX;
        imomt = UMY;
        imomtt = UMZ;

      } else if (idir == 1) {
        ivel = QV;
        ivelt = QU;
        iveltt = QW;

        imom = UMY;
        imomt = UMX;
        imomtt = UMZ;

      } else {
        ivel = QW;
        ivelt = QU;
        iveltt = QV;

        imom = UMZ;
        imomt = UMX;
        imomtt = UMY;
      }

      Real rhol_sqrt = std::sqrt(ql[QRHO]);
      Real rhor_sqrt = std::sqrt(qr[QRHO]);

      Real rhod = 1.0_rt/(rhol_sqrt + rhor_sqrt);


      // compute the average sound speed. This uses an approximation from
      // E88, eq. 5.6, 5.7 that assumes gamma falls between 1
      // and 5/3
      Real cavg = std::sqrt( (rhol_sqrt*cl*cl + rhor_sqrt*cr*cr)*rhod +
                             0.5_rt*rhol_sqrt*rhor_sqrt*rhod*rhod*std::pow(qr[ivel] - ql[ivel], 2));


      // Roe eigenvalues (E91, eq. 5.3b)
      Real uavg = (rhol_sqrt*ql[ivel] + rhor_sqrt*qr[ivel])*rhod;

      Real a1 = uavg - cavg;
      Real a4 = uavg + cavg;


      // signal speeds (E91, eq. 4.5)
      Real bl = std::min(a1, ql[ivel] - cl);
      Real br = std::max(a4, qr[ivel] + cr);

      Real bm = std::min(0.0_rt, bl);
      Real bp = std::max(0.0_rt, br);

      Real bd = bp - bm;

      if (std::abs(bd) < small_hll * std::max(std::abs(bm), std::abs(bp))) {
          return;
      }

      // we'll overwrite the passed in flux with the HLL flux

      bd = 1.0_rt/bd;

      // compute the fluxes according to E91, eq. 4.4b -- note that the
      // min/max above picks the correct flux if we are not in the star
      // region

      // density flux
      Real fl_tmp = ql[QRHO]*ql[ivel];
      Real fr_tmp = qr[QRHO]*qr[ivel];

      flux_hll[URHO] = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr[QRHO] - ql[QRHO]);

      // normal momentum flux.  Note for 1-d and 2-d non cartesian
      // r-coordinate, we leave off the pressure term and handle that
      // separately in the update, to accommodate different geometries
      fl_tmp = ql[QRHO]*ql[ivel]*ql[ivel];
      fr_tmp = qr[QRHO]*qr[ivel]*qr[ivel];
      if (mom_flux_has_p(idir, idir, coord)) {
        fl_tmp = fl_tmp + ql[QPRES];
        fr_tmp = fr_tmp + qr[QPRES];
      }

      flux_hll[imom] = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr[QRHO]*qr[ivel] - ql[QRHO]*ql[ivel]);

      // transverse momentum flux
      fl_tmp = ql[QRHO]*ql[ivel]*ql[ivelt];
      fr_tmp = qr[QRHO]*qr[ivel]*qr[ivelt];

      flux_hll[imomt] = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr[QRHO]*qr[ivelt] - ql[QRHO]*ql[ivelt]);


      fl_tmp = ql[QRHO]*ql[ivel]*ql[iveltt];
      fr_tmp = qr[QRHO]*qr[ivel]*qr[iveltt];

      flux_hll[imomtt] = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr[QRHO]*qr[iveltt] - ql[QRHO]*ql[iveltt]);

      // total energy flux
      Real rhoEl = ql[QREINT] + 0.5_rt*ql[QRHO]*(ql[ivel]*ql[ivel] + ql[ivelt]*ql[ivelt] + ql[iveltt]*ql[iveltt]);
      fl_tmp = ql[ivel]*(rhoEl + ql[QPRES]);

      Real rhoEr = qr[QREINT] + 0.5_rt*qr[QRHO]*(qr[ivel]*qr[ivel] + qr[ivelt]*qr[ivelt] + qr[iveltt]*qr[iveltt]);
      fr_tmp = qr[ivel]*(rhoEr + qr[QPRES]);

      flux_hll[UEDEN] = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(rhoEr - rhoEl);


      // eint flux
      fl_tmp = ql[QREINT]*ql[ivel];
      fr_tmp = qr[QREINT]*qr[ivel];

      flux_hll[UEINT] = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr[QREINT] - ql[QREINT]);


      // passively-advected scalar fluxes
      for (int ipassive = 0; ipassive < npassive; ipassive++) {
          const int n  = upassmap(ipassive);
          const int nqs = qpassmap(ipassive);

          fl_tmp = ql[QRHO]*ql[nqs]*ql[ivel];
          fr_tmp = qr[QRHO]*qr[nqs]*qr[ivel];

          flux_hll[n] = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr[QRHO]*qr[nqs] - ql[QRHO]*ql[nqs]);
      }
    }


    ///
    /// A HLLC Riemann solver for pure hydrodynamics
    ///
    /// @param ql          the left interface state
    /// @param qr          the right interface state
    /// @param qaux_arr    the auxiliary state
    /// @param uflx        the flux through the interface
    /// @param qint        an approximate Godunov state on the interface
    /// @param idir        coordinate direction for the solve (0 = x, 1 = y, 2 = z)
    ///
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void
    HLLC(const int i, const int j, const int k, const int idir,
         Array4<Real const> const& ql,
         Array4<Real const> const& qr,
         Array4<Real const> const& qaux_arr,
         Array4<Real> const& uflx,
         Array4<Real> const& qgdnv, const bool store_full_state,
         const GeometryData& geom,
         const bool special_bnd_lo, const bool special_bnd_hi,
         GpuArray<int, 3> const& domlo, GpuArray<int, 3> const& domhi) {

        // this is an implementation of the HLLC solver described in Toro's
        // book.  it uses the simplest estimate of the wave speeds, since
        // those should work for a general EOS.  We also initially do the
        // CGF Riemann construction to get pstar and ustar, since we'll
        // need to know the pressure and velocity on the interface for the
        // pdV term in the internal energy update.

        const Real small = 1.e-8_rt;

        int iu;
        int sx, sy, sz;

        if (idir == 0) {
            iu = QU;
            sx = 1;
            sy = 0;
            sz = 0;

        } else if (idir == 1) {
            iu = QV;
            sx = 0;
            sy = 1;
            sz = 0;

        } else {
            iu = QW;
            sx = 0;
            sy = 0;
            sz = 1;

        }


        int coord = geom.Coord();

        // deal with hard walls
        Real bnd_fac = 1.0_rt;

        if (idir == 0) {
            if ((i == domlo[0] && special_bnd_lo) ||
                (i == domhi[0]+1 && special_bnd_hi)) {
                bnd_fac = 0.0_rt;
            }

        } else if (idir == 1) {
            if ((j == domlo[1] && special_bnd_lo) ||
                (j == domhi[1]+1 && special_bnd_hi)) {
                bnd_fac = 0.0_rt;
            }

        } else {
            if ((k == domlo[2] && special_bnd_lo) ||
                (k == domhi[2]+1 && special_bnd_hi)) {
                bnd_fac = 0.0_rt;
            }
        }


        Real rl = std::max(ql(i,j,k,QRHO), small_dens);

        // pick left velocities based on direction
        Real ul  = ql(i,j,k,iu);

        Real pl = std::max(ql(i,j,k,QPRES), small_pres);

        Real rr = std::max(qr(i,j,k,QRHO), small_dens);

        // pick right velocities based on direction
        Real ur  = qr(i,j,k,iu);

        Real pr = std::max(qr(i,j,k,QPRES), small_pres);

        // now we essentially do the CGF solver to get p and u on the
        // interface, but we won't use these in any flux construction.
        Real csmall = amrex::max(small, small * qaux_arr(i,j,k,QC), small * qaux_arr(i-sx,j-sy,k-sz,QC));
        Real cavg = 0.5_rt*(qaux_arr(i,j,k,QC) + qaux_arr(i-sx,j-sy,k-sz,QC));

        Real gamcl = qaux_arr(i-sx,j-sy,k-sz,QGAMC);
        Real gamcr = qaux_arr(i,j,k,QGAMC);

    #ifdef TRUE_SDC
        if (use_reconstructed_gamma1 == 1) {
            gamcl = ql(i,j,k,QGC);
            gamcr = qr(i,j,k,QGC);
        }
    #endif

        Real wsmall = small_dens*csmall;
        Real wl = std::max(wsmall, std::sqrt(std::abs(gamcl*pl*rl)));
        Real wr = std::max(wsmall, std::sqrt(std::abs(gamcr*pr*rr)));

        Real wwinv = 1.0_rt/(wl + wr);
        Real pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))*wwinv;
        Real ustar = ((wl*ul + wr*ur) + (pl - pr))*wwinv;

        pstar = std::max(pstar, small_pres);
        // for symmetry preservation, if ustar is really small, then we
        // set it to zero
        if (std::abs(ustar) < riemann_constants::smallu*0.5_rt*(std::abs(ul) + std::abs(ur))){
            ustar = 0.0_rt;
        }

        Real ro;
        Real uo;
        Real po;
        Real gamco;

        if (ustar > 0.0_rt) {
            ro = rl;
            uo = ul;
            po = pl;
            gamco = gamcl;

        } else if (ustar < 0.0_rt) {
            ro = rr;
            uo = ur;
            po = pr;
            gamco = gamcr;

        } else {
            ro = 0.5_rt*(rl + rr);
            uo = 0.5_rt*(ul + ur);
            po = 0.5_rt*(pl + pr);
            gamco = 0.5_rt*(gamcl + gamcr);
        }

        ro = std::max(small_dens, ro);

        Real roinv = 1.0_rt/ro;
        Real co = std::sqrt(std::abs(gamco*po*roinv));
        co = std::max(csmall, co);
        Real co2inv = 1.0_rt/(co*co);

        Real rstar = ro + (pstar - po)*co2inv;
        rstar = std::max(small_dens, rstar);

        Real cstar = std::sqrt(std::abs(gamco*pstar/rstar));
        cstar = std::max(cstar, csmall);

        Real sgnm = std::copysign(1.0_rt, ustar);
        Real spout = co - sgnm*uo;
        Real spin = cstar - sgnm*ustar;
        Real ushock = 0.5_rt*(spin + spout);

        if (pstar-po > 0.0_rt) {
            spin = ushock;
            spout = ushock;
        }

        Real scr = spout-spin;
        if (spout-spin == 0.0_rt) {
            scr = small*cavg;
        }

        Real frac = (1.0_rt + (spout + spin)/scr)*0.5_rt;
        frac = std::clamp(frac, 0.0_rt, 1.0_rt);

        Real qint[NQ] = {0.0};

        qint[QRHO] = frac*rstar + (1.0_rt - frac)*ro;
        qint[iu] = frac*ustar + (1.0_rt - frac)*uo;
        qint[QPRES] = frac*pstar + (1.0_rt - frac)*po;


        // now we do the HLLC construction

        // use the simplest estimates of the wave speeds
        Real S_l = std::min(ul - std::sqrt(gamcl*pl/rl), ur - std::sqrt(gamcr*pr/rr));
        Real S_r = std::max(ul + std::sqrt(gamcl*pl/rl), ur + std::sqrt(gamcr*pr/rr));

        // estimate of the contact speed -- this is Toro Eq. 10.8
        Real S_c = (pr - pl + rl*ul*(S_l - ul) - rr*ur*(S_r - ur))/
            (rl*(S_l - ul) - rr*(S_r - ur));

        Real q_zone[NQ];
        Real U_state[NUM_STATE];
        Real U_hllc_state[NUM_STATE];
        Real F_state[NUM_STATE];

        if (S_r <= 0.0_rt) {
            // R region
            for (int n = 0; n < NQ; n++) {
                q_zone[n] = qr(i,j,k,n);
            }
            cons_state(q_zone, U_state);
            compute_flux(idir, bnd_fac, coord,
                         U_state, pr, F_state);

        } else if (S_r > 0.0_rt && S_c <= 0.0_rt) {
            // R* region
            for (int n = 0; n < NQ; n++) {
                q_zone[n] = qr(i,j,k,n);
            }
            cons_state(q_zone, U_state);
            compute_flux(idir, bnd_fac, coord,
                         U_state, pr, F_state);
            HLLC_state(idir, S_r, S_c, q_zone, U_hllc_state);

            // correct the flux
            for (int n = 0; n < NUM_STATE; n++) {
                F_state[n] = F_state[n] + S_r*(U_hllc_state[n] - U_state[n]);
            }

        } else if (S_c > 0.0_rt && S_l < 0.0_rt) {
            // L* region
            for (int n = 0; n < NQ; n++) {
                q_zone[n] = ql(i,j,k,n);
            }
            cons_state(q_zone, U_state);
            compute_flux(idir, bnd_fac, coord,
                         U_state, pl, F_state);
            HLLC_state(idir, S_l, S_c, q_zone, U_hllc_state);

            // correct the flux
            for (int n = 0; n < NUM_STATE; n++) {
                F_state[n] = F_state[n] + S_l*(U_hllc_state[n] - U_state[n]);
            }

        } else {
            // L region
            for (int n = 0; n < NQ; n++) {
                q_zone[n] = ql(i,j,k,n);
            }
            cons_state(q_zone, U_state);
            compute_flux(idir, bnd_fac, coord,
                         U_state, pl, F_state);
        }

        for (int n = 0; n < NUM_STATE; n++) {
            uflx(i,j,k,n) = F_state[n];
        }


        // now store the state

        if (store_full_state) {

            for (int n = 0; n < NQ; n++) {
                qgdnv(i,j,k,n) = qint[n];
            }

        } else {

    #ifdef HYBRID_MOMENTUM
            qgdnv(i,j,k,GDRHO) = qint[QRHO];
    #endif
            qgdnv(i,j,k,GDU) = qint[QU];
            qgdnv(i,j,k,GDV) = qint[QV];
            qgdnv(i,j,k,GDW) = qint[QW];
            qgdnv(i,j,k,GDPRES) = qint[QPRES];

        }


    }
}

#endif
