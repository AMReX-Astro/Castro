#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"
#include "Castro_util.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

#include <cmath>

#include <eos.H>
#include <riemann.H>

using namespace amrex;

void
Castro::HLLC(const Box& bx,
             Array4<Real const> const ql,
             Array4<Real const> const qr,
             Array4<Real const> const qaux_arr,
             Array4<Real> const uflx,
             Array4<Real> const qint,
             const int idir) {

  // this is an implementation of the HLLC solver described in Toro's
  // book.  it uses the simplest estimate of the wave speeds, since
  // those should work for a general EOS.  We also initially do the
  // CGF Riemann construction to get pstar and ustar, since we'll
  // need to know the pressure and velocity on the interface for the
  // pdV term in the internal energy update.

  const auto domlo = geom.Domain().loVect3d();
  const auto domhi = geom.Domain().hiVect3d();

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

  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();

  // do we want to force the flux to zero at the boundary?
  const bool special_bnd_lo = (lo_bc[idir] == Symmetry ||
                               lo_bc[idir] == SlipWall ||
                               lo_bc[idir] == NoSlipWall);
  const bool special_bnd_hi = (hi_bc[idir] == Symmetry ||
                               hi_bc[idir] == SlipWall ||
                               hi_bc[idir] == NoSlipWall);

  int coord = geom.Coord();

  AMREX_PARALLEL_FOR_3D(bx, i, j, k,
  {

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


    Real rl = amrex::max(ql(i,j,k,QRHO), castro::small_dens);

    // pick left velocities based on direction
    Real ul  = ql(i,j,k,iu);

    Real pl = amrex::max(ql(i,j,k,QPRES), castro::small_pres);
    Real rel = ql(i,j,k,QREINT);

    Real rr = amrex::max(qr(i,j,k,QRHO), castro::small_dens);

    // pick right velocities based on direction
    Real ur  = qr(i,j,k,iu);

    Real pr = amrex::max(qr(i,j,k,QPRES), castro::small_pres);
    Real rer = qr(i,j,k,QREINT);

    // now we essentially do the CGF solver to get p and u on the
    // interface, but we won't use these in any flux construction.
    Real csmall = amrex::max(SMALL, amrex::max(SMALL * qaux_arr(i,j,k,QC), SMALL * qaux_arr(i-sx,j-sy,k-sz,QC)));
    Real cavg = 0.5_rt*(qaux_arr(i,j,k,QC) + qaux_arr(i-sx,j-sy,k-sz,QC));

    Real gamcl = qaux_arr(i-sx,j-sy,k-sz,QGAMC);
    Real gamcr = qaux_arr(i,j,k,QGAMC);

#ifdef TRUE_SDC
    if (castro::use_reconstructed_gamma1 == 1) {
      gamcl = ql(i,j,k,QGC);
      gamcr = qr(i,j,k,QGC);
    }
#endif

    Real wsmall = castro::small_dens * csmall;
    Real wl = amrex::max(wsmall, std::sqrt(std::abs(gamcl*pl*rl)));
    Real wr = amrex::max(wsmall, std::sqrt(std::abs(gamcr*pr*rr)));

    Real wwinv = 1.0_rt/(wl + wr);
    Real pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))*wwinv;
    Real ustar = ((wl*ul + wr*ur) + (pl - pr))*wwinv;

    pstar = amrex::max(pstar, castro::small_pres);
    // for symmetry preservation, if ustar is really small, then we
    // set it to zero
    if (std::abs(ustar) < smallu*0.5_rt*(std::abs(ul) + std::abs(ur))){
      ustar = 0.0_rt;
    }

    Real ro;
    Real uo;
    Real po;
    Real reo;
    Real gamco;

    if (ustar > 0.0_rt) {
      ro = rl;
      uo = ul;
      po = pl;
      reo = rel;
      gamco = gamcl;

    } else if (ustar < 0.0_rt) {
      ro = rr;
      uo = ur;
      po = pr;
      reo = rer;
      gamco = gamcr;

    } else {
      ro = 0.5_rt*(rl + rr);
      uo = 0.5_rt*(ul + ur);
      po = 0.5_rt*(pl + pr);
      reo = 0.5_rt*(rel + rer);
      gamco = 0.5_rt*(gamcl + gamcr);
    }

    ro = amrex::max(castro::small_dens, ro);

    Real roinv = 1.0_rt/ro;
    Real co = std::sqrt(std::abs(gamco*po*roinv));
    co = amrex::max(csmall, co);
    Real co2inv = 1.0_rt/(co*co);

    Real rstar = ro + (pstar - po)*co2inv;
    rstar = amrex::max(castro::small_dens, rstar);

    Real cstar = std::sqrt(std::abs(gamco*pstar/rstar));
    cstar = max(cstar, csmall);

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
      scr = SMALL * cavg;
    }

    Real frac = (1.0_rt + (spout + spin)/scr)*0.5_rt;
    frac = amrex::max(0.0_rt, amrex::min(1.0_rt, frac));

    qint(i,j,k,iu) = frac*ustar + (1.0_rt - frac)*uo;
    qint(i,j,k,QPRES) = frac*pstar + (1.0_rt - frac)*po;


    // now we do the HLLC construction

    // use the simplest estimates of the wave speeds
    Real S_l = amrex::min(ul - std::sqrt(gamcl*pl/rl), ur - std::sqrt(gamcr*pr/rr));
    Real S_r = amrex::max(ul + std::sqrt(gamcl*pl/rl), ur + std::sqrt(gamcr*pr/rr));

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
  });
}


AMREX_GPU_HOST_DEVICE
void
Castro::HLL(const Real* ql, const Real* qr,
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
  Real bl = amrex::min(a1, ql[ivel] - cl);
  Real br = amrex::max(a4, qr[ivel] + cr);

  Real bm = amrex::min(0.0_rt, bl);
  Real bp = amrex::max(0.0_rt, br);

  Real bd = bp - bm;

  if (std::abs(bd) < small_hll*amrex::max(std::abs(bm), std::abs(bp))) return;

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
    int n  = upassmap(ipassive);
    int nqs = qpassmap(ipassive);

    fl_tmp = ql[QRHO]*ql[nqs]*ql[ivel];
    fr_tmp = qr[QRHO]*qr[nqs]*qr[ivel];

    flux_hll[n] = (bp*fl_tmp - bm*fr_tmp)*bd + bp*bm*bd*(qr[QRHO]*qr[nqs] - ql[QRHO]*ql[nqs]);
  }
}
