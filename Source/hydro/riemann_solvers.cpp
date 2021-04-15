#include <Castro.H>
#include <Castro_F.H>
#include <Castro_util.H>

#ifdef RADIATION
#include <Radiation.H>
#endif

#include <cmath>

#include <eos.H>
#include <riemann.H>

using namespace amrex;


void
Castro::HLLC(const Box& bx,
             Array4<Real const> const& ql,
             Array4<Real const> const& qr,
             Array4<Real const> const& qaux_arr,
             Array4<Real> const& uflx,
             Array4<Real> const& qint,
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

  const Real lsmall_dens = small_dens;
  const Real lsmall_pres = small_pres;
  const Real lsmall = riemann_constants::small;

  int coord = geom.Coord();

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
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


    Real rl = amrex::max(ql(i,j,k,QRHO), lsmall_dens);

    // pick left velocities based on direction
    Real ul  = ql(i,j,k,iu);

    Real pl = amrex::max(ql(i,j,k,QPRES), lsmall_pres);

    Real rr = amrex::max(qr(i,j,k,QRHO), lsmall_dens);

    // pick right velocities based on direction
    Real ur  = qr(i,j,k,iu);

    Real pr = amrex::max(qr(i,j,k,QPRES), lsmall_pres);

    // now we essentially do the CGF solver to get p and u on the
    // interface, but we won't use these in any flux construction.
    Real csmall = amrex::max(lsmall, amrex::max(lsmall * qaux_arr(i,j,k,QC), lsmall * qaux_arr(i-sx,j-sy,k-sz,QC)));
    Real cavg = 0.5_rt*(qaux_arr(i,j,k,QC) + qaux_arr(i-sx,j-sy,k-sz,QC));

    Real gamcl = qaux_arr(i-sx,j-sy,k-sz,QGAMC);
    Real gamcr = qaux_arr(i,j,k,QGAMC);

#ifdef TRUE_SDC
    if (use_reconstructed_gamma1 == 1) {
      gamcl = ql(i,j,k,QGC);
      gamcr = qr(i,j,k,QGC);
    }
#endif

    Real wsmall = lsmall_dens*csmall;
    Real wl = amrex::max(wsmall, std::sqrt(std::abs(gamcl*pl*rl)));
    Real wr = amrex::max(wsmall, std::sqrt(std::abs(gamcr*pr*rr)));

    Real wwinv = 1.0_rt/(wl + wr);
    Real pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))*wwinv;
    Real ustar = ((wl*ul + wr*ur) + (pl - pr))*wwinv;

    pstar = amrex::max(pstar, lsmall_pres);
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

    ro = amrex::max(lsmall_dens, ro);

    Real roinv = 1.0_rt/ro;
    Real co = std::sqrt(std::abs(gamco*po*roinv));
    co = amrex::max(csmall, co);
    Real co2inv = 1.0_rt/(co*co);

    Real rstar = ro + (pstar - po)*co2inv;
    rstar = amrex::max(lsmall_dens, rstar);

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
      scr = lsmall*cavg;
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
