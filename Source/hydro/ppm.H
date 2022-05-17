
#include <cmath>
#include <reconstruction.H>

#ifndef _ppm_H_
#define _ppm_H_

using namespace amrex;


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int
check_trace_source(Array4<Real const> const& srcQ, const int idir,
                   const int i, const int j, const int k, const int ncomp) {

    int do_trace = 0;
    if (idir == 0) {
        for (int b = i-2; b <= i+2; b++) {
            if (std::abs(srcQ(b,j,k,ncomp)) > 0.0_rt) {
                do_trace = 1;
                break;
            }
        }
    } else if (idir == 1) {
        for (int b = j-2; b <= j+2; b++) {
            if (std::abs(srcQ(i,b,k,ncomp)) > 0.0_rt) {
                do_trace = 1;
                break;
            }
        }
    } else {
        for (int b = k-2; b <= k+2; b++) {
            if (std::abs(srcQ(i,j,b,ncomp)) > 0.0_rt) {
                do_trace = 1;
                break;
            }
        }
    }

    return do_trace;
}



///
/// Compute the coefficients of a parabolic reconstruction of the data in a zone.
/// This uses the standard PPM limiters described in Colella & Woodward (1984)
///
/// @param s      Real[5] the state to be reconstructed in zones i-2, i-1, i, i+1, i+2
/// @param flatn  The flattening coefficient
/// @param sm     The value of the parabola on the left edge of the zone
/// @param sp     The value of the parabola on the right edge of the zone
///
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
ppm_reconstruct(const Real* s,
                const Real flatn, Real& sm, Real& sp) {


  // Compute van Leer slopes

  Real dsl = 2.0_rt * (s[im1] - s[im2]);
  Real dsr = 2.0_rt * (s[i0] - s[im1]);

  Real dsvl_l = 0.0_rt;
  if (dsl*dsr > 0.0_rt) {
    Real dsc = 0.5_rt * (s[i0] - s[im2]);
    dsvl_l = std::copysign(1.0_rt, dsc) * amrex::min(std::abs(dsc),amrex::min(std::abs(dsl),std::abs(dsr)));
  }

  dsl = 2.0_rt * (s[i0] - s[im1]);
  dsr = 2.0_rt * (s[ip1] - s[i0]);

  Real dsvl_r = 0.0_rt;
  if (dsl*dsr > 0.0_rt) {
    Real dsc = 0.5_rt * (s[ip1] - s[im1]);
    dsvl_r = std::copysign(1.0_rt, dsc) * amrex::min(std::abs(dsc),amrex::min(std::abs(dsl),std::abs(dsr)));
  }

  // Interpolate s to edges

  sm = 0.5_rt * (s[i0] + s[im1]) - (1.0_rt/6.0_rt) * (dsvl_r - dsvl_l);

  // Make sure sedge lies in between adjacent cell-centered values

  sm = amrex::max(sm, amrex::min(s[i0], s[im1]));
  sm = amrex::min(sm, amrex::max(s[i0], s[im1]));


  // Compute van Leer slopes

  dsl = 2.0_rt * (s[i0] - s[im1]);
  dsr = 2.0_rt * (s[ip1] - s[i0]);

  dsvl_l = 0.0_rt;
  if (dsl*dsr > 0.0_rt) {
    Real dsc = 0.5_rt * (s[ip1] - s[im1]);
    dsvl_l = std::copysign(1.0_rt, dsc) * amrex::min(std::abs(dsc), amrex::min(std::abs(dsl),std::abs(dsr)));
  }

  dsl = 2.0_rt * (s[ip1] - s[i0]);
  dsr = 2.0_rt * (s[ip2] - s[ip1]);

  dsvl_r = 0.0_rt;
  if (dsl*dsr > 0.0_rt) {
    Real dsc = 0.5_rt * (s[ip2] - s[i0]);
    dsvl_r = std::copysign(1.0_rt, dsc) * amrex::min(std::abs(dsc), amrex::min(std::abs(dsl),std::abs(dsr)));
  }

  // Interpolate s to edges

  sp = 0.5_rt * (s[ip1] + s[i0]) - (1.0_rt/6.0_rt) * (dsvl_r - dsvl_l);

  // Make sure sedge lies in between adjacent cell-centered values

  sp = amrex::max(sp, amrex::min(s[ip1], s[i0]));
  sp = amrex::min(sp, amrex::max(s[ip1], s[i0]));


  // Flatten the parabola

  sm = flatn * sm + (1.0_rt - flatn) * s[i0];
  sp = flatn * sp + (1.0_rt - flatn) * s[i0];

  // Modify using quadratic limiters -- note this version of the limiting comes
  // from Colella and Sekora (2008), not the original PPM paper.

  if ((sp - s[i0]) * (s[i0] - sm) <= 0.0_rt) {
    sp = s[i0];
    sm = s[i0];

  } else if (std::abs(sp - s[i0]) >= 2.0_rt * std::abs(sm - s[i0])) {
    sp = 3.0_rt * s[i0] - 2.0_rt * sm;

  } else if (std::abs(sm - s[i0]) >= 2.0_rt * std::abs(sp - s[i0])) {
    sm = 3.0_rt * s[i0] - 2.0_rt * sp;
  }

}


///
/// Compute the coefficients of a parabolic reconstruction of the pressure by
/// first subtracting off the hydrostatic pressure.  This should better balance
/// the gravitational force by only limiting on the wave-generating pressure.
/// This uses the standard PPM limiters described in Colella & Woodward (1984)
///
/// @param rho          Real[5] giving the density in zones i-2, i-1, i, i+1, i+2
/// @param p           Real[5] giving the pressure in zones i-2, i-1, i, i+1, i+2
/// @param src         Real[5] the source in the velocity equation (e.g. g) in zones
///                    i-2, i-1, i, i+2, i+2
/// @param flatn       flattening coefficient
/// @param lo_bc_test  true if a physical boundary is just to the left of this zone
/// @param hi_bc_test  true if a physical boundary is just to the right of this zone
///
/// @param sm     The value of the parabola on the left edge of the zone
/// @param sp     The value of the parabola on the right edge of the zone
///
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
ppm_reconstruct_pressure(const Real* rho, const Real* p, const Real* src, const Real flatn,
                         bool lo_bc_test, bool hi_bc_test,
                         const Real dx,
                         Real& sm, Real& sp) {

    Real p0;
    Real pp1;
    Real pm1;
    Real pp2;
    Real pm2;

    if (rho[i0] >= pslope_cutoff_density) {

        Real p0_hse;
        Real pp1_hse;
        Real pp2_hse;
        Real pm1_hse;
        Real pm2_hse;

        // compute the hydrostatic pressure in each zone center starting with i0

        p0_hse = p[i0];

        pp1_hse = p0_hse + 0.25_rt*dx * (rho[i0] + rho[ip1]) * (src[i0] + src[ip1]);
        pp2_hse = pp1_hse + 0.25_rt*dx * (rho[ip1] + rho[ip2]) * (src[ip1] + src[ip2]);

        pm1_hse = p0_hse - 0.25_rt*dx * (rho[i0] + rho[im1]) * (src[i0] + src[im1]);
        pm2_hse = pm1_hse - 0.25_rt*dx * (rho[im1] + rho[im2]) * (src[im1] + src[im2]);

        // subtract off the hydrostatic pressure

        p0 = 0.0_rt;
        pp1 = p[ip1] - pp1_hse;
        pp2 = p[ip2] - pp2_hse;

        pm1 = p[im1] - pm1_hse;
        pm2 = p[im2] - pm2_hse;

        if (lo_bc_test) {
            pm1 = 0.0_rt;  // HSE is perfectly satisfied
            pm2 = 0.0_rt;
        }

        if (hi_bc_test) {
            pp1 = 0.0_rt;
            pp2 = 0.0_rt;
        }
    } else {

        // don't subtract off HSE

        pm2 = p[im2];
        pm1 = p[im1];
        p0 = p[i0];
        pp1 = p[ip1];
        pp2 = p[ip2];
    }

  // Compute van Leer slopes

  Real dsl = 2.0_rt * (pm1 - pm2);
  Real dsr = 2.0_rt * (p0 - pm1);

  Real dsvl_l = 0.0_rt;
  if (dsl*dsr > 0.0_rt) {
    Real dsc = 0.5_rt * (p0 - pm2);
    dsvl_l = std::copysign(1.0_rt, dsc) * amrex::min(std::abs(dsc),amrex::min(std::abs(dsl),std::abs(dsr)));
  }

  dsl = 2.0_rt * (p0 - pm1);
  dsr = 2.0_rt * (pp1 - p0);

  Real dsvl_r = 0.0_rt;
  if (dsl*dsr > 0.0_rt) {
    Real dsc = 0.5_rt * (pp1 - pm1);
    dsvl_r = std::copysign(1.0_rt, dsc) * amrex::min(std::abs(dsc),amrex::min(std::abs(dsl),std::abs(dsr)));
  }

  // Interpolate s to edges

  sm = 0.5_rt * (p0 + pm1) - (1.0_rt/6.0_rt) * (dsvl_r - dsvl_l);

  // Make sure sedge lies in between adjacent cell-centered values

  sm = amrex::max(sm, amrex::min(p0, pm1));
  sm = amrex::min(sm, amrex::max(p0, pm1));


  // Compute van Leer slopes

  dsl = 2.0_rt * (p0 - pm1);
  dsr = 2.0_rt * (pp1 - p0);

  dsvl_l = 0.0_rt;
  if (dsl*dsr > 0.0_rt) {
    Real dsc = 0.5_rt * (pp1 - pm1);
    dsvl_l = std::copysign(1.0_rt, dsc) * amrex::min(std::abs(dsc), amrex::min(std::abs(dsl),std::abs(dsr)));
  }

  dsl = 2.0_rt * (pp1 - p0);
  dsr = 2.0_rt * (pp2 - pp1);

  dsvl_r = 0.0_rt;
  if (dsl*dsr > 0.0_rt) {
    Real dsc = 0.5_rt * (pp2 - p0);
    dsvl_r = std::copysign(1.0_rt, dsc) * amrex::min(std::abs(dsc), amrex::min(std::abs(dsl),std::abs(dsr)));
  }

  // Interpolate s to edges

  sp = 0.5_rt * (pp1 + p0) - (1.0_rt/6.0_rt) * (dsvl_r - dsvl_l);

  // Make sure sedge lies in between adjacent cell-centered values

  sp = amrex::max(sp, amrex::min(pp1, p0));
  sp = amrex::min(sp, amrex::max(pp1, p0));


  // Flatten the parabola

  sm = flatn * sm + (1.0_rt - flatn) * p0;
  sp = flatn * sp + (1.0_rt - flatn) * p0;

  // Modify using quadratic limiters -- note this version of the limiting comes
  // from Colella and Sekora (2008), not the original PPM paper.

  if ((sp - p0) * (p0 - sm) <= 0.0_rt) {
    sp = p0;
    sm = p0;

  } else if (std::abs(sp - p0) >= 2.0_rt * std::abs(sm - p0)) {
    sp = 3.0_rt * p0 - 2.0_rt * sm;

  } else if (std::abs(sm - p0) >= 2.0_rt * std::abs(sp - p0)) {
    sm = 3.0_rt * p0 - 2.0_rt * sp;
  }

  // now correct sm and sp to be back to the full pressure by
  // adding in the hydrostatic pressure at the interface

  if (rho[i0] >= pslope_cutoff_density) {
      sp += 0.5 * dx * rho[i0] * src[i0];
      sm -= 0.5 * dx * rho[i0] * src[i0];
  }

}


///
/// Integrate under the parabola using from the left and right edges
/// with the wave speeds u-c, u, u+c
///
/// @param sm    The parabola value on the left edge of the zone
/// @param sp    The parabola value on the right edge of the zone
/// @param sc    The state value at the zone center
/// @param u     The fluid velocity
/// @param c     The sound speed
/// @param dtdx  dt/dx (timestep / cell width)
/// @param Ip    Real[3], the integrals under the parabola over
///              the 3 waves from the right zone edge
/// @param Im    Real[3], the integrals under the parabola over
///              the 3 waves from the left zone edge
///
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
ppm_int_profile(const Real sm, const Real sp, const Real sc,
                const Real u, const Real c, const Real dtdx,
                Real* Ip, Real* Im) {

  // Integrate the parabolic profile to the edge of the cell.

  // compute x-component of Ip and Im
  Real s6 = 6.0_rt * sc - 3.0_rt * (sm + sp);

  // Ip/m is the integral under the parabola for the extent
  // that a wave can travel over a timestep
  //
  // Ip integrates to the right edge of a cell
  // Im integrates to the left edge of a cell

  // u-c wave
  Real speed = u - c;
  Real sigma = std::abs(speed) * dtdx;

  // if speed == ZERO, then either branch is the same
  if (speed <= 0.0_rt) {
    Ip[0] = sp;
    Im[0] = sm + 0.5_rt * sigma * (sp - sm + (1.0_rt - (2.0_rt/3.0_rt) * sigma) * s6);
  } else {
    Ip[0] = sp - 0.5_rt * sigma * (sp - sm - (1.0_rt - (2.0_rt/3.0_rt) * sigma) * s6);
    Im[0] = sm;
  }

  // u wave
  speed = u;
  sigma = std::abs(speed) * dtdx;

  if (speed <= 0.0_rt) {
    Ip[1] = sp;
    Im[1] = sm + 0.5_rt * sigma * (sp - sm + (1.0_rt - (2.0_rt/3.0_rt) * sigma) * s6);
  } else {
    Ip[1] = sp - 0.5_rt * sigma * (sp - sm - (1.0_rt - (2.0_rt/3.0_rt) * sigma) * s6);
    Im[1] = sm;
  }

  // u+c wave
  speed = u + c;
  sigma = std::abs(speed) * dtdx;

  if (speed <= 0.0_rt) {
    Ip[2] = sp;
    Im[2] = sm + 0.5_rt * sigma * (sp - sm + (1.0_rt - (2.0_rt/3.0_rt) * sigma) * s6);
  } else {
    Ip[2] = sp - 0.5_rt * sigma * (sp - sm - (1.0_rt - (2.0_rt/3.0_rt) * sigma) * s6);
    Im[2] = sm;
  }

}

///
/// Integrate under the parabola from the left and right edges for a
/// single wave speed
///
/// @param sm    The parabola value on the left edge of the zone
/// @param sp    The parabola value on the right edge of the zone
/// @param sc    The state value at the zone center
/// @param lam   The wave speed
/// @param dtdx  dt/dx (timestep / cell width)
/// @param Ip    integral under the parabola from the right zone edge
/// @param Im    integral under the parabola from the left zone edge
///
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
ppm_int_profile_single(const Real sm, const Real sp, const Real sc,
                       const Real lam, const Real dtdx,
                       Real& Ip, Real& Im) {

  // Integrate the parabolic profile to the edge of the cell.
  // This is the MHD version.  We come in already with the eigenvalues.

  // compute x-component of Ip and Im
  Real s6 = 6.0_rt * sc - 3.0_rt * (sm + sp);

  // Ip/m is the integral under the parabola for the extent
  // that a wave can travel over a timestep
  //
  // Ip integrates to the right edge of a cell
  // Im integrates to the left edge of a cell

  Real sigma = std::abs(lam) * dtdx;

  if (lam <= 0.0_rt) {
    Ip = sp;
    Im = sm + 0.5_rt * sigma * (sp - sm + (1.0_rt - (2.0_rt/3.0_rt) * sigma) * s6);
  } else {
    Ip = sp - 0.5_rt * sigma * (sp - sm - (1.0_rt - (2.0_rt/3.0_rt) * sigma) * s6);
    Im = sm;
  }
}

#endif
