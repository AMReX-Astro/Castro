#ifndef SLOPE_H
#define SLOPE_H

#include <Castro.H>

#ifdef RADIATION
#include <Radiation.H>
#endif

#ifdef GRAVITY
#include <Gravity.H>
#endif

#include <cmath>

#include <reconstruction.H>

using namespace amrex;
using namespace reconstruction;

///
/// Construct the piecewise linear slopes for an input state q
///
/// @param q               Real[nslp] giving q_{i-2} through q_{i+2}, used in the reconstruction
/// @param flatn           flattening coefficient
/// @param bnd_lo_reflect  true if a symmetry boundary is just to the left of this zone
/// @param bnd_hi_reflect  true if a symmetry boundary is just to the right of this zone
///
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real
uslope(const Real* q, const Real flatn,
       bool bnd_lo_reflect, bool bnd_hi_reflect) {

  Real dq;

  if (plm_iorder == 1) {

    // first order -- piecewise constant slopes
    dq = 0.0;

  } else {

    // second-order -- piecewise linear slopes
    if (plm_limiter == 1) {
      // the 2nd order MC limiter

      Real dlft = 2.0_rt*(q[i0] - q[im1]);
      Real drgt = 2.0_rt*(q[ip1] - q[i0]);
      Real dcen = 0.25_rt * (dlft + drgt);
      Real dsgn = std::copysign(1.0_rt, dcen);
      Real slop = amrex::min(std::abs(dlft), std::abs(drgt));
      Real dlim = dlft*drgt >= 0.0_rt ? slop : 0.0_rt;

      dq = flatn * dsgn * amrex::min(dlim, std::abs(dcen));

    } else {
      // the 4th order MC limiter

      // special consideration for reflecting BCs -- see
      // Saltzman p. 162 (but note that Saltzman has a
      // sign error)

      Real qm2 = q[im2];
      Real qm1 = q[im1];
      Real q0 = q[i0];
      Real qp1 = q[ip1];
      Real qp2 = q[ip2];

      if (bnd_lo_reflect) {
        qm2 = -qp1;
        qm1 = -3.0_rt*q0 + qp1 - 0.125_rt*(qp2 + qp1);
      }

      if (bnd_hi_reflect) {
        qp2 = -qm1;
        qp1 = -3.0_rt*q0 + qm1 - 0.125_rt*(qm2 + qm1);
      }

      // First compute Fromm slopes

      // df at i+1
      Real dlftp1 = 2.0_rt*(qp1 - q0);
      Real drgtp1 = 2.0_rt*(qp2 - qp1);
      Real dcen = 0.25_rt * (dlftp1 + drgtp1);
      Real dsgn = std::copysign(1.0_rt, dcen);
      Real slop = amrex::min(std::abs(dlftp1), std::abs(drgtp1));
      Real dlim = dlftp1*drgtp1 >= 0.0_rt ? slop : 0.0_rt;
      Real dfp1 = dsgn*amrex::min(dlim, std::abs(dcen));

      // df at i-1
      Real dlftm1 = 2.0_rt*(qm1 - qm2);
      Real drgtm1 = 2.0_rt*(q0 - qm1);
      dcen = 0.25_rt * (dlftm1 + drgtm1);
      dsgn = std::copysign(1.0_rt, dcen);
      slop = amrex::min(std::abs(dlftm1), std::abs(drgtm1));
      dlim = dlftm1*drgtm1 >= 0.0_rt ? slop : 0.0_rt;
      Real dfm1 = dsgn*amrex::min(dlim, std::abs(dcen));

      // Now compute limited fourth order slopes at i
      Real dlft = drgtm1;
      Real drgt = dlftp1;
      dcen = 0.25_rt * (dlft + drgt);
      dsgn = std::copysign(1.0_rt, dcen);
      slop = amrex::min(std::abs(dlft), std::abs(drgt));
      dlim = dlft*drgt >= 0.0_rt ? slop : 0.0_rt;

      Real dq1 = (4.0_rt/3.0_rt)*dcen - (1.0_rt/6.0_rt)*(dfp1 + dfm1);
      dq = flatn * dsgn * amrex::min(dlim, std::abs(dq1));
    }
  }

  return dq;
}

///
/// Compute the pressure slope using a well-balanced procedure where we subtract
/// HSE off of the pressure before limiting
///
/// @param rho          Real[nslp] giving the density in zones i-2, i-1, i, i+1, i+2
/// @param p           Real[nslp] giving the pressure in zones i-2, i-1, i, i+1, i+2
/// @param src         Real[nslp] the source in the velocity equation (e.g. g) in zones
///                    i-2, i-1, i, i+2, i+2
/// @param flatn       flattening coefficient
/// @param lo_bc_test  true if a physical boundary is just to the left of this zone
/// @param hi_bc_test  true if a physical boundary is just to the right of this zone
///
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
pslope(const Real* rho, const Real* p, const Real* src, const Real flatn,
       bool lo_bc_test, bool hi_bc_test,
       const Real dx, Real& dp) {

  // dp should come in as the current estimate of the pressure slope
  // and we overwrite it if needed

  if (plm_iorder == 1) {

    dp = 0.0_rt;

  } else {

    if (rho[i0] < pslope_cutoff_density) {
      return;
    }

    Real p0;
    Real pp1;
    Real pm1;
    Real pp2;
    Real pm2;

    Real p0_hse;
    Real pp1_hse;
    Real pp2_hse;
    Real pm1_hse;
    Real pm2_hse;

    p0_hse = p[i0];

    pp1_hse = p0_hse + 0.25_rt*dx * (rho[i0] + rho[ip1]) * (src[i0] + src[ip1]);
    pp2_hse = pp1_hse + 0.25_rt*dx * (rho[ip1] + rho[ip2]) * (src[ip1] + src[ip2]);

    pm1_hse = p0_hse - 0.25_rt*dx * (rho[i0] + rho[im1]) * (src[i0] + src[im1]);
    pm2_hse = pm1_hse - 0.25_rt*dx * (rho[im1] + rho[im2]) * (src[im1] + src[im2]);

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


    // First compute Fromm slopes

    // Dp at i+1

    // we need dp_{i+1} = p_{i+1} - p_i
    // and     dp_{i+2} = p_{i+2} - p_{i-1}
    //
    // then we compute
    //   Dp_{i+1} = min{ 2 |dp_{i+1}|, 2 |dp_{i+2}|, 1/2 |dp_{i+1} + dp_{i+2}|}
    //
    // we use the value above with HSE removed

    Real dlftp1 = pp1 - p0;
    Real drgtp1 = pp2 - pp1;

    Real dcen = 0.5_rt*(dlftp1 + drgtp1);
    Real dsgn = std::copysign(1.0_rt, dcen);
    Real dlim = dlftp1*drgtp1 >= 0.0_rt ? 2.0_rt * amrex::min(std::abs(dlftp1), std::abs(drgtp1)) : 0.0_rt;
    Real dfp1 = dsgn*amrex::min(dlim, std::abs(dcen));

    // Dp at i-1

    // we need dp_{i-1} = p_{i-1} - p_{i-2}
    // and     dp_i = p_i - p_{i-1}
    //
    // then we compute
    //   Dp_{i-1} = min{ 2 |dp_{ii1}|, 2 |dp_i|, 1/2 |dp_{i-1} + dp_i|}
    //
    // again use the values with HSE removed

    Real dlftm1 = pm1 - pm2;
    Real drgtm1 = p0 - pm1;

    dcen = 0.5_rt*(dlftm1 + drgtm1);
    dsgn = std::copysign(1.0_rt, dcen);
    dlim = dlftm1*drgtm1 >= 0.0_rt ? 2.0_rt * amrex::min(std::abs(dlftm1), std::abs(drgtm1)) : 0.0_rt;

    Real dfm1 = dsgn*amrex::min(dlim, std::abs(dcen));

    // Now limited fourth order slopes at i

    // this uses Dp_{i+1} and Dp_{i-1} computed above as well as
    // dp_{i+1} and dp_i

    // dp_i
    Real dlft = drgtm1;

    // dp_{i+1}
    Real drgt = dlftp1;

    dcen = 0.5_rt*(dlft + drgt);
    dsgn = std::copysign(1.0_rt, dcen);
    dlim = dlft*drgt >= 0.0_rt ? 2.0_rt * amrex::min(std::abs(dlft), std::abs(drgt)) : 0.0_rt;

    Real dp1 = (4.0_rt/3.0_rt)*dcen - (1.0_rt/6.0_rt)*(dfp1 + dfm1);
    dp = flatn * dsgn*amrex::min(dlim, std::abs(dp1));

    dp += rho[i0] * src[i0] * dx;
  }

}
#endif
