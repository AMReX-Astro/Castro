
#include <cmath>

#ifndef _ppm_H_
#define _ppm_H_

using namespace amrex;

constexpr int im3 = 0;
constexpr int im2 = 1;
constexpr int im1 = 2;
constexpr int i0 = 3;
constexpr int ip1 = 4;
constexpr int ip2 = 5;
constexpr int ip3 = 6;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
ppm_reconstruct(const Real* s,
                const Real flatn, Real& sm, Real& sp) {

  // This routine does the reconstruction of the zone data into a parabola.

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


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
ppm_reconstruct_pslope(const Real* p, const Real* rho, const Real* src,
                       const Real dx,
                       const int i, const int j, const int k, const int idir,
                       const bool lo_bc_test, const bool hi_bc_test,
                       GpuArray<int, 3> const domlo, GpuArray<int, 3> const domhi,
                       const Real flatn, Real& pm, Real& pp) {

  // This routine does the reconstruction of the zone data into a parabola.

  // first integrate HSE outward from zone i and subtract it off from
  // the neighboring zones.
  Real p0_hse = 0.0_rt;
  Real pp1_hse = p[ip1] - (p0_hse + 0.25_rt * dx * (rho[i0] + rho[ip1]) * (src[i0] + src[ip1]));
  Real pp2_hse = p[ip2] - (p0_hse + 0.25_rt * dx * (rho[i0] + rho[ip1]) * (src[i0] + src[ip1])
                                  + 0.25_rt * dx * (rho[ip1] + rho[ip2]) * (src[ip1] + src[ip2]));
  Real pm1_hse = p[im1] - (p0_hse - 0.25_rt * dx * (rho[i0] + rho[im1]) * (src[i0] + src[im1]));
  Real pm2_hse = p[im2] - (p0_hse - 0.25_rt * dx * (rho[i0] + rho[im1]) * (src[i0] + src[im1])
                                  - 0.25_rt * dx * (rho[im1] + rho[im2]) * (src[im1] + src[im2]));

  Real p0;
  Real pp1;
  Real pp2;
  Real pm1;
  Real pm2;

  int idx;
  if (idir == 0) {
    idx = i;
  } else if (idir == 1) {
    idx = j;
  } else {
    idx = k;
  }

  // i-1/2 interface

  bool left_is_hse = false;

  if (lo_bc_test && (idx == domlo[idir])) {

    // we are computing p_{domlo-1/2} and two zones are to the left of
    // the interface are ghost zones.
    left_is_hse = true;

    pm2 = 0.0_rt;     // perfect HSE
    pm1 = 0.0_rt;     // perfect HSE
    p0 = p0_hse;
    pp1 = pp1_hse;
    pp2 = pp2_hse;

  } else if (lo_bc_test && (idx == domlo[idir]+1)) {

    // we are computing p_{domlo+1/2} and one zone in the cubic
    // interpolant is to the left of the domain boundary.
    left_is_hse = true;

    pm2 = 0.0_rt;
    pm1 = pm1_hse;
    p0 = p0_hse;
    pp1 = pp1_hse;
    pp2 = pp2_hse;

  } else if (hi_bc_test && (idx == domhi[idir])) {

    // we are computing p_{domhi-1/2} and one zone in the cubic
    // interpolant is to the right of the domain boundary.
    left_is_hse = true;

    pm2 = pm2_hse;
    pm1 = pm1_hse;
    p0 = p0_hse;
    pp1 = pp1_hse;
    pp2 = 0.0_rt;

  } else {

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
    dsvl_l = std::copysign(1.0_rt, dsc) * amrex::min(std::abs(dsc),
                                                     amrex::min(std::abs(dsl), std::abs(dsr)));
  }

  dsl = 2.0_rt * (p0 - pm1);
  dsr = 2.0_rt * (pp1 - p0);

  Real dsvl_r = 0.0_rt;
  if (dsl*dsr > 0.0_rt) {
    Real dsc = 0.5_rt * (pp1 - pm1);
    dsvl_r = std::copysign(1.0_rt, dsc) * amrex::min(std::abs(dsc),amrex::min(std::abs(dsl),std::abs(dsr)));
  }

  // Interpolate s to edges

  pm = 0.5_rt * (p0 + pm1) - (1.0_rt/6.0_rt) * (dsvl_r - dsvl_l);

  // Make sure sedge lies in between adjacent cell-centered values

  pm = amrex::max(pm, amrex::min(p0, pm1));
  pm = amrex::min(pm, amrex::max(p0, pm1));

  // correct for HSE if needed
  if (left_is_hse) {
    pm += p[i0] - 0.5_rt * dx * rho[i0] * src[i0];
  }

  // i+1/2 interface

  bool right_is_hse = false;

  if (lo_bc_test && (idx == domlo[idir])) {

    // we are computing p_{domlo+1/2} and one of the zones in the
    // cubic interpolant will be ghost zones
    right_is_hse = true;

    pm2 = 0.0_rt;
    pm1 = pm1_hse;
    p0 = p0_hse;
    pp1 = pp1_hse;
    pp2 = pp2_hse;

  } else if (hi_bc_test && (idx == domhi[idir]-1)) {

    // we are computing p_{domhi-1/2} and one zone in the cubic
    // interpolant is a ghost zone.
    right_is_hse = true;

    pm2 = pm2_hse;
    pm1 = pm1_hse;
    p0 = p0_hse;
    pp1 = pp1_hse;
    pp2 = 0.0_rt;

  } else if (hi_bc_test && (idx == domhi[idir])) {

    // we are computing p_{domhi+1/2} and two zones in the cubic
    // interpolant are ghost cells
    right_is_hse = true;

    pm2 = pm2_hse;
    pm1 = pm1_hse;
    p0 = p0_hse;
    pp1 = 0.0_rt;
    pp2 = 0.0_rt;

  } else {

    pm2 = p[im2];
    pm1 = p[im1];
    p0 = p[i0];
    pp1 = p[ip1];
    pp2 = p[ip2];

  }

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

  pp = 0.5_rt * (pp1 + p0) - (1.0_rt/6.0_rt) * (dsvl_r - dsvl_l);

  // Make sure sedge lies in between adjacent cell-centered values

  pp = amrex::max(pp, amrex::min(pp1, p0));
  pp = amrex::min(pp, amrex::max(pp1, p0));

  // correct for HSE if needed
  if (right_is_hse) {
    pp += p[i0] + 0.5_rt * dx * rho[i0] * src[i0];
  }

  // Flatten the parabola

  pm = flatn * pm + (1.0_rt - flatn) * p0;
  pp = flatn * pp + (1.0_rt - flatn) * p0;

  // Modify using quadratic limiters -- note this version of the limiting comes
  // from Colella and Sekora (2008), not the original PPM paper.

  if ((pp - p0) * (p0 - pm) <= 0.0_rt) {
    pp = p0;
    pm = p0;

  } else if (std::abs(pp - p0) >= 2.0_rt * std::abs(pm - p0)) {
    pp = 3.0_rt * p0 - 2.0_rt * pm;

  } else if (std::abs(pm - p0) >= 2.0_rt * std::abs(pp - p0)) {
    pm = 3.0_rt * p0 - 2.0_rt * pp;
  }

}


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


#endif
