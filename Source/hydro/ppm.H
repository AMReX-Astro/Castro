
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
                const int i, const int j, const int k, const int idir,
                const bool lo_bc_test, const bool hi_bc_test,
                GpuArray<int, 3> const domlo, GpuArray<int, 3> const domhi,
                const Real flatn, Real& sm, Real& sp) {

  // This routine does the reconstruction of the zone data into a parabola.
  // Note the i, j, k here refer to zone centers, not interfaces

  // first we compute s_{i-1/2} -- the left interface value for zone i

  if (lo_bc_test && ((idir == 0 && i == domlo[0]) ||
                     (idir == 1 && j == domlo[1]) ||
                     (idir == 2 && k == domlo[2]))) {

    // use a stencil for when the current zone is on the left physical
    // boundary.  Then the left interface is on the physical boundary,
    // MC Eq. 21
    sm = (1.0_rt/12.0_rt)*(25.0_rt*s[i0] - 23.0_rt*s[ip1] + 13.0_rt*s[ip2] - 3.0_rt*s[ip3]);

  } else if (lo_bc_test && ((idir == 0 && i == domlo[0]+1) ||
                            (idir == 1 && j == domlo[1]+1) ||
                            (idir == 2 && k == domlo[2]+1))) {

    // use a stencil for when the current zone is one zone away from
    // the left physical boundary, and then the left interface is one
    // zone away from the boundary, MC Eq. 22
    sm  = (1.0_rt/12.0_rt)*(3.0_rt*s[im1] + 13.0_rt*s[i0] - 5.0_rt*s[ip1] + s[ip2]);

  } else if (hi_bc_test && ((idir == 0 && i == domhi[0]) ||
                            (idir == 1 && j == domhi[1]) ||
                            (idir == 2 && j == domhi[2]))) {

    // use a stencil for when the current zone is on the right physical boundary
    // then the left interface is one zone away from the physical boundary
    sm = (1.0_rt/12.0_rt)*(3.0_rt*s[i0] + 13.0_rt*s[im1] - 5.0_rt*s[im2] + s[im3]);

  } else {

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
  }

  // Make sure sedge lies in between adjacent cell-centered values

  sm = amrex::max(sm, amrex::min(s[i0], s[im1]));
  sm = amrex::min(sm, amrex::max(s[i0], s[im1]));


  // now we compute s_{i+1/2} -- the right interface value for zone i

  if (lo_bc_test && ((idir == 0 && i == domlo[0]) ||
                     (idir == 1 && j == domlo[1]) ||
                     (idir == 2 && k == domlo[2]))) {

    // use a stencil for when the current zone is on the left physical
    // boundary.  Then the right interface is one zone away from the
    // physical boundary,
    sp = (1.0_rt/12.0_rt)*(3.0_rt*s[i0] + 13.0_rt*s[ip1] - 5.0_rt*s[ip2] + s[ip3]);

  } else if (hi_bc_test && ((idir == 0 && i == domhi[0]) ||
                            (idir == 1 && j == domhi[1]) ||
                            (idir == 2 && j == domhi[2]))) {

    // use a stencil for when the current zone is on the right physical boundary
    // then the right interface is on the physical boundary
    sp = (1.0_rt/12.0_rt)*(25.0_rt*s[i0] - 23.0_rt*s[im1] + 13.0_rt*s[im2] - 3.0_rt*s[im3]);

  } else if (hi_bc_test && ((idir == 0 && i == domhi[0]-1) ||
                            (idir == 1 && j == domhi[1]-1) ||
                            (idir == 2 && k == domhi[2]-1))) {

    // use a stencil for when the current zone is one zone away from
    // the right physical boundary, and then the right interface is one
    // zone away from the boundary
    sp  = (1.0_rt/12.0_rt)*(3.0_rt*s[ip1] + 13.0_rt*s[i0] - 5.0_rt*s[im1] + s[im2]);

  } else {

    // Compute van Leer slopes

    Real dsl = 2.0_rt * (s[i0] - s[im1]);
    Real dsr = 2.0_rt * (s[ip1] - s[i0]);

    Real dsvl_l = 0.0_rt;
    if (dsl*dsr > 0.0_rt) {
      Real dsc = 0.5_rt * (s[ip1] - s[im1]);
      dsvl_l = std::copysign(1.0_rt, dsc) * amrex::min(std::abs(dsc), amrex::min(std::abs(dsl),std::abs(dsr)));
    }

    dsl = 2.0_rt * (s[ip1] - s[i0]);
    dsr = 2.0_rt * (s[ip2] - s[ip1]);

    Real dsvl_r = 0.0_rt;
    if (dsl*dsr > 0.0_rt) {
      Real dsc = 0.5_rt * (s[ip2] - s[i0]);
      dsvl_r = std::copysign(1.0_rt, dsc) * amrex::min(std::abs(dsc), amrex::min(std::abs(dsl),std::abs(dsr)));
    }

    // Interpolate s to edges

    sp = 0.5_rt * (s[ip1] + s[i0]) - (1.0_rt/6.0_rt) * (dsvl_r - dsvl_l);

  }

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
                       const Real flatn, Real& pm, Real& pp) {

  // This routine does the reconstruction of the zone data into a parabola.

  // first integrate HSE outward from zone i and subtract it off from
  // the neighboring zones.
  Real p0 = 0.0_rt;
  Real pp1 = p[ip1] - (p0 + 0.25_rt * dx * (rho[i0] + rho[ip1]) * (src[i0] + src[ip1]));
  Real pp2 = p[ip2] - (p0 + 0.25_rt * dx * (rho[i0] + rho[ip1]) * (src[i0] + src[ip1])
                          + 0.25_rt * dx * (rho[ip1] + rho[ip2]) * (src[ip1] + src[ip2]));
  Real pm1 = p[im1] - (p0 - 0.25_rt * dx * (rho[i0] + rho[im1]) * (src[i0] + src[im1]));
  Real pm2 = p[im2] - (p0 - 0.25_rt * dx * (rho[i0] + rho[im1]) * (src[i0] + src[im1])
                          - 0.25_rt * dx * (rho[im1] + rho[im2]) * (src[im1] + src[im2]));

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

  pm = 0.5_rt * (p0 + pm1) - (1.0_rt/6.0_rt) * (dsvl_r - dsvl_l);

  // Make sure sedge lies in between adjacent cell-centered values

  pm = amrex::max(pm, amrex::min(p0, pm1));
  pm = amrex::min(pm, amrex::max(p0, pm1));


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

  // now correct the state to have the HSE contribution
  pm += p[i0] - 0.5_rt * dx * rho[i0] * src[i0];
  pp += p[i0] + 0.5_rt * dx * rho[i0] * src[i0];

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
