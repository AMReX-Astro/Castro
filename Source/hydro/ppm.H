
#include <cmath>
#include <reconstruction.H>

#ifndef PPM_H
#define PPM_H

using namespace amrex;
using namespace reconstruction;

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
/// @param s      Real[nslp] the state to be reconstructed in zones i-3, i-2, i-1, i, i+1, i+2, i+3
/// @param flatn  The flattening coefficient
/// @param sm     The value of the parabola on the left edge of the zone
/// @param sp     The value of the parabola on the right edge of the zone
///
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
ppm_reconstruct(const Real* s,
                const int i, const int j, const int k, const int idir,
                const Real dx,
                const bool lo_bc_test, const bool hi_bc_test, const bool is_axisymmetric,
                const GpuArray<int, 3>& domlo, const GpuArray<int, 3>& domhi,
                const Real flatn, Real& sm, Real& sp) {

    // first we compute s_{i-1/2} -- the left interface value for zone i

    if (ppm_do_limiting) {


        if (lo_bc_test && ((idir == 0 && i == domlo[0]) ||
                           (idir == 1 && j == domlo[1]) ||
                           (idir == 2 && k == domlo[2]))) {

            // use a stencil for when the current zone is on the left physical
            // boundary.  Then the left interface is on the physical boundary,
            // MC Eq. 21

            if (idir == 0 && is_axisymmetric) {
                sm = (1.0_rt/144.0_rt)*(415.0_rt*s[i0] - 483.0_rt*s[ip1] + 275.0_rt*s[ip2] - 63.0_rt*s[ip3]);
            } else {
          // Cartesian
                sm = (1.0_rt/12.0_rt)*(25.0_rt*s[i0] - 23.0_rt*s[ip1] + 13.0_rt*s[ip2] - 3.0_rt*s[ip3]);
            }

        } else if (lo_bc_test && ((idir == 0 && i == domlo[0]+1) ||
                                  (idir == 1 && j == domlo[1]+1) ||
                                  (idir == 2 && k == domlo[2]+1))) {

            // use a stencil for when the current zone is one zone away from
            // the left physical boundary, and then the left interface is one
            // zone away from the boundary, MC Eq. 22

            if (idir == 0 && is_axisymmetric) {
                sm  = (1.0_rt/96.0_rt)*(37.0_rt*s[im1] + 87.0_rt*s[i0] - 35.0_rt*s[ip1] + 7.0_rt*s[ip2]);
            } else {
                sm  = (1.0_rt/12.0_rt)*(3.0_rt*s[im1] + 13.0_rt*s[i0] - 5.0_rt*s[ip1] + s[ip2]);
            }

            // Make sure sedge lies in between adjacent cell-centered values

            sm = amrex::max(sm, amrex::min(s[i0], s[im1]));
            sm = amrex::min(sm, amrex::max(s[i0], s[im1]));

        } else if (hi_bc_test && ((idir == 0 && i == domhi[0]) ||
                                  (idir == 1 && j == domhi[1]) ||
                                  (idir == 2 && k == domhi[2]))) {

            // use a stencil for when the current zone is on the right physical boundary
            // then the left interface is one zone away from the physical boundary
            sm = (1.0_rt/12.0_rt)*(3.0_rt*s[i0] + 13.0_rt*s[im1] - 5.0_rt*s[im2] + s[im3]);

            // Make sure sedge lies in between adjacent cell-centered values

            sm = amrex::max(sm, amrex::min(s[i0], s[im1]));
            sm = amrex::min(sm, amrex::max(s[i0], s[im1]));

        } else {

            if (idir == 0 && is_axisymmetric) {

                const Real r_m12 = (static_cast<int>(i) - 0.5_rt) * dx;

                Real dr4_term = std::pow(dx, 4) * (-12.0_rt * s[im2] + 60.0_rt * s[im1]  + 60.0_rt * s[i0] - 12.0_rt * s[ip1]);
                Real dr3_term = r_m12 * std::pow(dx, 3) * (-s[im2]  - 27.0_rt * s[im1] + 27.0_rt * s[i0] + s[ip1]);
                Real dr2_term = std::pow(r_m12, 2) * std::pow(dx, 2) * (30.0_rt * s[im2] - 210.0_rt * s[im1] - 210.0_rt * s[i0] + 30.0_rt * s[ip1]);
                Real dr1_term = std::pow(r_m12, 3) * dx * (-s[im2] - 13.0_rt * s[im1] - 13.0_rt * s[i0] + s[ip1]);
                Real dr0_term = std::pow(r_m12, 4) * (-10.0_rt * s[im2] + 70.0_rt * s[im1] + 70.0_rt * s[i0] - 10.0_rt * s[ip1]);

                Real denom = 24.0_rt * (4.0_rt * std::pow(dx, 4) - 15.0_rt * std::pow(dx, 2) * std::pow(r_m12, 2) + 5.0_rt * std::pow(r_m12, 4));

                sm = (dr4_term + dr3_term + dr2_term + dr1_term + dr0_term) / denom;

            } else {

                // Compute van Leer slopes

                Real dsl = 2.0_rt * (s[im1] - s[im2]);
                Real dsr = 2.0_rt * (s[i0] - s[im1]);

                Real dsvl_l = 0.0_rt;
                if (dsl*dsr > 0.0_rt) {
                    Real dsc = 0.5_rt * (s[i0] - s[im2]);
                    dsvl_l = std::copysign(1.0_rt, dsc) * amrex::min(std::abs(dsc), std::abs(dsl), std::abs(dsr));
                }

                dsl = 2.0_rt * (s[i0] - s[im1]);
                dsr = 2.0_rt * (s[ip1] - s[i0]);

                Real dsvl_r = 0.0_rt;
                if (dsl*dsr > 0.0_rt) {
                    Real dsc = 0.5_rt * (s[ip1] - s[im1]);
                    dsvl_r = std::copysign(1.0_rt, dsc) * amrex::min(std::abs(dsc), std::abs(dsl), std::abs(dsr));
                }

                // Interpolate s to edges

                sm = 0.5_rt * (s[i0] + s[im1]) - (1.0_rt/6.0_rt) * (dsvl_r - dsvl_l);

            }

            // Make sure sedge lies in between adjacent cell-centered values

            sm = std::clamp(sm, std::min(s[i0], s[im1]), std::max(s[i0], s[im1]));

        }

        // now we compute s_{i+1/2} -- the right interface value for zone i

        if (lo_bc_test && ((idir == 0 && i == domlo[0]) ||
                           (idir == 1 && j == domlo[1]) ||
                           (idir == 2 && k == domlo[2]))) {

            // use a stencil for when the current zone is on the left physical
            // boundary.  Then the right interface is one zone away from the
            // physical boundary,

            if (idir == 0 && is_axisymmetric) {
                sp = (1.0_rt/96.0_rt)*(37.0_rt*s[i0] + 87.0_rt*s[ip1] - 35.0_rt*s[ip2] + 7.0_rt*s[ip3]);
            } else {
                sp = (1.0_rt/12.0_rt)*(3.0_rt*s[i0] + 13.0_rt*s[ip1] - 5.0_rt*s[ip2] + s[ip3]);
            }

            // Make sure sedge lies in between adjacent cell-centered values

            sp = amrex::max(sp, amrex::min(s[ip1], s[i0]));
            sp = amrex::min(sp, amrex::max(s[ip1], s[i0]));

        } else if (hi_bc_test && ((idir == 0 && i == domhi[0]) ||
                                  (idir == 1 && j == domhi[1]) ||
                                  (idir == 2 && k == domhi[2]))) {

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

            // Make sure sedge lies in between adjacent cell-centered values

            sp = amrex::max(sp, amrex::min(s[ip1], s[i0]));
            sp = amrex::min(sp, amrex::max(s[ip1], s[i0]));

        } else {

            if (idir == 0 && is_axisymmetric) {

                const Real r_p12 = (static_cast<int>(i) + 0.5_rt) * dx;

                Real dr4_term = std::pow(dx, 4) * (-12.0_rt * s[im1] + 60.0_rt * s[i0]  + 60.0_rt * s[ip1] - 12.0_rt * s[ip2]);
                Real dr3_term = r_p12 * std::pow(dx, 3) * (-s[im1]  - 27.0_rt * s[i0] + 27.0_rt * s[ip1] + s[ip2]);
                Real dr2_term = std::pow(r_p12, 2) * std::pow(dx, 2) * (30.0_rt * s[im1] - 210.0_rt * s[i0] - 210.0_rt * s[ip1] + 30.0_rt * s[ip2]);
                Real dr1_term = std::pow(r_p12, 3) * dx * (-s[im1] - 13.0_rt * s[i0] - 13.0_rt * s[ip1] + s[ip2]);
                Real dr0_term = std::pow(r_p12, 4) * (-10.0_rt * s[im1] + 70.0_rt * s[i0] + 70.0_rt * s[ip1] - 10.0_rt * s[ip2]);

                Real denom = 24.0_rt * (4.0_rt * std::pow(dx, 4) - 15.0_rt * std::pow(dx, 2) * std::pow(r_p12, 2) + 5.0_rt * std::pow(r_p12, 4));

                sp = (dr4_term + dr3_term + dr2_term + dr1_term + dr0_term) / denom;

            } else {

                // Compute van Leer slopes

                Real dsl = 2.0_rt * (s[i0] - s[im1]);
                Real dsr = 2.0_rt * (s[ip1] - s[i0]);

                Real dsvl_l = 0.0_rt;
                if (dsl*dsr > 0.0_rt) {
                    Real dsc = 0.5_rt * (s[ip1] - s[im1]);
                    dsvl_l = std::copysign(1.0_rt, dsc) * amrex::min(std::abs(dsc), std::abs(dsl), std::abs(dsr));
                }

                dsl = 2.0_rt * (s[ip1] - s[i0]);
                dsr = 2.0_rt * (s[ip2] - s[ip1]);

                Real dsvl_r = 0.0_rt;
                if (dsl*dsr > 0.0_rt) {
                    Real dsc = 0.5_rt * (s[ip2] - s[i0]);
                    dsvl_r = std::copysign(1.0_rt, dsc) * amrex::min(std::abs(dsc), std::abs(dsl), std::abs(dsr));
                }

                // Interpolate s to edges

                sp = 0.5_rt * (s[ip1] + s[i0]) - (1.0_rt/6.0_rt) * (dsvl_r - dsvl_l);

            }

            // Make sure sedge lies in between adjacent cell-centered values

            sp = std::clamp(sp, std::min(s[ip1], s[i0]), std::max(s[ip1], s[i0]));

        }

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

    } else {

        // unlimited PPM reconstruction  (Eq. 1.9 in the PPM paper)

        sm = (7.0_rt/12.0_rt) * (s[i0] + s[im1]) - (1.0_rt/12.0_rt) * (s[im2] + s[ip1]);
        sp = (7.0_rt/12.0_rt) * (s[ip1] + s[i0]) - (1.0_rt/12.0_rt) * (s[im1] + s[ip2]);

    }

}


///
/// Compute the coefficients of a parabolic reconstruction of the pressure by
/// first subtracting off the hydrostatic pressure.  This should better balance
/// the gravitational force by only limiting on the wave-generating pressure.
/// This uses the standard PPM limiters described in Colella & Woodward (1984)
///
/// @param rho          Real[nslp] giving the density in zones i-2, i-1, i, i+1, i+2
/// @param p           Real[nslp] giving the pressure in zones i-2, i-1, i, i+1, i+2
/// @param src         Real[nslp] the source in the velocity equation (e.g. g) in zones
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
ppm_reconstruct_pslope(const Real* rho, const Real* p, const Real* src,
                       const Real dx,
                       const int i, const int j, const int k, const int idir,
                       bool lo_bc_test, bool hi_bc_test, const bool is_axisymmetric,
                       const GpuArray<int, 3>& domlo, const GpuArray<int, 3>& domhi,
                       const Real flatn,
                       Real& sm, Real& sp) {

    Real tp[nslp];

    // compute the hydrostatic pressure in each zone center starting with i0

    Real p0_hse = p[i0];

    Real pp1_hse = p0_hse + 0.5_rt * dx * (rho[i0] * src[i0] + rho[ip1] * src[ip1]);
    Real pp2_hse = pp1_hse + 0.5_rt * dx * (rho[ip1] * src[ip1] + rho[ip2] * src[ip2]);
    Real pp3_hse = pp2_hse + 0.5_rt * dx * (rho[ip2] * src[ip2] + rho[ip3] * src[ip3]);

    Real pm1_hse = p0_hse - 0.5_rt * dx * (rho[i0] * src[i0] + rho[im1] * src[im1]);
    Real pm2_hse = pm1_hse - 0.5_rt * dx * (rho[im1] * src[im1] + rho[im2] * src[im2]);
    Real pm3_hse = pm2_hse - 0.5_rt * dx * (rho[im2] * src[im2] + rho[im3] * src[im3]);

    // this only makes sense if the pressures are positive
    bool ptest = p0_hse < 0.0 ||
                 pp1_hse < 0.0 ||
                 pp2_hse < 0.0 ||
                 pp3_hse < 0.0 ||
                 pm1_hse < 0.0 ||
                 pm2_hse < 0.0 ||
                 pm3_hse < 0.0;


    if (!ptest && rho[i0] >= pslope_cutoff_density) {

        // subtract off the hydrostatic pressure

        tp[i0] = p[i0] - p0_hse;

        tp[ip1] = p[ip1] - pp1_hse;
        tp[ip2] = p[ip2] - pp2_hse;
        tp[ip3] = p[ip3] - pp3_hse;

        tp[im1] = p[im1] - pm1_hse;
        tp[im2] = p[im2] - pm2_hse;
        tp[im3] = p[im3] - pm3_hse;

    } else {

        // don't subtract off HSE

        tp[im3] = p[im3];
        tp[im2] = p[im2];
        tp[im1] = p[im1];
        tp[i0] = p[i0];
        tp[ip1]= p[ip1];
        tp[ip2] = p[ip2];
        tp[ip3] = p[ip3];
    }

    ppm_reconstruct(tp, i, j, k, idir, dx, lo_bc_test, hi_bc_test, is_axisymmetric, domlo, domhi, flatn, sm, sp);


    // now correct sm and sp to be back to the full pressure by
    // adding in the hydrostatic pressure at the interface

    if (!ptest && rho[i0] >= pslope_cutoff_density) {
        sp += p[i0] + 0.5_rt * dx * rho[i0] * src[i0];
        sm += p[i0] - 0.5_rt * dx * rho[i0] * src[i0];
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
