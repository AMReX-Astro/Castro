
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
/// @param s      Real[nslp] the state to be reconstructed in zones i-2, i-1, i, i+1, i+2
/// @param flatn  The flattening coefficient
/// @param sm     The value of the parabola on the left edge of the zone
/// @param sp     The value of the parabola on the right edge of the zone
///
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
ppm_reconstruct(const Real* s,
                const Real flatn, Real& sm, Real& sp) {


    if (ppm_do_limiting) {

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
          dsvl_r = std::copysign(1.0_rt, dsc) * amrex::min(std::abs(dsc), std::abs(dsl),std::abs(dsr));
        }

        // Interpolate s to edges

        sm = 0.5_rt * (s[i0] + s[im1]) - (1.0_rt/6.0_rt) * (dsvl_r - dsvl_l);

        // Make sure sedge lies in between adjacent cell-centered values

        sm = std::clamp(sm, std::min(s[i0], s[im1]), std::max(s[i0], s[im1]));

        // Compute van Leer slopes

        dsl = 2.0_rt * (s[i0] - s[im1]);
        dsr = 2.0_rt * (s[ip1] - s[i0]);

        dsvl_l = 0.0_rt;
        if (dsl*dsr > 0.0_rt) {
          Real dsc = 0.5_rt * (s[ip1] - s[im1]);
          dsvl_l = std::copysign(1.0_rt, dsc) * amrex::min(std::abs(dsc), std::abs(dsl), std::abs(dsr));
        }

        dsl = 2.0_rt * (s[ip1] - s[i0]);
        dsr = 2.0_rt * (s[ip2] - s[ip1]);

        dsvl_r = 0.0_rt;
        if (dsl*dsr > 0.0_rt) {
          Real dsc = 0.5_rt * (s[ip2] - s[i0]);
          dsvl_r = std::copysign(1.0_rt, dsc) * amrex::min(std::abs(dsc), std::abs(dsl), std::abs(dsr));
        }

        // Interpolate s to edges

        sp = 0.5_rt * (s[ip1] + s[i0]) - (1.0_rt/6.0_rt) * (dsvl_r - dsvl_l);

        // Make sure sedge lies in between adjacent cell-centered values

        sp = std::clamp(sp, std::min(s[ip1], s[i0]), std::max(s[ip1], s[i0]));


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
ppm_reconstruct_pslope(const Real* rho, const Real* p, const Real* src, const Real flatn,
                       const Real dx,
                       Real& sm, Real& sp) {

    Real tp[nslp];

    // compute the hydrostatic pressure in each zone center starting with i0

    Real p0_hse = p[i0];

    Real pp1_hse = p0_hse + 0.25_rt*dx * (rho[i0] + rho[ip1]) * (src[i0] + src[ip1]);
    Real pp2_hse = pp1_hse + 0.25_rt*dx * (rho[ip1] + rho[ip2]) * (src[ip1] + src[ip2]);

    Real pm1_hse = p0_hse - 0.25_rt*dx * (rho[i0] + rho[im1]) * (src[i0] + src[im1]);
    Real pm2_hse = pm1_hse - 0.25_rt*dx * (rho[im1] + rho[im2]) * (src[im1] + src[im2]);

    // this only makes sense if the pressures are positive
    bool ptest = p0_hse < 0.0 ||
                 pp1_hse < 0.0 ||
                 pp2_hse < 0.0 ||
                 pm1_hse < 0.0 ||
                 pm2_hse < 0.0;


    if (!ptest && rho[i0] >= pslope_cutoff_density) {

        // subtract off the hydrostatic pressure

        tp[i0] = p[i0] - p0_hse;

        tp[ip1] = p[ip1] - pp1_hse;
        tp[ip2] = p[ip2] - pp2_hse;

        tp[im1] = p[im1] - pm1_hse;
        tp[im2] = p[im2] - pm2_hse;

    } else {

        // don't subtract off HSE

        tp[im2] = p[im2];
        tp[im1] = p[im1];
        tp[i0] = p[i0];
        tp[ip1]= p[ip1];
        tp[ip2] = p[ip2];
    }

    ppm_reconstruct(tp, flatn, sm, sp);


    // now correct sm and sp to be back to the full pressure by
    // adding in the hydrostatic pressure at the interface

    if (!ptest && rho[i0] >= pslope_cutoff_density) {
        sp += p[i0] + 0.5 * dx * rho[i0] * src[i0];
        sm += p[i0] - 0.5 * dx * rho[i0] * src[i0];
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
