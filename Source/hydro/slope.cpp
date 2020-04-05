#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

#ifdef GRAVITY
#include "Gravity.H"
#endif

#include <cmath>

#include <ppm.H>

using namespace amrex;

void
Castro::uslope(const Box& bx, const int idir,
               Array4<Real const> const q_arr, const int n,
               Array4<Real const> const flatn_arr,
               Array4<Real> const dq) {

  const auto dx = geom.CellSizeArray();

  const auto domlo = geom.Domain().loVect3d();
  const auto domhi = geom.Domain().hiVect3d();

  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();

  bool lo_bc_test = lo_bc[idir] == Symmetry;
  bool hi_bc_test = hi_bc[idir] == Symmetry;

#ifdef GRAVITY
  Real lconst_grav =  gravity->get_const_grav();
#else
  Real lconst_grav = 0.0_rt;
#endif

  if (plm_iorder == 1) {

    // first order -- piecewise constant slopes
    AMREX_PARALLEL_FOR_3D(bx, i, j, k,
    {
      dq(i,j,k,n) = 0.0;
    });

  } else {

    // second-order -- piecewise linear slopes

    const int lplm_limiter = plm_limiter;

    AMREX_PARALLEL_FOR_3D(bx, i, j, k,
    {

      if (lplm_limiter == 1) {
        // the 2nd order MC limiter

        Real qm1;
        Real q0;
        Real qp1;

        if (idir == 0) {
          qm1 = q_arr(i-1,j,k,n);
          q0 = q_arr(i,j,k,n);
          qp1 = q_arr(i+1,j,k,n);

        } else if (idir == 1) {
          qm1 = q_arr(i,j-1,k,n);
          q0 = q_arr(i,j,k,n);
          qp1 = q_arr(i,j+1,k,n);

        } else {
          qm1 = q_arr(i,j,k-1,n);
          q0 = q_arr(i,j,k,n);
          qp1 = q_arr(i,j,k+1,n);
        }

        Real dlft = 2.0_rt*(q0 - qm1);
        Real drgt = 2.0_rt*(qp1 - q0);
        Real dcen = 0.25_rt * (dlft + drgt);
        Real dsgn = std::copysign(1.0_rt, dcen);
        Real slop = amrex::min(std::abs(dlft), std::abs(drgt));
        Real dlim = dlft*drgt >= 0.0_rt ? slop : 0.0_rt;

        dq(i,j,k,n) = flatn_arr(i,j,k)*dsgn*amrex::min(dlim, std::abs(dcen));

      } else {
        // the 4th order MC limiter

        Real qm2;
        Real qm1;
        Real q0;
        Real qp1;
        Real qp2;

        if (idir == 0) {
          qm2 = q_arr(i-2,j,k,n);
          qm1 = q_arr(i-1,j,k,n);
          q0 = q_arr(i,j,k,n);
          qp1 = q_arr(i+1,j,k,n);
          qp2 = q_arr(i+2,j,k,n);

          // special consideration for reflecting BCs -- see
          // Saltzman p. 162 (but note that Saltzman has a
          // sign error)
          if (i == domlo[0] && n == QU && lo_bc_test) {
            qm2 = -qp1;
            qm1 = -3.0_rt*q0 + qp1 - 0.125_rt*(qp2 + qp1);
          }

          if (i == domhi[0] && n == QU && hi_bc_test) {
            qp2 = -qm1;
            qp1 = -3.0_rt*q0 + qm1 - 0.125_rt*(qm2 + qm1);
          }

        } else if (idir == 1) {
          qm2 = q_arr(i,j-2,k,n);
          qm1 = q_arr(i,j-1,k,n);
          q0 = q_arr(i,j,k,n);
          qp1 = q_arr(i,j+1,k,n);
          qp2 = q_arr(i,j+2,k,n);

          // special consideration for reflecting BCs -- see
          // Saltzmann p. 162 (but note that Saltzmann has a
          // sign error)
          if (j == domlo[1] && n == QV && lo_bc_test) {
            qm2 = -qp1;
            qm1 = -3.0_rt*q0 + qp1 - 0.125_rt*(qp2 + qp1);
          }

          if (j == domhi[1] && n == QV && hi_bc_test) {
            qp2 = -qm1;
            qp1 = -3.0_rt*q0 + qm1 - 0.125_rt*(qm2 + qm1);
          }

        } else {
          qm2 = q_arr(i,j,k-2,n);
          qm1 = q_arr(i,j,k-1,n);
          q0 = q_arr(i,j,k,n);
          qp1 = q_arr(i,j,k+1,n);
          qp2 = q_arr(i,j,k+2,n);

          // special consideration for reflecting BCs -- see
          // Saltzmann p. 162 (but note that Saltzmann has a
          // sign error)
          if (k == domlo[2] && n == QW && lo_bc_test) {
            qm2 = -qp1;
            qm1 = -3.0_rt*q0 + qp1 - 0.125_rt*(qp2 + qp1);
          }

          if (k == domhi[2] && n == QW && hi_bc_test) {
            qp2 = -qm1;
            qp1 = -3.0_rt*q0 + qm1 - 0.125_rt*(qm2 + qm1);
          }
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
        dq(i,j,k,n) = flatn_arr(i,j,k)*dsgn*amrex::min(dlim, std::abs(dq1));
      }
    });

  }
}


void
Castro::pslope(const Box& bx, const int idir,
               Array4<Real const> const q_arr,
               Array4<Real const> const flatn_arr,
               Array4<Real> const dq,
               Array4<Real const> const src) {

  const auto dx = geom.CellSizeArray();

  const auto domlo = geom.Domain().loVect3d();
  const auto domhi = geom.Domain().hiVect3d();

  const int* lo_bc = phys_bc.lo();
  const int* hi_bc = phys_bc.hi();

  bool lo_bc_test = lo_bc[idir] == Symmetry;
  bool hi_bc_test = hi_bc[idir] == Symmetry;


  if (plm_iorder == 1) {

    // first order -- piecewise constant slopes
    AMREX_PARALLEL_FOR_3D(bx, i, j, k,
    {
     dq(i,j,k,QPRES) = 0.0_rt;
    });

  } else {

    AMREX_PARALLEL_FOR_3D(bx, i, j, k,
    {

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

      if (idir == 0) {

        p0_hse = q_arr(i,j,k,QPRES);

        pp1_hse = p0_hse + 0.25_rt*dx[0] *
          (q_arr(i,j,k,QRHO) + q_arr(i+1,j,k,QRHO)) *
          (src(i,j,k,QU) + src(i+1,j,k,QU));
        pp2_hse = pp1_hse + 0.25_rt*dx[0] *
          (q_arr(i+1,j,k,QRHO) + q_arr(i+2,j,k,QRHO)) *
          (src(i+1,j,k,QU) + src(i+2,j,k,QU));

        pm1_hse = p0_hse - 0.25_rt*dx[0] *
          (q_arr(i,j,k,QRHO) + q_arr(i-1,j,k,QRHO)) *
          (src(i,j,k,QU) + src(i-1,j,k,QU));
        pm2_hse = pm1_hse - 0.25_rt*dx[0] *
          (q_arr(i-1,j,k,QRHO) + q_arr(i-2,j,k,QRHO)) *
          (src(i-1,j,k,QU) + src(i-2,j,k,QU));

        p0 = 0.0_rt;
        pp1 = q_arr(i+1,j,k,QPRES) - pp1_hse;
        pp2 = q_arr(i+2,j,k,QPRES) - pp2_hse;

        pm1 = q_arr(i-1,j,k,QPRES) - pm1_hse;
        pm2 = q_arr(i-2,j,k,QPRES) - pm2_hse;

        if (i == domlo[0] && lo_bc_test) {
          pm1 = 0.0_rt;  // HSE is perfectly satisfied
          pm2 = 0.0_rt;
        }

        if (i == domhi[0] && hi_bc_test) {
          pp1 = 0.0_rt;
          pp2 = 0.0_rt;
        }

      } else if (idir == 1) {

        p0_hse = q_arr(i,j,k,QPRES);

        pp1_hse = p0_hse + 0.25_rt*dx[1] *
          (q_arr(i,j,k,QRHO) + q_arr(i,j+1,k,QRHO)) *
          (src(i,j,k,QV) + src(i,j+1,k,QV));
        pp2_hse = pp1_hse + 0.25_rt*dx[1] *
          (q_arr(i,j+1,k,QRHO) + q_arr(i,j+2,k,QRHO)) *
          (src(i,j+1,k,QV) + src(i,j+2,k,QV));

        pm1_hse = p0_hse - 0.25_rt*dx[1] *
          (q_arr(i,j,k,QRHO) + q_arr(i,j-1,k,QRHO)) *
          (src(i,j,k,QV) + src(i,j-1,k,QV));
        pm2_hse = pm1_hse - 0.25_rt*dx[1] *
          (q_arr(i,j-1,k,QRHO) + q_arr(i,j-2,k,QRHO)) *
          (src(i,j-1,k,QV) + src(i,j-2,k,QV));

        p0 = 0.0_rt;
        pp1 = q_arr(i,j+1,k,QPRES) - pp1_hse;
        pp2 = q_arr(i,j+2,k,QPRES) - pp2_hse;

        pm1 = q_arr(i,j-1,k,QPRES) - pm1_hse;
        pm2 = q_arr(i,j-2,k,QPRES) - pm2_hse;

      } else {

        p0_hse = q_arr(i,j,k,QPRES);

        pp1_hse = p0_hse + 0.25_rt*dx[2] *
          (q_arr(i,j,k,QRHO) + q_arr(i,j,k+1,QRHO)) *
          (src(i,j,k,QW) + src(i,j,k+1,QW));
        pp2_hse = pp1_hse + 0.25_rt*dx[2] *
          (q_arr(i,j,k+1,QRHO) + q_arr(i,j,k+2,QRHO)) *
          (src(i,j,k+1,QW) + src(i,j,k+2,QW));

        pm1_hse = p0_hse - 0.25_rt*dx[2] *
          (q_arr(i,j,k,QRHO) + q_arr(i,j,k-1,QRHO)) *
          (src(i,j,k,QW) + src(i,j,k-1,QW));
        pm2_hse = pm1_hse - 0.25_rt*dx[2] *
          (q_arr(i,j,k-1,QRHO) + q_arr(i,j,k-2,QRHO)) *
          (src(i,j,k-1,QW) + src(i,j,k-2,QW));

        p0 = 0.0_rt;
        pp1 = q_arr(i,j,k+1,QPRES) - pp1_hse;
        pp2 = q_arr(i,j,k+2,QPRES) - pp2_hse;

        pm1 = q_arr(i,j,k-1,QPRES) - pm1_hse;
        pm2 = q_arr(i,j,k-2,QPRES) - pm2_hse;

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
      dq(i,j,k,QPRES) = flatn_arr(i,j,k)*dsgn*amrex::min(dlim, std::abs(dp1));

      if (idir == 0) {
        dq(i,j,k,QPRES) += q_arr(i,j,k,QRHO)*src(i,j,k,QU)*dx[0];
      } else if (idir == 1) {
        dq(i,j,k,QPRES) += q_arr(i,j,k,QRHO)*src(i,j,k,QV)*dx[1];
      } else {
        dq(i,j,k,QPRES) += q_arr(i,j,k,QRHO)*src(i,j,k,QW)*dx[2];
      }
    });
  }
}
