#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

#ifdef RADIATION
#include "Radiation.H"
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
  Real lconst_grav = gravity->const_grav;
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

    AMREX_PARALLEL_FOR_3D(bx, i, j, k,
    {

      if (plm_well_balanced == 1 && n == QPRES && idir == AMREX_SPACEDIM-1) {
        // we'll only do a second-order pressure slope,
        // but we'll follow the well-balanced scheme of
        // Kappeli.  Note at the moment we are assuming
        // constant gravity.

        Real p0;
        Real pp1;
        Real pm1;

        if (idir == 0) {
          p0 = 0.0_rt;
          pp1 = q_arr(i+1,j,k,QPRES) - (p0 + 0.5_rt*dx[0]*(q_arr(i,j,k,QRHO) + q_arr(i+1,j,k,QRHO))*lconst_grav);
          pm1 = q_arr(i-1,j,k,QPRES) - (p0 - 0.5_rt*dx[0]*(q_arr(i,j,k,QRHO) + q_arr(i-1,j,k,QRHO))*lconst_grav);

          if (i == domlo[0] && lo_bc_test) {
            pm1 = 0.0_rt;  // HSE is perfectly satisfied
          }

          if (i == domhi[0] && hi_bc_test) {
            pp1 = 0.0_rt;
          }

        } else if (idir == 1) {
          p0 = 0.0_rt;
          pp1 = q_arr(i,j+1,k,QPRES) - (p0 + 0.5_rt*dx[1]*(q_arr(i,j,k,QRHO) + q_arr(i,j+1,k,QRHO))*lconst_grav);
          pm1 = q_arr(i,j-1,k,QPRES) - (p0 - 0.5_rt*dx[1]*(q_arr(i,j,k,QRHO) + q_arr(i,j-1,k,QRHO))*lconst_grav);

          if (j == domlo[1] && lo_bc_test) {
            pm1 = 0.0_rt;  // HSE is perfectly satisfied
          }

          if (j == domhi[1] && hi_bc_test) {
            pp1 = 0.0_rt;
          }

        } else {
          p0 = 0.0_rt;
          pp1 = q_arr(i,j,k+1,QPRES) - (p0 + 0.5_rt*dx[2]*(q_arr(i,j,k,QRHO) + q_arr(i,j,k+1,QRHO))*lconst_grav);
          pm1 = q_arr(i,j,k-1,QPRES) - (p0 - 0.5_rt*dx[2]*(q_arr(i,j,k,QRHO) + q_arr(i,j,k-1,QRHO))*lconst_grav);

          if (k == domlo[2] && lo_bc_test) {
            pm1 = 0.0_rt;  // HSE is perfectly satisfied
          }

          if (k == domhi[2] && hi_bc_test) {
            pp1 = 0.0_rt;
          }
        }

        Real dlft = 2.0_rt*(p0 - pm1);
        Real drgt = 2.0_rt*(pp1 - p0);
        Real dcen = 0.25_rt * (dlft + drgt);
        Real dsgn = std::copysign(1.0_rt, dcen);
        Real slop = amrex::min(std::abs(dlft), std::abs(drgt));
        Real dlim = dlft*drgt >= 0.0_rt ? slop : 0.0_rt;

        dq(i,j,k,n) = flatn_arr(i,j,k)*dsgn*amrex::min(dlim, std::abs(dcen));

      } else if (plm_limiter == 1) {
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
          // Saltzmann p. 162 (but note that Saltzmann has a
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

  if (plm_iorder == 1) {

    // first order -- piecewise constant slopes
    AMREX_PARALLEL_FOR_3D(bx, i, j, k,
    {
     dq(i,j,k,QPRES) = 0.0_rt;
    });

  } else {

    AMREX_PARALLEL_FOR_3D(bx, i, j, k,
    {

      // First compute Fromm slopes

      // df at i+1
      Real dlftp1;
      Real drgtp1;

      if (idir == 0) {
        dlftp1 = q_arr(i+1,j,k,QPRES) - q_arr(i,j,k,QPRES);
        drgtp1 = q_arr(i+2,j,k,QPRES) - q_arr(i+1,j,k,QPRES);

        // Subtract off (rho * acceleration) so as not to limit that part of the slope
        dlftp1 += -0.25_rt * (q_arr(i+1,j,k,QRHO) + q_arr(i,j,k,QRHO)) *
          (src(i+1,j,k,QU)+src(i,j,k,QU))*dx[0];
        drgtp1 += -0.25_rt * (q_arr(i+1,j,k,QRHO) + q_arr(i+2,j,k,QRHO)) *
          (src(i+1,j,k,QU)+src(i+2,j,k,QU))*dx[0];

      } else if (idir == 1) {
        dlftp1 = q_arr(i,j+1,k,QPRES) - q_arr(i,j,k,QPRES);
        drgtp1 = q_arr(i,j+2,k,QPRES) - q_arr(i,j+1,k,QPRES);

        // Subtract off (rho * acceleration) so as not to limit that part of the slope
        dlftp1 += -0.25_rt * (q_arr(i,j+1,k,QRHO) + q_arr(i,j,k,QRHO)) *
          (src(i,j+1,k,QV)+src(i,j,k,QV))*dx[1];
        drgtp1 += -0.25_rt * (q_arr(i,j+1,k,QRHO) + q_arr(i,j+2,k,QRHO)) *
          (src(i,j+1,k,QV)+src(i,j+2,k,QV))*dx[1];

      } else {
        dlftp1 = q_arr(i,j,k+1,QPRES) - q_arr(i,j,k,QPRES);
        drgtp1 = q_arr(i,j,k+2,QPRES) - q_arr(i,j,k+1,QPRES);

        // Subtract off (rho * acceleration) so as not to limit that part of the slope
        dlftp1 += -0.25_rt * (q_arr(i,j,k+1,QRHO) + q_arr(i,j,k,QRHO)) *
          (src(i,j,k+1,QW)+src(i,j,k,QW))*dx[2];
        drgtp1 += -0.25_rt * (q_arr(i,j,k+1,QRHO) + q_arr(i,j,k+2,QRHO)) *
          (src(i,j,k+1,QW)+src(i,j,k+2,QW))*dx[2];

      }


      Real dcen = 0.5_rt*(dlftp1 + drgtp1);
      Real dsgn = std::copysign(1.0_rt, dcen);
      Real dlim = dlftp1*drgtp1 >= 0.0_rt ? 2.0_rt * std::min(std::abs(dlftp1), std::abs(drgtp1)) : 0.0_rt;
      Real dfp1 = dsgn*amrex::min(dlim, std::abs(dcen));

      // df at i-1
      Real dlftm1;
      Real drgtm1;

      if (idir == 0) {
        dlftm1 = q_arr(i-1,j,k,QPRES) - q_arr(i-2,j,k,QPRES);
        drgtm1 = q_arr(i,j,k,QPRES) - q_arr(i-1,j,k,QPRES);

        // Subtract off (rho * acceleration) so as not to limit that part of the slope
        dlftm1 += -0.25_rt * (q_arr(i-1,j,k,QRHO) + q_arr(i-2,j,k,QRHO)) *
          (src(i-1,j,k,QU)+src(i-2,j,k,QU))*dx[0];
        drgtm1 += -0.25_rt * (q_arr(i-1,j,k,QRHO) + q_arr(i,j,k,QRHO)) *
          (src(i-1,j,k,QU)+src(i,j,k,QU))*dx[0];

      } else if (idir == 1) {
        dlftm1 = q_arr(i,j-1,k,QPRES) - q_arr(i,j-2,k,QPRES);
        drgtm1 = q_arr(i,j,k,QPRES) - q_arr(i,j-1,k,QPRES);

        // Subtract off (rho * acceleration) so as not to limit that part of the slope
        dlftm1 += -0.25_rt * (q_arr(i,j-1,k,QRHO) + q_arr(i,j-2,k,QRHO)) *
          (src(i,j-1,k,QV)+src(i,j-2,k,QV))*dx[1];
        drgtm1 += -0.25_rt * (q_arr(i,j-1,k,QRHO) + q_arr(i,j,k,QRHO)) *
          (src(i,j-1,k,QV)+src(i,j,k,QV))*dx[1];

      } else {
        dlftm1 = q_arr(i,j,k-1,QPRES) - q_arr(i,j,k-2,QPRES);
        drgtm1 = q_arr(i,j,k,QPRES) - q_arr(i,j,k-1,QPRES);

        // Subtract off (rho * acceleration) so as not to limit that part of the slope
        dlftm1 += -0.25_rt * (q_arr(i,j,k-1,QRHO) + q_arr(i,j,k-2,QRHO)) *
          (src(i,j,k-1,QW)+src(i,j,k-2,QW))*dx[2];
        drgtm1 += -0.25_rt * (q_arr(i,j,k-1,QRHO) + q_arr(i,j,k,QRHO)) *
          (src(i,j,k-1,QW)+src(i,j,k,QW))*dx[2];

      }

      dcen = 0.5_rt*(dlftm1 + drgtm1);
      dsgn = std::copysign(1.0_rt, dcen);
      dlim = dlftm1*drgtm1 >= 0.0_rt ? 2.0_rt * amrex::min(std::abs(dlftm1), std::abs(drgtm1)) : 0.0_rt;
      Real dfm1 = dsgn*amrex::min(dlim, std::abs(dcen));

      // Now limited fourth order slopes at i
      Real dlft = drgtm1;
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
