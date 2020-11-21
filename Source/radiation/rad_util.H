#ifndef CASTRO_RAD_UTIL_H
#define CASTRO_RAD_UTIL_H

#include <cmath>

#include <Castro_util.H>

#include <fluxlimiter.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void cell_center_metric(int i, int j, int k, const GeometryData& geomdata,
                        Real& r, Real& s)
{
    if (geomdata.Coord() == 0)
    {
        r = 1.0_rt;
        s = 1.0_rt;
    }
    else if (geomdata.Coord() == 1)
    {
        GpuArray<Real, 3> loc;
        position(i, j, k, geomdata, loc);

        r = loc[0];
        s = 1.0_rt;
    }
    else if (geomdata.Coord() == 2)
    {
       GpuArray<Real, 3> loc;
       position(i, j, k, geomdata, loc);

       Real h1 = 0.5e0_rt * geomdata.CellSize(0);
       Real d1 = 1.e0_rt / (3.e0_rt * geomdata.CellSize(0));

       r = loc[0];
       r = d1 * (std::pow(r + h1, 3) - std::pow(r - h1, 3));

#if AMREX_SPACEDIM >= 2
       Real h2 = 0.5e0_rt * geomdata.CellSize(1);
       Real d2 = 1.e0_rt / geomdata.CellSize(1);

       s = loc[1];
       s = d2 * (std::cos(s - h2) - std::cos(s + h2));
#else
       s = 1.0_rt;
#endif
    }
}



AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
compute_ptot_ctot(const Real* lam, const Real* q,
                  int is_comoving, int limiter, int closure,
                  Real cg,
                  Real& ptot, Real& ctot, Real& gamc_tot) {

  Real csrad2 = 0.0_rt;
  Real prad = 0.0_rt;

  for (int g = 0; g < NGROUPS; g++) {
    Real gamr = 0.0;
    if (is_comoving) {
      Real Eddf = Edd_factor(lam[g], limiter, closure);
      gamr = (3.0_rt - Eddf)/2.0_rt;
    } else {
      gamr = lam[g] + 1.0_rt;
    }

    prad += lam[g]*q[QRAD+g];
    csrad2 += gamr * (lam[g]*q[QRAD+g]) / q[QRHO];
  }

  ptot = q[QPRES] + prad;

  ctot = cg*cg + csrad2;
  gamc_tot = ctot * q[QRHO] / ptot;

  ctot = std::sqrt(ctot);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real FLDlambda (Real r, int limiter)
{
    Real lambda;

    if (limiter == 0) {
        // no limiter
        lambda = 1.e0_rt / 3.e0_rt;
    }
    else if (limiter < 10) {
        // approximate LP
        lambda = (2.e0_rt + r) / (6.e0_rt + r * (3.e0_rt + r));
    }
    else if (limiter < 20) {
        // Bruenn
        lambda = 1.e0_rt / (3.e0_rt + r);
    }
    else if (limiter < 30) {
        // Larsen's square root
        lambda = 1.e0_rt / std::sqrt(9.e0_rt + r * r);
    }
    else if (limiter < 40) {
        // Minerbo
        if (r < 1.5e0_rt) {
            lambda = 2.e0_rt / (3.e0_rt + std::sqrt(9.e0_rt + 12.e0_rt * r * r));
        } else {
            lambda = 1.e0_rt / (1.e0_rt + r + std::sqrt(1.e0_rt + 2.e0_rt * r));
        }
    }
    else {
#ifndef AMREX_USE_GPU
        amrex::Error("Unknown limiter type");
#endif
    }

    return lambda;
}

#endif
