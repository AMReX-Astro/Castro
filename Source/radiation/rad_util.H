#ifndef CASTRO_RAD_UTIL_H
#define CASTRO_RAD_UTIL_H

#include <cmath>

#include <Castro_util.H>
#include <radiation_params.H>
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



AMREX_GPU_HOST_DEVICE AMREX_INLINE
void edge_center_metric (int i, int j, int k, int idir, const GeometryData& geomdata,
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

        if (idir == 0)
        {
            bool ccx = false;
            position(i, j, k, geomdata, loc, ccx);
        }
        else
        {
            position(i, j, k, geomdata, loc);
        }

        r = loc[0];
        s = 1.0_rt;
    }
    else if (geomdata.Coord() == 2)
    {
        GpuArray<Real, 3> loc;

        if (idir == 0)
        {
            bool ccx = false;
            position(i, j, k, geomdata, loc, ccx);

            r = loc[0] * loc[0];

#if AMREX_SPACEDIM >= 2
            position(i, j, k, geomdata, loc);
            s = loc[1];

            Real h2 = 0.5e0_rt * geomdata.CellSize(1);
            Real d2 = 1.e0_rt / geomdata.CellSize(1);

            s = d2 * (std::cos(s - h2) - std::cos(s + h2));
#else
            s = 1.0_rt;
#endif
        }
        else
        {
            GpuArray<Real, 3> loc;
            position(i, j, k, geomdata, loc);

            r = loc[0];

            if (AMREX_SPACEDIM >= 2)
            {
                bool ccx = true;
                bool ccy = false;
                position(i, j, k, geomdata, loc, ccx, ccy);
            }
            else
            {
                bool ccx = false;
                position(i, j, k, geomdata, loc, ccx);
            }

            s = loc[1];

            Real h1 = 0.5e0_rt * geomdata.CellSize(0);
            Real d1 = 1.e0_rt / (3.e0_rt * geomdata.CellSize(0));

            r = d1 * (std::pow(r + h1, 3) - std::pow(r - h1, 3));
            s = std::sin(s);
        }
    }
}



AMREX_GPU_HOST_DEVICE AMREX_INLINE
void face_metric (int i, int j, int k, int lo, int hi,
                  const GeometryData& geomdata, int ori_dir, int ori_lo, Real& r)
{
    const Real* problo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();

    if (ori_dir == 0) {

        if (geomdata.Coord() == 0) {

            r = 1.0_rt;

        }
        else {

            if (ori_lo == 1) {

                r = problo[0] + lo * dx[0];

            }
            else {

                r = problo[0] + (hi + 1) * dx[0];

            }

            if (geomdata.Coord() == 2) {

                r = r * r;

            }

        }

    }
    else {

        if (geomdata.Coord() == 0) {

            r = 1.0_rt;

        }
        else {

            // This must be RZ since we only support spherical coordinates in 1D.

            GpuArray<Real, 3> loc;
            position(i, j, k, geomdata, loc);

            r = loc[0];

        }

    }
}



AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
compute_ptot_ctot(const Real* lam, const Real* q,
                  Real cg,
                  Real& ptot, Real& ctot, Real& gamc_tot) {

  Real csrad2 = 0.0_rt;
  Real prad = 0.0_rt;

  for (int g = 0; g < NGROUPS; g++) {
    Real gamr = 0.0;
    if (radiation::comoving) {
      Real Eddf = Edd_factor(lam[g]);
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
Real FLDlambda (Real r)
{
    Real lambda;

    if (radiation::limiter == 0) {
        // no limiter
        lambda = 1.e0_rt / 3.e0_rt;
    }
    else if (radiation::limiter < 10) {
        // approximate LP
        lambda = (2.e0_rt + r) / (6.e0_rt + r * (3.e0_rt + r));
    }
    else if (radiation::limiter < 20) {
        // Bruenn
        lambda = 1.e0_rt / (3.e0_rt + r);
    }
    else if (radiation::limiter < 30) {
        // Larsen's square root
        lambda = 1.e0_rt / std::sqrt(9.e0_rt + r * r);
    }
    else if (radiation::limiter < 40) {
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

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real kavg(Real a, Real b, Real d, int opt)
{
#ifndef AMREX_USE_GPU
    if (opt > 2) {
        amrex::Error("kavg: invalid averaging option");
    }
#endif

    Real k;
    const Real tiny = 1.e-50_rt;

    if (opt == 0)
    {
        // arithmetic average, geometrically correct(?) but underestimates surface flux
        k = 0.5e0_rt * (a + b + tiny);

    }
    else if (opt == 1)
    {
        // harmonic average, overestimates surface flux
        k = (2.e0_rt * a * b) / (a + b + tiny) + tiny;
    }
    else
    {
        // chooses arithmetic where optically thin, harmonic where optically thick,
        // surface flux approximation at a thick/thin boundary
        k = amrex::min(0.5e0_rt * (a + b + tiny),
                       amrex::max((2.e0_rt * a * b) / (a + b + tiny) + tiny,
                                  4.e0_rt / (3.e0_rt * d)));
    }

    return k;
}

AMREX_INLINE
void rfface (Array4<Real> const fine,
             Array4<Real const> const crse,
             int idim, const IntVect& irat)
{
    const Dim3 flo = amrex::lbound(fine);
    const Dim3 fhi = amrex::ubound(fine);

    const Dim3 clo = amrex::lbound(crse);
    const Dim3 chi = amrex::ubound(crse);

    Real ifac = 1.0_rt;
    Real jfac = 1.0_rt;
    Real kfac = 1.0_rt;
    Real rfac = 1.0_rt;

    if (idim == 0) {
#if AMREX_SPACEDIM >= 2
        jfac = irat[1];
#endif
#if AMREX_SPACEDIM == 3
        kfac = irat[2];
#endif

        Real rfac = jfac * kfac;

        int i = flo.x;
        for (int k = flo.z; k <= fhi.z; ++k) {
            for (int j = flo.y; j <= fhi.y; ++j) {
                fine(i,j,k) = crse(clo.x, j / jfac, k / kfac) / rfac;
            }
        }
    }
    else if (idim == 1) {
#if AMREX_SPACEDIM >= 2
        ifac = irat[0];
#endif
#if AMREX_SPACEDIM == 3
        kfac = irat[2];
#endif

        rfac = ifac * kfac;

        int j = flo.y;
        for (int k = flo.z; k <= fhi.z; ++k) {
            for (int i = flo.x; i <= fhi.x; ++i) {
                fine(i,j,k) = crse(i / ifac, clo.y, k / kfac) / rfac;
            }
        }
    }
    else {
#if AMREX_SPACEDIM >= 2
        ifac = irat[0];
#endif
#if AMREX_SPACEDIM == 3
        jfac = irat[1];
#endif

        rfac = ifac * jfac;

        int k = flo.z;
        for (int j = flo.y; j <= fhi.y; ++j) {
            for (int i = flo.x; i <= fhi.x; ++i) {
                fine(i,j,k) = crse(i / ifac, j / jfac, clo.z) / rfac;
            }
        }
    }
}

#endif
