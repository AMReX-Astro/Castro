#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

using namespace amrex;


void
Castro::src_to_prim(const Box& bx,
                    Array4<Real const> const q,
                    Array4<Real const> const qaux,
                    Array4<Real const> const src,
                    Array4<Real> const srcQ)
{

  AMREX_PARALLEL_FOR_3D(bx, i, j, k,
  {


      for (int n = 0; n < NQSRC; ++n) {
        srcQ(i,j,k,n) = 0.0;
      }

      Real rhoinv = 1.0 / q(i,j,k,QRHO);

      srcQ(i,j,k,QRHO) = src(i,j,k,URHO);
      srcQ(i,j,k,QU) = (src(i,j,k,UMX) - q(i,j,k,QU) * srcQ(i,j,k,QRHO)) * rhoinv;
      srcQ(i,j,k,QV) = (src(i,j,k,UMY) - q(i,j,k,QV) * srcQ(i,j,k,QRHO)) * rhoinv;
      srcQ(i,j,k,QW) = (src(i,j,k,UMZ) - q(i,j,k,QW) * srcQ(i,j,k,QRHO)) * rhoinv;
      srcQ(i,j,k,QREINT) = src(i,j,k,UEINT);
      srcQ(i,j,k,QPRES ) = qaux(i,j,k,QDPDE) *
        (srcQ(i,j,k,QREINT) - q(i,j,k,QREINT)*srcQ(i,j,k,QRHO)*rhoinv) *
        rhoinv + qaux(i,j,k,QDPDR)*srcQ(i,j,k,QRHO);

#ifdef PRIM_SPECIES_HAVE_SOURCES
      for (int ipassive = 0; ipassive < npassive; ++ipassive) {
        int n = upass_map[ipassive];
        int iq = qpass_map[ipassive];

       // we may not be including the ability to have species sources,
       //  so check to make sure that we are < NQSRC
        srcQ(i,j,k,iq) = (src(i,j,k,n) - q(i,j,k,iq) * srcQ(i,j,k,QRHO) ) /
          q(i,j,k,QRHO);
      }
#endif
  });

}


void
Castro::shock(const Box& bx,
              Array4<Real const> const q,
              Array4<Real> const shk) {

  // This is a basic multi-dimensional shock detection algorithm.
  // This implementation follows Flash, which in turn follows
  // AMRA and a Woodward (1995) (supposedly -- couldn't locate that).
  //
  // The spirit of this follows the shock detection in Colella &
  // Woodward (1984)
  //

  constexpr Real small = 1.e-10;
  constexpr Real eps = 0.33e0;

  const auto dx = geom.CellSizeArray();
  const int coord_type = geom.Coord();

  Real dxinv = 1.0/dx[0];
  Real dyinv = 1.0/dx[1];
  Real dzinv = 1.0/dx[2];

  AMREX_PARALLEL_FOR_3D(bx, i, j, k,
  {
    Real div_u = 0;
    // construct div{U}
    if (coord_type == 0) {

      // Cartesian
      div_u += 0.5*(q(i+1,j,k,QU) - q(i-1,j,k,QU))*dxinv;
#if (AMREX_SPACEDIM >= 2)
      div_u += 0.5*(q(i,j+1,k,QV) - q(i,j-1,k,QV))*dyinv;
#endif
#if (AMREX_SPACEDIM == 3)
     div_u += 0.5*(q(i,j,k+1,QW) - q(i,j,k-1,QW))*dzinv;
#endif

   } else if (coord_type == 1) {

     // r-z
     Real rc = (i + 0.5)*dx[0];
     Real rm = (i - 1 + 0.5)*dx[0];
     Real rp = (i + 1 + 0.5)*dx[0];

#if (AMREX_SPACEDIM == 1)
     div_u += 0.5*(rp*q(i+1,j,k,QU) - rm*q(i-1,j,k,QU))/(rc*dx[0]);
#endif
#if (AMREX_SPACEDIM == 2)
     div_u += 0.5*(rp*q(i+1,j,k,QU) - rm*q(i-1,j,k,QU))/(rc*dx[0]) +
              0.5*(q(i,j+1,k,QV) - q(i,j-1,k,QV)) * dyinv;
#endif

    } else if (coord_type == 2) {

      // 1-d spherical
      Real rc = (i + 0.5)*dx[0];
      Real rm = (i - 1 + 0.5)*dx[0];
      Real rp = (i + 1 + 0.5)*dx[0];

      div_u += 0.5*(rp*rp*q(i+1,j,k,QU) - rm*rm*q(i-1,j,k,QU))/(rc*rc*dx[0]);

#ifndef AMREX_USE_CUDA

    } else {
      amrex::Error("ERROR: invalid coord_type in shock");
#endif
    }

    // find the pre- and post-shock pressures in each direction
    Real px_pre;
    Real px_post;
    Real e_x;

    if (q(i+1,j,k,QPRES) - q(i-1,j,k,QPRES) < 0.0) {
      px_pre = q(i+1,j,k,QPRES);
      px_post = q(i-1,j,k,QPRES);
    } else {
      px_pre = q(i-1,j,k,QPRES);
      px_post = q(i+1,j,k,QPRES);
    }

    // use compression to create unit vectors for the shock direction
    e_x = std::pow(q(i+1,j,k,QU) - q(i-1,j,k,QU), 2);

    Real py_pre;
    Real py_post;
    Real e_y;

#if (AMREX_SPACEDIM >= 2)
    if (q(i,j+1,k,QPRES) - q(i,j-1,k,QPRES) < 0.0) {
      py_pre = q(i,j+1,k,QPRES);
      py_post = q(i,j-1,k,QPRES);
    } else {
      py_pre = q(i,j-1,k,QPRES);
      py_post = q(i,j+1,k,QPRES);
    }

    e_y = std::pow(q(i,j+1,k,QV) - q(i,j-1,k,QV), 2);

#else
    py_pre = 0.0;
    py_post = 0.0;

    e_y = 0.0;
#endif

    Real pz_pre;
    Real pz_post;
    Real e_z;

#if (AMREX_SPACEDIM == 3)
    if (q(i,j,k+1,QPRES) - q(i,j,k-1,QPRES) < 0.0) {
      pz_pre  = q(i,j,k+1,QPRES);
      pz_post = q(i,j,k-1,QPRES);
    } else {
      pz_pre  = q(i,j,k-1,QPRES);
      pz_post = q(i,j,k+1,QPRES);
    }

    e_z = std::pow(q(i,j,k+1,QW) - q(i,j,k-1,QW), 2);

#else
    pz_pre = 0.0;
    pz_post = 0.0;

    e_z = 0.0;
#endif

    Real d = 1.0/(e_x + e_y + e_z + small);

    e_x = e_x*d;
    e_y = e_y*d;
    e_z = e_z*d;

    // project the pressures onto the shock direction
    Real p_pre  = e_x*px_pre + e_y*py_pre + e_z*pz_pre;
    Real p_post = e_x*px_post + e_y*py_post + e_z*pz_post;

    // test for compression + pressure jump to flag a shock
    // this avoid U = 0, so e_x, ... = 0
    Real pjump = p_pre == 0 ? 0.0 : eps - (p_post - p_pre)/p_pre;

    if (pjump < 0.0 && div_u < 0.0) {
      shk(i,j,k) = 1.0;
    } else {
      shk(i,j,k) = 0.0;
    }
  });

}


void
Castro::divu(const Box& bx,
             Array4<Real const> const q,
             Array4<Real> const div) {
  // this computes the *node-centered* divergence

  const auto dx = geom.CellSizeArray();
  const int coord_type = geom.Coord();

  const auto problo = geom.ProbLoArray();

  Real dxinv = 1.0/dx[0];
#if AMREX_SPACEDIM >= 2
  Real dyinv = 1.0/dx[1];
#else
  Real dyinv = 0.0;
#endif
#if AMREX_SPACEDIM == 3
  Real dzinv = 1.0/dx[2];
#else
  Real dzinv = 0.0;
#endif

  AMREX_PARALLEL_FOR_3D(bx, i, j, k,
  {

#if AMREX_SPACEDIM == 1
    if (coord_type == 0) {
      div(i,j,k) = (q(i,j,k,QU) - q(i-1,j,k,QU)) * dxinv;

    } else if (coord_type == 1) {
      // axisymmetric
      if (i == 0) {
        div(i,j,k) = 0.0;
      } else {
        Real rl = (i - 0.5) * dx[0] + problo[0];
        Real rr = (i + 0.5) * dx[0] + problo[0];
        Real rc = (i) * dx[0] + problo[0];

        div(i,j,k) = (rr*q(i,j,k,QU) - rl*q(i-1,j,k,QU)) * dxinv / rc;
      }
    } else {
      // spherical
      if (i == 0) {
        div(i,j,k) = 0.0;
      } else {
        Real rl = (i - 0.5) * dx[0] + problo[0];
        Real rr = (i + 0.5) * dx[0] + problo[0];
        Real rc = (i) * dx[0] + problo[0];

        div(i,j,k) = (rr*rr*q(i,j,k,QU) - rl*rl*q(i-1,j,k,QU)) * dxinv / (rc*rc);
      }
    }
#endif

#if AMREX_SPACEDIM == 2
    Real ux = 0.0;
    Real vy = 0.0;

    if (coord_type == 0) {
      ux = 0.5*(q(i,j,k,QU) - q(i-1,j,k,QU) + q(i,j-1,k,QU) - q(i-1,j-1,k,QU)) * dxinv;
      vy = 0.5*(q(i,j,k,QV) - q(i,j-1,k,QV) + q(i-1,j,k,QV) - q(i-1,j-1,k,QV)) * dyinv;

    } else {
      if (i == 0) {
        ux = 0.0;
        vy = 0.0;  // is this part correct?
      } else {
        Real rl = (i - 0.5) * dx[0] + problo[0];
        Real rr = (i + 0.5) * dx[0] + problo[0];
        Real rc = (i) * dx[0] + problo[0];

        // These are transverse averages in the y-direction
        Real ul = 0.5 * (q(i-1,j,k,QU) + q(i-1,j-1,k,QU));
        Real ur = 0.5 * (q(i,j,k,QU) + q(i,j-1,k,QU));

        // Take 1/r d/dr(r*u)
        ux = (rr*ur - rl*ul) * dxinv / rc;

        // These are transverse averages in the x-direction
        Real vb = 0.5 * (q(i,j-1,k,QV) + q(i-1,j-1,k,QV));
        Real vt = 0.5 * (q(i,j,k,QV) + q(i-1,j,k,QV));

        vy = (vt - vb) * dyinv;
      }
    }

    div(i,j,k) = ux + vy;
#endif

#if AMREX_SPACEDIM == 3
    Real ux = 0.25 * (q(i,j,k,QU) - q(i-1,j,k,QU) +
                      q(i,j,k-1,QU) - q(i-1,j,k-1,QU) +
                      q(i,j-1,k,QU) - q(i-1,j-1,k,QU) +
                      q(i,j-1,k-1,QU) - q(i-1,j-1,k-1,QU)) * dxinv;

    Real vy = 0.25 * (q(i,j,k,QV) - q(i,j-1,k,QV) +
                      q(i,j,k-1,QV) - q(i,j-1,k-1,QV) +
                      q(i-1,j,k,QV) - q(i-1,j-1,k,QV) +
                      q(i-1,j,k-1,QV) - q(i-1,j-1,k-1,QV)) * dyinv;

    Real wz = 0.25 * (q(i,j,k,QW) - q(i,j,k-1,QW) +
                      q(i,j-1,k,QW) - q(i,j-1,k-1,QW) +
                      q(i-1,j,k,QW) - q(i-1,j,k-1,QW) +
                      q(i-1,j-1,k,QW) - q(i-1,j-1,k-1,QW)) * dzinv;

    div(i,j,k) = ux + vy + wz;
#endif

  });

}


void
Castro::apply_av(const Box& bx,
                 const int idir,
                 Array4<Real const> const div,
                 Array4<Real const> const uin,
                 Array4<Real> const flux) {

  const auto dx = geom.CellSizeArray();

  AMREX_PARALLEL_FOR_4D(bx, NUM_STATE, i, j, k, n,
  {

    if (n == UTEMP) continue;
#ifdef SHOCK_VAR
    if (n == USHK) continue;
#endif

    Real div1;
    if (idir == 0) {

      div1 = 0.25 * (div(i,j,k) + div(i,j+dg[1],k) +
                     div(i,j,k+dg[2]) + div(i,j+dg[1],k+dg[2]));
      div1 = difmag * std::min(0.0, div1);
      div1 = div1 * (uin(i,j,k,n) - uin(i-1,j,k,n));

    } else if (idir == 1) {

      div1 = 0.25 * (div(i,j,k) + div(i+1,j,k) +
                     div(i,j,k+dg[2]) + div(i+1,j,k+dg[2]));
      div1 = difmag * std::min(0.0, div1);
      div1 = div1 * (uin(i,j,k,n) - uin(i,j-dg[1],k,n));

    } else {

      div1 = 0.25 * (div(i,j,k) + div(i+1,j,k) +
                     div(i,j+dg[1],k) + div(i+1,j+dg[1],k));
      div1 = difmag * std::min(0.0, div1);
      div1 = div1 * (uin(i,j,k,n) - uin(i,j,k-dg[2],n));

    }

    flux(i,j,k,n) += dx[idir] * div1;
  });
}


#ifdef RADIATION
void
Castro::apply_av_rad(const Box& bx,
                     const int idir,
                     Array4<Real const> const div,
                     Array4<Real const> const Erin,
                     Array4<Real> const radflux) {

  const auto dx = geom.CellSizeArray();

  AMREX_PARALLEL_FOR_4D(bx, Radiation::nGroups, i, j, k, n,
  {

    Real div1;
    if (idir == 0) {

      div1 = 0.25 * (div(i,j,k) + div(i,j+dg[1],k) +
                     div(i,j,k+dg[2]) + div(i,j+dg[1],k+dg[2]));
      div1 = difmag * std::min(0.0, div1);
      div1 = div1 * (Erin(i,j,k,n) - Erin(i-1,j,k,n));

    } else if (idir == 1) {

      div1 = 0.25 * (div(i,j,k) + div(i+1,j,k) +
                     div(i,j,k+dg[2]) + div(i+1,j,k+dg[2]));
      div1 = difmag * std::min(0.0, div1);
      div1 = div1 * (Erin(i,j,k,n) - Erin(i,j-dg[1],k,n));

    } else {

      div1 = 0.25 * (div(i,j,k) + div(i+1,j,k) +
                     div(i,j+dg[1],k) + div(i+1,j+dg[1],k));
      div1 = difmag * std::min(0.0, div1);
      div1 = div1 * (Erin(i,j,k,n) - Erin(i,j,k-dg[2],n));

    }

    radflux(i,j,k,n) += dx[idir] * div1;
  });
}
#endif
