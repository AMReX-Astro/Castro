#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

using namespace amrex;

void
Castro::consup_hydro(const Box& bx,
                     Array4<Real const> const shk,
                     Array4<Real> const update,
                     Array4<Real> const flux0,
                     Array4<Real const> const qx,
                     Array4<Real const> const area0,
#if AMREX_SPACEDIM >= 2
                     Array4<Real> const flux1,
                     Array4<Real const> const qy,
                     Array4<Real const> const area1,
#endif
#if AMREX_SPACEDIM == 3
                     Array4<Real> const flux2,
                     Array4<Real const> const qz,
                     Array4<Real const> const area2,
#endif
                     Array4<Real const> const vol,
                     const Real dt)
{


  const auto dx = geom.CellSizeArray();


  // For hydro, we will create an update source term that is
  // essentially the flux divergence.  This can be added with dt to
  // get the update

  const int flux_has_p = momx_flux_has_p[0];

  GpuArray<Real, 3> center;
  ca_get_center(center.begin());

  AMREX_PARALLEL_FOR_4D(bx, NUM_STATE, i, j, k, n,
  {

    Real volinv = 1.0 / vol(i,j,k);

    update(i,j,k,n) = update(i,j,k,n) +
      ( flux0(i,j,k,n) * area0(i,j,k) - flux0(i+1,j,k,n) * area0(i+1,j,k)
#if AMREX_SPACEDIM >= 2
      + flux1(i,j,k,n) * area1(i,j,k) - flux1(i,j+1,k,n) * area1(i,j+1,k)
#endif
#if AMREX_SPACEDIM == 3
      + flux2(i,j,k,n) * area2(i,j,k) - flux2(i,j,k+1,n) * area2(i,j,k+1)
#endif
        ) * volinv;

   // Add the p div(u) source term to (rho e).
    if (n == UEINT) {

      Real pdu = (qx(i+1,j,k,GDPRES) + qx(i,j,k,GDPRES)) *
                 (qx(i+1,j,k,GDU) * area0(i+1,j,k) - qx(i,j,k,GDU) * area0(i,j,k));

#if AMREX_SPACEDIM >= 2
      pdu += (qy(i,j+1,k,GDPRES) + qy(i,j,k,GDPRES)) *
             (qy(i,j+1,k,GDV) * area1(i,j+1,k) - qy(i,j,k,GDV) * area1(i,j,k));
#endif

#if AMREX_SPACEDIM == 3
      pdu += (qz(i,j,k+1,GDPRES) + qz(i,j,k,GDPRES)) *
             (qz(i,j,k+1,GDW) * area2(i,j,k+1) - qz(i,j,k,GDW) * area2(i,j,k));
#endif

      pdu = 0.5 * pdu * volinv;

      update(i,j,k,n) = update(i,j,k,n) - pdu;

#ifdef SHOCK_VAR
    } else if (n == USHK) {
      update(i,j,k,USHK) = shk(i,j,k,1) / dt;
#endif

#ifndef RADIATION
    } else if (n == UMX) {
      // Add gradp term to momentum equation -- only for axisymmetric
      // coords (and only for the radial flux).

      if (! flux_has_p) {
        update(i,j,k,UMX) += - (qx(i+1,j,k,GDPRES) - qx(i,j,k,GDPRES)) / dx[0];
      }
#endif

#ifdef HYBRID_MOMENTUM
    } else if (n == UMR) {
      // update the radial momentum with the hybrid advection source
      add_hybrid_advection_source_c(i, j, k,
                                    dt, dx,
                                    center,
                                    update,
                                    qx, qy, qz);
#endif

    }



  });

}
