#ifndef advection_util_H
#define advection_util_H

#include "Castro_util.H"

#ifdef HYBRID_MOMENTUM
#include "hybrid.H"
#endif

AMREX_GPU_HOST_DEVICE
void
dflux(const GpuArray<Real, NUM_STATE>& u,
      const GpuArray<Real, NQ>& q,
      int dir, int coord,
      const GeometryData& geomdata,
      const GpuArray<Real, 3>& center,
      const GpuArray<int, 3>& idx,
      GpuArray<Real, NUM_STATE>& flux)
{
    // Given a conservative state and its corresponding primitive state, calculate the
    // corresponding flux in a given direction.

    // Set everything to zero; this default matters because some
    // quantities like temperature are not updated through fluxes.

    for (int n = 0; n < NUM_STATE; ++n) {
        flux[n] = 0.0_rt;
    }

    // Determine the advection speed based on the flux direction.

    Real v_adv = q[QU + dir];

    // Core quantities (density, momentum, energy).

    flux[URHO] = u[URHO] * v_adv;
    flux[UMX] = u[UMX] * v_adv;
    flux[UMY] = u[UMY] * v_adv;
    flux[UMZ] = u[UMZ] * v_adv;
    flux[UEDEN] = (u[UEDEN] + q[QPRES]) * v_adv;
    flux[UEINT] = u[UEINT] * v_adv;

    // Optionally include the pressure term in the momentum flux.
    // It is optional because for some geometries we cannot write
    // the pressure term in a conservative form.

    if (mom_flux_has_p(dir, dir, coord)) {
        flux[UMX + dir] = flux[UMX + dir] + q[QPRES];
    }

    // Hybrid flux.

#ifdef HYBRID_MOMENTUM
    // Create a temporary edge-based q for this routine.
    GpuArray<Real, NGDNV> qgdnv;
    for (int n = 0; n < NGDNV; ++n) {
        qgdnv[n] = 0.0_rt;
    }
    qgdnv[GDRHO] = q[QRHO];
    qgdnv[GDU] = q[QU];
    qgdnv[GDV] = q[QV];
    qgdnv[GDW] = q[QW];
    qgdnv[GDPRES] = q[QPRES];
    bool cell_centered = true;
    compute_hybrid_flux(qgdnv, geomdata, center, dir,
                        idx[0], idx[1], idx[2],
                        flux, cell_centered);
#endif

    // Passively advected quantities.

    for (int ipassive = 0; ipassive < npassive; ++ipassive) {

        int n = upassmap(ipassive);
        flux[n] = u[n] * v_adv;

    }

}

#endif
