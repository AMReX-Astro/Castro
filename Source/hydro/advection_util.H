#ifndef advection_util_H
#define advection_util_H

#include <Castro_util.H>

#ifdef HYBRID_MOMENTUM
#include <hybrid.H>
#endif

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
dflux(const GpuArray<Real, NUM_STATE>& u,
      const GpuArray<Real, NQ>& q,
      int dir, int coord,
      const GeometryData& geomdata,
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
    compute_hybrid_flux(qgdnv, geomdata, dir,
                        idx[0], idx[1], idx[2],
                        flux, cell_centered);
#endif

    // Passively advected quantities.

    for (int ipassive = 0; ipassive < npassive; ++ipassive) {

        int n = upassmap(ipassive);
        flux[n] = u[n] * v_adv;

    }

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
normalize_species_flux (int i, int j, int k, Array4<Real> const& flux)
{
    // Normalize the fluxes of the mass fractions so that
    // they sum to 0. This is essentially the CMA procedure that is
    // defined in Plewa & Muller, 1999, A&A, 342, 179.

    Real sum = 0.0_rt;

    for (int n = UFS; n < UFS+NumSpec; n++) {
      sum += flux(i,j,k,n);
    }

    Real fac = 1.0_rt;

    // We skip the normalization if the sum is zero or within epsilon.
    // There can be numerical problems here if the density flux is
    // approximately zero at the interface but not exactly, resulting in
    // division by a small number and/or resulting in one of the species
    // fluxes being negative because of roundoff error. There are also other
    // terms like artificial viscosity which can cause these problems.
    // So checking that sum is sufficiently large helps avoid this.

    if (std::abs(sum) > std::numeric_limits<Real>::epsilon() * std::abs(flux(i,j,k,URHO))) {
      fac = flux(i,j,k,URHO) / sum;
    }

    for (int n = UFS; n < UFS+NumSpec; n++) {
      flux(i,j,k,n) = flux(i,j,k,n) * fac;
    }
}
#endif
