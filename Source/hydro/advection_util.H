#ifndef advection_util_H
#define advection_util_H

#include <Castro_util.H>

#ifdef HYBRID_MOMENTUM
#include <hybrid.H>
#endif

namespace hydro {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void src_to_prim(int i, int j, int k, const Real dt,
                 Array4<Real const> const& q_arr,
                 Array4<Real const> const& old_src,
#ifndef TRUE_SDC
                 Array4<Real const> const& src_corr,
#endif
                 Array4<Real> const& srcQ)
{
    for (int n = 0; n < NQSRC; ++n) {
      srcQ(i,j,k,n) = 0.0_rt;
    }

    // the conserved source may have a predictor that time-centers it

    Real srcU[NSRC] = {0.0_rt};

    for (int n = 0; n < NSRC; n++) {

#ifndef TRUE_SDC
        if (time_integration_method == CornerTransportUpwind && source_term_predictor == 1) {
            if (n == UMX || n == UMY || n == UMZ) {
                srcU[n] += 0.5 * dt * src_corr(i,j,k,n);
            }
        } else if (time_integration_method == SimplifiedSpectralDeferredCorrections  && source_term_predictor == 1) {
            srcU[n] += src_corr(i,j,k,n);
        }
#endif

        srcU[n] += old_src(i,j,k,n);
    }

    Real rhoinv = 1.0_rt / q_arr(i,j,k,QRHO);

    // get the needed derivatives
    eos_rep_t eos_state;
    eos_state.T = q_arr(i,j,k,QTEMP);
    eos_state.rho = q_arr(i,j,k,QRHO);
    eos_state.e = q_arr(i,j,k,QREINT) * rhoinv;
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n]  = q_arr(i,j,k,QFS+n);
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; n++) {
      eos_state.aux[n] = q_arr(i,j,k,QFX+n);
    }
#endif

    eos(eos_input_re, eos_state);

    srcQ(i,j,k,QRHO) = srcU[URHO];
    srcQ(i,j,k,QU) = (srcU[UMX] - q_arr(i,j,k,QU) * srcQ(i,j,k,QRHO)) * rhoinv;
    srcQ(i,j,k,QV) = (srcU[UMY] - q_arr(i,j,k,QV) * srcQ(i,j,k,QRHO)) * rhoinv;
    srcQ(i,j,k,QW) = (srcU[UMZ] - q_arr(i,j,k,QW) * srcQ(i,j,k,QRHO)) * rhoinv;
    srcQ(i,j,k,QREINT) = srcU[UEINT];
    srcQ(i,j,k,QPRES ) = eos_state.dpde *
      (srcQ(i,j,k,QREINT) - q_arr(i,j,k,QREINT) * srcQ(i,j,k,QRHO)*rhoinv) *
      rhoinv + eos_state.dpdr_e * srcQ(i,j,k,QRHO);

#ifdef PRIM_SPECIES_HAVE_SOURCES
    for (int ipassive = 0; ipassive < npassive; ++ipassive) {
      int n = upassmap(ipassive);
      int iq = qpassmap(ipassive);

      // note: this does not include any SDC sources -- those are
      // computed and stored separately
      srcQ(i,j,k,iq) = (srcU[n] - q_arr(i,j,k,iq) * srcQ(i,j,k,QRHO) ) /
        q_arr(i,j,k,QRHO);
    }
#endif
}

}

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

#endif
