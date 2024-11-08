#ifndef CASTRO_HYBRID_H
#define CASTRO_HYBRID_H

#include <Castro_util.H>

using namespace amrex;

///
/// For hybrid momentum advection, convert the linear momentum into the
/// hybrid / angular momentum components
///
/// @param loc         spatial location
/// @param linear_mom  linear momentum
/// @param hybrid_mom  hybrid / angular momentum
///
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void linear_to_hybrid(const GpuArray<Real, 3>& loc,
                      const GpuArray<Real, 3>& linear_mom,
                      GpuArray<Real, 3>& hybrid_mom)
{
    // Convert a linear momentum into the "hybrid" scheme
    // that has radial and angular components.

    Real R = amrex::max(std::sqrt(loc[0] * loc[0] + loc[1] * loc[1]),
                        std::numeric_limits<Real>::min());
    Real RInv = 1.0_rt / R;

    // This conversion is Eqs. 25 and 26 in Byerly et al. 2014.
    // Note that we expect the linear momentum to be consistent
    // with which frame we're measuring the fluid quantities in.
    // So we're effectively always using the first form of those
    // equalities, not the second.

    hybrid_mom[0] = linear_mom[0] * (loc[0] * RInv) + linear_mom[1] * (loc[1] * RInv);
    hybrid_mom[1] = linear_mom[1] * loc[0] - linear_mom[0] * loc[1];
    hybrid_mom[2] = linear_mom[2];
}


///
/// Convert the hybrid / angular momentum back to linear momentum
///
/// @param loc          spatial location
/// @param hybrid_mom   hybrid / angular momentum
/// @param linear_mom   linear momentum
///
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void hybrid_to_linear(const GpuArray<Real, 3>& loc,
                      const GpuArray<Real, 3>& hybrid_mom,
                      GpuArray<Real, 3>& linear_mom)
{
    // Convert a "hybrid" momentum into a linear one.

    Real R = amrex::max(std::sqrt(loc[0] * loc[0] + loc[1] * loc[1]),
                        std::numeric_limits<Real>::min());
    Real RInv = 1.0_rt / R;

    // This is the inverse of Byerly et al., Equations 25 and 26.

    linear_mom[0] = hybrid_mom[0] * (loc[0] * RInv) - hybrid_mom[1] * (loc[1] * RInv * RInv);
    linear_mom[1] = hybrid_mom[1] * (loc[0] * RInv * RInv) + hybrid_mom[0] * (loc[1] * RInv);
    linear_mom[2] = hybrid_mom[2];
}


///
/// Create the hydrodynamic source term for the hybrid momenta based on the
/// input linear momentum source term
///
/// @param loc             spatial location
/// @param linear_source   linear momentum source term
/// @param hybrid_source   hybrid / angular momentum source term
///
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void set_hybrid_momentum_source(const GpuArray<Real, 3>& loc,
                                const GpuArray<Real, 3>& linear_source,
                                GpuArray<Real, 3>& hybrid_source)
{
    Real R = amrex::max(std::sqrt(loc[0] * loc[0] + loc[1] * loc[1]),
                        std::numeric_limits<Real>::min());
    Real RInv = 1.0_rt / R;

    // This is analogous to the conversion of linear momentum to hybrid momentum.

    hybrid_source[0] = linear_source[0] * (loc[0] * RInv) + linear_source[1] * (loc[1] * RInv);
    hybrid_source[1] = linear_source[1] * loc[0] - linear_source[0] * loc[1];
    hybrid_source[2] = linear_source[2];
}


///
/// Compute and store the flux of the hybrid momentum components
///
/// @param state          the Godunov state
/// @param geomdata       the Castro geometry data object
/// @param idir           coordinate direction for flux computation (0 = x, 1 = y, 2 = z)
/// @param i              x zone index
/// @param j              y zone index
/// @param k              z zone indxe
/// @param flux           hydrodynamic flux in direction idir
/// @param cell_centered  is this a cell-centered flux or an interface flux?
///
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compute_hybrid_flux(const GpuArray<Real, NGDNV>& state, const GeometryData& geomdata,
                         int idir, int i, int j, int k,
                         GpuArray<Real, NUM_STATE>& flux, bool cell_centered = false)
{
    bool ccx = true;
    bool ccy = true;
    bool ccz = true;

    if (idir == 0) {
        ccx = cell_centered;
    }
    else if (idir == 1) {
        ccy = cell_centered;
    }
    else {
        ccz = cell_centered;
    }

    GpuArray<Real, 3> loc;

    position(i, j, k, geomdata, loc, ccx, ccy, ccz);

    GpuArray<Real, 3> linear_mom;

    for (int dir = 0; dir < 3; ++dir) {
        linear_mom[dir] = state[GDRHO] * state[GDU + dir];
    }

    GpuArray<Real, 3> hybrid_mom;

    linear_to_hybrid(loc, linear_mom, hybrid_mom);

    Real u_adv = state[GDU + idir];

    if (idir == 0) {
        flux[UMR] = hybrid_mom[0] * u_adv;
        flux[UML] = hybrid_mom[1] * u_adv - loc[1] * state[GDPRES];
        flux[UMP] = hybrid_mom[2] * u_adv;
    }
    else if (idir == 1) {
        flux[UMR] = hybrid_mom[0] * u_adv;
        flux[UML] = hybrid_mom[1] * u_adv + loc[0] * state[GDPRES];
        flux[UMP] = hybrid_mom[2] * u_adv;
    }
    else {
        flux[UMR] = hybrid_mom[0] * u_adv;
        flux[UML] = hybrid_mom[1] * u_adv;
        flux[UMP] = hybrid_mom[2] * u_adv + state[GDPRES];
    }

}

#endif
