#ifndef wdmerger_util_H
#define wdmerger_util_H

#include <AMReX.H>

#include <math.H>
#include <fundamental_constants.H>

#ifdef ROTATION
#include <Rotation.H>
#endif

void freefall_velocity (Real mass, Real distance, Real& vel);

void kepler_third_law (Real radius_1, Real mass_1, Real radius_2, Real mass_2,
                       Real& period, Real eccentricity, Real phi, Real& a,
                       Real& r_1, Real& r_2, Real& v_1r, Real& v_2r, Real& v_1p, Real& v_2p);

void set_wd_composition (Real mass, Real& envelope_mass, Real core_comp[NumSpec], Real envelope_comp[NumSpec], const std::string& star_type);

void ensure_primary_mass_larger ();

void initialize_problem ();

void finalize_probdata ();

void binary_setup ();

void set_small ();

void update_roche_radii ();

// Define a mask for all zones considered to be within either the
// primary or secondary star. It uses the effective potential of
// the star (phiRot + phiGrav) and the stellar density threshold
// input parameter to determine what parts of the domain should
// be considered stellar material. The convention will be that
// the mask is positive (1) for zones inside the star and negative
// (-1) for zones outside the star.

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real stellar_mask (int i, int j, int k,
                   const amrex::GeometryData& geomdata,
                   Real rho, bool is_primary = true)
{
    // By default, assume we're not inside the star.

    Real mask = -1.0_rt;

    // Don't do anything here if the star no longer exists,
    // or if it never existed.

    if (is_primary) {
        if (problem::mass_P <= 0.0_rt) return mask;
    }
    else {
        if (problem::mass_S <= 0.0_rt) return mask;
    }

    GpuArray<Real, 3> loc;
    loc[0] = geomdata.ProbLo(0) + (static_cast<Real>(i) + 0.5_rt) * geomdata.CellSize(0);

#if AMREX_SPACEDIM >= 2
    loc[1] = geomdata.ProbLo(1) + (static_cast<Real>(j) + 0.5_rt) * geomdata.CellSize(1);
#else
    loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
    loc[2] = geomdata.ProbLo(2) + (static_cast<Real>(k) + 0.5_rt) * geomdata.CellSize(2);
#else
    loc[2] = 0.0_rt;
#endif

    // Ignore zones whose density is too low.

    if (rho < problem::stellar_density_threshold) return mask;

    Real r_P = std::sqrt((loc[0] - problem::com_P[0]) * (loc[0] - problem::com_P[0]) +
                         (loc[1] - problem::com_P[1]) * (loc[1] - problem::com_P[1]) +
                         (loc[2] - problem::com_P[2]) * (loc[2] - problem::com_P[2]));

    Real r_S = std::sqrt((loc[0] - problem::com_S[0]) * (loc[0] - problem::com_S[0]) +
                         (loc[1] - problem::com_S[1]) * (loc[1] - problem::com_S[1]) +
                         (loc[2] - problem::com_S[2]) * (loc[2] - problem::com_S[2]));


    for (int iloc = 0; iloc < AMREX_SPACEDIM; ++iloc) {
        loc[iloc] -= problem::center[iloc];
    }

    auto omega = get_omega_vec(geomdata, j);

    Real phi_rot = rotational_potential(loc, omega, geomdata.Coord());

    Real phi_p = -C::Gconst * problem::mass_P / r_P + phi_rot;
    Real phi_s = -C::Gconst * problem::mass_S / r_S + phi_rot;

    if (is_primary) {
        if (phi_p < 0.0_rt && phi_p < phi_s) {
            mask = 1.0_rt;
        }
    }
    else {
        if (phi_s < 0.0_rt && phi_s < phi_p) {
            mask = 1.0_rt;
        }
    }

    return mask;
}

#endif
