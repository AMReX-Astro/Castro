#ifndef CASTRO_FLUX_UTIL_H
#define CASTRO_FLUX_UTIL_H

#include <AMReX_Array4.H>
#include <AMReX_GpuContainers.H>
#include <Castro_util.H>

#include <algorithm>
#include <cmath>
#include <limits>

namespace hydro {

// These face-local operations can be composed in one kernel. Keep their
// order: artificial viscosity, density limiting, then species normalization.

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void apply_av(int i, int j, int k, int idir,
               amrex::Array4<amrex::Real const> const& div,
               amrex::Array4<amrex::Real const> const& uin,
               amrex::Array4<amrex::Real> const& flux,
               amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
               int coord,
               amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& problo,
               amrex::Real diff_coeff)
{
    Real dL = dx[idir];
    if (coord == 2 && idir == 1) {
        Real r = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0];
        dL *= r;
    }

    Real div1;
    if (idir == 0) {
        div1 = 0.25_rt * (div(i,j,k) + div(i,j+dg1,k) +
                          div(i,j,k+dg2) + div(i,j+dg1,k+dg2));
    } else if (idir == 1) {
        div1 = 0.25_rt * (div(i,j,k) + div(i+1,j,k) +
                          div(i,j,k+dg2) + div(i+1,j,k+dg2));
    } else {
        div1 = 0.25_rt * (div(i,j,k) + div(i+1,j,k) +
                          div(i,j+dg1,k) + div(i+1,j+dg1,k));
    }

    div1 = diff_coeff * std::min(0.0_rt, div1);

    for (int n = 0; n < NUM_STATE; ++n) {

        if (n == UTEMP) {
            continue;
        }
#ifdef SHOCK_VAR
        if (n == USHK) {
            continue;
        }
#endif
#ifdef NSE_NET
        if (n == UMUP || n == UMUN) {
            continue;
        }
#endif

        Real div_var{};

        if (idir == 0) {
            div_var = div1 * (uin(i,j,k,n) - uin(i-1,j,k,n));
        } else if (idir == 1) {
            div_var = div1 * (uin(i,j,k,n) - uin(i,j-dg1,k,n));
        } else {
            div_var = div1 * (uin(i,j,k,n) - uin(i,j,k-dg2,n));
        }

        flux(i,j,k,n) += dL * div_var;
    }
}

#ifdef RADIATION
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void apply_av_rad(int i, int j, int k, int idir,
                   amrex::Array4<amrex::Real const> const& div,
                   amrex::Array4<amrex::Real const> const& Erin,
                   amrex::Array4<amrex::Real> const& radflux,
                   amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
                   amrex::Real diff_coeff, int ngroups)
{

    Real div1;
    if (idir == 0) {
        div1 = 0.25_rt * (div(i,j,k) + div(i,j+dg1,k) +
                          div(i,j,k+dg2) + div(i,j+dg1,k+dg2));
    } else if (idir == 1) {
        div1 = 0.25_rt * (div(i,j,k) + div(i+1,j,k) +
                          div(i,j,k+dg2) + div(i+1,j,k+dg2));
    } else {
        div1 = 0.25_rt * (div(i,j,k) + div(i+1,j,k) +
                          div(i,j+dg1,k) + div(i+1,j+dg1,k));
    }

    div1 = diff_coeff * std::min(0.0_rt, div1);

    for (int n = 0; n < ngroups; ++n) {

        Real div_var{};

        if (idir == 0) {
            div_var = div1 * (Erin(i,j,k,n) - Erin(i-1,j,k,n));
        } else if (idir == 1) {
            div_var = div1 * (Erin(i,j,k,n) - Erin(i,j-dg1,k,n));
        } else {
            div_var = div1 * (Erin(i,j,k,n) - Erin(i,j,k-dg2,n));
        }

        radflux(i,j,k,n) += dx[idir] * div_var;
    }
}
#endif

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void normalize_species_fluxes(int i, int j, int k, amrex::Array4<amrex::Real> const& flux)
{

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

#ifndef MHD
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void limit_hydro_fluxes_on_small_dens(int i, int j, int k, int idir,
                                    amrex::Array4<amrex::Real const> const& u,
                                    amrex::Array4<amrex::Real const> const& vol,
                                    amrex::Array4<amrex::Real> const& flux,
                                    amrex::Array4<amrex::Real const> const& area_arr,
                                    amrex::Real dt, amrex::Real density_floor,
                                    bool scale_by_dAdt)
{
    // Grab the states on either side of the interface we are working with,
    // depending on which dimension we're currently calling this with.

    Real rhoR = u(i,j,k,URHO);
    Real volR = vol(i,j,k);

    Real rhoL, volL;

    if (idir == 0) {
        rhoL = u(i-1,j,k,URHO);
        volL = vol(i-1,j,k);
    }
    else if (idir == 1) {
        rhoL = u(i,j-1,k,URHO);
        volL = vol(i,j-1,k);
    }
    else {
        rhoL = u(i,j,k-1,URHO);
        volL = vol(i,j,k-1);
    }

    // Coefficients of fluxes on either side of the interface.

    Real flux_coefR = 1.0_rt / volR;
    Real flux_coefL = 1.0_rt / volL;

    if (scale_by_dAdt) {
        flux_coefR *= dt * area_arr(i,j,k);
        flux_coefL *= dt * area_arr(i,j,k);
    }

    // Updates to the zones on either side of the interface.

    Real drhoR = flux_coefR * flux(i,j,k,URHO);
    Real drhoL = flux_coefL * flux(i,j,k,URHO);

    // Limit all fluxes such that the zone does not go negative in density.
    // We apply this flux limiter on a per-edge basis. So we can guarantee
    // that any individual flux cannot cause a small density in one step,
    // but with the above floor we cannot guarantee that the sum of the
    // fluxes will enforce this constraint. A simple way to guarantee that
    // is to note that in the worst case, all six faces adjacent to a zone
    // could be evacuating that zone (more generally, 2 * AMREX_SPACEDIM).
    // So instead of checking whether rho + drho is less than the floor,
    // we can check whether rho + 6 * drho is less than the floor, and
    // strengthen the limiter by that same factor of 6, so that even in the
    // worst case, the zone will remain positive. This is simple to implement
    // at the risk of being overly aggressive on flux limiting. A more
    // sophisticated version of this might check the ratio of the fluxes
    // (or, similarly, the ratio of the edge velocities coming from the
    // Riemann solver) and scale the larger fluxes more strongly than the
    // weaker fluxes, and perhaps only consider fluxes that subtract from
    // the zone, but this would also be more complicated to implement.

    if (rhoR >= density_floor && std::abs(drhoR) > 0.0_rt && rhoR + 2 * AMREX_SPACEDIM * drhoR < density_floor) {
        Real limiting_factor = std::abs((density_floor - rhoR) / (2 * AMREX_SPACEDIM * drhoR));
        for (int n = 0; n < NUM_STATE; ++n) {
            flux(i,j,k,n) = flux(i,j,k,n) * limiting_factor;
        }
    }
    else if (rhoL >= density_floor && std::abs(drhoL) > 0.0_rt && rhoL - 2 * AMREX_SPACEDIM * drhoL < density_floor) {
        Real limiting_factor = std::abs((density_floor - rhoL) / (2 * AMREX_SPACEDIM * drhoL));
        for (int n = 0; n < NUM_STATE; ++n) {
            flux(i,j,k,n) = flux(i,j,k,n) * limiting_factor;
        }
    }
}
#endif

}

#endif
