#ifndef wdmerger_data_H
#define wdmerger_data_H

#include <AMReX.H>
#include <AMReX_REAL.H>

#include <prob_parameters.H>

namespace wdmerger
{
    // Data we want to store over time

    // Rate of change of stellar masses

    extern AMREX_GPU_MANAGED amrex::Real mdot_P;
    extern AMREX_GPU_MANAGED amrex::Real mdot_S;

    // Radii of the WDs at various density thresholds.

    extern AMREX_GPU_MANAGED amrex::Real rad_P[7];
    extern AMREX_GPU_MANAGED amrex::Real rad_S[7];

    // Effective volume of the stars at various density cutoffs.

    extern AMREX_GPU_MANAGED amrex::Real vol_P[7];
    extern AMREX_GPU_MANAGED amrex::Real vol_S[7];

    // Average density of the stars.

    extern AMREX_GPU_MANAGED amrex::Real rho_avg_P;
    extern AMREX_GPU_MANAGED amrex::Real rho_avg_S;
}

#endif
