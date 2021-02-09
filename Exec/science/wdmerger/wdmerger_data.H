#ifndef wdmerger_data_H
#define wdmerger_data_H

#include <AMReX.H>
#include <AMReX_REAL.H>

#include <prob_parameters.H>
#include <prob_parameters_F.H>

namespace wdmerger
{
    // Input parameters

    extern AMREX_GPU_MANAGED int use_stopping_criterion;
    extern AMREX_GPU_MANAGED int use_energy_stopping_criterion;
    extern AMREX_GPU_MANAGED amrex::Real ts_te_stopping_criterion;
    extern AMREX_GPU_MANAGED amrex::Real T_stopping_criterion;

    // Data we want to store, for interfacing with the Fortran

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

    // Current values of the above extrema

    extern AMREX_GPU_MANAGED amrex::Real T_curr_max;
    extern AMREX_GPU_MANAGED amrex::Real rho_curr_max;
    extern AMREX_GPU_MANAGED amrex::Real ts_te_curr_max;
}

#endif
