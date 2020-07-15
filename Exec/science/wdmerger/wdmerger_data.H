#ifndef wdmerger_data_H
#define wdmerger_data_H

#include <AMReX.H>
#include <AMReX_REAL.H>

namespace wdmerger
{
    // Input parameters

    extern AMREX_GPU_MANAGED int use_stopping_criterion;
    extern AMREX_GPU_MANAGED int use_energy_stopping_criterion;
    extern AMREX_GPU_MANAGED amrex::Real ts_te_stopping_criterion;
    extern AMREX_GPU_MANAGED amrex::Real T_stopping_criterion;

    // Data we want to store, for interfacing with the Fortran

    extern AMREX_GPU_MANAGED int relaxation_is_done;
    extern AMREX_GPU_MANAGED int problem;

    // Stellar masses

    extern AMREX_GPU_MANAGED amrex::Real mass_p;
    extern AMREX_GPU_MANAGED amrex::Real mass_s;

    // Rate of change of stellar masses

    extern AMREX_GPU_MANAGED amrex::Real mdot_p;
    extern AMREX_GPU_MANAGED amrex::Real mdot_s;

    // Stellar centers of mass and velocities

    extern AMREX_GPU_MANAGED amrex::Real com_p[3];
    extern AMREX_GPU_MANAGED amrex::Real com_s[3];

    extern AMREX_GPU_MANAGED amrex::Real vel_p[3];
    extern AMREX_GPU_MANAGED amrex::Real vel_s[3];

    // Radii of the WDs at various density thresholds.

    extern AMREX_GPU_MANAGED amrex::Real rad_p[7];
    extern AMREX_GPU_MANAGED amrex::Real rad_s[7];

    // Effective volume of the stars at various density cutoffs.

    extern AMREX_GPU_MANAGED amrex::Real vol_p[7];
    extern AMREX_GPU_MANAGED amrex::Real vol_s[7];

    // Average density of the stars.

    extern AMREX_GPU_MANAGED amrex::Real rho_avg_p;
    extern AMREX_GPU_MANAGED amrex::Real rho_avg_s;

    // Gravitational free-fall timescale of the stars.

    extern AMREX_GPU_MANAGED amrex::Real t_ff_p;
    extern AMREX_GPU_MANAGED amrex::Real t_ff_s;

    // Global extrema of various quantities over the whole simulation

    extern AMREX_GPU_MANAGED amrex::Real T_global_max;
    extern AMREX_GPU_MANAGED amrex::Real rho_global_max;
    extern AMREX_GPU_MANAGED amrex::Real ts_te_global_max;

    // Current values of the above extrema

    extern AMREX_GPU_MANAGED amrex::Real T_curr_max;
    extern AMREX_GPU_MANAGED amrex::Real rho_curr_max;
    extern AMREX_GPU_MANAGED amrex::Real ts_te_curr_max;

    // Value of the total energy on the domain over the last several timesteps

    const int num_previous_ener_timesteps = 5;
    extern AMREX_GPU_MANAGED amrex::Real total_ener_array[num_previous_ener_timesteps];
}

#endif
