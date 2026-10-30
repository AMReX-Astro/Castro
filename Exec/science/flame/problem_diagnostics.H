#ifndef problem_diagnostics_H
#define problem_diagnostics_H

#include <AMReX_REAL.H>

AMREX_INLINE
void
Castro::problem_diagnostics ()
{
    int finest_level = parent->finestLevel();
    amrex::Real time = state[State_Type].curTime();

    // for flame width
    amrex::Real T_max = -1.0e200;
    amrex::Real T_min =  1.0e200;
    amrex::Real grad_T_max = -1.0e200;

    // for flame speed
    amrex::Real rho_fuel_dot = 0.0;
    amrex::Real rho_fuel_initial = 1.0;

    // for peak energy generation location
    amrex::Real rho_enuc_max = -1.e200;
    amrex::Real peak_x = -1.e200;

    // derived quantities
    amrex::Real flame_width = 0.0;
    amrex::Real flame_speed = 0.0;

    for (int lev = 0; lev <= finest_level; lev++)
    {
        Castro& ca_lev = getLevel(lev);

        ca_lev.flame_width_properties(time, T_max, T_min, grad_T_max);

        ca_lev.flame_peak_properties(time, rho_enuc_max, peak_x);

        // For now we only do the flame speed calculation on the coarse grid since
        // it doesn't share the same timestep as the fine grids.

        if (lev == 0) {
            ca_lev.flame_speed_properties(time, rho_fuel_dot);
        }

    }

    ParallelDescriptor::ReduceRealMax(T_max);
    ParallelDescriptor::ReduceRealMin(T_min);
    ParallelDescriptor::ReduceRealMax(grad_T_max);

    flame_width = (T_max - T_min) / grad_T_max;

    ParallelDescriptor::ReduceRealMax(rho_enuc_max);
    ParallelDescriptor::ReduceRealMax(peak_x);

    ParallelDescriptor::ReduceRealSum(rho_fuel_dot);

    // Note that rho_fuel_dot has already been multiplied by dx, so
    // the dimensionality here checks out.

    flame_speed = -rho_fuel_dot / rho_fuel_initial;

    if (ParallelDescriptor::IOProcessor()) {

        if (verbose > 0) {
            std::cout << '\n';
            std::cout << "TIME= " << time << " FLAME WIDTH = "   << flame_width << '\n';
            std::cout << "TIME= " << time << " FLAME SPEED = "   << flame_speed << '\n';
        }

        std::ostream& datalog = *Castro::problem_data_logs[0];

        if (time == 0.0) {
            datalog << std::setw(14) <<  "      time    ";
            datalog << std::setw(14) <<  "  flame width ";
            datalog << std::setw(14) <<  "  flame speed " << std::endl;
        }

        // Write the quantities at this time
        datalog << std::setw(14) <<  time;
        datalog << std::setw(14) <<  std::setprecision(6) << flame_width;
        datalog << std::setw(14) <<  std::setprecision(6) << flame_speed << std::endl;

    }

    // check if the simulation is over -- have we crossed enough of
    // the domain?

    const auto prob_lo = geom.ProbLoArray();
    const auto prob_hi = geom.ProbHiArray();
    const amrex::Real abort_x = prob_lo[0] +
        problem::flame_domain_abort * (prob_hi[0] - prob_lo[0]);

    if (rho_enuc_max > 0.0_rt && peak_x > abort_x) {
        amrex::Abort("flame has reached the edge of the domain");
    }
}

#endif
