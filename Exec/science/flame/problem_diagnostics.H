#ifndef problem_diagnostics_H
#define problem_diagnostics_H

AMREX_INLINE
void
Castro::problem_diagnostics ()
{
    int finest_level = parent->finestLevel();
    Real time        = state[State_Type].curTime();

    Real T_max = -1.0e200;
    Real T_min =  1.0e200;
    Real grad_T_max = -1.0e200;
    Real rho_fuel_dot = 0.0;
    Real rho_fuel_initial = 1.0;

    Real flame_width = 0.0;
    Real flame_speed = 0.0;

    for (int lev = 0; lev <= finest_level; lev++)
    {
        Castro& ca_lev = getLevel(lev);

        ca_lev.flame_width_properties(time, T_max, T_min, grad_T_max);

        // For now we only do the flame speed calculation on the coarse grid since
        // it doesn't share the same timestep as the fine grids.

        if (lev == 0)
            ca_lev.flame_speed_properties(time, rho_fuel_dot);
    }

    ParallelDescriptor::ReduceRealMax(T_max);
    ParallelDescriptor::ReduceRealMin(T_min);
    ParallelDescriptor::ReduceRealMax(grad_T_max);

    flame_width = (T_max - T_min) / grad_T_max;

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
}

#endif
