#ifndef problem_diagnostics_H
#define problem_diagnostics_H

AMREX_INLINE
void
Castro::problem_diagnostics ()
{
    int finest_level = parent->finestLevel();
    Real time        = state[State_Type].curTime();
    Real mass        = 0.0;
    Real rho_E       = 0.0;
    Real Tmax        = -std::numeric_limits<Real>::max();
    Real MachMax     = -std::numeric_limits<Real>::max();
    Real Tmax_level;
    Real MachMax_level;

    for (int lev = 0; lev <= finest_level; lev++)
    {
        Castro& ca_lev = getLevel(lev);

        mass     += ca_lev.volWgtSum("density", time);

        rho_E    += ca_lev.volWgtSum("rho_E", time);

        auto temp_mf = ca_lev.derive("Temp", time, 0);
        Tmax_level = temp_mf->max(0, 0);
        if (Tmax_level > Tmax)
        {
            Tmax = Tmax_level;
        }

        auto mach_mf = ca_lev.derive("MachNumber", time, 0);
        MachMax_level = mach_mf->max(0, 0);
        if (MachMax_level > MachMax)
        {
            MachMax = MachMax_level;
        }

    }


    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        std::cout << '\n';
        std::cout << "TIME= " << time << " MASS  = "   << mass  << '\n';
        std::cout << "TIME= " << time << " RHO*E = "   << rho_E     << '\n';

        // get output files
        std::ostream& data_log1 = *Castro::problem_data_logs[0];

        if (time == 0.0) {
           data_log1 << std::setw(14) <<  "#     time    ";
           data_log1 << std::setw(14) <<  "       max(T) ";
           data_log1 << std::setw(14) <<  "       max(M) ";
           data_log1 << std::setw(14) <<  "        rho_E ";
           data_log1 << std::endl;
        }

        // Write the quantities at this time
        data_log1 << std::setw(14) <<  time;
        data_log1 << std::setw(14) <<  std::setprecision(6) << Tmax;
        data_log1 << std::setw(14) <<  std::setprecision(6) << MachMax;
        data_log1 << std::setw(14) <<  std::setprecision(6) << rho_E;
        data_log1 << std::endl;

        std::cout<<'\n';
    }
}

#endif
