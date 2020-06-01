#include <iomanip>

#include <Castro.H>
#include <Castro_F.H>

using namespace amrex;

void
Castro::sum_integrated_quantities ()
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

	Tmax_level = ca_lev.maxVal("Temp", time);
	if (Tmax_level > Tmax)
	  {
	    Tmax = Tmax_level;
	  }

	MachMax_level = ca_lev.maxVal("MachNumber", time);
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
	std::ostream& data_log1 = parent->DataLog(0);

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


Real
Castro::maxVal (const std::string& name,
                Real               time)
{
  Real        maxval  = 0.0;
  const Real* dx      = geom.CellSize();
  auto        mf      = derive(name,time,0);
  BL_ASSERT(mf);

  maxval = (*mf).max(0,0);

  return maxval;
}

