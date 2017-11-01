#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::just_the_mhd(Real time, Real dt)
{
      if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... mhd ...!!! " << std::endl << std::endl;

      const int finest_level = parent->finestLevel();
      MultiFab& S_old        = get_old_data(State_Type);
      MultiFab& Bx_old       = get_old_data(Mag_Type_x);
      
      ca_advance_mhd();

}
