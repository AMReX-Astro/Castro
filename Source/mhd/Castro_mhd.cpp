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
      MultiFab& By_old       = get_old_data(Mag_Type_y);
      MultiFab& Bz_old       = get_old_data(Mag_Type_z);
      MultiFab& S_new        = get_old_data(State_Type);
      MultiFab& Bx_new       = get_new_data(Mag_Type_x);
      MultiFab& By_new       = get_new_data(Mag_Type_y);
      MultiFab& Bz_new       = get_new_data(Mag_Type_z);
      ca_advance_mhd();

}
