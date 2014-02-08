#include <winstd.H>

#include "Castro.H"
#include "Castro_F.H"

using std::string;

#ifdef RADIATION
#include "Radiation.H"

void
Castro::final_radiation_call (MultiFab& S_new, int iteration, int ncycle) 
{
  if (do_radiation) {

      if (Radiation::pure_hydro) {
	  radiation->computeTemp(S_new, 0);
	  return;
      }

    if (add_ext_src<=0 
#ifdef GRAVITY
	&& do_grav<=0
#endif
	) 
      {
	if (radiation->do_real_eos > 0) {
	  radiation->computeTemp(S_new, 1);
	} 
      }

    if (Radiation::filter_prim_int > 0 && Radiation::filter_prim_T>0 
	&& iteration==ncycle) {
      int nstep = parent->levelSteps(0);
      if (nstep % Radiation::filter_prim_int == 0) {
	radiation->filter_prim(level, S_new);
      }
    }

    if (Radiation::SolverType == Radiation::MGFLDSolver) {
      if (! Radiation::rad_hydro_combined) {
	MultiFab& Er_old = get_old_data(Rad_Type);
	MultiFab& Er_new = get_new_data(Rad_Type);
	Er_new.copy(Er_old);	
      }
      radiation->MGFLD_implicit_update(level, iteration, ncycle);
    }
    else if (Radiation::SolverType == Radiation::SGFLDSolver) {
      if (! Radiation::rad_hydro_combined) {
	MultiFab& Er_old = get_old_data(Rad_Type);
	MultiFab& Er_new = get_new_data(Rad_Type);
	Er_new.copy(Er_old);	
      }
      radiation->single_group_update(level, iteration, ncycle);
    }
    else {
      MultiFab& Er_old = get_old_data(Rad_Type);
      MultiFab& Er_new = get_new_data(Rad_Type);
      Er_new.copy(Er_old);
      radiation->single_group_update(level, iteration, ncycle);
    }
    
    // Recompute temperatures after radiation update
    if (Radiation::SolverType != Radiation::MGFLDSolver) {
      if (radiation->do_real_eos > 0) {
	radiation->computeTemp(S_new, 0);
      }
    }
  }
}
#endif

