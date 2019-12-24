#include "Castro.H"
#include "Castro_F.H"

#include "Gravity.H"
#include <Gravity_F.H>

#include "AMReX_ParmParse.H"

#include "AMReX_buildInfo.H"

using namespace amrex;

int Castro::problem = -1;
Real Castro::diameter = 0.0;
Real Castro::density = 0.0;

void Castro::problem_post_init() {

  // Read in inputs.

  ParmParse pp("castro");

  // Get the problem number fom Fortran.

  get_problem_number(&problem);

  // Get the diameter.

  get_diameter(&diameter);

  // Get the density

  get_density(&density);

  // If we're doing problem 2, the normalized sphere,
  // add up the mass on the domain and then update
  // the density in the sphere so that it has the 'correct'
  // amount of total mass.

  if (problem == 2) {

      Real actual_mass = 0.0;

      bool local_flag = true;
      Real time       = state[State_Type].curTime();

      for (int lev = 0; lev <= parent->finestLevel(); lev++)
	  actual_mass += getLevel(lev).volWgtSum("density", time, local_flag);

      ParallelDescriptor::ReduceRealSum(actual_mass);

      // The correct amount of mass is the mass of a sphere
      // with the given diameter and density.

      Real target_mass = density * (1.0e0 / 6.0e0) * M_PI * std::pow(diameter, 3);

      Real update_factor = target_mass / actual_mass;

      // Now update the density given this factor.

      amrex::Print() << "\n";
      amrex::Print() << "  Updating density by the factor " << update_factor << " to ensure total mass matches target mass.\n";
      amrex::Print() << "\n";

      for (int lev = 0; lev <= parent->finestLevel(); lev++) {

	  MultiFab& state = getLevel(lev).get_new_data(State_Type);

	  const Real* dx = getLevel(lev).geom.CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
	  for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

	    const Box& box = mfi.tilebox();

	    const int* lo  = box.loVect();
	    const int* hi  = box.hiVect();

#pragma gpu box(box)
	    update_density(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
			   BL_TO_FORTRAN_ANYD(state[mfi]),
			   dx, update_factor);

	  }

      }

      // Do a final check, for sanity purposes.

      actual_mass = 0.0;

      for (int lev = 0; lev <= parent->finestLevel(); lev++)
	  actual_mass += getLevel(lev).volWgtSum("density", time, local_flag);

      ParallelDescriptor::ReduceRealSum(actual_mass);

      if (std::abs( (actual_mass - target_mass) / target_mass ) > 1.0e-6) {
          amrex::Print() << "\n";
          amrex::Print() << "Actual mass: " << actual_mass << "\n";
          amrex::Print() << "Target mass: " << target_mass << "\n";
          amrex::Print() << "\n";
	  amrex::Abort("Sphere does not have the right amount of mass.");
      }

  }

}
