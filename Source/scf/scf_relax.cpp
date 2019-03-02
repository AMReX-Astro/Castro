#include "Castro.H"
#include "Castro_F.H"

#include "Gravity.H"

using namespace amrex;

#ifdef GRAVITY
#ifdef ROTATION
void Castro::scf_relaxation() {

  int finest_level = parent->finestLevel();

  MultiFab& S_new = get_new_data(State_Type);

  int j = 1;
  int relax_max_iterations = 30;

  const Real* dx = parent->Geom(level).CellSize();

  // First do some sanity checks.

  // Maximum density and temperature must be set.

  if (scf_maximum_density <= 0) {
      amrex::Error("castro.scf_maximum_density must be set for SCF relaxation.");
  }

  if (scf_temperature <= 0) {
      amrex::Error("castro.scf_temperature must be set for SCF relaxation.");
  }

  // Equatorial radius and polar radius must both
  // be positive, and the polar radius cannot be
  // larger than the equatorial radius.

  if (scf_equatorial_radius <= 0) {
      amrex::Error("Equatorial radius must be positive for SCF relaxation.");
  }

  if (scf_polar_radius <= 0) {
      amrex::Error("Polar radius must be positive for SCF relaxation.");
  }

  if (scf_polar_radius > scf_equatorial_radius) {
      amrex::Error("Polar radius must not be larger than equatorial radius for SCF relaxation.");
  }

  // Do the initial relaxation setup.

  scf_setup_relaxation();

  // Get the phi MultiFab.

  MultiFab& phi = get_new_data(PhiGrav_Type);

  Real time = state[State_Type].curTime();

  // Iterate until the system is relaxed by filling the level data
  // and then doing a multilevel gravity solve.

  int is_relaxed = 0;

  while (j <= relax_max_iterations) {

     // First step is to find the rotational frequency.

     Real phi_A = 0.0;
     Real psi_A = 0.0;
     Real phi_B = 0.0;
     Real psi_B = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:phi_A, psi_A, phi_B, psi_B)
#endif
     for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi) {

       const Box& bx = mfi.tilebox();

       scf_update_for_omegasq(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                              BL_TO_FORTRAN_ANYD(S_new[mfi]),
                              BL_TO_FORTRAN_ANYD(phi[mfi]),
                              AMREX_ZFILL(dx),
                              &phi_A, &psi_A, &phi_B, &psi_B);

     }

     ParallelDescriptor::ReduceRealSum(phi_A);
     ParallelDescriptor::ReduceRealSum(psi_A);
     ParallelDescriptor::ReduceRealSum(phi_B);
     ParallelDescriptor::ReduceRealSum(psi_B);

     // Now update the square of the rotation frequency, following Hachisu (Equation 16).
     // Deal carefully with the special case where phi_A and phi_B are equal -- we assume
     // this is the non-rotating case, and we want to handle that gracefully.

     Real omegasq;
     Real omega;

     if (std::abs(phi_A - phi_B) / std::abs(phi_A) < 1.e-6) {

         omegasq = 0.0;
         omega = 0.0;
         rotational_period = 0.0;
         do_rotation = 0;

     } else {

         omegasq = -(phi_A - phi_B) / (psi_A - psi_B);

         if (omegasq < 0.0 && ParallelDescriptor::IOProcessor()) {
             std::cerr << "Omega squared = " << omegasq << " is negative in the relaxation step; aborting." << std::endl;
             amrex::Error();
         }

         omega = sqrt(omegasq);

         // Rotational period is 2 pi / omega.

         rotational_period = 2.0 * M_PI / omega;

     }

     // Now save the updated rotational frequency in the Fortran module.

     set_rot_period(&rotational_period);



     // Second step is to evaluate the Bernoulli constant.

     Real bernoulli = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:bernoulli)
#endif
     for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {

       const Box& bx = mfi.tilebox();

       scf_get_bernoulli_const(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
			       BL_TO_FORTRAN_ANYD(S_new[mfi]),
			       BL_TO_FORTRAN_ANYD(phi[mfi]),
			       AMREX_ZFILL(dx), omega, &bernoulli);

     }

     ParallelDescriptor::ReduceRealSum(bernoulli);



     // Third step is to construct the enthalpy field and
     // find the maximum enthalpy for the star.

     Real h_max = 0.0;

     MultiFab enthalpy(grids, dmap, 1, 0);
     enthalpy.setVal(0.0);

#ifdef _OPENMP
#pragma omp parallel reduction(max:h_max)
#endif
     for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {

       const Box& bx = mfi.tilebox();

       scf_construct_enthalpy(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
			      BL_TO_FORTRAN_ANYD(S_new[mfi]),
			      BL_TO_FORTRAN_ANYD(phi[mfi]),
			      BL_TO_FORTRAN_ANYD(enthalpy[mfi]),
			      AMREX_ZFILL(dx), omega,
			      &bernoulli, &h_max);

     }

     ParallelDescriptor::ReduceRealMax(h_max);

     Real kin_eng = 0.0;
     Real pot_eng = 0.0;
     Real int_eng = 0.0;
     Real mass = 0.0;
     Real Linf_norm = 0.0;

     // Finally, update the density using the enthalpy field.

#ifdef _OPENMP
#pragma omp parallel reduction(+:kin_eng,pot_eng,int_eng, mass) \
                     reduction(max:Linf_norm)
#endif
     for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {

       const Box& bx = mfi.tilebox();

       scf_update_density(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
			  BL_TO_FORTRAN_ANYD(S_new[mfi]),
			  BL_TO_FORTRAN_ANYD(phi[mfi]),
			  BL_TO_FORTRAN_ANYD(enthalpy[mfi]),
			  AMREX_ZFILL(dx), omega, h_max,
			  &kin_eng, &pot_eng, &int_eng, &mass,
                          &Linf_norm);

     }

     ParallelDescriptor::ReduceRealSum(kin_eng);
     ParallelDescriptor::ReduceRealSum(pot_eng);
     ParallelDescriptor::ReduceRealSum(int_eng);
     ParallelDescriptor::ReduceRealSum(mass);
     ParallelDescriptor::ReduceRealMax(Linf_norm);



     // Now check to see if we're converged.

     scf_check_convergence(kin_eng, pot_eng, int_eng, mass,
			   Linf_norm,
			   &is_relaxed, j);

     //	for (int k = finest_level-1; k >= 0; k--)
     //	  getLevel(k).avgDown();

    gravity->multilevel_solve_for_new_phi(level,finest_level);

    if (is_relaxed == 1) break;

    j++;

  }

  for (int k = 0; k <= parent->finestLevel(); k++)
  {
     MultiFab& grav_new = getLevel(k).get_new_data(Gravity_Type);
     gravity->get_new_grav_vector(k,grav_new,time);
  }

}
#endif
#endif
