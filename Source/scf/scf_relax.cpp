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

  const Real* dx   = parent->Geom(level).CellSize();
  
  const int* domlo = geom.Domain().loVect();
  const int* domhi = geom.Domain().hiVect();

  scf_setup_relaxation(dx);

  // Get the phi MultiFab.

  MultiFab& phi = get_new_data(PhiGrav_Type);

  Real time = state[State_Type].curTime();

  // Iterate until the system is relaxed by filling the level data 
  // and then doing a multilevel gravity solve.

  int is_relaxed = 0;

  while ( j <= relax_max_iterations ) {

     // First step is to find the rotational frequency.

     Real phi_A = 0.0;
     Real psi_A = 0.0;
     Real phi_B = 0.0;
     Real psi_B = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:phi_A, psi_A, phi_B, psi_B)
#endif    	
     for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi) {

       const Box& box  = mfi.tilebox();

       scf_update_for_omegasq(AMREX_ARLIM_ANYD(box.loVect()), AMREX_ARLIM_ANYD(box.hiVect()),
                              AMREX_ARLIM_ANYD(domlo), AMREX_ARLIM_ANYD(domhi),
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

     Real omegasq = -(phi_A - phi_B) / (psi_A - psi_B);

     if (omegasq < 0.0 && ParallelDescriptor::IOProcessor()) {
	 std::cerr << "Omega squared = " << omegasq << " is negative in the relaxation step; aborting." << std::endl;
	 amrex::Error();
     }

     // Rotational period is 2 pi / omega.

     rotational_period = 2.0 * M_PI / sqrt(omegasq);

     // Now save the updated rotational frequency in the Fortran module.

     set_rot_period(&rotational_period);



     // Second step is to evaluate the Bernoulli constant.

     Real bernoulli = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:bernoulli)
#endif
     for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {

       const Box& box  = mfi.tilebox();
       const int* lo   = box.loVect();
       const int* hi   = box.hiVect();

       Real b = 0.0;

       scf_get_bernoulli_const(ARLIM_3D(lo), ARLIM_3D(hi),
			       ARLIM_3D(domlo), ARLIM_3D(domhi),
			       BL_TO_FORTRAN_3D(S_new[mfi]),
			       BL_TO_FORTRAN_3D(phi[mfi]),
			       ZFILL(dx), &time, &b);

       bernoulli += b;

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

       const Box& box  = mfi.tilebox();
       const int* lo   = box.loVect();
       const int* hi   = box.hiVect();

       Real h = 0.0;

       scf_construct_enthalpy(ARLIM_3D(lo), ARLIM_3D(hi),
			      ARLIM_3D(domlo), ARLIM_3D(domhi),
			      BL_TO_FORTRAN_3D(S_new[mfi]),
			      BL_TO_FORTRAN_3D(phi[mfi]),
			      BL_TO_FORTRAN_3D(enthalpy[mfi]),
			      ZFILL(dx), &time,
			      &bernoulli, &h);

       if (h > h_max) h_max = h;

     }

     ParallelDescriptor::ReduceRealMax(h_max);

     Real kin_eng = 0.0;
     Real pot_eng = 0.0;
     Real int_eng = 0.0;
     Real l2_norm_resid = 0.0;
     Real l2_norm_source = 0.0;
     Real mass = 0.0;
     Real delta_rho = 0.0;

     // Finally, update the density using the enthalpy field.

#ifdef _OPENMP
#pragma omp parallel reduction(+:kin_eng,pot_eng,int_eng)      \
                     reduction(+:l2_norm_resid,l2_norm_source) \
                     reduction(+:mass)                         \
                     reduction(max:delta_rho)
#endif
     for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {

       const Box& box  = mfi.tilebox();
       const int* lo   = box.loVect();
       const int* hi   = box.hiVect();

       Real ke = 0.0;
       Real pe = 0.0;
       Real ie = 0.0;
       Real m  = 0.0;
       Real dr = 0.0;
       Real nr = 0.0;
       Real ns = 0.0;

       scf_update_density(ARLIM_3D(lo), ARLIM_3D(hi),
			  ARLIM_3D(domlo), ARLIM_3D(domhi),
			  BL_TO_FORTRAN_3D(S_new[mfi]),
			  BL_TO_FORTRAN_3D(phi[mfi]),
			  BL_TO_FORTRAN_3D(enthalpy[mfi]),
			  ZFILL(dx), &time, 
			  &h_max,
			  &ke, &pe, &ie,
			  &m,
			  &dr, &nr, &ns);

       kin_eng += ke;
       pot_eng += pe;
       int_eng += ie;
       l2_norm_resid += nr;
       l2_norm_source += ns;
       mass += m;
       if (dr > delta_rho) delta_rho = dr;

     }

     ParallelDescriptor::ReduceRealSum(kin_eng);
     ParallelDescriptor::ReduceRealSum(pot_eng);
     ParallelDescriptor::ReduceRealSum(int_eng);

     ParallelDescriptor::ReduceRealSum(l2_norm_resid);
     ParallelDescriptor::ReduceRealSum(l2_norm_source);

     Real l2_norm = l2_norm_resid / l2_norm_source;

     ParallelDescriptor::ReduceRealMax(delta_rho);

     ParallelDescriptor::ReduceRealSum(mass);



     // Now check to see if we're converged.

     scf_check_convergence(&kin_eng, &pot_eng, &int_eng, 
			   &mass,
			   &delta_rho, &l2_norm,
			   &is_relaxed, &j);

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
