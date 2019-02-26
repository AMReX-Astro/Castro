#include "Castro.H"
#include "Castro_F.H"

#include "Gravity.H"

#ifdef GRAVITY
#ifdef ROTATION
#ifdef do_problem_post_init
void Castro::scf_relaxation() {

  int finest_level = parent->finestLevel();

  MultiFab& S_new = get_new_data(State_Type);

  int j = 1;
  int relax_max_iterations = 30;

  const Real* dx   = parent->Geom(level).CellSize();
  
  const int* domlo = geom.Domain().loVect();
  const int* domhi = geom.Domain().hiVect();

  MultiFab gcoeff(grids,1,0,Fab_allocate);
  gcoeff.setVal(0.0);

  Box box_A( IntVect(), IntVect(2, 2, 2) );
  Box box_B( IntVect(), IntVect(2, 2, 2) );
  Box box_C( IntVect(), IntVect(2, 2, 2) );

  FArrayBox cA(box_A);
  FArrayBox cB(box_B);
  FArrayBox cC(box_C);

  cA.setVal(0.0);
  cB.setVal(0.0);
  cC.setVal(0.0);

  int loc_A[3] = {0};
  int loc_B[3] = {0};
  int loc_C[3] = {0};

  scf_setup_relaxation(dx);

  scf_get_coeff_info(loc_A, loc_B, loc_C, cA.dataPtr(), cB.dataPtr(), cC.dataPtr());

  // Get the phi MultiFab.

  MultiFab& phi = get_new_data(PhiGrav_Type);

  Real time = state[State_Type].curTime();

  // Iterate until the system is relaxed by filling the level data 
  // and then doing a multilevel gravity solve.

  int is_relaxed = 0;

  while ( j <= relax_max_iterations ) {

     // First step is to find the rotational frequency.

     Real omegasq = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:omegasq)
#endif    	
     for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi) {

       const Box& box  = mfi.tilebox();
       const int* lo   = box.loVect();
       const int* hi   = box.hiVect();

       Real osq = 0.0;

       scf_get_omegasq(ARLIM_3D(lo), ARLIM_3D(hi),
		       ARLIM_3D(domlo), ARLIM_3D(domhi),
		       BL_TO_FORTRAN_3D(S_new[mfi]),
		       BL_TO_FORTRAN_3D(phi[mfi]),
		       ZFILL(dx), &time, &osq);

       omegasq += osq;

     }

     ParallelDescriptor::ReduceRealSum(omegasq);

     if (omegasq < 0.0 && ParallelDescriptor::IOProcessor()) {
	 std::cerr << "Omega squared = " << omegasq << " is negative in the relaxation step; aborting." << std::endl;
	 BoxLib::Error();
     }

     // Rotational period is 2 pi / omega.

     rotational_period = 2.0 * M_PI / sqrt(omegasq);

     // Now save the updated rotational frequency in the Fortran module.

     set_period(&rotational_period);



     // Second step is to evaluate the Bernoulli constants.

     Real bernoulli_1 = 0.0;
     Real bernoulli_2 = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:bernoulli_1,bernoulli_2)
#endif
     for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi) {

       const Box& box  = mfi.tilebox();
       const int* lo   = box.loVect();
       const int* hi   = box.hiVect();

       Real b1 = 0.0;
       Real b2 = 0.0;

       scf_get_bernoulli_const(ARLIM_3D(lo), ARLIM_3D(hi),
			       ARLIM_3D(domlo), ARLIM_3D(domhi),
			       BL_TO_FORTRAN_3D(S_new[mfi]),
			       BL_TO_FORTRAN_3D(phi[mfi]),
			       ZFILL(dx), &time, &b1, &b2);

       bernoulli_1 += b1;
       bernoulli_2 += b2;

     }

     ParallelDescriptor::ReduceRealSum(bernoulli_1);
     ParallelDescriptor::ReduceRealSum(bernoulli_2);



     // Third step is to construct the enthalpy field and 
     // find the maximum enthalpy for each star.

     Real h_max_1 = 0.0;
     Real h_max_2 = 0.0;

     // Define the enthalpy MultiFAB to have one component and zero ghost cells.

     MultiFab enthalpy(grids,1,0,Fab_allocate);
     enthalpy.setVal(0.0);

#ifdef _OPENMP
#pragma omp parallel reduction(max:h_max_1,h_max_2)
#endif    	
     for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi) {

       const Box& box  = mfi.tilebox();
       const int* lo   = box.loVect();
       const int* hi   = box.hiVect();

       Real h1 = 0.0;
       Real h2 = 0.0;

       scf_construct_enthalpy(ARLIM_3D(lo), ARLIM_3D(hi),
			      ARLIM_3D(domlo), ARLIM_3D(domhi),
			      BL_TO_FORTRAN_3D(S_new[mfi]),
			      BL_TO_FORTRAN_3D(phi[mfi]),
			      BL_TO_FORTRAN_3D(enthalpy[mfi]),
			      ZFILL(dx), &time,
			      &bernoulli_1, &bernoulli_2, &h1, &h2);

       if (h1 > h_max_1) h_max_1 = h1;
       if (h2 > h_max_2) h_max_2 = h2;

     }

     ParallelDescriptor::ReduceRealMax(h_max_1);
     ParallelDescriptor::ReduceRealMax(h_max_2);

     Real kin_eng = 0.0;
     Real pot_eng = 0.0;
     Real int_eng = 0.0;
     Real l2_norm_resid = 0.0;
     Real l2_norm_source = 0.0;
     Real left_mass = 0.0;
     Real right_mass = 0.0;
     Real delta_rho = 0.0;

     // Finally, update the density using the enthalpy field.

#ifdef _OPENMP
#pragma omp parallel reduction(+:kin_eng,pot_eng,int_eng)      \
		 reduction(+:l2_norm_resid,l2_norm_source) \
		 reduction(+:left_mass,right_mass)         \
		 reduction(max:delta_rho)
#endif
     for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi) {

       const Box& box  = mfi.tilebox();
       const int* lo   = box.loVect();
       const int* hi   = box.hiVect();

       Real ke = 0.0;
       Real pe = 0.0;
       Real ie = 0.0;
       Real lm = 0.0;
       Real rm = 0.0;
       Real dr = 0.0;
       Real nr = 0.0;
       Real ns = 0.0;

       scf_update_density(ARLIM_3D(lo), ARLIM_3D(hi),
			  ARLIM_3D(domlo), ARLIM_3D(domhi),
			  BL_TO_FORTRAN_3D(S_new[mfi]),
			  BL_TO_FORTRAN_3D(phi[mfi]),
			  BL_TO_FORTRAN_3D(enthalpy[mfi]),
			  ZFILL(dx), &time, 
			  &h_max_1, &h_max_2,
			  &ke, &pe, &ie,
			  &lm, &rm,
			  &dr, &nr, &ns);

       kin_eng += ke;
       pot_eng += pe;
       int_eng += ie;
       l2_norm_resid += nr;
       l2_norm_source += ns;
       left_mass += lm;
       right_mass += rm;
       if (dr > delta_rho) delta_rho = dr;

     }

     ParallelDescriptor::ReduceRealSum(kin_eng);
     ParallelDescriptor::ReduceRealSum(pot_eng);
     ParallelDescriptor::ReduceRealSum(int_eng);

     ParallelDescriptor::ReduceRealSum(l2_norm_resid);
     ParallelDescriptor::ReduceRealSum(l2_norm_source);

     Real l2_norm = l2_norm_resid / l2_norm_source;

     ParallelDescriptor::ReduceRealMax(delta_rho);

     ParallelDescriptor::ReduceRealSum(left_mass);
     ParallelDescriptor::ReduceRealSum(right_mass);



     // Now check to see if we're converged.

     scf_check_convergence(&kin_eng, &pot_eng, &int_eng, 
			   &left_mass, &right_mass,
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
#endif
