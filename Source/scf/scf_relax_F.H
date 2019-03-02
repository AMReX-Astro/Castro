
  void scf_setup_relaxation();

  void scf_update_for_omegasq(const int* lo, const int* hi,
                              BL_FORT_FAB_ARG_3D(state),
                              BL_FORT_FAB_ARG_3D(phi),
                              const amrex::Real* dx,  
                              amrex::Real* phi_A, amrex::Real* psi_A,
                              amrex::Real* phi_B, amrex::Real* psi_B);

  void scf_get_bernoulli_const(const int* lo, const int* hi,
			       BL_FORT_FAB_ARG_3D(state),
			       BL_FORT_FAB_ARG_3D(phi),
			       const amrex::Real* dx, const amrex::Real omega,
			       amrex::Real* bernoulli);

  void scf_construct_enthalpy(const int* lo, const int* hi,
			      BL_FORT_FAB_ARG_3D(state),
			      BL_FORT_FAB_ARG_3D(phi),
			      BL_FORT_FAB_ARG_3D(enthalpy),
			      const amrex::Real* dx, const amrex::Real omega,
			      amrex::Real* bernoulli, amrex::Real* h_max);

  void scf_update_density(const int* lo, const int* hi,
			  BL_FORT_FAB_ARG_3D(state),
			  BL_FORT_FAB_ARG_3D(phi),
			  BL_FORT_FAB_ARG_3D(enthalpy),
			  const amrex::Real* dx, amrex::Real omega, amrex::Real h_max,
			  amrex::Real* kin_eng, amrex::Real* pot_eng, amrex::Real* int_eng, amrex::Real* mass,
			  amrex::Real* Linf_norm);

  void scf_check_convergence(amrex::Real kin_eng, amrex::Real pot_eng, amrex::Real int_eng, amrex::Real mass,
			     amrex::Real Linf_norm,
			     int* is_relaxed, int num_iterations);
