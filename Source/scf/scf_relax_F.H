
  void scf_setup_relaxation(const amrex::Real* dx);

  void scf_get_omegasq(const int* lo, const int* hi,
		       const int* domlo, const int* domhi,
		       BL_FORT_FAB_ARG_3D(state),
		       BL_FORT_FAB_ARG_3D(phi),
		       const amrex::Real* dx,  
		       const amrex::Real* time, amrex::Real* omegasq);

  void scf_get_bernoulli_const(const int* lo, const int* hi,
			       const int* domlo, const int* domhi,
			       BL_FORT_FAB_ARG_3D(state),
			       BL_FORT_FAB_ARG_3D(phi),
			       const amrex::Real* dx, const amrex::Real* time,
			       amrex::Real* bernoulli);

  void scf_construct_enthalpy(const int* lo, const int* hi,
			      const int* domlo, const int* domhi,
			      BL_FORT_FAB_ARG_3D(state),
			      BL_FORT_FAB_ARG_3D(phi),
			      BL_FORT_FAB_ARG_3D(enthalpy),
			      const amrex::Real* dx,
			      const amrex::Real* time,
			      amrex::Real* bernoulli, amrex::Real* h_max);

  void scf_update_density(const int* lo, const int* hi,
			  const int* domlo, const int* domhi,
			  BL_FORT_FAB_ARG_3D(state),
			  BL_FORT_FAB_ARG_3D(phi),
			  BL_FORT_FAB_ARG_3D(enthalpy),
			  const amrex::Real* dx,
			  const amrex::Real* time,
			  amrex::Real* h_max,
			  amrex::Real* kin_eng, amrex::Real* pot_eng, amrex::Real* int_eng,
			  amrex::Real* mass,
			  amrex::Real* delta_rho, amrex::Real* l2_norm_resid, amrex::Real* l2_norm_source);

  void scf_check_convergence(amrex::Real* kin_eng, amrex::Real* pot_eng, amrex::Real* int_eng,
			     amrex::Real* mass,
			     amrex::Real* delta_rho, amrex::Real* l2_norm,
			     int* is_relaxed, int* num_iterations);
