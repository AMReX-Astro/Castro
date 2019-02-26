extern "C"
{

  void scf_setup_relaxation(const Real* dx);

  void scf_get_coeff_info(int* loc_A, int* loc_B, int* loc_C,
			  Real* c_A, Real* c_B, Real* c_C);

  void scf_get_omegasq(const int* lo, const int* hi,
		       const int* domlo, const int* domhi,
		       BL_FORT_FAB_ARG_3D(state),
		       BL_FORT_FAB_ARG_3D(phi),
		       const Real* dx,  
		       const Real* time, Real* omegasq);

  void scf_get_bernoulli_const(const int* lo, const int* hi,
			       const int* domlo, const int* domhi,
			       BL_FORT_FAB_ARG_3D(state),
			       BL_FORT_FAB_ARG_3D(phi),
			       const Real* dx, const Real* time,
			       Real* bernoulli_1, Real* bernoulli_2);

  void scf_construct_enthalpy(const int* lo, const int* hi,
			      const int* domlo, const int* domhi,
			      BL_FORT_FAB_ARG_3D(state),
			      BL_FORT_FAB_ARG_3D(phi),
			      BL_FORT_FAB_ARG_3D(enthalpy),
			      const Real* dx,
			      const Real* time,
			      Real* bernoulli_1, Real* bernoulli_2,
			      Real* h_max_1, Real* h_max_2);

  void scf_update_density(const int* lo, const int* hi,
			  const int* domlo, const int* domhi,
			  BL_FORT_FAB_ARG_3D(state),
			  BL_FORT_FAB_ARG_3D(phi),
			  BL_FORT_FAB_ARG_3D(enthalpy),
			  const Real* dx,
			  const Real* time,
			  Real* h_max_1, Real* h_max_2,
			  Real* kin_eng, Real* pot_eng, Real* int_eng,
			  Real* left_mass, Real* right_mass,
			  Real* delta_rho, Real* l2_norm_resid, Real* l2_norm_source);

  void scf_check_convergence(Real* kin_eng, Real* pot_eng, Real* int_eng,
			     Real* left_mass, Real* right_mass,
			     Real* delta_rho, Real* l2_norm,
			     int* is_relaxed, int* num_iterations);

}
