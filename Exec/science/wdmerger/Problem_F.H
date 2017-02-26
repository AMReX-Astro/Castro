extern "C"
{

  void problem_checkpoint(int * int_dir_name, int * len);

  void problem_restart(int * int_dir_name, int * len);

  void wdcom(BL_FORT_FAB_ARG_3D(rho), 
	     BL_FORT_FAB_ARG_3D(xmom), 
	     BL_FORT_FAB_ARG_3D(ymom), 
	     BL_FORT_FAB_ARG_3D(zmom),
	     BL_FORT_FAB_ARG_3D(pmask),
	     BL_FORT_FAB_ARG_3D(smask),
	     BL_FORT_FAB_ARG_3D(vol),
	     const int* lo, const int* hi,
	     const Real* dx, const Real* time,
	     Real* com_p_x, Real* com_p_y, Real* com_p_z,
	     Real* com_s_x, Real* com_s_y, Real* com_s_z,
	     Real* vel_p_x, Real* vel_p_y, Real* vel_p_z,
	     Real* vel_s_x, Real* vel_s_y, Real* vel_s_z,
	     Real* m_p, Real* m_s);

  void ca_volumeindensityboundary(BL_FORT_FAB_ARG_3D(rho),
				  BL_FORT_FAB_ARG_3D(pmask),
				  BL_FORT_FAB_ARG_3D(smask),
				  BL_FORT_FAB_ARG_3D(vol),
				  const int* lo, const int* hi,
				  const Real* dx,  
				  Real* volp, Real* vols,
				  Real* rho_cutoff);

  void get_single_star(int& flag);

  void set_star_data(const Real com_p[], const Real com_s[], 
		     const Real vel_p[], const Real vel_s[], 
		     const Real* mass_p, const Real* mass_s,
		     const Real* t_ff_p, const Real* t_ff_s);

  void get_star_data(Real com_p[], Real com_s[], 
		     Real vel_p[], Real vel_s[], 
		     Real* mass_p, Real* mass_s,
		     Real* t_ff_p, Real* t_ff_s);

  void quadrupole_tensor_double_dot(BL_FORT_FAB_ARG_3D(rho),
				    BL_FORT_FAB_ARG_3D(xmom),
				    BL_FORT_FAB_ARG_3D(ymom),
				    BL_FORT_FAB_ARG_3D(zmom),
				    BL_FORT_FAB_ARG_3D(gx),
				    BL_FORT_FAB_ARG_3D(gy),
				    BL_FORT_FAB_ARG_3D(gz),
				    BL_FORT_FAB_ARG_3D(vol),
				    const int* lo, const int* hi,
				    const Real* dx, const Real* time,
				    Real* Qtt);

  void gw_strain_tensor(Real* h_plus_1, Real* h_cross_1, 
			Real* h_plus_2, Real* h_cross_2,
			Real* h_plus_3, Real* h_cross_3,
			Real* Qtt, const Real* time);

  void get_period(Real* period);
  void set_period(const Real* period);

  void update_center(const Real* time);

  void get_omega_vec(const Real omega[], const Real& time);

  void get_lagrange_points(const Real& mass_1, const Real& mass_2,
	                   const Real* com_1, const Real* com_2,
			   Real* L1, Real* L2, Real* L3);

  void get_critical_roche_potential(BL_FORT_FAB_ARG_3D(phieff),
				    const int* lo, const int* hi,
				    Real* L1, Real* potential);

  void transform_to_inertial_frame(BL_FORT_FAB_ARG_3D(state),
				   const int* lo, const int* hi,
				   const Real* time);

  void check_relaxation(BL_FORT_FAB_ARG_3D(state),
			BL_FORT_FAB_ARG_3D(phieff),
			const int* lo, const int* hi,
			const Real* potential, int* is_done);

  void turn_off_relaxation(const Real* time);

  void sum_force_on_stars(const int* lo, const int* hi,
                          BL_FORT_FAB_ARG_3D(force),
                          BL_FORT_FAB_ARG_3D(state),
                          BL_FORT_FAB_ARG_3D(vol),
                          BL_FORT_FAB_ARG_3D(pmask),
                          BL_FORT_FAB_ARG_3D(smask),
                          Real* fpx, Real* fpy, Real* fpz,
                          Real* fsx, Real* fsy, Real* fsz);

  void get_axes(int& axis_1, int& axis_2, int& axis_3);

  void get_problem_number(int* problem);

  void get_frame_choice(int* accurate_IC_frame);

  void set_extrema(const Real* T_max, const Real* rho_max, const Real* ts_te_max);

  void get_job_status(int* job_status);

  void set_job_status(int* job_status);

  void get_relaxation_status(int* relaxation_status);

  void set_relaxation_status(int* relaxation_status);

  void get_total_ener_array(const Real* ener_array);

  void set_total_ener_array(const Real* const ener_array);

}
