#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void problem_checkpoint(int* int_dir_name, int* len);

  void problem_restart(int* int_dir_name, int* len);

  void wdcom(BL_FORT_FAB_ARG_3D(rho),
	     BL_FORT_FAB_ARG_3D(xmom),
	     BL_FORT_FAB_ARG_3D(ymom),
	     BL_FORT_FAB_ARG_3D(zmom),
	     BL_FORT_FAB_ARG_3D(pmask),
	     BL_FORT_FAB_ARG_3D(smask),
	     BL_FORT_FAB_ARG_3D(vol),
	     const int* lo, const int* hi,
	     const amrex::Real* dx, const amrex::Real time,
	     amrex::Real* com_p_x, amrex::Real* com_p_y, amrex::Real* com_p_z,
	     amrex::Real* com_s_x, amrex::Real* com_s_y, amrex::Real* com_s_z,
	     amrex::Real* vel_p_x, amrex::Real* vel_p_y, amrex::Real* vel_p_z,
	     amrex::Real* vel_s_x, amrex::Real* vel_s_y, amrex::Real* vel_s_z,
	     amrex::Real* m_p, amrex::Real* m_s);

  void ca_volumeindensityboundary(BL_FORT_FAB_ARG_3D(rho),
				  BL_FORT_FAB_ARG_3D(pmask),
				  BL_FORT_FAB_ARG_3D(smask),
				  BL_FORT_FAB_ARG_3D(vol),
				  const int* lo, const int* hi,
				  const amrex::Real* dx,
				  amrex::Real* volp, amrex::Real* vols,
				  amrex::Real rho_cutoff);

  void get_single_star(int* flag);

  void set_star_data(const amrex::Real* com_p, const amrex::Real* com_s,
		     const amrex::Real* vel_p, const amrex::Real* vel_s,
		     const amrex::Real* mass_p, const amrex::Real* mass_s,
		     const amrex::Real* t_ff_p, const amrex::Real* t_ff_s);

  void get_star_data(amrex::Real* com_p, amrex::Real* com_s,
		     amrex::Real* vel_p, amrex::Real* vel_s,
		     amrex::Real* mass_p, amrex::Real* mass_s,
		     amrex::Real* t_ff_p, amrex::Real* t_ff_s);

  void quadrupole_tensor_double_dot(BL_FORT_FAB_ARG_3D(rho),
				    BL_FORT_FAB_ARG_3D(xmom),
				    BL_FORT_FAB_ARG_3D(ymom),
				    BL_FORT_FAB_ARG_3D(zmom),
				    BL_FORT_FAB_ARG_3D(gx),
				    BL_FORT_FAB_ARG_3D(gy),
				    BL_FORT_FAB_ARG_3D(gz),
				    BL_FORT_FAB_ARG_3D(vol),
				    const int* lo, const int* hi,
				    const amrex::Real* dx, const amrex::Real time,
				    amrex::Real* Qtt);

  void gw_strain_tensor(amrex::Real* h_plus_1, amrex::Real* h_cross_1,
			amrex::Real* h_plus_2, amrex::Real* h_cross_2,
			amrex::Real* h_plus_3, amrex::Real* h_cross_3,
			amrex::Real* Qtt, const amrex::Real* time);

  void get_period(amrex::Real* period);
  void set_period(const amrex::Real* period);

  void update_center(const amrex::Real* time);

  void get_omega_vec(const amrex::Real* omega, const amrex::Real time);

  void get_lagrange_points(const amrex::Real mass_1, const amrex::Real mass_2,
	                   const amrex::Real* com_1, const amrex::Real* com_2,
			   amrex::Real* L1, amrex::Real* L2, amrex::Real* L3);

  void get_critical_roche_potential(BL_FORT_FAB_ARG_3D(phieff),
				    const int* lo, const int* hi,
				    amrex::Real* L1, amrex::Real* potential);

  void transform_to_inertial_frame(BL_FORT_FAB_ARG_3D(state),
				   const int* lo, const int* hi,
				   const amrex::Real* time);

  void get_relaxation_cutoff_time(amrex::Real* relaxation_cutoff_time);

  void check_relaxation(BL_FORT_FAB_ARG_3D(state),
			BL_FORT_FAB_ARG_3D(phieff),
			const int* lo, const int* hi,
			const amrex::Real* potential, int* is_done);

  void set_relaxation_damping_factor(const amrex::Real factor);

  void sum_force_on_stars(const int* lo, const int* hi,
                          BL_FORT_FAB_ARG_3D(force),
                          BL_FORT_FAB_ARG_3D(state),
                          BL_FORT_FAB_ARG_3D(vol),
                          BL_FORT_FAB_ARG_3D(pmask),
                          BL_FORT_FAB_ARG_3D(smask),
                          amrex::Real* fpx, amrex::Real* fpy, amrex::Real* fpz,
                          amrex::Real* fsx, amrex::Real* fsy, amrex::Real* fsz);

  void get_inertial_velocity(const amrex::Real* loc, const amrex::Real* vel, const amrex::Real* time, const amrex::Real* inertial_vel);

  void get_axes(int* axis_1, int* axis_2, int* axis_3);

  void get_problem_number(int* problem);

  void get_frame_choice(int* accurate_IC_frame);

  void set_extrema(const amrex::Real* T_max, const amrex::Real* rho_max, const amrex::Real* ts_te_max);

  void get_extrema(amrex::Real* T_max, amrex::Real* rho_max, amrex::Real* ts_te_max);

  void get_job_status(int* job_status);

  void set_job_status(int* job_status);

  void get_relaxation_status(int* relaxation_status);

  void set_relaxation_status(int* relaxation_status);

  void get_total_ener_array(const amrex::Real* ener_array);

  void set_total_ener_array(const amrex::Real* const ener_array);

#ifdef __cplusplus
}
#endif
