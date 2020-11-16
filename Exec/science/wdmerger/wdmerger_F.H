#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void set_star_data(const amrex::Real* com_p, const amrex::Real* com_s,
		     const amrex::Real* vel_p, const amrex::Real* vel_s,
		     const amrex::Real* mass_p, const amrex::Real* mass_s,
		     const amrex::Real* t_ff_p, const amrex::Real* t_ff_s);

  void get_period(amrex::Real* period);

  void update_center(const amrex::Real* time);

  void get_relaxation_density_cutoff(amrex::Real* relaxation_density_cutoff);

  void get_relaxation_cutoff_time(amrex::Real* relaxation_cutoff_time);

  void set_relaxation_damping_factor(const amrex::Real factor);

  void get_axes(int* axis_1, int* axis_2, int* axis_3);

  void get_relaxation_status(int* relaxation_status);

  void set_relaxation_status(int* relaxation_status);

#ifdef __cplusplus
}
#endif
