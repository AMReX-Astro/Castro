/* problem-specific stuff goes here */

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

void problem_checkpoint(int* int_dir_name, int* len);

void problem_restart(int* int_dir_name, int* len);

void flame_width_temp(BL_FORT_FAB_ARG_3D(temp), const int* lo, const int* hi,
                      const amrex::Real* dx, const amrex::Real* time,
                      amrex::Real* T_max, amrex::Real* T_min,
                      amrex::Real* grad_T_max);

void flame_speed_data(BL_FORT_FAB_ARG_3D(omegadot), const int* lo,
                      const int* hi, const amrex::Real* dx,
                      amrex::Real* rho_X_dot);

#ifdef __cplusplus
}
#endif
