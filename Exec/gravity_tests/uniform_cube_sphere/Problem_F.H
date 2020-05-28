#include <AMReX_BLFort.H>

extern "C" {

void problem_checkpoint(int* int_dir_name, int* len);

void problem_restart(int* int_dir_name, int* len);

void get_problem_number(int* problem);

void get_diameter(amrex::Real* diameter);

void get_density(amrex::Real* density);

void set_density(amrex::Real density);

void update_density(const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state),
                    const amrex::Real* dx, const amrex::Real update_factor);

void compute_norm(const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(phi),
                  BL_FORT_FAB_ARG_3D(vol), const amrex::Real* dx,
                  const int norm_power, amrex::Real* norm_diff,
                  amrex::Real* norm_exact);
}
