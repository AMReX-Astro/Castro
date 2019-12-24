#include <AMReX_BLFort.H>

extern "C"
{

    void problem_checkpoint(int* int_dir_name, int* len);

    void problem_restart(int* int_dir_name, int* len);

    void get_problem_number(int* problem);

    void get_diameter(amrex::Real* diameter);

    void get_density(amrex::Real* density);

    void update_density(const int* lo, const int* hi,
                        BL_FORT_FAB_ARG_3D(state),
                        const amrex::Real* dx, const amrex::Real update_factor);

}
