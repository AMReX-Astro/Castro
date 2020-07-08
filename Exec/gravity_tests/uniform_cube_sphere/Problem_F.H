#include <AMReX_BLFort.H>

extern "C"
{

    void problem_checkpoint(int* int_dir_name, int* len);

    void problem_restart(int* int_dir_name, int* len);

}
