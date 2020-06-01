#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

  void problem_checkpoint(int * int_dir_name, int * len);

  void problem_restart(int * int_dir_name, int * len);

  void check_stopping_criteria(const int* lo, const int* hi,
                               const amrex::Real* v, const int* v_lo, const int* v_hi,
                               const amrex::Real* T, const int* t_lo, const int* t_hi,
                               const amrex::Real* ts_te, const int* s_lo, const int* s_hi,
                               const amrex::Real T_stopping_criterion,
                               const amrex::Real ts_te_stopping_criterion,
                               int* to_stop);

#ifdef __cplusplus
}
#endif

