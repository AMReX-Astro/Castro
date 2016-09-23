
BL_FORT_PROC_DECL(PROBLEM_CHECKPOINT,problem_checkpoint)
  (int * int_dir_name, int * len);

BL_FORT_PROC_DECL(PROBLEM_RESTART,problem_restart)
  (int * int_dir_name, int * len);

extern "C"
{

    void get_problem_number(int* problem);

    void get_diameter(Real* diameter);

    void get_density(Real* density);

    void update_density(const int* lo, const int* hi, const Real* dx,
			Real* state, const int* s_lo, const int* s_hi,
			const Real* update_factor);

}
