
extern "C"
{
    void problem_checkpoint(int * int_dir_name, int * len);

    void problem_restart(int * int_dir_name, int * len);

    void get_problem_number(int* problem);

    void get_diameter(Real* diameter);

    void get_density(Real* density);

    void update_density(const int* lo, const int* hi, const Real* dx,
			Real* state, const int* s_lo, const int* s_hi,
			const Real* update_factor);

}
