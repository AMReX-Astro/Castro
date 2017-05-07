
#ifdef __cplusplus
extern "C"
{
#endif

    void problem_checkpoint(int * int_dir_name, int * len);

    void problem_restart(int * int_dir_name, int * len);

    void get_problem_number(int* problem);

    void get_diameter(amrex::Real* diameter);

    void get_density(amrex::Real* density);

    void update_density(const int* lo, const int* hi, const amrex::Real* dx,
			amrex::Real* state, const int* s_lo, const int* s_hi,
			const amrex::Real* update_factor);

#ifdef __cplusplus
}
#endif
