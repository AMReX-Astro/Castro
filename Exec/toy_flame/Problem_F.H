/* problem-specific stuff goes here */

BL_FORT_PROC_DECL(PROBLEM_CHECKPOINT,problem_checkpoint)(int * int_dir_name, int * len);

BL_FORT_PROC_DECL(PROBLEM_RESTART,problem_restart)(int * int_dir_name, int * len);

BL_FORT_PROC_DECL(FLAME_WIDTH_TEMP,flame_width_temp)
    (BL_FORT_FAB_ARG_3D(temp), 
     const int* lo, const int* hi,
     const Real* dx, const Real* time,
     Real* T_max, Real* T_min, Real* grad_T_max);
