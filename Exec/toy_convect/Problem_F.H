BL_FORT_PROC_DECL(PROBLEM_CHECKPOINT,problem_checkpoint)(int * int_dir_name, int * len);

BL_FORT_PROC_DECL(PROBLEM_RESTART,problem_restart)(int * int_dir_name, int * len);



BL_FORT_PROC_DECL(CA_STATE_ERROR,ca_state_error)
  (int* tag, const int* tag_lo, const int* tag_hi, 
   const int* tagval, const int* clearval,
   BL_FORT_FAB_ARG_3D(state),
   const int* lo, const int* hi,
   const int* ncomp,
   const int* domlo, const int* domhi,
   const Real* dx, const Real* xlo, const Real* problo,
   const Real* time, const int* level);

