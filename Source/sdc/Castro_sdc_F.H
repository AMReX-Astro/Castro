#ifndef AMREX_USE_CUDA
  void ca_sdc_update_advection_o2_lobatto
    (const int* lo, const int* hi,
     const amrex::Real* dt_m,
     const amrex::Real* dt,
     const BL_FORT_FAB_ARG_3D(k_m),
     BL_FORT_FAB_ARG_3D(k_n),
     BL_FORT_FAB_ARG_3D(A_m),
     BL_FORT_FAB_ARG_3D(A_0_old),
     BL_FORT_FAB_ARG_3D(A_1_old),
     const int* m_start);

  void ca_sdc_update_advection_o2_radau
    (const int* lo, const int* hi,
     const amrex::Real* dt_m,
     const amrex::Real* dt,
     const BL_FORT_FAB_ARG_3D(k_m),
     BL_FORT_FAB_ARG_3D(k_n),
     BL_FORT_FAB_ARG_3D(A_m),
     BL_FORT_FAB_ARG_3D(A_0_old),
     BL_FORT_FAB_ARG_3D(A_1_old),
     BL_FORT_FAB_ARG_3D(A_2_old),
     const int* m_start);

  void ca_sdc_update_advection_o4_lobatto
    (const int* lo, const int* hi,
     const amrex::Real* dt_m,
     const amrex::Real* dt,
     const BL_FORT_FAB_ARG_3D(k_m),
     BL_FORT_FAB_ARG_3D(k_n),
     BL_FORT_FAB_ARG_3D(A_m),
     BL_FORT_FAB_ARG_3D(A_0_old),
     BL_FORT_FAB_ARG_3D(A_1_old),
     BL_FORT_FAB_ARG_3D(A_2_old),
     const int* m_start);

  void ca_sdc_update_advection_o4_radau
    (const int* lo, const int* hi,
     const amrex::Real* dt_m,
     const amrex::Real* dt,
     const BL_FORT_FAB_ARG_3D(k_m),
     BL_FORT_FAB_ARG_3D(k_n),
     BL_FORT_FAB_ARG_3D(A_m),
     BL_FORT_FAB_ARG_3D(A_0_old),
     BL_FORT_FAB_ARG_3D(A_1_old),
     BL_FORT_FAB_ARG_3D(A_2_old),
     BL_FORT_FAB_ARG_3D(A_3_old),
     const int* m_start);

#ifdef REACTIONS
  void ca_sdc_compute_C2_lobatto
    (const int* lo, const int* hi,
     const amrex::Real* dt_m,
     const amrex::Real* dt,
     const BL_FORT_FAB_ARG_3D(A_m),
     const BL_FORT_FAB_ARG_3D(A_0_old),
     const BL_FORT_FAB_ARG_3D(A_1_old),
     const BL_FORT_FAB_ARG_3D(R_0_old),
     const BL_FORT_FAB_ARG_3D(R_1_old),
     BL_FORT_FAB_ARG_3D(C2),
     const int* m_start);

  void ca_sdc_compute_C2_radau
    (const int* lo, const int* hi,
     const amrex::Real* dt_m,
     const amrex::Real* dt,
     const BL_FORT_FAB_ARG_3D(A_m),
     const BL_FORT_FAB_ARG_3D(A_0_old),
     const BL_FORT_FAB_ARG_3D(A_1_old),
     const BL_FORT_FAB_ARG_3D(A_2_old),
     const BL_FORT_FAB_ARG_3D(R_0_old),
     const BL_FORT_FAB_ARG_3D(R_1_old),
     const BL_FORT_FAB_ARG_3D(R_2_old),
     BL_FORT_FAB_ARG_3D(C2),
     const int* m_start);

  void ca_sdc_update_o2
    (const int* lo, const int* hi,
     const amrex::Real* dt_m,
     const BL_FORT_FAB_ARG_3D(k_m),
     BL_FORT_FAB_ARG_3D(k_n),
     const BL_FORT_FAB_ARG_3D(A_m),
     const BL_FORT_FAB_ARG_3D(R_m_old),
     const BL_FORT_FAB_ARG_3D(C2),
     const int* sdc_iteration,
     const int* m_start);

  void ca_sdc_update_centers_o4
    (const int* lo, const int* hi,
     const amrex::Real* time,
     const BL_FORT_FAB_ARG_3D(U_center),
     BL_FORT_FAB_ARG_3D(U_new_center),
     const BL_FORT_FAB_ARG_3D(C_center),
     const int* sdc_iteration);

  void ca_sdc_conservative_update
    (const int* lo, const int* hi,
     const amrex::Real* time,
     const BL_FORT_FAB_ARG_3D(U_old),
     BL_FORT_FAB_ARG_3D(U_new),
     const BL_FORT_FAB_ARG_3D(C),
     const BL_FORT_FAB_ARG_3D(R_new));

  void ca_sdc_compute_C4_lobatto(const int* lo, const int* hi,
                                 const amrex::Real* dt_m,
                                 const amrex::Real* dt,
                                 const BL_FORT_FAB_ARG_3D(A_m),
                                 const BL_FORT_FAB_ARG_3D(A_0),
                                 const BL_FORT_FAB_ARG_3D(A_1),
                                 const BL_FORT_FAB_ARG_3D(A_2),
                                 const BL_FORT_FAB_ARG_3D(R_0),
                                 const BL_FORT_FAB_ARG_3D(R_1),
                                 const BL_FORT_FAB_ARG_3D(R_2),
                                 BL_FORT_FAB_ARG_3D(C),
                                 const int* m_start);

  void ca_sdc_compute_C4_radau(const int* lo, const int* hi,
                               const amrex::Real* dt_m,
                               const amrex::Real* dt,
                               const BL_FORT_FAB_ARG_3D(A_m),
                               const BL_FORT_FAB_ARG_3D(A_0),
                               const BL_FORT_FAB_ARG_3D(A_1),
                               const BL_FORT_FAB_ARG_3D(A_2),
                               const BL_FORT_FAB_ARG_3D(A_3),
                               const BL_FORT_FAB_ARG_3D(R_0),
                               const BL_FORT_FAB_ARG_3D(R_1),
                               const BL_FORT_FAB_ARG_3D(R_2),
                               const BL_FORT_FAB_ARG_3D(R_3),
                               BL_FORT_FAB_ARG_3D(C),
                               const int* m_start);

  void ca_sdc_compute_initial_guess(const int* lo, const int* hi,
                                    const BL_FORT_FAB_ARG_3D(U_old),
                                    const BL_FORT_FAB_ARG_3D(U_new),
                                    const BL_FORT_FAB_ARG_3D(A_old),
                                    const BL_FORT_FAB_ARG_3D(R_old),
                                    BL_FORT_FAB_ARG_3D(U_guess),
                                    const amrex::Real* dt_m, const int* sdc_iteration);
#endif
#endif // CUDA
