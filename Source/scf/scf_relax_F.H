#ifndef _scf_relax_F_H_
#define _scf_relax_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void scf_setup_relaxation(amrex::Real d_A, int axis_A, amrex::Real d_B, int axis_B);

  void scf_calculate_target_h_max(const int* lo, const int* hi,
                                  BL_FORT_FAB_ARG_3D(state),
                                  amrex::Real maximum_density,
                                  amrex::Real* target_h_max);

  void scf_update_for_omegasq(const int* lo, const int* hi,
                              BL_FORT_FAB_ARG_3D(phi),
                              BL_FORT_FAB_ARG_3D(psi),
                              const amrex::Real* dx,  
                              amrex::Real* phi_A, amrex::Real* psi_A,
                              amrex::Real* phi_B, amrex::Real* psi_B);

  void scf_get_bernoulli_const(const int* lo, const int* hi,
                               BL_FORT_FAB_ARG_3D(phi),
                               BL_FORT_FAB_ARG_3D(psi),
                               const amrex::Real* dx,
                               amrex::Real* bernoulli);

  void scf_construct_enthalpy(const int* lo, const int* hi,
                              BL_FORT_FAB_ARG_3D(phi),
                              BL_FORT_FAB_ARG_3D(psi),
                              BL_FORT_FAB_ARG_3D(enthalpy),
                              const amrex::Real* dx, amrex::Real bernoulli);

  void scf_update_density(const int* lo, const int* hi,
                          BL_FORT_FAB_ARG_3D(state),
                          BL_FORT_FAB_ARG_3D(enthalpy),
                          const amrex::Real* dx, amrex::Real actual_rho_max,
                          amrex::Real actual_h_max, amrex::Real target_h_max,
                          amrex::Real* Linf_norm);

  void scf_diagnostics(const int* lo, const int* hi,
                       BL_FORT_FAB_ARG_3D(state),
                       BL_FORT_FAB_ARG_3D(phi),
                       BL_FORT_FAB_ARG_3D(psi),
                       const amrex::Real* dx,
                       amrex::Real* kin_eng, amrex::Real* pot_eng,
                       amrex::Real* int_eng, amrex::Real* mass);

  void scf_get_solar_mass(amrex::Real* M_solar);

#ifdef __cplusplus
}
#endif

#endif
