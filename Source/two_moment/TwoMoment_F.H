#ifndef _TwoMoment_F_H_
#define _TwoMoment_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
  void InitThornado(const int* nDimsX, const int* nDimsE);

  void InitThornado_Patch(const int* len, int* swX,
                          const amrex::Real* grid_lo, const amrex::Real* grid_hi, 
                          const int* n_energy, 
                          const int* swE, 
                          const amrex::Real* eL, 
                          const amrex::Real* eR, 
                          const int* n_species); 

  void FreeThornado_Patch();

  void call_to_thornado(const int lo[], const int hi[], amrex::Real* dt_sub,
                        amrex::Real* U_F, 
                        BL_FORT_FAB_ARG_3D(dU_F), const int* nc,
                        amrex::Real* U_R_o, 
                        BL_FORT_FAB_ARG_3D( U_R_n), const int* nr,
                        int* n_fluid_dof, int* n_energy, int* n_species, 
                        int* n_rad_dof, int* n_moments);

#ifdef __cplusplus
}
#endif

#endif
