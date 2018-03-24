#ifndef _TwoMoment_F_H_
#define _TwoMoment_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
  void init_thornado(const int lo[], const int hi[],
                     BL_FORT_FAB_ARG_3D( U_F), const int* nc_uf,
                     BL_FORT_FAB_ARG_3D( U_R), const int* nc ur);

  void call_to_thornado(const int lo[], const int hi[],
                        BL_FORT_FAB_ARG_3D( U_F), const int* nc_uf,
                        amrex::Real* U_R_o, 
                        BL_FORT_FAB_ARG_3D( U_R_n), const int* nc ur,
                        BL_FORT_FAB_ARG_3D(dU_F), const int* nc_duf,
                        amrex::Real* dt_sub);

#ifdef __cplusplus
}
#endif

#endif
