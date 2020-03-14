#ifndef _Castro_diffusion_F_H_
#define _Castro_diffusion_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void ca_fill_temp_cond
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(coeff_cc));

#ifdef __cplusplus
}
#endif

#endif
