#ifndef CASTRO_MHD_F_H
#define CASTRO_MHD_F_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void hlld(const int lo[], const int hi[],
            const BL_FORT_FAB_ARG_3D(qleft),
            const BL_FORT_FAB_ARG_3D(qright),
            BL_FORT_FAB_ARG_3D(flx),
            const int dir);


#ifdef __cplusplus
}
#endif

#endif

