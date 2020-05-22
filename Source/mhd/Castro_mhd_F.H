#ifndef CASTRO_MHD_F_H
#define CASTRO_MHD_F_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void plm
    (const int lo[], const int hi[],
     const int idir,
     const BL_FORT_FAB_ARG_3D(s),
     const BL_FORT_FAB_ARG_3D(qaux),
     const BL_FORT_FAB_ARG_3D(flatn),
     const BL_FORT_FAB_ARG_3D(bx),
     const BL_FORT_FAB_ARG_3D(by),
     const BL_FORT_FAB_ARG_3D(bz),
     BL_FORT_FAB_ARG_3D(Ip),
     BL_FORT_FAB_ARG_3D(Im),
     const BL_FORT_FAB_ARG_3D(srcq),
     const amrex::Real dx[], const amrex::Real dt);

  void hlld(const int lo[], const int hi[],
            const BL_FORT_FAB_ARG_3D(qleft),
            const BL_FORT_FAB_ARG_3D(qright),
            BL_FORT_FAB_ARG_3D(flx),
            const int dir);


#ifdef __cplusplus
}
#endif

#endif

