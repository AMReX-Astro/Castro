#ifndef CASTRO_MHD_F_H
#define CASTRO_MHD_F_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void ctoprim_mhd
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG_3D(state),
     const BL_FORT_FAB_ARG_3D(Bx),
     const BL_FORT_FAB_ARG_3D(By),
     const BL_FORT_FAB_ARG_3D(Bz),
     BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(cx),
     BL_FORT_FAB_ARG_3D(cy),
     BL_FORT_FAB_ARG_3D(cz));

  void mhd_flatten
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(flatn));

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

  void corner_transport
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG_3D(q),
     const BL_FORT_FAB_ARG_3D(qm),
     const BL_FORT_FAB_ARG_3D(qp),
     BL_FORT_FAB_ARG_3D(flxx),
     BL_FORT_FAB_ARG_3D(flxy),
     BL_FORT_FAB_ARG_3D(flxz),
     BL_FORT_FAB_ARG_3D(Ex),
     BL_FORT_FAB_ARG_3D(Ey),
     BL_FORT_FAB_ARG_3D(Ez),
     const amrex::Real dx[], const amrex::Real dt);

#ifdef __cplusplus
}
#endif

#endif

