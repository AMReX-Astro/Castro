#ifndef _Castro_mhd_F_H_
#define _Castro_mhd_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void ctoprim_mhd
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(bcc),
     const BL_FORT_FAB_ARG_3D(Bx),
     const BL_FORT_FAB_ARG_3D(By),
     const BL_FORT_FAB_ARG_3D(Bz),
     BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(cx),
     BL_FORT_FAB_ARG_3D(cy),
     BL_FORT_FAB_ARG_3D(cz));

  void check_for_mhd_cfl_violation
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG_3D(q),
     const BL_FORT_FAB_ARG_3D(cx),
     const BL_FORT_FAB_ARG_3D(cy),
     const BL_FORT_FAB_ARG_3D(cz),
     amrex::Real& courno, const amrex::Real dx[], const amrex::Real dt);

  void mhd_flatten
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(flatn));

  void plm
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG_3D(s),
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

  void consup
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG_3D(uin),
     BL_FORT_FAB_ARG_3D(uout),
     BL_FORT_FAB_ARG_3D(bcc),
     const BL_FORT_FAB_ARG_3D(fluxx),
     const BL_FORT_FAB_ARG_3D(fluxy),
     const BL_FORT_FAB_ARG_3D(fluxz),
     const amrex::Real dx[], const amrex::Real dt);

  void magup
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG_3D(bxin),
     const BL_FORT_FAB_ARG_3D(byin),
     const BL_FORT_FAB_ARG_3D(bzin),
     BL_FORT_FAB_ARG_3D(bxout),
     BL_FORT_FAB_ARG_3D(byout),
     BL_FORT_FAB_ARG_3D(bzout),
     const BL_FORT_FAB_ARG_3D(Ex),
     const BL_FORT_FAB_ARG_3D(Ey),
     const BL_FORT_FAB_ARG_3D(Ez),
     const amrex::Real* dx, const amrex::Real dt);

#ifdef __cplusplus
}
#endif

#endif

