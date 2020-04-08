#ifndef _Castro_mhd_F_H_
#define _Castro_mhd_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void ca_advance_mhd
    (const amrex::Real* time,
     const int lo[], const int hi[],
     const BL_FORT_FAB_ARG_3D(state), BL_FORT_FAB_ARG_3D(stateout),
     const BL_FORT_FAB_ARG_3D(Bx),
     const BL_FORT_FAB_ARG_3D(By),
     const BL_FORT_FAB_ARG_3D(Bz),
     BL_FORT_FAB_ARG_3D(Bxout),
     BL_FORT_FAB_ARG_3D(Byout),
     BL_FORT_FAB_ARG_3D(Bzout),
     const BL_FORT_FAB_ARG_3D(src),
     const amrex::Real dx[], const amrex::Real* dt,
     D_DECL(BL_FORT_FAB_ARG_3D(xflux),
            BL_FORT_FAB_ARG_3D(yflux),
            BL_FORT_FAB_ARG_3D(zflux)),
     D_DECL(BL_FORT_FAB_ARG_3D(Ex),
            BL_FORT_FAB_ARG_3D(Ey),
            BL_FORT_FAB_ARG_3D(Ez)),
     const amrex::Real* cflLoc,
     const amrex::Real* e_added, const amrex::Real* ke_added,
     const int* print_fortran_warnings);

#ifdef __cplusplus
}
#endif

#endif

