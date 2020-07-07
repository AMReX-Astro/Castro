#ifndef CASTRO_HYDRO_F_H
#define CASTRO_HYDRO_F_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

#ifdef RADIATION
  void ctu_rad_consup
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(hydro_source),
     const BL_FORT_FAB_ARG_3D(Erborder),
     BL_FORT_FAB_ARG_3D(S_new),
     BL_FORT_FAB_ARG_3D(Er_out),
     const BL_FORT_FAB_ARG_3D(rad_flux0),
     const BL_FORT_FAB_ARG_3D(qe0),
     const BL_FORT_FAB_ARG_3D(area0),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(rad_flux1),
     const BL_FORT_FAB_ARG_3D(qe1),
     const BL_FORT_FAB_ARG_3D(area1),
#endif
#if AMREX_SPACEDIM == 3
     const BL_FORT_FAB_ARG_3D(rad_flux2),
     const BL_FORT_FAB_ARG_3D(qe2),
     const BL_FORT_FAB_ARG_3D(area2),
#endif
     int* nstep_fsp,
     const BL_FORT_FAB_ARG_3D(volume),
     const amrex::Real* dx, const amrex::Real dt);
#endif


#ifdef __cplusplus
}
#endif

#endif
