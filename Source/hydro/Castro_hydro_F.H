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

  void ca_track_grid_losses
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(flux0),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(flux1),
#endif
#if AMREX_SPACEDIM == 3
     const BL_FORT_FAB_ARG_3D(flux2),
#endif
     amrex::Real* mass_lost, amrex::Real* xmom_lost, amrex::Real* ymom_lost, amrex::Real* zmom_lost,
     amrex::Real* eden_lost, amrex::Real* xang_lost, amrex::Real* yang_lost, amrex::Real* zang_lost);

  void edge_state_temp_to_pres
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(qm),
     BL_FORT_FAB_ARG_3D(qp));

  void ca_mol_consup
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(shk),
     const BL_FORT_FAB_ARG_3D(statein),
     const BL_FORT_FAB_ARG_3D(srcU),
     BL_FORT_FAB_ARG_3D(update),
     const amrex::Real* dx, const amrex::Real dt,
     BL_FORT_FAB_ARG_3D(flux1),
#if AMREX_SPACEDIM >= 2
     BL_FORT_FAB_ARG_3D(flux2),
#endif
#if AMREX_SPACEDIM == 3
     BL_FORT_FAB_ARG_3D(flux3),
#endif
     const BL_FORT_FAB_ARG_3D(area1),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(area2),
#endif
#if AMREX_SPACEDIM == 3
     const BL_FORT_FAB_ARG_3D(area3),
#endif
     const BL_FORT_FAB_ARG_3D(q1),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(q2),
#endif
#if AMREX_SPACEDIM == 3
     const BL_FORT_FAB_ARG_3D(q3),
#endif
     const BL_FORT_FAB_ARG_3D(volume));

  void ca_mol_diffusive_flux
    (const int* lo, const int* hi,
     const int idir,
     const BL_FORT_FAB_ARG_3D(U),
     const BL_FORT_FAB_ARG_3D(cond),
     BL_FORT_FAB_ARG_3D(flux),
     const amrex::Real* dx);


#ifdef __cplusplus
}
#endif

#endif
