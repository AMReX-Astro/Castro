#ifndef CASTRO_SOURCES_F_H
#define CASTRO_SOURCES_F_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

#ifdef SPONGE
  void update_sponge_params(const amrex::Real* time);
#endif

  void ca_ext_src
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(old_state),
     const BL_FORT_FAB_ARG_3D(new_state),
     BL_FORT_FAB_ARG_3D(ext_src),
     const amrex::Real* prob_lo, const amrex::Real* dx,
     const amrex::Real time, const amrex::Real dt);

#ifdef __cplusplus
}
#endif

#endif
