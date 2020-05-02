#ifndef _Castro_bc_fill_nd_H_
#define _Castro_bc_fill_nd_H_


#include <AMReX_BLFort.H>


#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

void ca_statefill(amrex::Box const& bx, amrex::FArrayBox& data,
                  const int dcomp, const int numcomp,
                  amrex::Geometry const& geom, const amrex::Real time,
                  const amrex::Vector<amrex::BCRec>& bcr, const int bcomp,
                  const int scomp);

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef MHD
  void ca_face_fillx
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

  void ca_face_filly
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

  void ca_face_fillz
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);
#endif
#ifdef __cplusplus
}
#endif

#endif
