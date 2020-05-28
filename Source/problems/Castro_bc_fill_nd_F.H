#ifndef CASTRO_BC_FILL_ND_F_H
#define CASTRO_BC_FILL_ND_F_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

void hypfill(const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state),
             const int* dlo, const int* dhi, const amrex::Real* dx,
             const amrex::Real* glo, const amrex::Real time, const int* bc);

void denfill(const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state),
             const int* dlo, const int* dhi, const amrex::Real* dx,
             const amrex::Real* glo, const amrex::Real time, const int* bc);

#ifdef MHD
void face_fillx(const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state),
                const int* dlo, const int* dhi, const amrex::Real* dx,
                const amrex::Real* glo, const amrex::Real time, const int* bc);

void face_filly(const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state),
                const int* dlo, const int* dhi, const amrex::Real* dx,
                const amrex::Real* glo, const amrex::Real time, const int* bc);

void face_fillz(const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state),
                const int* dlo, const int* dhi, const amrex::Real* dx,
                const amrex::Real* glo, const amrex::Real time, const int* bc);
#endif

inline void ca_nullfill(BL_FORT_FAB_ARG_3D(state), const int* dlo,
                        const int* dhi, const amrex::Real* dx,
                        const amrex::Real* glo, const amrex::Real* time,
                        const int* bc)
#ifdef __cplusplus
{
    amrex::Abort("Not supposed to call ca_nullfill");
}
#else
    ;
#endif

#ifdef __cplusplus
}
#endif

#endif
