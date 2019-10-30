#ifndef _DustCollapse_F_H_
#define _DustCollapse_F_H_
#include <AMReX_BLFort.H>


#ifdef __cplusplus
extern "C"
{
#endif

void fdustcollapse1d(const int* lo, const int* hi,
                     BL_FORT_FAB_ARG_3D(p),const int* nc_p,
                     const int nbins, amrex::Real* dens,
                     int* imask, const int mask_size, const int r1,
                     const int dens_comp, int* cnt);

void fdustcollapse2d(const int* lo, const int* hi,
                     BL_FORT_FAB_ARG_3D(p),const int* nc_p,
                     const int nbins, amrex::Real* dens, amrex::Real* volcount,
                     int* imask, const int mask_size, const int r1,
                     const amrex::Real* dx, const amrex::Real dx_fine,
                     const amrex::Real yctr, const int dens_comp);

void fdustcollapse3d(const int* lo, const int* hi,
                     BL_FORT_FAB_ARG_3D(p),const int* nc_p,
                     const int nbins, amrex::Real* dens, amrex::Real* ncount,
                     int* imask, const int mask_size, const int r1,
                     const amrex::Real* dx, const amrex::Real dx_fine,
                     const amrex::Real xctr,
                     const amrex::Real yctr, const amrex::Real zctr,
                     const int dens_comp);

#ifdef __cplusplus
}
#endif

#endif
