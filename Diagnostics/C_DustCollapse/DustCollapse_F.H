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
                     const int nbins,
                     amrex::Real* dens_bin, amrex::Real* vel_bin,
                     amrex::Real* pres_bin, amrex::Real* e_bin,
                     int* ncount,
                     int* imask, const int mask_size, const int r1,
                     const int dens_comp, const int xmom_comp,
                     const int ymom_comp,
                     const int pres_comp, const int rhoe_comp,
                     const amrex::Real dx_fine, const amrex::Real* dx,
                     const amrex::Real xctr, const amrex::Real yctr);

void fdustcollapse3d(const int* lo, const int* hi,
                     BL_FORT_FAB_ARG_3D(p),const int* nc_p,
                     const int nbins,
                     amrex::Real* dens_bin, amrex::Real* vel_bin,
                     amrex::Real* pres_bin,
                     int* ncount,
                     int* imask, const int mask_size, const int r1,
                     const int dens_comp, const int xmom_comp,
                     const int ymom_comp, const int zmom_comp,
                     const int pres_comp,
                     const amrex::Real dx_fine, const amrex::Real* dx,
                     const amrex::Real xctr, const amrex::Real yctr);

#ifdef __cplusplus
}
#endif

#endif
