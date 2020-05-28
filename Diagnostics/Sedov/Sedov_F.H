#ifndef _Sedov_F_H_
#define _Sedov_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

void fextract1d(const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(p),
                const int* nc_p, const int nbins, amrex::Real* dens_bin,
                amrex::Real* vel_bin, amrex::Real* pres_bin, amrex::Real* e_bin,
                int* imask, const int mask_size, const int r1,
                const int dens_comp, const int xmom_comp, const int pres_comp,
                const int rhoe_comp);

void fextract2d_cyl(const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(p),
                    const int* nc_p, const int nbins, amrex::Real* dens_bin,
                    amrex::Real* vel_bin, amrex::Real* pres_bin,
                    amrex::Real* e_bin, int* ncount, int* imask,
                    const int mask_size, const int r1, const int dens_comp,
                    const int xmom_comp, const int ymom_comp,
                    const int pres_comp, const int rhoe_comp,
                    const amrex::Real dx_fine, const amrex::Real* dx,
                    const amrex::Real xctr, const amrex::Real yctr);

void fextract2d_sph(const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(p),
                    const int* nc_p, const int nbins, amrex::Real* dens_bin,
                    amrex::Real* vel_bin, amrex::Real* pres_bin,
                    amrex::Real* e_bin, amrex::Real* volcount, int* imask,
                    const int mask_size, const int r1, const int dens_comp,
                    const int xmom_comp, const int ymom_comp,
                    const int pres_comp, const int rhoe_comp,
                    const amrex::Real dx_fine, const amrex::Real* dx,
                    const amrex::Real xctr, const amrex::Real yctr);

void fextract3d_cyl(const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(p),
                    const int* nc_p, const int nbins, amrex::Real* dens_bin,
                    amrex::Real* vel_bin, amrex::Real* pres_bin, int* ncount,
                    int* imask, const int mask_size, const int r1,
                    const int dens_comp, const int xmom_comp,
                    const int ymom_comp, const int zmom_comp,
                    const int pres_comp, const amrex::Real dx_fine,
                    const amrex::Real* dx, const amrex::Real xctr,
                    const amrex::Real yctr);

void fextract3d_sph(const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(p),
                    const int* nc_p, const int nbins, amrex::Real* dens_bin,
                    amrex::Real* vel_bin, amrex::Real* pres_bin,
                    amrex::Real* e_bin, int* ncount,

                    int* imask, const int mask_size, const int r1,
                    const int dens_comp, const int xmom_comp,
                    const int ymom_comp, const int zmom_comp,
                    const int pres_comp, const int rhoe_comp,
                    const amrex::Real dx_fine, const amrex::Real* dx,
                    const amrex::Real xctr, const amrex::Real yctr,
                    const amrex::Real zctr);

#ifdef __cplusplus
}
#endif

#endif
