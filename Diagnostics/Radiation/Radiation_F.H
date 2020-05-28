#ifndef _Radiation_F_H_
#define _Radiation_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

void fgaussian_pulse(const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(p),
                     const int* nc_p, const int nbins, amrex::Real* rad_bin,
                     int* ncount, int* imask, const int mask_size, const int r1,
                     const int rad_comp, const amrex::Real* dx,
                     const amrex::Real dx_fine, const amrex::Real xctr,
                     const amrex::Real yctr);

void flgt_frnt1d(const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(p),
                 const int* nc_p, const int nbins, amrex::Real* dens_bin,
                 amrex::Real* vel_bin, amrex::Real* pres_bin,
                 amrex::Real* rad_bin, int* imask, const int mask_size,
                 const int r1, const int dens_comp, const int xmom_comp,
                 const int pres_comp, const int rad_comp, const amrex::Real* dx,
                 const amrex::Real dx_fine);

void fradshock(const int* lo, const int* hi, const amrex::Real* problo,
               const amrex::Real* probhi, BL_FORT_FAB_ARG_3D(p),
               const int* nc_p, const int nbins, amrex::Real* vars_bin,
               int* imask, const int mask_size, const int r1, const int rr,
               const amrex::Real* dx, const int idir, int* cnt);

void fradsource(const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(p),
                const int* nc_p, amrex::Real* rhoe, amrex::Real* rad,
                const int rhoe_comp, const int rad_comp);

void fradsphere(const int* lo, const int* hi, const amrex::Real* problo,
                const amrex::Real* probhi, BL_FORT_FAB_ARG_3D(p),
                const int* nc_p, const int nbins, amrex::Real* vars_bin,
                int* imask, const int mask_size, const int r1,
                const amrex::Real* dx, int* cnt);

void frhdshocktube(const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(p),
                   const int* nc_p, const int nbins, amrex::Real* dens_bin,
                   amrex::Real* vel_bin, amrex::Real* pres_bin,
                   amrex::Real* rad_bin, const int dens_comp,
                   const int xmom_comp, const int pres_comp, const int rad_comp,
                   const int ngroups);

#ifdef __cplusplus
}
#endif

#endif
