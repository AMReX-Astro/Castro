#ifndef _Radiation_F_H_
#define _Radiation_F_H_
#include <AMReX_BLFort.H>


#ifdef __cplusplus
extern "C"
{
#endif

void fgaussian_pulse(const int* lo, const int* hi,
                BL_FORT_FAB_ARG_3D(p),const int* nc_p,
                const int nbins,
                amrex::Real* rad_bin,
                int* ncount,
                int* imask, const int mask_size, const int r1,
                const int rad_comp, const amrex::Real* dx,
                const amrex::Real dx_fine,
                const amrex::Real xctr, const amrex::Real yctr);


#ifdef __cplusplus
}
#endif

#endif
