#ifndef _Limiter_F_H_
#define _Limiter_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

void microphysics_initialize(const int* name, const int* namlen);
void microphysics_finalize();

void find_timestep_limiter(const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(p),
                           const int* nc_p, const int dens_comp,
                           const int xmom_comp, const int ymom_comp,
                           const int zmom_comp, const int pres_comp,
                           const int rhoe_comp, const int spec_comp,
                           const int time_integration_method,
                           const amrex::Real* dx, amrex::Real* dt,
                           amrex::Real* dt_loc);

void find_timestep_limiter_burning(const int* lo, const int* hi,
                                   BL_FORT_FAB_ARG_3D(p), const int* nc_p,
                                   const int dens_comp, const int temp_comp,
                                   const int rhoe_comp, const int spec_comp,
                                   const amrex::Real* dx, amrex::Real* dt,
                                   amrex::Real* dt_loc);

void find_timestep_limiter_diffusion(const int* lo, const int* hi,
                                     BL_FORT_FAB_ARG_3D(p), const int* nc_p,
                                     const int dens_comp, const int temp_comp,
                                     const int rhoe_comp, const int spec_comp,
                                     const amrex::Real* dx, amrex::Real* dt,
                                     amrex::Real* dt_loc);

#ifdef __cplusplus
}
#endif

#endif
