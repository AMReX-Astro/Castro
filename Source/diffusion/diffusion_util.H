#ifndef _DIFFUSION_UTIL_H_
#define _DIFFUSION_UTIL_H_

void
fill_temp_cond(const amrex::Box& bx,
               amrex::Array4<amrex::Real const> const U_arr,
               amrex::Array4<amrex::Real> const coeff_arr);

void
fill_temp_diff_coeff(const amrex::Box& bx,
                     amrex::Array4<amrex::Real const> const U_arr,
                     amrex::Array4<amrex::Real> const coeff_arr);

#endif
