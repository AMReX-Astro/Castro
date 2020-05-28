#ifndef CASTRO_DIFFUSION_UTIL_H
#define CASTRO_DIFFUSION_UTIL_H

void fill_temp_cond(const amrex::Box& bx,
                    amrex::Array4<amrex::Real const> const& U_arr,
                    amrex::Array4<amrex::Real> const& coeff_arr);

void fill_temp_diff_coeff(const amrex::Box& bx,
                          amrex::Array4<amrex::Real const> const& U_arr,
                          amrex::Array4<amrex::Real> const& coeff_arr);

#endif
