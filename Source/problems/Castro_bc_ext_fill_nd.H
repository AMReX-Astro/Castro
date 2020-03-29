#ifndef _Castro_bc_ext_fill_nd_H_
#define _Castro_bc_ext_fill_nd_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

void ca_ext_fill(amrex::Box const& bx, amrex::FArrayBox& data,
                 const int dcomp, const int numcomp,
                 amrex::Geometry const& geom, const amrex::Real time,
                 const int* bc_f);

void ca_ext_denfill(amrex::Box const& bx, amrex::FArrayBox& data,
                    const int dcomp, const int numcomp,
                    amrex::Geometry const& geom, const amrex::Real time,
                    const int* bc_f);

#endif
