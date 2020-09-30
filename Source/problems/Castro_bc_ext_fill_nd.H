#ifndef CASTRO_BC_EXT_FILL_ND_H
#define CASTRO_BC_EXT_FILL_ND_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

void ca_ext_fill(amrex::Box const& bx, amrex::FArrayBox& data,
                 const int dcomp, const int numcomp,
                 amrex::Geometry const& geom, const amrex::Real time,
                 const int* bc_f);

#ifndef CXX_MODEL_PARSER
void ca_ext_denfill(amrex::Box const& bx, amrex::FArrayBox& data,
                    const int dcomp, const int numcomp,
                    amrex::Geometry const& geom, const amrex::Real time,
                    const int* bc_f);
#endif

#endif
