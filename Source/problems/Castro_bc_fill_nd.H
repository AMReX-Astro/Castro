#ifndef _Castro_bc_fill_nd_H_
#define _Castro_bc_fill_nd_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

void ca_statefill(amrex::Box const& bx, amrex::FArrayBox& data,
                  const int dcomp, const int numcomp,
                  amrex::Geometry const& geom, const amrex::Real time,
                  const amrex::Vector<amrex::BCRec>& bcr, const int bcomp,
                  const int scomp);

#endif
