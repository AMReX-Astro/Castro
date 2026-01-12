#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>

#ifndef CASTRO_MGUTILS_H
#define CASTRO_MGUTILS_H

using namespace amrex;

void
apply_metric(const Box& bx,
             Array4<Real> const& rhs, const Box& rbx,
             Array4<Real> const& ecx, const Box& xbx,

#if AMREX_SPACEDIM >= 2
             Array4<Real> const& ecy, const Box& ybx,
#endif
             GpuArray<Real, AMREX_SPACEDIM> dx,
             const int coord_type);

#endif
