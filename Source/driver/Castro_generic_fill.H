#ifndef CASTRO_GENERIC_FILL_H
#define CASTRO_GENERIC_FILL_H

#include <AMReX_BLFort.H>
#include <Castro.H>

#ifdef AMREX_USE_CUDA
#include <cuda_runtime_api.h>
#include <AMReX_Arena.H>
#endif

struct CastroGenericFill
{
    AMREX_GPU_DEVICE
    void operator() (const IntVect& /*iv*/, Array4<Real> const& /*data*/,
                     const int /*dcomp*/, const int /*numcomp*/,
                     GeometryData const& /*geom*/, const Real /*time*/,
                     const BCRec* /*bcr*/, const int /*bcomp*/,
                     const int /*orig_comp*/) const
        {
            // Dummy routine that does nothing for inflow boundaries.
            // We assume that there are no inflow boundaries for a
            // generic fill, since we overwrote them with first-order
            // extrapolation in variableSetUp().
        }
};

void ca_generic_fill(amrex::Box const& bx, amrex::FArrayBox& data,
                     const int dcomp, const int numcomp,
                     amrex::Geometry const& geom, const amrex::Real time,
                     const amrex::Vector<amrex::BCRec>& bcr, const int bcomp,
                     const int scomp);

#endif
