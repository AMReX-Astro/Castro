#include <AMReX_BLFort.H>
#include <Castro.H>
#include <Castro_bc_fill_nd.H>
#include <Castro_bc_ext_fill_nd_F.H>

using namespace amrex;

void ca_ext_fill(Box const& bx, FArrayBox& data,
                 const int /*dcomp*/, const int /*numcomp*/,
                 Geometry const& geom, const Real time,
                 const int* bc_f)
{
#pragma gpu box(bx)
    ext_fill(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
             BL_TO_FORTRAN_ANYD(data),
             AMREX_INT_ANYD(geom.Domain().loVect()), AMREX_INT_ANYD(geom.Domain().hiVect()),
             AMREX_REAL_ANYD(geom.CellSize()), AMREX_REAL_ANYD(geom.ProbLo()), time, bc_f);
}

void ca_ext_denfill(Box const& bx, FArrayBox& data,
                    const int dcomp, const int /*numcomp*/,
                    Geometry const& geom, const Real time,
                    const int* bc_f)
{
#pragma gpu box(bx)
    ext_denfill(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                BL_TO_FORTRAN_N_ANYD(data, dcomp),
                AMREX_INT_ANYD(geom.Domain().loVect()), AMREX_INT_ANYD(geom.Domain().hiVect()),
                AMREX_REAL_ANYD(geom.CellSize()), AMREX_REAL_ANYD(geom.ProbLo()), time, bc_f);
}
