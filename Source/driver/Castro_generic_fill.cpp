
#include <AMReX_BLFort.H>
#include <Castro.H>
#include <Castro_generic_fill.H>
#include <Castro_generic_fill_F.H>

using namespace amrex;

void ca_generic_fill(Box const& bx, FArrayBox& data,
                     const int dcomp, const int numcomp,
                     Geometry const& geom, const Real time,
                     const Vector<BCRec>& bcr, const int bcomp,
                     const int scomp)
{
    // Make a copy of the raw BCRec data in the format
    // our BC routines can handle (a contiguous array
    // of integers).

    Vector<int> bcrs(2 * AMREX_SPACEDIM * numcomp);

    for (int n = 0; n < numcomp; ++n)
        for (int k = 0; k < 2 * AMREX_SPACEDIM; ++k)
            bcrs[2 * AMREX_SPACEDIM * n + k] = bcr[n].vect()[k];

#ifdef AMREX_USE_CUDA
    int* bc_f = prepare_bc(bcrs.data(), numcomp);
    set_bc_launch_config();
#else
    const int* bc_f = bcrs.data();
#endif

#pragma gpu box(bx)
    generic_fill(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                 BL_TO_FORTRAN_N_ANYD(data, dcomp), numcomp,
                 AMREX_INT_ANYD(geom.Domain().loVect()), AMREX_INT_ANYD(geom.Domain().hiVect()),
                 AMREX_REAL_ANYD(geom.CellSize()), AMREX_REAL_ANYD(geom.ProbLo()), bc_f);

#ifdef AMREX_USE_CUDA
    clean_bc_launch_config();
    clean_bc(bc_f);
#endif
}
