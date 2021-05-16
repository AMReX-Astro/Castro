#ifndef filt_prim_H
#define filt_prim_H

#include <AMReX.H>
#include <filter.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void filt_prim (int i, int j, int k,
                Array4<Real> const& Stmp,
                Array4<Real> const& Snew,
                int filt_T, int S,
                const GeometryData& geomdata,
                Real time)
{
    // this is a stub -- a problem can override this in its own directory
    // to implement filtering
}

#endif
