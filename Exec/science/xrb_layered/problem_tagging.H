#ifndef problem_tagging_H
#define problem_tagging_H

#include <Castro_util.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_tagging(int i, int j, int k,
                     Array4<char> const& tag,
                     Array4<Real const> const& state,
                     int level, const GeometryData& geomdata)
{
    // Tag on regions of with H > H_min and rho < cutoff_density.
    // Also allow us to not refine the very base.
    // Note that H is the first species variable and so is in index
    // UFS of the state array.

    GpuArray<Real, 3> loc;
    position(i, j, k, geomdata, loc);


    if (state(i,j,k,URHO) > problem::cutoff_density &&
        state(i,j,k,UFS) > problem::H_min &&
        loc[AMREX_SPACEDIM-1] > problem::refine_height) {
        tag(i,j,k) = TagBox::SET;
    }
}

#endif
