#ifndef problem_tagging_H
#define problem_tagging_H

#include <state_indices.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_tagging(int i, int j, int k,
                     Array4<char> const& tag,
                     Array4<Real const> const& state,
                     int level, const GeometryData& geomdata)
{
    // Tag on regions of with H > H_min and rho < cutoff_density.
    // Note that H is the first species variable and so is in index UFS of the state array.

    if (state(i,j,k,URHO) > problem::cutoff_density &&
        (state(i,j,k,UFS) > state(i,j,k,URHO) * problem::H_min ||
         state(i,j,k,UFS+1) > state(i,j,k,URHO) * problem::He_min)) {
        tag(i,j,k) = TagBox::SET;
    }
}

#endif
