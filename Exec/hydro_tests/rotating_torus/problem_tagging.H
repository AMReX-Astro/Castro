#ifndef problem_tagging_H
#define problem_tagging_H

#include <Castro_util.H>
#include <prob_parameters.H>
#include <AMReX_TagBox.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_tagging(int i, int j, int k,
                     Array4<char> const& tag,
                     Array4<Real> const& state,
                     const GeometryData& geomdata)
{
    GpuArray<Real, 3> loc;
    position(i, j, k, geomdata, loc);
    for (int i = 0; i < 3; ++i) {
        loc[i] -= problem::center[i];
    }

    Real R = std::sqrt(loc[0] * loc[0] + loc[1] * loc[1]);

    Real Z = loc[2];

    Real dist2 = (problem::torus_center - R) * (problem::torus_center - R) + Z * Z;

    if (dist2 < problem::torus_width * problem::torus_width) {
        tag(i,j,k) = TagBox::SET;
    } else {
        tag(i,j,k) = TagBox::CLEAR;
    }
}

#endif
