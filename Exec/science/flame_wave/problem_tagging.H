#ifndef problem_tagging_H
#define problem_tagging_H
#include <AMReX_TagBox.H>
#include <prob_parameters.H>
#include <Castro_util.H>

///
/// Define problem-specific tagging criteria
///
/// @param i         x-index
/// @param j         y-index
/// @param k         z-index
/// @param tag       tag array (TagBox)
/// @param state     simulation state (Fab)
/// @param level     AMR level
/// @param geomdata  geometry data
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_tagging(int i, int j, int k,
                     Array4<char> const& tag,
                     Array4<Real const> const& state,
                     int level, const GeometryData& geomdata)
{

    GpuArray<Real, 3> loc;
    position(i, j, k, geomdata, loc);

    if (state(i,j,k,URHO) > problem::cutoff_density &&
        state(i,j,k,UFS) / state(i,j,k,URHO) > problem::X_min) {

#if AMREX_SPACEDIM == 2
        Real dist = std::abs(loc[0] - problem::center[0]);
#else
        Real dist = std::sqrt((loc[0] - problem::center[0]) * (loc[0] - problem::center[0]) +
                              (loc[1] - problem::center[1]) * (loc[1] - problem::center[1]));
#endif
        if (level < problem::max_hse_tagging_level && dist < problem::x_refine_distance) {
            tag(i,j,k) = TagBox::SET;
        }
    }

    if (state(i,j,k,URHO) > problem::cutoff_density) {
        if (level < problem::max_base_tagging_level) {
            tag(i,j,k) = TagBox::SET;
        }
    }
}
#endif
