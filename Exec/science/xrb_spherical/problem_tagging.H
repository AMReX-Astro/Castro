#ifndef problem_tagging_H
#define problem_tagging_H
#include <AMReX_TagBox.H>
#include <prob_parameters.H>
#include <Castro_util.H>
#include <network.H>

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

    if (problem::tag_by_thermodynamics) {

        // Consider refining when density and temp are above some threshold
        if (state(i,j,k,URHO)  > problem::cutoff_density &&
            state(i,j,k,UTEMP) > problem::cutoff_temp) {

            if (level < problem::max_base_tagging_level) {
                tag(i,j,k) = TagBox::SET;
            }

            // find the distance in cylindrical radial direction
            Real dist = std::abs(loc[0] * std::sin(loc[1]));

            // Allow higher level refinement to better resolve
            // the initial burning phase which is at the pole of the NS.
            if (dist < problem::refine_distance &&
                state(i,j,k,UFS+1) >  problem::X_min * state(i,j,k,URHO) &&
                level < problem::max_hse_tagging_level) {
                tag(i,j,k) = TagBox::SET;
            }
        }

    } else {

        // tag everything below a certain height
        if (loc[0] < geomdata.ProbLo(0) + problem::refine_height) {
            if (level < problem::max_base_tagging_level) {
                tag(i,j,k) = TagBox::SET;
            }
        }
    }

}
#endif
