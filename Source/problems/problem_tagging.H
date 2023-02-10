#ifndef problem_tagging_H
#define problem_tagging_H

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
void problem_tagging([[maybe_unused]] int i,
                     [[maybe_unused]] int j,
                     [[maybe_unused]] int k,
                     [[maybe_unused]] Array4<char> const& tag,
                     [[maybe_unused]] Array4<Real const> const& state,
                     [[maybe_unused]] int level,
                     [[maybe_unused]] const GeometryData& geomdata) {}

#endif
