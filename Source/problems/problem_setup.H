#ifndef problem_setup_H
#define problem_setup_H

///
/// Initialize problem-specific variables
///
AMREX_INLINE
void problem_initialize () {}

///
/// Initialize the simulation state data
///
/// @param i         x-index
/// @param j         y-index
/// @param k         z-index
/// @param state     simulation state (Fab)
/// @param geomdata  geometry data
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k, Array4<Real> const& state, const GeometryData& geomdata) {}

#endif
