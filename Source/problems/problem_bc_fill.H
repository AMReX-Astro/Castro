#ifndef problem_bc_fill_H
#define problem_bc_fill_H

///
/// Problem-specific boundary condition fill
///
/// @param i         x-index
/// @param j         y-index
/// @param k         z-index
/// @param state     simulation state (Fab)
/// @param time      simulation time
/// @param bcs       boundary conditions
/// @param geomdata  geometry data
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_bc_fill([[maybe_unused]] int i,
                     [[maybe_unused]] int j,
                     [[maybe_unused]] int k,
                     [[maybe_unused]] Array4<Real> const& state,
                     [[maybe_unused]] Real time,
                     [[maybe_unused]] Array1D<BCRec, 0, NUM_STATE-1> bcs,
                     [[maybe_unused]] GeometryData geomdata) {}

#endif
