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
void problem_bc_fill(int i, int j, int k,
                      Array4<Real> const& state,
                      Real time, const Array1D<BCRec, 0, NUM_STATE-1>& bcs,
                      const GeometryData& geomdata) {

    amrex::ignore_unused(i);
    amrex::ignore_unused(j);
    amrex::ignore_unused(k);
    amrex::ignore_unused(state);
    amrex::ignore_unused(time);
    amrex::ignore_unused(bcs);
    amrex::ignore_unused(geomdata);
}

#endif
