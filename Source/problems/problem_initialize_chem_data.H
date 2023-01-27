#ifndef problem_initialize_chem_data_H
#define problem_initialize_chem_data_H

///
/// Initialize the chemical potential data for NSE
///
/// @param i         x-index
/// @param j         y-index
/// @param k         z-index
/// @param state     chemical potential data (Fab)
/// @param geomdata  geometry data
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_chem_data (int i, int j, int k, Array4<Real> const& state, const GeometryData& geomdata) {}

#endif
