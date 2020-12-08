#ifndef problem_initialize_rad_data_H
#define problem_initialize_rad_data_H

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_rad_data (int i, int j, int k,
                                  Array4<Real> const& rad,
                                  GpuArray<Real, NGROUPS+1> const& xnu,
                                  const GeometryData& geomdata)
{

    for (int g = 0; g < NGROUPS; g++) {
        rad(i,j,k,g) = 0.0_rt;
    }
}
#endif
