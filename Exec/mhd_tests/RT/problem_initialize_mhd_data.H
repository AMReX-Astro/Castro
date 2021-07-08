#ifndef problem_initialize_mhd_data_H
#define problem_initialize_mhd_data_H

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_mhd_data (int i, int j, int k,
                                  Array4<Real> const& B, const int ibdir,
                                  const GeometryData& geomdata)
{


    if (ibdir == 0) {
        B(i,j,k) = problem::B_0;
    } else {
        B(i,j,k) = 0.0_rt;
    }

}
#endif
