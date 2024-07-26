#ifndef problem_initialize_mhd_data_H
#define problem_initialize_mhd_data_H

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_mhd_data (int i, int j, int k,
                                  Array4<Real> const& B, const int ibdir,
                                  const GeometryData& geomdata)
{
    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    // these are cell-center coordinates

    Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

    Real y = 0.0;
#if AMREX_SPACEDIM >= 2
    y = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
#endif

    if (ibdir == 0) {
        B(i,j,k) = -problem::B_0 * std::sin(2.0_rt * M_PI * y);
    } else if (ibdir == 1) {
        B(i,j,k) = problem::B_0 * std::sin(4.0_rt * M_PI * x);
    } else {
        B(i,j,k) = 0.0_rt;
    }
}
#endif

