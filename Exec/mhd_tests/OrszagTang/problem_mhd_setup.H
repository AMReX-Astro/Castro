#ifndef problem_mhd_setup_H
#define problem_mhd_setup_H

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_mhd_data (int i, int j, int k,
                                  Array4<Real> const& B, const int ibdir,
                                  const GeometryData& geomdata)
{
    int coord_type = geomdata.Coord();

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    // these are cell-center coordinates

    Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

    Real y = 0.0;
#if AMREX_SPACEDIM >= 2
    y = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
#endif

    Real z = 0.0;
#if AMREX_SPACEDIM == 3
    z = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
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

