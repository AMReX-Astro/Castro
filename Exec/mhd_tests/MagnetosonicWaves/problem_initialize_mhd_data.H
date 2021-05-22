#ifndef problem_initialize_mhd_data_H
#define problem_initialize_mhd_data_H

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_mhd_data (int i, int j, int k,
                                  Array4<Real> const& B, const int ibdir,
                                  const GeometryData& geomdata)
{
    int coord_type = geomdata.Coord();

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

    Real y = 0.0;
#if AMREX_SPACEDIM >= 2
    y = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
#endif

    Real z = 0.0;
#if AMREX_SPACEDIM == 3
    z = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
#endif

    Real pert = 1.0e-5_rt * std::sin(2.0_rt * M_PI * (problem::k_x * x + problem::k_y * y));

    if (ibdir == 0) {
        B(i,j,k) = (problem::B_0 * problem::u_bx) -
            std::sqrt(2.0_rt) * problem::B_0 * (problem::cw * problem::cw -
                                                problem::c_s) * problem::u_ky * pert;

    } else if (ibdir == 1) {
        B(i,j,k) = (problem::B_0 * problem::u_by) +
            std::sqrt(2.0_rt) * problem::B_0 * (problem::cw * problem::cw -
                                                problem::c_s) * problem::u_kx * pert;

    } else {
        B(i,j,k) = 0.0_rt;
    }

}
#endif
