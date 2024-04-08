#ifndef problem_initialize_mhd_data_H
#define problem_initialize_mhd_data_H

///
/// x component of dipole vector potential
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real compute_A_x(const Real x, const Real y, const Real z)
{

    // this should be called with a node-centered y and z and a cell-centered x

    Real theta = std::atan2(std::sqrt(x * x + y * y), z);
    Real phi = std::atan2(y, x);

    Real dist = std::sqrt(x * x + y * y + z * z);

    Real A_x;

    if (dist <= 1.0e8 ) {
        A_x = -problem::m_0 * dist * std::sin(theta) * std::sin(phi) / (1.0e24);

    } else {
        // outside the star radius threshold
        A_x = -problem::m_0 * std::sin(theta) * std::sin(phi) / (dist*dist);
    }

    return A_x;
}

///
/// y component of dipole vector potential
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real compute_A_y(const Real x, const Real y, const Real z)
{

    // this should be called with a node-centered x and z and a cell-centered y

    Real theta = std::atan2(std::sqrt(x * x + y * y), z);
    Real phi = std::atan2(y, x);

    Real dist = std::sqrt(x * x + y * y + z * z);

    Real A_y;

    if (dist <= 1.0e8 ) {
        A_y = problem::m_0 * dist * std::sin(theta) * std::cos(phi) / (1.0e24);

    } else {
        // outside the star radius threshold
        A_y = problem::m_0 * std::sin(theta) * std::cos(phi) / (dist*dist);
    }

    return A_y;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_mhd_data (int i, int j, int k,
                                  Array4<Real> const& B, const int ibdir,
                                  const GeometryData& geomdata)
{
    int coord_type = geomdata.Coord();

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real xcen = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt) - problem::center[0];
    Real xl = problo[0] + dx[0] * (static_cast<Real>(i)) - problem::center[0];
    Real xr = problo[0] + dx[0] * (static_cast<Real>(i) + 1.0_rt) - problem::center[0];

    Real ycen = 0.0;
    Real yl = 0.0;
    Real yr = 0.0;
#if AMREX_SPACEDIM >= 2
    ycen = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt) - problem::center[1];
    yl = problo[1] + dx[1] * (static_cast<Real>(j)) - problem::center[1];
    yr = problo[1] + dx[1] * (static_cast<Real>(j) + 1.0_rt) - problem::center[1];
#endif

    Real zcen = 0.0;
    Real zl = 0.0;
    Real zr = 0.0;
#if AMREX_SPACEDIM == 3
    zcen = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt) - problem::center[2];
    zl = problo[2] + dx[2] * (static_cast<Real>(k)) - problem::center[2];
    zr = problo[2] + dx[2] * (static_cast<Real>(k) + 1.0_rt) - problem::center[2];
#endif

    // We'll compute B from the vector potential A
    // A is of the form A = (A_x, A_y, 0)

    // inside the star radius threshold (R) the vector potential
    // will be of the form A_x = constant * x
    //                         = constant * rsin(theta) * cos(phi)
    // the constant is then obtained by matching B.C with
    // the field outside at r = R
    //
    // with an analogous form for A_y

    // centerings:
    //   A_x, i, j-1/2, k-1/2
    //   A_y, i-1/2, j, k-1/2

    if (ibdir == 0) {

        // initialize Bx on x-faces
        // B_x = -dA_y/dz

        // we need A_y_{i-1/2,j,k+1/2} and A_y_{i-1/2, j, k-1/2}
        // differencing across the x-face (i-1/2) in the z-direction
        // then gives us B_x_{i-1/2,j,k}

        Real A_y_zp = compute_A_y(xl, ycen, zr);
        Real A_y_zm = compute_A_y(xl, ycen, zl);

        B(i,j,k) = -(A_y_zp - A_y_zm) / dx[2];

    } else if (ibdir == 1) {

        // initialize By on y-faces
        // B_y = dA_x/dz

        Real A_x_zp = compute_A_x(xcen, yl, zr);
        Real A_x_zm = compute_A_x(xcen, yl, zl);

        B(i,j,k) = (A_x_zp - A_x_zm) / dx[2];

    } else {

        // initialize Bz on z-faces
        // B_z = dA_y/dx - dA_x/dy

        Real A_y_xp = compute_A_y(xr, ycen, zl);
        Real A_y_xm = compute_A_y(xl, ycen, zl);

        Real A_x_yp = compute_A_x(xcen, yr, zl);
        Real A_x_ym = compute_A_x(xcen, yl, zl);

        B(i,j,k) = (A_y_xp - A_y_xm) / dx[0] - (A_x_yp - A_x_ym) / dx[1];

    }

}
#endif
