#ifndef problem_initialize_mhd_data_H
#define problem_initialize_mhd_data_H


AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real A_z (int i, int j, int k,
          const GeometryData& geomdata)
{

    amrex::ignore_unused(k);

    // Compute A_z.  This lives on edges, e.g., {A_z}_{i-1/2,j-1/2,k}
    // so it is centered only in the z-direction.

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    // x_{i-1/2}
    Real x = problo[0] + dx[0] * (static_cast<Real>(i));

    // y_{j-1/2}
    Real y = problo[1] + dx[1] * (static_cast<Real>(j));

    Real r = std::sqrt(x * x + y * y);

    Real A;

    if (r <= 0.3_rt ) {
        A = problem::B_0 * (0.3_rt - r);
    } else {
        A = 0.0_rt;
    }

    return A;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_mhd_data (int i, int j, int k,
                                  Array4<Real> const& B, const int ibdir,
                                  const GeometryData& geomdata)
{

    const Real* dx = geomdata.CellSize();

    if (ibdir == 0) {

        Real A_z_jp1 = A_z(i, j+1, k, geomdata);
        Real A_z_0 = A_z(i, j, k, geomdata);

        B(i,j,k,0) = (A_z_jp1 - A_z_0) / dx[1];

    } else if (ibdir == 1) {

        Real A_z_ip1 = A_z(i+1, j, k, geomdata);
        Real A_z_0 = A_z(i, j, k, geomdata);

        B(i,j,k,0) = (A_z_0 - A_z_ip1) / dx[0];

    } else {

        B(i,j,k,0) = 0.0_rt;
    }
}
#endif
