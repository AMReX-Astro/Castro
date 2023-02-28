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

    if (problem::idir == 1) {
        // initialize along the x-direction

        if (ibdir == 0) {
            if (x <= problem::split[0]+1) {
                B(i,j,k) = problem::B_x_l;
            } else {
                B(i,j,k) = problem::B_x_r;
            }
        } else if (ibdir == 1) {
            if (x <= problem::split[0]) {
                B(i,j,k) = problem::B_y_l;
            } else {
                B(i,j,k) = problem::B_y_r;
            }
        } else if (ibdir == 2) {
            if (x <= problem::split[0]) {
                B(i,j,k) = problem::B_z_l;
            } else {
                B(i,j,k) = problem::B_z_r;
            }
        }

    } else if (problem::idir == 2) {
        // initialize along the y-direction

        if (ibdir == 0) {
            if (y <= problem::split[1]) {
                B(i,j,k) = problem::B_x_l;
            } else {
                B(i,j,k) = problem::B_x_r;
            }
        } else if (ibdir == 1) {
            if (y <= problem::split[1]+1) {
                B(i,j,k) = problem::B_y_l;
            } else {
                B(i,j,k) = problem::B_y_r;
            }
        } else if (ibdir == 2) {
            if (y <= problem::split[1]) {
                B(i,j,k) = problem::B_z_l;
            } else {
                B(i,j,k) = problem::B_z_r;
            }
        }

    } else {
        // initialize along the z-direction

        if (ibdir == 0) {
            if (z <= problem::split[2]) {
                B(i,j,k) = problem::B_x_l;
            } else {
                B(i,j,k) = problem::B_x_r;
            }
        } else if (ibdir == 1) {
            if (z <= problem::split[2]) {
                B(i,j,k) = problem::B_y_l;
            } else {
                B(i,j,k) = problem::B_y_r;
            }
        } else if (ibdir == 2) {
            if (z <= problem::split[2]+1) {
                B(i,j,k) = problem::B_z_l;
            } else {
                B(i,j,k) = problem::B_z_r;
            }
        }

    }
}
#endif
