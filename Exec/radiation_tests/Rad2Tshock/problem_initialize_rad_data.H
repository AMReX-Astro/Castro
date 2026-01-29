#ifndef problem_initialize_rad_data_H
#define problem_initialize_rad_data_H

#include <blackbody.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_rad_data (int i, int j, int k,
                                  Array4<Real> const& rad,
                                  GpuArray<Real, NGROUPS+1> const& xnu,
                                  GpuArray<Real, NGROUPS> const& nugroup,
                                  GpuArray<Real, NGROUPS> const& dnugroup,
                                  const GeometryData& geomdata)
{

    amrex::ignore_unused(nugroup);
    amrex::ignore_unused(dnugroup);

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    // Provides the simulation to be run in the x,y,or z direction
    // where length direction is the length side in a square prism

    Real length_cell;
    if (problem::idir == 1) {
        length_cell = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);
    } else if (problem::idir == 2) {
        length_cell = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
    } else if (problem::idir == 3) {
        length_cell = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
    } else {
#ifndef AMREX_USE_GPU
        amrex::Error("Invalid direction please input idir = [1,3]");
#endif
    }

    Real T;
    if (length_cell < 0.0_rt) {
        T = problem::T0;
    } else {
        T = problem::T1;
    }

    if (NGROUPS == 1) {
        rad(i,j,k,0) = C::a_rad * std::pow(T, 4);

    } else {
        for (int g = 0; g < NGROUPS; g++) {
            rad(i,j,k,g) = BGroup(T, xnu[g], xnu[g+1]);
        }
    }
}
#endif
