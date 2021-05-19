#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

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

    state(i,j,k,URHO) = problem::rho;
    state(i,j,k,UMX) = problem::rho * problem::u_x;
    state(i,j,k,UMY) = problem::rho * problem::u_y;
    state(i,j,k,UMZ) = problem::rho * problem::u_z;
    state(i,j,k,UEDEN) = problem::rhoe +
        0.5_rt * problem::rho * (problem::u_x * problem::u_x +
                                 problem::u_y * problem::u_y +
                                 problem::u_z * problem::u_z);
    state(i,j,k,UEINT) = problem::rhoe;
    state(i,j,k,UTEMP) = problem::T;

    Real r2 = std::pow(x - 0.5_rt, 2);
#if AMREX_SPACEDIM >= 2
    r2 += std::pow(y - 0.5_rt, 2);
#endif
#if AMREX_SPACEDIM == 3
    r2 += std::pow(z - 0.5_rt, 2);
#endif
    r2 /= 0.01_rt;

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = 0.0_rt;
    }
    state(i,j,k,UFS) = std::exp(-r2);
    state(i,j,k,UFS+1)= 1.0_rt - std::exp(-r2);

}
#endif
