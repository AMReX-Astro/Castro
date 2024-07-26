#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{
    int coord_type = geomdata.Coord();

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real xx = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

    Real yy = 0.0;
#if AMREX_SPACEDIM >= 2
    yy = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
#endif

    Real zz = 0.0;
#if AMREX_SPACEDIM == 3
    zz = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
#endif

    Real dist = 0;
    if (problem::init_as_1d == 0) {
        dist = std::sqrt((problem::center[0] - xx) * (problem::center[0] - xx) +
                         (problem::center[1] - yy) * (problem::center[1] - yy) +
                         (problem::center[2] - zz) * (problem::center[2] - zz));

    } else if (problem::init_as_1d == 1) {
        dist = std::abs(problem::center[0] - xx);

    } else if (problem::init_as_1d == 2) {
        dist = std::abs(problem::center[1] - yy);

    } else if (problem::init_as_1d == 3) {
        dist = std::abs(problem::center[2] - zz);

    } else {
#ifndef AMREX_USE_GPU
        amrex::Error("invalid init_as_1d");
#endif
    }

    if (dist <= 0.5_rt) {
        state(i,j,k,URHO) = problem::rho0 + problem::drho0 * std::exp(-16.0_rt * dist * dist) *
            std::pow(std::cos(M_PI*dist), 6);

    } else {
        state(i,j,k,URHO) = problem::rho0;
    }

    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

    // we are isentropic, so p = (dens/rho0)**Gamma_1
    Real p = std::pow(state(i,j,k,URHO) / problem::rho0, eos_gamma);
    Real eint = p / (eos_gamma - 1.0_rt);

    state(i,j,k,UEDEN) = eint;
    state(i,j,k,UEINT) = eint;

    Real xn_zone[NumSpec] = {0.0};
    xn_zone[0] = 1.0_rt;

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * xn_zone[n];
    }
}

#endif
