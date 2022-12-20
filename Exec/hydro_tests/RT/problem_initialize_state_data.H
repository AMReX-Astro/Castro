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

    Real r[3] = {0.0_rt};

    r[0] = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

#if AMREX_SPACEDIM >= 2
    r[1] = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
#endif

#if AMREX_SPACEDIM == 3
    r[2] = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
#endif

    Real presmid  = problem::p0_base - problem::rho_1 *
        problem::split[AMREX_SPACEDIM-1];

    if (r[AMREX_SPACEDIM-1] < problem::split[AMREX_SPACEDIM-1]) {
        Real pres = problem::p0_base - problem::rho_1 * r[AMREX_SPACEDIM-1];
        state(i,j,k,UEDEN) = pres / (eos_gamma - 1.0_rt);
        state(i,j,k,UEINT) = pres / (eos_gamma - 1.0_rt);
    } else {
        Real pres = presmid - problem::rho_2 *
            (r[AMREX_SPACEDIM-1] - problem::split[AMREX_SPACEDIM-1]);
        state(i,j,k,UEDEN) = pres / (eos_gamma - 1.0_rt);
        state(i,j,k,UEINT) = pres / (eos_gamma - 1.0_rt);
    }

    Real pertheight;

#if AMREX_SPACEDIM == 2
    // We explicitly make the perturbation symmetric here.
    // This prevents the RT from bending.
    pertheight = 0.01_rt * 0.5_rt * (std::cos(2.0_rt * M_PI * r[0] / problem::L_x) +
                                     std::cos(2.0_rt * M_PI * (problem::L_x - r[0]) / problem::L_x)) + 0.5_rt;
#else
    // In 3D, the perturbation is in both x and y.
    Real r2d = amrex::min(0.5e0_rt * L_x,
                          std::sqrt((r[0] - problem::split[0]) * (r[0] - problem::split[0]) +
                                    (r[1] - problem::split[1]) * (r[1] - problem::split[1])));
    pertheight = 0.5_rt - 0.01_rt * std::cos(2.0_rt * M_PI * r2d / problem::L_x);
#endif

    state(i,j,k,URHO) = problem::rho_1 + ((problem::rho_2 - problem::rho_1) / 2.0_rt) *
        (1 + std::tanh((r[AMREX_SPACEDIM-1] - pertheight) / 0.005_rt));

    state(i,j,k,UFS) = state(i,j,k,URHO);

    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;
    state(i,j,k,UTEMP) = small_temp;

}

#endif
