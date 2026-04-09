#ifndef problem_source_H
#define problem_source_H

#include <prob_parameters.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_source (int i, int j, int k,
                     GeometryData const& geomdata,
                     Array4<const Real> const& state,
                     Array4<Real> const& src,
                     const Real dt, const Real time)
{

    amrex::ignore_unused(dt);

    using namespace problem;

    const Real H_0 = heating_peak;
    const Real W_0 = heating_sigma;
    const Real r_0 = heating_rad;
    const Real t_stop = heating_time;

    if (time > t_stop) return;

    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();

    if (prob_type == 1) {

#ifndef AMREX_USE_GPU
        if (i == 0 && j == 0 && k == 0) {
            AllPrint() << "TIME vs TSTOP " << time << " " << t_stop << std::endl;
        }
#endif

        // For heating at the center
        auto x = prob_lo[0] + (Real(i) + 0.5_rt)*dx[0] - center[0];
#if AMREX_SPACEDIM == 1
        auto y = 0.0_rt;
#else
        auto y = prob_lo[1] + (Real(j) + 0.5_rt)*dx[1] - center[1];
#endif
#if AMREX_SPACEDIM < 3
        auto z = 0.0_rt;
#else
        auto z = prob_lo[2] + (Real(k) + 0.5_rt)*dx[2] - center[2];
#endif

        auto dist = std::sqrt(x*x + y*y + z*z);

        auto Hext = H_0 * std::exp(-((dist - r_0)*(dist - r_0))/(W_0*W_0));

        src(i,j,k,UEINT) = state(i,j,k,URHO) * Hext;
        src(i,j,k,UEDEN) = state(i,j,k,URHO) * Hext;

    } else if (prob_type == 3) {

#ifndef AMREX_USE_GPU
        if (i == 0 && j == 0 && k == 0) {
            AllPrint() << "TIME vs TSTOP " << time << " " << t_stop << std::endl;
        }
#endif

        // sub-chandra heating -- modulate by He

        auto x = prob_lo[0] + (Real(i) + 0.5_rt)*dx[0] - center[0];
#if AMREX_SPACEDIM == 1
        auto y = 0.0_rt;
#else
        auto y = prob_lo[1] + (Real(j) + 0.5_rt)*dx[1] - center[1];
#endif
#if AMREX_SPACEDIM < 3
        auto z = 0.0_rt;
#else
        auto z = prob_lo[2] + (Real(k) + 0.5_rt)*dx[2] - center[2];
#endif

        auto dist = std::sqrt(x*x + y*y + z*z);

        auto Hext = H_0 * std::exp(-((dist - r_0)*(dist - r_0))/(W_0*W_0)) *
            state(i,j,k,UFS+problem::ifuel) / state(i,j,k,URHO);

        src(i,j,k,UEINT) = state(i,j,k,URHO) * Hext;
        src(i,j,k,UEDEN) = state(i,j,k,URHO) * Hext;
    }
}

#endif
