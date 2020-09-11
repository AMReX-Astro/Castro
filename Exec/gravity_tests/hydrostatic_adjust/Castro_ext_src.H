#ifndef CASTRO_EXT_SRC_H
#define CASTRO_EXT_SRC_H

#include <prob_parameters.H>

using namespace amrex;

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void do_ext_src(int i, int j, int k,
                GeometryData const& geomdata, 
                Array4<const Real> const& state,
                Array4<Real> const& src,
                GpuArray<Real, 3> const& center,
                const Real dt, const Real time)
{
    const Real H_0 = heating_peak;
    const Real W_0 = heating_sigma;
    const Real r_0 = heating_rad;
    const Real t_stop = heating_time;

    if (time > t_stop) return;

    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();

    if (prob_type == 1) {

        if (i == 0 && j == 0 && k == 0) {
            AllPrint() << "TIME vs TSTOP " << time << " " << t_stop << std::endl;
        }

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

        if (i == 0 && j == 0 && k == 0) {
            AllPrint() << "TIME vs TSTOP " << time << " " << t_stop << std::endl;
        }

        // sub-chandra heating -- modulate by He

        const auto ihe4 = network_spec_index("helium-4");

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
            state(i,j,k,UFS+ihe4) / state(i,j,k,URHO);

        src(i,j,k,UEINT) = state(i,j,k,URHO) * Hext;
        src(i,j,k,UEDEN) = state(i,j,k,URHO) * Hext;
    }
}

#endif