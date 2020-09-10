#ifndef CASTRO_EXT_SRC_H
#define CASTRO_EXT_SRC_H

#include <prob_parameters.H>

using namespace amrex;

const Real H_0 = heating_peak;
const Real W_0 = heating_sigma;
const Real r_0 = heating_rad;
const Real t_stop = heating_time;

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void do_ext_src(int i, int j, int k,
                GeometryData const& geomdata, 
                Array4<const Real> const& state,
                Array4<Real> const& src,
                GpuArray<Real, 3> center,
                const Real dt)
{
    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();

    if (prob_type == 1) {

        // For heating at the center
        auto x = prob_lo[0] + (Real(i) + 0.5_rt)*dx[0] - center[0];
        auto y = prob_lo[1] + (Real(j) + 0.5_rt)*dx[1] - center[1];
        auto z = prob_lo[2] + (Real(k) + 0.5_rt)*dx[2] - center[2];

        auto dist = std::sqrt(x*x + y*y + z*z);

        auto Hext = H_0 * std::exp(-((dist - r_0)*(dist - r_0))/(W_0*W_0));

        src(i,j,k,UEINT) = state(i,j,k,URHO) * Hext;
        src(i,j,k,UEDEN) = state(i,j,k,URHO) * Hext;

    } else if (prob_type == 3) {

        // sub-chandra heating -- modulate by He

        const auto ihe4 = network_spec_index("helium-4");

        auto x = prob_lo[0] + (Real(i) + 0.5_rt)*dx[0] - center[0];
        auto y = prob_lo[1] + (Real(j) + 0.5_rt)*dx[1] - center[1];
        auto z = prob_lo[2] + (Real(k) + 0.5_rt)*dx[2] - center[2];

        auto dist = std::sqrt(x*x + y*y + z*z);

        auto Hext = H_0 * std::exp(-((dist - r_0)*(dist - r_0))/(W_0*W_0)) * 
            state(i,j,k,UFS-1+ihe4) / state(i,j,k,URHO);

        src(i,j,k,UEINT) = state(i,j,k,URHO) * Hext;
        src(i,j,k,UEDEN) = state(i,j,k,URHO) * Hext;
    }
}

#endif