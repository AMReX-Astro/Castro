#ifndef problem_emissivity_H
#define problem_emissivity_H

#include <fundamental_constants.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_emissivity (int i, int j, int k, int g,
                         const GpuArray<Real, NGROUPS>& nu,
                         const GpuArray<Real, NGROUPS+1>& xnu,
                         Real T, Real kg, Real dkdT,
                         Real& jg, Real& djdT)
{
    static_assert(NGROUPS == 2, "NGROUPS must be 2 for this problem");

    Real pfc[2] = {0.5_rt, 0.5_rt};

    Real Bg = C::a_rad * std::pow(T, 4);
    Real dBdT = 4.e0_rt * C::a_rad * std::pow(T, 3);

    jg = pfc[g] * Bg * kg;
    djdT = pfc[g] * (dkdT * Bg + dBdT * kg);
}

#endif
