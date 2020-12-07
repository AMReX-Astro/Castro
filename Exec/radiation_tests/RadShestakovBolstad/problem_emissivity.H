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
    const Real Tf = 1.1604518737857157e6_rt; // 1 (dimensionless)

    Real hoverk = C::hplanck / C::k_B;
    Real cB = 8.0_rt * M_PI * C::k_B / std::pow(C::c_light, 3);

    Real num = xnu[g];
    Real nup = xnu[g+1];

    Real dBdT = cB * std::pow(nu[g], 3) * (std::exp(-hoverk * num / Tf) -
                                           std::exp(-hoverk * nup / Tf));
    Real Bg = dBdT * T;

    jg = Bg * kg;
    djdT = dkdT * Bg + dBdT * kg;
}

#endif
