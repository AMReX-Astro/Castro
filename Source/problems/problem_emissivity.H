#ifndef problem_emissivity_H
#define problem_emissivity_H

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_emissivity (int i, int j, int k, int g,
                         const GpuArray<Real, NGROUPS>& nu,
                         const GpuArray<Real, NGROUPS+1>& xnu,
                         Real T, Real kg, Real dkdT,
                         Real& jg, Real& djdT) {}

#endif
