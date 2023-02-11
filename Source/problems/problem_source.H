#ifndef problem_source_H
#define problem_source_H

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_source ([[maybe_unused]] int i,
                     [[maybe_unused]] int j,
                     [[maybe_unused]] int k,
                     [[maybe_unused]] GeometryData const& geomdata,
                     [[maybe_unused]] Array4<const Real> const& state,
                     [[maybe_unused]] Array4<Real> const& src,
                     [[maybe_unused]] const Real dt,
                     [[maybe_unused]] const Real time) {}

#endif
