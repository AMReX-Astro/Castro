#ifndef CASTRO_EXT_SRC_H
#define CASTRO_EXT_SRC_H

#include <prob_parameters.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void do_ext_src(int i, int j, int k,
                GeometryData const& geomdata, 
                Array4<const Real> const& state,
                Array4<Real> const& src,
                GpuArray<Real, 3> const& center,
                const Real dt, const Real time)
{
    const Real* problo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();

    for (int n = 0; n < NSRC; ++n) {
        src(i,j,k,n) = 0.0_rt;
    }

    auto y = problo[1] + (Real(j) + 0.5_rt) * dx[1];

    if (y < 1.125_rt * 4.8e8_rt) {
        auto fheat = std::sin(8.0_rt * M_PI * (y / 4.e8_rt - 1.0_rt));

        src(i,j,k,UEDEN) = heating_factor * fheat;
        src(i,j,k,UEINT) = src(i,j,k,UEDEN);
    }
}

#endif
