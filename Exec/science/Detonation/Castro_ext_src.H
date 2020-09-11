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

    // Add a mock gravitational acceleration which points to the center
    // with uniform magnitude on either side of the center.

    const Real* problo = geomdata.ProbLo();
    const Real* probhi = geomdata.ProbHi();
    const Real* dx     = geomdata.CellSize();

    auto c_T = problo[0] + center_T * (probhi[0] - problo[0]);

    auto x = problo[0] + (Real(i) + 0.5_rt) * dx[0] - center[0];

    auto g = (x > c_T) ? -grav_acceleration : grav_acceleration;

    src(i,j,k,UMX) = state(i,j,k,URHO) * g;

    // Energy source is v . momentum source

    src(i,j,k,UEDEN) = state(i,j,k,UMX) * g;
}

#endif