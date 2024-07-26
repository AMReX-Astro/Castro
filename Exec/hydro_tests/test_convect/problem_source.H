#ifndef problem_source_H
#define problem_source_H

#include <prob_parameters.H>

using namespace amrex;

constexpr Real y_layer = 1.25e8_rt;
constexpr Real L_x = 2.5e8_rt;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_source (int i, int j, int k,
                     GeometryData const& geomdata,
                     Array4<const Real> const& state,
                     Array4<Real> const& src,
                     const Real dt, const Real time)
{
    using namespace problem;

    const Real* problo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();

    auto x = (Real(i)+0.5e0_rt)*dx[0] + problo[0];
    auto y = (Real(j)+0.5e0_rt)*dx[1] + problo[1];
    auto ey = std::exp(-(y-y_layer)*(y-y_layer)/1.e14_rt);

    auto H = ey * (1.e0_rt +
        0.00625_rt * std::sin( 2*M_PI*x/L_x)
        + 0.01875_rt * std::sin((6*M_PI*x/L_x) + M_PI/3.e0_rt)
        + 0.01250_rt * std::sin((8*M_PI*x/L_x) + M_PI/5.e0_rt));

    // Source terms
    src(i,j,k,UEDEN) = state(i,j,k,URHO) * H * 2.5e16_rt;
    src(i,j,k,UEINT) = src(i,j,k,UEDEN);
}

#endif
