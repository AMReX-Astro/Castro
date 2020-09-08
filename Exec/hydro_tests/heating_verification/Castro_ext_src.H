#ifndef CASTRO_EXT_SRC_H
#define CASTRO_EXT_SRC_H

#include "prob_parameters.H"

using namespace amrex;

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void do_ext_src(int i, int j, int k,
                GeometryData const& geomdata, 
                Array4<const Real> const& sold,
                Array4<const Real> const& snew,
                Array4<Real> const& src,
                const Real dt)
{
    const Real* problo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();

    auto x = problo[0] + (Real(i) + 0.5_rt) * dx[0];

    src(i,j,k,UEINT) = snew(i,j,k,URHO) * 
        (A1*std::exp(-(x - 0.5_rt)*(x - 0.5_rt)/(sigma1*sigma1)) + 
         A2*std::exp(-(x - 1.0_rt)*(x - 1.0_rt)/(sigma2*sigma2)));

    src(i,j,k,UEDEN) = src(i,j,k,UEINT); 
}

#endif