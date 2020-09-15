#ifndef CASTRO_EXT_SRC_H
#define CASTRO_EXT_SRC_H

// #include <prob_parameters.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void do_ext_src(int i, int j, int k,
                GeometryData const& geomdata, 
                Array4<const Real> const& state,
                Array4<Real> const& src,
                GpuArray<Real, 3> const& center,
                const Real dt, const Real time)
{
    // const Real* problo = geomdata.ProbLo();
    // const Real* probhi = geomdata.ProbHi();
    // const Real* dx     = geomdata.CellSize();
    // const int*  domlo  = geomdata.domain.loVect();
    // const int*  domhi  = geomdata.domain.hiVect();
}

#endif
