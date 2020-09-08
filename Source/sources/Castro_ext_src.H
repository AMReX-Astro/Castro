#ifndef CASTRO_EXT_SRC_H
#define CASTRO_EXT_SRC_H

using namespace amrex;

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void do_ext_src(int i, int j, int k,
                GeometryData const& geomdata, 
                Array4<const Real> const& sold,
                Array4<const Real> const& snew,
                Array4<Real> const& src,
                const Real dt)
{
    // const Real* problo = geomdata.ProbLo();
    // const Real* probhi = geomdata.ProbHi();
    // const Real* dx     = geomdata.CellSize();
    // const int*  domlo  = geomdata.domain.loVect();
    // const int*  domhi  = geomdata.domain.hiVect();
}

#endif