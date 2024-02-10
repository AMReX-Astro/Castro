#ifndef problem_initialize_H
#define problem_initialize_H

AMREX_INLINE
void problem_initialize()
{
    const Geometry& dgeom = DefaultGeometry();

    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    for (int n=0; n < AMREX_SPACEDIM; ++n)
    {
        problem::center[n] = 0.5_rt * (problo[n] + probhi[n]);
    }



}

#endif
