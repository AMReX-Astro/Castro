#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>
#include <network.H>

AMREX_INLINE
void problem_initialize ()
{
    const Geometry& dgeom = DefaultGeometry();

    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
        problem::center[idir] = 0.5_rt * (problo[idir] + probhi[idir]);
    }

    for (int idir = AMREX_SPACEDIM; idir < 3; idir++) {
        problem::center[idir] = 0.0_rt;
    }

    // we only work in NSE mode, so the mass fractions are not used by
    // the EOS

#ifndef NSE_THERMO
    amrex::Error("Error: this problem requires USE_NSE=TRUE");
#endif

}

#endif
