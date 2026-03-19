#ifndef problem_diagnostics_H
#define problem_diagnostics_H

#include <AMReX_REAL.H>

AMREX_INLINE
void
Castro::problem_diagnostics ()
{
    int finest_level = parent->finestLevel();
    amrex::Real time = state[State_Type].curTime();

    amrex::Real rho_enuc_max = -1.e200;
    amrex::Real peak_x = -1.e200;

    for (int lev = 0; lev <= finest_level; lev++) {
        Castro& ca_lev = getLevel(lev);
        ca_lev.det_peak_properties(time, rho_enuc_max, peak_x);
    }

    ParallelDescriptor::ReduceRealMax(rho_enuc_max);
    ParallelDescriptor::ReduceRealMax(peak_x);

    const auto prob_lo = geom.ProbLoArray();
    const auto prob_hi = geom.ProbHiArray();
    const amrex::Real abort_x = prob_lo[0] +
        problem::flame_domain_abort * (prob_hi[0] - prob_lo[0]);

    if (rho_enuc_max > 0.0_rt && peak_x > abort_x) {
        amrex::Abort("flame has reached the edge of the domain");
    }
}

#endif
