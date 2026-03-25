// Implementations of functions in Problem.H go here

#include <AMReX_REAL.H>

#include <Castro.H>

using namespace amrex;

void
Castro::det_peak_properties (amrex::Real time,
                             amrex::Real& rho_enuc_max, amrex::Real& peak_x)
{
    BL_PROFILE("Castro::det_peak_properties()");

    const auto dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();

    auto mf = derive("rho_enuc", time, 0);
    BL_ASSERT(mf != nullptr);

    // find the peak energy generation rate and its index
    // these functions work on the entire MultiFab and do a reduction so
    // the information is available on all processors

    const amrex::Real level_rho_enuc_max = mf->max(0);
    const amrex::IntVect level_max_index = mf->maxIndex(0);

    // now get the physical coordinate of the maximum -- we are assuming
    // the detonation is in the x-direction

    const amrex::Real level_peak_x = prob_lo[0] +
        (static_cast<amrex::Real>(level_max_index[0]) + 0.5_rt) * dx[0];

    // find the max over all previous levels

    if (level_rho_enuc_max > rho_enuc_max) {
        rho_enuc_max = level_rho_enuc_max;
        peak_x = level_peak_x;
    }
}
