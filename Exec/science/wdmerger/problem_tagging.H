#ifndef problem_tagging_H
#define problem_tagging_H

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_tagging(int i, int j, int k,
                     Array4<char> const& tag,
                     Array4<Real const> const& state,
                     int level, const GeometryData& geomdata)
{
    const Real* problo = geomdata.ProbLo();
    const Real* probhi = geomdata.ProbHi();
    const Real* dx = geomdata.CellSize();

    Real loc[3] = {0.0};

    loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0];
#if AMREX_SPACEDIM >= 2
    loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1];
#endif
#if AMREX_SPACEDIM == 3
    loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2];
#endif

    // Clear all tagging that occurs outside the radius set by max_tagging_radius.

    Real r = std::sqrt((loc[0] - problem::center[0]) * (loc[0] - problem::center[0]) +
                       (loc[1] - problem::center[1]) * (loc[1] - problem::center[1]) +
                       (loc[2] - problem::center[2]) * (loc[2] - problem::center[2]));

    Real max_dist_lo = 0.0;
    Real max_dist_hi = 0.0;

    for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
        max_dist_lo = amrex::max(max_dist_lo, std::abs(problo[dim] - problem::center[dim]));
        max_dist_hi = amrex::max(max_dist_hi, std::abs(probhi[dim] - problem::center[dim]));
    }

    if (r > problem::max_tagging_radius * amrex::max(max_dist_lo, max_dist_hi)) {
        tag(i,j,k) = TagBox::CLEAR;
    }
}

#endif
