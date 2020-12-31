#ifndef problem_tagging_H
#define problem_tagging_H

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_tagging(int i, int j, int k,
                     Array4<char> const& tag,
                     Array4<Real const> const& state,
                     int level, const GeometryData& geomdata)
{
    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real y = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1];

#if AMREX_SPACEDIM == 3
    Real z = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2];
#else
    Real z = 0.0;
#endif

    Real height;

    if (AMREX_SPACEDIM == 2) {
        height = y;
    } else {
        height = z;
    }

    if (height > problem::refine_cutoff_height) {
        tag(i,j,k) = TagBox::CLEAR;
    }
}

#endif
