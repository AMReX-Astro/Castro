#ifndef problem_rad_source_H
#define problem_rad_source_H

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_rad_source (int i, int j, int k,
                         Array4<Real> const& rhs,
                         const GeometryData& geomdata, Real time,
                         Real dt, int igroup)
{
    const Real x0 = 0.5e0_rt;
    const Real t0 = 3.3356409519815202e-10_rt;
    const Real qn = 1.134074546528399e20_rt;

    GpuArray<Real, 3> loc;
    position(i, j, k, geomdata, loc);

    Real x = loc[0];

    if (time <= t0 && std::abs(x) <= x0) {
        rhs(i,j,k) = rhs(i,j,k) + qn; // (qn / dt) * dt
    }
}

#endif
