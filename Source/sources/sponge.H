#ifndef SPONGE_H
#define SPONGE_H

#ifdef SPONGE

#include <Castro_util.H>

namespace sponge
{
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    int in_sponge_region (int i, int j, int k,
                          const GeometryData& geomdata,
                          Real density, Real pressure)
    {
        GpuArray<Real, 3> loc;
        position(i, j, k, geomdata, loc);

        for (int d = 0; d < 3; ++d) {
            loc[d] -= problem::center[d];
        }

        Real r = std::sqrt(loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2]);

        if ((sponge_lower_radius > 0.0 && r >= sponge_lower_radius) ||
            (sponge_upper_density > 0.0 && density <= sponge_upper_density) ||
            (sponge_upper_pressure > 0.0 && pressure <= sponge_upper_pressure)) {
            return 1;
        }
        else {
            return 0;
        }
    }
}

#endif

#endif
