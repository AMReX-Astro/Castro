#ifndef GRAVITY_UTIL_H
#define GRAVITY_UTIL_H

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real direct_sum_symmetric_add(const GpuArray<Real, 3>& loc, const GpuArray<Real, 3>& locb,
                              const GpuArray<Real, 3>& problo, const GpuArray<Real, 3>& probhi,
                              Real rho, Real dV,
                              const GpuArray<bool, 3>& doSymmetricAddLo, const GpuArray<bool, 3>& doSymmetricAddHi)
{
    // Add contributions from any symmetric boundaries.

    Real bcTerm = 0.0_rt;

    Real x, y, z, r;

    if (doSymmetricAddLo[0]) {

        x = 2.0_rt * problo[0] - loc[0];
        y = loc[1];
        z = loc[2];

        r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

        bcTerm -= C::Gconst * rho * dV / r;

        if (doSymmetricAddLo[1]) {

            x = 2.0_rt * problo[0] - loc[0];
            y = 2.0_rt * problo[1] - loc[1];
            z = loc[2];

            r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

            bcTerm -= C::Gconst * rho * dV / r;

        }

        if (doSymmetricAddLo[2]) {

            x = 2.0_rt * problo[0] - loc[0];
            y = loc[1];
            z = 2.0_rt * problo[2] - loc[2];

            r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

            bcTerm -= C::Gconst * rho * dV / r;

        }

        if (doSymmetricAddLo[1] && doSymmetricAddLo[2]) {

            x = 2.0_rt * problo[0] - loc[0];
            y = 2.0_rt * problo[1] - loc[1];
            z = 2.0_rt * problo[2] - loc[2];

            r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

            bcTerm -= C::Gconst * rho * dV / r;

        }

    }

    if (doSymmetricAddLo[1]) {

        x = loc[0];
        y = 2.0_rt * problo[1] - loc[1];
        z = loc[2];

        r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

        bcTerm -= C::Gconst * rho * dV / r;

        if (doSymmetricAddLo[2]) {

            x = loc[0];
            y = 2.0_rt * problo[1] - loc[1];
            z = 2.0_rt * problo[2] - loc[2];

            r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

            bcTerm -= C::Gconst * rho * dV / r;

        }

    }

    if (doSymmetricAddLo[2]) {

        x = loc[0];
        y = loc[1];
        z = 2.0_rt * problo[2] - loc[2];

        r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

        bcTerm -= C::Gconst * rho * dV / r;

    }



    if (doSymmetricAddHi[0]) {

        x = 2.0_rt * probhi[0] - loc[0];
        y = loc[1];
        z = loc[2];

        r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

        bcTerm -= C::Gconst * rho * dV / r;

        if (doSymmetricAddHi[1]) {

            x = 2.0_rt * probhi[0] - loc[0];
            y = 2.0_rt * probhi[1] - loc[1];
            z = loc[2];

            r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

            bcTerm -= C::Gconst * rho * dV / r;

        }

       if (doSymmetricAddHi[2]) {

           x = 2.0_rt * probhi[0] - loc[0];
           y = loc[1];
           z = 2.0_rt * probhi[2] - loc[2];

           r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

           bcTerm -= C::Gconst * rho * dV / r;

       }

       if (doSymmetricAddHi[1] && doSymmetricAddHi[2]) {

           x = 2.0_rt * probhi[0] - loc[0];
           y = 2.0_rt * probhi[1] - loc[1];
           z = 2.0_rt * probhi[2] - loc[2];

           r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

           bcTerm -= C::Gconst * rho * dV / r;

       }

    }

    if (doSymmetricAddHi[1]) {

        x = loc[0];
        y = 2.0_rt * probhi[1] - loc[1];
        z = loc[2];

        r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

        bcTerm -= C::Gconst * rho * dV / r;

        if (doSymmetricAddHi[2]) {

            x = loc[0];
            y = 2.0_rt * probhi[1] - loc[1];
            z = 2.0_rt * probhi[2] - loc[2];

            r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

            bcTerm -= C::Gconst * rho * dV / r;

        }

    }

    if (doSymmetricAddHi[2]) {

        x = loc[0];
        y = loc[1];
        z = 2.0_rt * probhi[2] - loc[2];

        r = std::sqrt((x - locb[0]) * (x - locb[0]) + (y - locb[1]) * (y - locb[1]) + (z - locb[2]) * (z - locb[2]));

        bcTerm -= C::Gconst * rho * dV / r;

    }

    return bcTerm;

}

#endif
