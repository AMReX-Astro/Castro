#include <Castro.H>

using namespace amrex;

void
Castro::pointmass_update(Real time, Real dt)
{

    amrex::ignore_unused(time);
    amrex::ignore_unused(dt);

    int finest_level = parent->finestLevel();

    if (level == finest_level && point_mass_fix_solution)
    {

        MultiFab& S_old = get_old_data(State_Type);
        MultiFab& S_new = get_new_data(State_Type);

        const auto dx = geom.CellSizeArray();
        const auto problo = geom.ProbLoArray();

        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.tilebox();

            Array4<Real const> const uin  = S_old.array(mfi);
            Array4<Real const> const uout = S_new.array(mfi);
            Array4<Real const> const vol  = volume.array(mfi);

            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                // This is just a small number to keep precision issues from making
                // icen, jcen, kcen one cell too low.
                const Real eps = 1.e-8_rt;

                // This should be the cell whose lower left corner is at center
                int icen = static_cast<int>(std::floor((problem::center[0] - problo[0]) / dx[0] + eps));
#if AMREX_SPACEDIM >= 2
                int jcen = static_cast<int>(std::floor((problem::center[1] - problo[1]) / dx[1] + eps));
#else
                int jcen = 0;
#endif
#if AMREX_SPACEDIM == 3
                int kcen = static_cast<int>(std::floor((problem::center[2] - problo[2]) / dx[2] + eps));
#endif

                // Make sure we only count contributions from this grid

                const int box_size = 2;

                int istart = amrex::max(icen - box_size, bx.smallEnd(0));
                int iend = amrex::min(icen + box_size - 1, bx.bigEnd(0));
#if AMREX_SPACEDIM >= 2
                int jstart = amrex::max(jcen - box_size, bx.smallEnd(1));
                int jend = amrex::min(jcen + box_size - 1, bx.bigEnd(1));
#else
                int jstart = 0;
                int jend = 0;
#endif
#if AMREX_SPACEDIM == 3
                int kstart = amrex::max(kcen - box_size, bx.smallEnd(2));
                int kend = amrex::min(kcen + box_size - 1, bx.bigEnd(2));
#else
                int kstart = 0;
                int kend = 0;
#endif

                Real delta_mass_tmp = 0.0_rt;

                if (i >= istart && i <= iend &&
                    j >= jstart && j <= jend &&
                    k >= kstart && k <= kend) {

                    delta_mass_tmp = vol(i,j,k) * (uout(i,j,k,URHO) - uin(i,j,k,URHO));

                }

                return delta_mass_tmp;

            });

        }

        ReduceTuple hv = reduce_data.value();
        Real mass_change_at_center = amrex::get<0>(hv);

        ParallelDescriptor::ReduceRealSum(mass_change_at_center);

        if (mass_change_at_center > 0.0)
        {

            if (verbose > 1) {
                amrex::Print() << "  Updating point mass from " << point_mass << ", by " << mass_change_at_center
                               << ", to " << point_mass + mass_change_at_center << std::endl << std::endl;
            }

            point_mass += mass_change_at_center;

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();

                Array4<Real const> const uin = S_old.array(mfi);
                Array4<Real> const uout = S_new.array(mfi);

                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // This is just a small number to keep precision issues from making
                    // icen, jcen, kcen one cell too low.

                    const Real eps = 1.e-8_rt;

                    // This should be the cell whose lower left corner is at center
                    int icen = static_cast<int>(std::floor((problem::center[0] - problo[0]) / dx[0] + eps));
#if AMREX_SPACEDIM >= 2
                    int jcen = static_cast<int>(std::floor((problem::center[1] - problo[1]) / dx[1] + eps));
#else
                    int jcen = 0;
#endif
#if AMREX_SPACEDIM == 3
                    int kcen = static_cast<int>(std::floor((problem::center[2] - problo[2]) / dx[2] + eps));
#endif

                    // Make sure we only count contributions from this grid

                    const int box_size = 2;

                    int istart = amrex::max(icen - box_size, bx.smallEnd(0));
                    int iend = amrex::min(icen + box_size - 1, bx.bigEnd(0));

#if AMREX_SPACEDIM >= 2
                    int jstart = amrex::max(jcen - box_size, bx.smallEnd(1));
                    int jend = amrex::min(jcen + box_size - 1, bx.bigEnd(1));
#else
                    int jstart = 0;
                    int jend = 0;
#endif

# if AMREX_SPACEDIM == 3
                    int kstart = amrex::max(kcen - box_size, bx.smallEnd(2));
                    int kend = amrex::min(kcen + box_size - 1, bx.bigEnd(2));
#else
                    int kstart = 0;
                    int kend = 0;
#endif

                    if (i >= istart && i <= iend &&
                        j >= jstart && j <= jend &&
                        k >= kstart && k <= kend) {

                        for (int n = 0; n < NUM_STATE; ++n) {
                            uout(i,j,k,n) = uin(i,j,k,n);
                        }

                    }

                });
            }
        }
    }
}
