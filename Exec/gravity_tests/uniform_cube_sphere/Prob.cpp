#include <Castro.H>
#include <Castro_F.H>

#include <Gravity.H>
#include <Gravity_F.H>

#include <AMReX_ParmParse.H>

#include <AMReX_buildInfo.H>

#include <prob_parameters.H>
#include <prob_parameters_F.H>

#include <fundamental_constants.H>

using namespace amrex;

void Castro::problem_post_init()
{
    BL_ASSERT(level == 0);

    // If we're doing problem 2, the normalized sphere,
    // add up the mass on the domain and then update
    // the density in the sphere so that it has the 'correct'
    // amount of total mass.

    if (problem::problem == 2) {

        Real actual_mass = 0.0;

        bool local_flag = true;
        Real time       = state[State_Type].curTime();

        for (int lev = 0; lev <= parent->finestLevel(); ++lev)
            actual_mass += getLevel(lev).volWgtSum("density", time, local_flag);

        ParallelDescriptor::ReduceRealSum(actual_mass);

        // The correct amount of mass is the mass of a sphere
        // with the given diameter and density.

        Real target_mass = problem::density * (1.0e0 / 6.0e0) * M_PI * std::pow(problem::diameter, 3);

        Real update_factor = target_mass / actual_mass;

        // Now update the density given this factor.

        amrex::Print() << "\n";
        amrex::Print() << "  Updating density by the factor " << update_factor << " to ensure total mass matches target mass.\n";
        amrex::Print() << "\n";

        problem::density = problem::density * update_factor;

        for (int lev = 0; lev <= parent->finestLevel(); lev++)
        {
            MultiFab& state = getLevel(lev).get_new_data(State_Type);

            const auto dx = getLevel(lev).geom.CellSizeArray();
            const auto problo = getLevel(lev).geom.ProbLoArray();

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

                const Box& box = mfi.tilebox();

                auto u = state[mfi].array();

                // Update the density field. This ensures that the sum of the
                // mass on the domain is what we intend it to be.

                amrex::ParallelFor(box,
                [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
                {
                    if (problem::problem != 2) return;

                    Real xx = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt) - problem::center[0];

#if AMREX_SPACEDIM >= 2
                    Real yy = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt) - problem::center[1];
#else
                    Real yy = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
                    Real zz = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt) - problem::center[2];
#else
                    Real zz = 0.0_rt;
#endif

                    if (std::sqrt(xx * xx + yy * yy + zz * zz) < problem::diameter / 2) {
                        u(i,j,k,URHO) *= update_factor;
                        for (int n = 0; n < NumSpec; ++n) {
                            u(i,j,k,UFS+n) *= update_factor;
                        }
                    }
                });
            }

        }

        // Do a final check, for sanity purposes.

        actual_mass = 0.0;

        for (int lev = 0; lev <= parent->finestLevel(); lev++)
            actual_mass += getLevel(lev).volWgtSum("density", time, local_flag);

        ParallelDescriptor::ReduceRealSum(actual_mass);

        if (std::abs( (actual_mass - target_mass) / target_mass ) > 1.0e-6) {
            amrex::Print() << "\n";
            amrex::Print() << "Actual mass: " << actual_mass << "\n";
            amrex::Print() << "Target mass: " << target_mass << "\n";
            amrex::Print() << "\n";
            amrex::Abort("Sphere does not have the right amount of mass.");
        }

    }

    gravity->multilevel_solve_for_new_phi(0, parent->finestLevel());

    const int norm_power = 2;

    ReduceOps<ReduceOpSum, ReduceOpSum> reduce_op;
    ReduceData<Real, Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (int lev = 0; lev <= parent->finestLevel(); lev++)
    {
        const auto dx = getLevel(lev).geom.CellSizeArray();
        const auto problo = getLevel(lev).geom.ProbLoArray();

        const Real time = getLevel(lev).state[State_Type].curTime();

        auto phiGrav = getLevel(lev).derive("phiGrav", time, 0);

        if (lev < parent->finestLevel())
        {
            const MultiFab& mask = getLevel(lev+1).build_fine_mask();
            MultiFab::Multiply(*phiGrav, mask, 0, 0, 1, 0);
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(*phiGrav, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.tilebox();

            auto phi = (*phiGrav)[mfi].array();
            auto vol = getLevel(lev).Volume()[mfi].array();

            // Compute the norm of the difference between the calculated potential
            // and the analytical solution.

            reduce_op.eval(box, reduce_data,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                Real radius = 0.5_rt * problem::diameter;
                Real mass = (4.0_rt / 3.0_rt) * M_PI * radius * radius * radius * problem::density;

                Real xx = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt) - problem::center[0];

#if AMREX_SPACEDIM >= 2
                Real yy = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt) - problem::center[1];
#else
                Real yy = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
                Real zz = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt) - problem::center[2];
#else
                Real zz = 0.0_rt;
#endif

                Real rr = std::sqrt(xx * xx + yy * yy + zz * zz);

                Real phiExact = 0.0_rt;

                if (problem::problem == 1 || problem::problem == 2) {

                    if (rr <= radius) {
                        phiExact = -C::Gconst * mass * (3 * radius * radius - rr * rr) / (2 * radius * radius * radius);
                    } else {
                        phiExact = -C::Gconst * mass / rr;
                    }

                }
                else if (problem::problem == 3) {

                    Real x[2];
                    Real y[2];
                    Real z[2];

                    x[0] = radius + xx;
                    x[1] = radius - xx;
                    y[0] = radius + yy;
                    y[1] = radius - yy;
                    z[0] = radius + zz;
                    z[1] = radius - zz;

                    for (int ii = 0; ii <= 1; ++ii) {
                        for (int jj = 0; jj <= 1; ++jj) {
                            for (int kk = 0; kk <= 1; ++kk) {

                                Real r = std::sqrt(x[ii] * x[ii] + y[jj] * y[jj] + z[kk] * z[kk]);

                                // This is Equation 20 in Katz et al. (2016). There is a special case
                                // where, e.g., x[ii] = y[jj] = 0, in which case z[kk] = r and the
                                // atanh will evaluate to infinity, making the product ill-defined.
                                // Handle this case by only doing the atanh evaluation away from those
                                // points, since the product should be zero anyway. We also want to
                                // avoid a divide by zero, so we'll skip all the cases where r = 0.

                                if (r / radius > 1.e-6_rt) {

                                    if (std::abs(x[ii]) / radius > 1.e-6_rt && std::abs(y[jj]) / radius > 1.e-6_rt) {
                                        phiExact -= C::Gconst * problem::density * (x[ii] * y[jj] * std::atanh(z[kk] / r));
                                    }

                                    if (std::abs(y[jj]) / radius > 1.e-6_rt && std::abs(z[kk]) / radius > 1.e-6_rt) {
                                        phiExact -= C::Gconst * problem::density * (y[jj] * z[kk] * std::atanh(x[ii] / r));
                                    }

                                    if (std::abs(z[kk]) / radius > 1.e-6_rt && std::abs(x[ii]) / radius > 1.e-6_rt) {
                                        phiExact -= C::Gconst * problem::density * (z[kk] * x[ii] * std::atanh(y[jj] / r));
                                    }

                                    // Also, avoid a divide-by-zero for the atan terms.

                                    if (std::abs(x[ii]) / radius > 1.e-6_rt) {
                                        phiExact += C::Gconst * problem::density * (x[ii] * x[ii] / 2.0_rt * std::atan(y[jj] * z[kk] / (x[ii] * r)));
                                    }

                                    if (std::abs(y[jj]) / radius > 1.e-6_rt) {
                                        phiExact += C::Gconst * problem::density * (y[jj] * y[jj] / 2.0_rt * std::atan(z[kk] * x[ii] / (y[jj] * r)));
                                    }

                                    if (std::abs(z[kk]) / radius > 1.e-6_rt) {
                                        phiExact += C::Gconst * problem::density * (z[kk] * z[kk] / 2.0_rt * std::atan(x[ii] * y[jj] / (z[kk] * r)));
                                    }

                                }
                            }
                        }
                    }

                }

                Real norm_diff = vol(i,j,k) * std::pow(phi(i,j,k) - phiExact, norm_power);
                Real norm_exact = vol(i,j,k) * std::pow(phiExact, norm_power);

                return {norm_diff, norm_exact};
            });
        }

    }

    ReduceTuple hv = reduce_data.value();
    Real norm_diff = amrex::get<0>(hv);
    Real norm_exact = amrex::get<1>(hv);

    ParallelDescriptor::ReduceRealSum(norm_diff);
    ParallelDescriptor::ReduceRealSum(norm_exact);

    norm_diff = std::pow(norm_diff, 1.0 / norm_power);
    norm_exact = std::pow(norm_exact, 1.0 / norm_power);

    const Real error = norm_diff / norm_exact;

    amrex::Print() << std::endl;
    amrex::Print() << "Error = " << error << std::endl;
    amrex::Print() << std::endl;

}
