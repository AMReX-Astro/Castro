#include <Castro.H>

#include <AMReX_ParmParse.H>

#include <Problem_F.H>

#include <fstream>

using namespace amrex;

#ifdef DO_PROBLEM_POST_TIMESTEP
void
Castro::problem_post_timestep()
{

    if (level != 0) return;

    int finest_level = parent->finestLevel();
    Real time = state[State_Type].curTime();
    Real dt = parent->dtLevel(0);

    int jobDoneStatus = 0;

    // Note that due to the velocity criterion,
    // this stopping criterion does not make any
    // sense for the non-collision inputs.

#ifdef REACTIONS
    if (problem::use_early_stopping) {

        int to_stop = 0;

        for (int lev = 0; lev <= finest_level; lev++) {

            std::unique_ptr<MultiFab> v = parent->getLevel(lev).derive("x_velocity", time, 0);
            std::unique_ptr<MultiFab> T = parent->getLevel(lev).derive("Temp", time, 0);
            std::unique_ptr<MultiFab> ts_te = parent->getLevel(lev).derive("t_sound_t_enuc", time, 0);

            if (lev < finest_level) {
                const MultiFab& mask = getLevel(lev+1).build_fine_mask();
                MultiFab::Multiply(*v, mask, 0, 0, 1, 0);
                MultiFab::Multiply(*T, mask, 0, 0, 1, 0);
                MultiFab::Multiply(*ts_te, mask, 0, 0, 1, 0);
            }

            ReduceOps<ReduceOpMax> reduce_op;
            ReduceData<int> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(*v, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

                const Box& bx = mfi.tilebox();

                auto xvel = (*v)[mfi].array();
                auto temp = (*T)[mfi].array();
                auto tste = (*ts_te)[mfi].array();

                // Note that this stopping criterion only makes sense for the collision problem.

                reduce_op.eval(bx, reduce_data,
                [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
                {
                    int stop = 0;

                    if (problem::vel > 0.0 && xvel(i,j,k) > 0.0) {
                        if (temp(i,j,k) >= problem::T_stop_cutoff || tste(i,j,k) >= problem::ts_te_stop_cutoff) {
                            stop = 1;
                        }
                    }

                    return {stop};
                });

            }

            ReduceTuple hv = reduce_data.value();
            to_stop = amrex::max(to_stop, amrex::get<0>(hv));

            v.reset();
            T.reset();
            ts_te.reset();

        }

        // Reduction

        amrex::ParallelDescriptor::ReduceIntMax(to_stop);

        if (to_stop > 0) {

            jobDoneStatus = 1;

            amrex::Print() << std::endl
                           << "Ending simulation because we are above the set threshold."
                           << std::endl;

        }

    }

    if (jobDoneStatus == 1) {

        // Write out a checkpoint. Note that this will
        // only happen if you have amr.message_int = 1.

        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::ofstream dump_file;
            dump_file.open("dump_and_stop", std::ofstream::out);
            dump_file.close();

            // Also write out a file signifying that we're done with the simulation.

            std::ofstream jobDoneFile;
            jobDoneFile.open("jobIsDone", std::ofstream::out);
            jobDoneFile.close();
        }

    }
#endif
}
#endif
