#include "Castro.H"

#include "AMReX_ParmParse.H"

#include "Problem_F.H"

#include <fstream>

using namespace amrex;

int Castro::use_stopping_criterion = 0;
Real Castro::ts_te_stopping_criterion = 1.e200;
Real Castro::T_stopping_criterion = 1.e200;

#ifdef DO_PROBLEM_POST_TIMESTEP
void Castro::problem_post_timestep() {

    if (level != 0) return;

    int finest_level = parent->finestLevel();
    Real time = state[State_Type].curTime();
    Real dt = parent->dtLevel(0);

    int jobDoneStatus = 0;

    // Note that due to the velocity criterion,
    // this stopping criterion does not make any
    // sense for the non-collision inputs.

#ifdef REACTIONS
    if (use_stopping_criterion) {

        int to_stop = 0;

        for (int lev = 0; lev <= finest_level; lev++) {

            std::unique_ptr<MultiFab> v =
                parent->getLevel(lev).derive("x_velocity", time, 0);
            std::unique_ptr<MultiFab> T =
                parent->getLevel(lev).derive("Temp", time, 0);
            std::unique_ptr<MultiFab> ts_te =
                parent->getLevel(lev).derive("t_sound_t_enuc", time, 0);

            if (lev < finest_level) {
                const MultiFab& mask = getLevel(lev + 1).build_fine_mask();
                MultiFab::Multiply(*v, mask, 0, 0, 1, 0);
                MultiFab::Multiply(*T, mask, 0, 0, 1, 0);
                MultiFab::Multiply(*ts_te, mask, 0, 0, 1, 0);
            }

            for (MFIter mfi(*v, true); mfi.isValid(); ++mfi) {

                const Box& bx = mfi.tilebox();

                check_stopping_criteria(AMREX_ARLIM_ANYD(bx.loVect()),
                                        AMREX_ARLIM_ANYD(bx.hiVect()),
                                        BL_TO_FORTRAN_ANYD((*v)[mfi]),
                                        BL_TO_FORTRAN_ANYD((*T)[mfi]),
                                        BL_TO_FORTRAN_ANYD((*ts_te)[mfi]),
                                        T_stopping_criterion,
                                        ts_te_stopping_criterion, &to_stop);
            }

            v.reset();
            T.reset();
            ts_te.reset();
        }

        // Reduction

        amrex::ParallelDescriptor::ReduceIntMax(to_stop);

        if (to_stop > 0) {

            jobDoneStatus = 1;

            amrex::Print()
                << std::endl
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

#ifdef DO_PROBLEM_POST_INIT

void Castro::problem_post_init() {

    // Read in inputs.

    ParmParse pp("castro");

    pp.query("use_stopping_criterion", use_stopping_criterion);
    pp.query("ts_te_stopping_criterion", ts_te_stopping_criterion);
    pp.query("T_stopping_criterion", T_stopping_criterion);
}

#endif

#ifdef DO_PROBLEM_POST_RESTART

void Castro::problem_post_restart() {

    // Read in inputs.

    ParmParse pp("castro");

    pp.query("use_stopping_criterion", use_stopping_criterion);
    pp.query("ts_te_stopping_criterion", ts_te_stopping_criterion);
    pp.query("T_stopping_criterion", T_stopping_criterion);
}

#endif
