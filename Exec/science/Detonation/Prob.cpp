#include "Castro.H"

#include "AMReX_ParmParse.H"

#include <fstream>

using namespace amrex;

int Castro::use_stopping_criterion = 0;
Real Castro::ts_te_stopping_criterion = 1.e200;
Real Castro::T_stopping_criterion = 1.e200;

#ifdef DO_PROBLEM_POST_TIMESTEP
void
Castro::problem_post_timestep()
{

    if (level != 0) return;

    int finest_level = parent->finestLevel();
    Real time = state[State_Type].curTime();
    Real dt = parent->dtLevel(0);

    int jobDoneStatus = 0;

    if (use_stopping_criterion) {

        // Compute extrema

        bool local_flag = true;

        Real T_max     = 0.0;
        Real ts_te_max = 0.0;

        for (int lev = 0; lev <= finest_level; lev++) {

            MultiFab& S_new = parent->getLevel(lev).get_new_data(State_Type);

            T_max = std::max(T_max, S_new.max(Temp, 0, local_flag));

#ifdef REACTIONS
            if (lev == finest_level) {

                auto ts_te_MF = parent->getLevel(lev).derive("t_sound_t_enuc", time, 0);
                ts_te_max = std::max(ts_te_max, ts_te_MF->max(0,0,local_flag));

            }
#endif

        }

        // Max reductions

        amrex::ParallelDescriptor::ReduceRealMax({T_max, ts_te_max});

        if (ts_te_max >= ts_te_stopping_criterion) {

            jobDoneStatus = 1;

            amrex::Print() << std::endl
                           << "Ending simulation because we are above the threshold for unstable burning."
                           << std::endl;

        }

        if (T_max >= T_stopping_criterion) {

            jobDoneStatus = 1;

            amrex::Print() << std::endl
                           << "Ending simulation because we are above the temperature threshold."
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
