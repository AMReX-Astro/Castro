#ifndef problem_restart_H
#define problem_restart_H

#include <iostream>
#include <fstream>

AMREX_INLINE
void problem_restart (std::string checkpoint_dir)
{
    std::ifstream com;
    com.open(checkpoint_dir + "/COM", std::ios::in);

    com >> problem::center[0] >> problem::center[1] >> problem::center[2];
    com >> problem::mass_P >> problem::mass_S;
    com >> problem::com_P[0] >> problem::com_S[0];
    com >> problem::com_P[1] >> problem::com_S[1];
    com >> problem::com_P[2] >> problem::com_S[2];
    com >> problem::vel_P[0] >> problem::vel_S[0];
    com >> problem::vel_P[1] >> problem::vel_S[1];
    com >> problem::vel_P[2] >> problem::vel_S[2];
    com >> problem::t_ff_P >> problem::t_ff_S;

    com.close();



    std::ifstream relaxation;
    relaxation.open(checkpoint_dir + "/Relaxation", std::ios::in);

    if (relaxation.is_open()) {
        relaxation >> problem::relaxation_is_done;
        relaxation.close();
    }
    else {
        if (problem::problem == 1) {
            amrex::Error("No Relaxation file found in the checkpoint.");
        }
    }

    if (problem::relaxation_is_done == 1) {
        problem::relaxation_damping_factor = -1.0_rt;
    }
}

#endif
