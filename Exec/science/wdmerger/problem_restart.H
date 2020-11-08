#ifndef problem_restart_H
#define problem_restart_H

#include <iostream>
#include <fstream>

#include <wdmerger_F.H>

AMREX_INLINE
void problem_restart(std::string checkpoint_dir)
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



    if (problem::problem == 1) {

        std::ifstream rotation;
        rotation.open(checkpoint_dir + "/Rotation", std::ios::in);

        if (rotation.is_open()) {
            rotation >> castro::rotational_period;
            amrex::Print() << "  Based on the checkpoint, updating the rotational period to "
                           << std::setprecision(7) << std::fixed << castro::rotational_period << " s.\n";
            rotation.close();
        }
        else {
            castro::rotational_period = -1.0;
        }

    }



    std::ifstream extrema;
    extrema.open(checkpoint_dir + "/Extrema", std::ios::in);

    if (extrema.is_open()) {
        extrema >> problem::T_global_max;
        extrema >> problem::rho_global_max;
        extrema >> problem::ts_te_global_max;
        extrema.close();
    }
    else {
        problem::T_global_max = 0.0;
        problem::rho_global_max = 0.0;
        problem::ts_te_global_max = 0.0;
    }



    std::ifstream energy;
    energy.open(checkpoint_dir + "/Energy", std::ios::in);

    if (energy.is_open()) {
        for (int i = 0; i < problem::num_previous_ener_timesteps; ++i) {
            energy >> problem::total_ener_array[i];
        }
        energy.close();
    }
    else {
        for (int i = 0; i < problem::num_previous_ener_timesteps; ++i) {
            problem::total_ener_array[i] = -1.e200_rt;
        }
    }


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
        set_relaxation_damping_factor(-1.0_rt);
    }



    // Pick up whether the job has been completed.

    std::ifstream jobisdone;
    jobisdone.open(checkpoint_dir + "/jobIsDone", std::ios::in);

    if (jobisdone.is_open()) {
        problem::jobIsDone = 1;
        jobisdone.close();
    }
}

#endif
