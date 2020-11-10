#ifndef problem_checkpoint_H
#define problem_checkpoint_H

#include <iostream>
#include <fstream>

#include <prob_parameters.H>

AMREX_INLINE
void problem_checkpoint (std::string checkpoint_dir)
{
    std::ofstream com;
    com.open(checkpoint_dir + "/COM", std::ios::out);

    com << std::scientific;
    com.precision(19);

    com << std::setw(30) << problem::center[0] << " "
        << std::setw(30) << problem::center[1] << " "
        << std::setw(30) << problem::center[2] << std::endl;
    com << std::setw(30) << problem::mass_P << " "
        << std::setw(30) << problem::mass_S << std::endl;
    com << std::setw(30) << problem::com_P[0] << " "
        << std::setw(30) << problem::com_S[0] << std::endl;
    com << std::setw(30) << problem::com_P[1] << " "
        << std::setw(30) << problem::com_S[1] << std::endl;
    com << std::setw(30) << problem::com_P[2] << " "
        << std::setw(30) << problem::com_S[2] << std::endl;
    com << std::setw(30) << problem::vel_P[0] << " "
        << std::setw(30) << problem::vel_S[0] << std::endl;
    com << std::setw(30) << problem::vel_P[1] << " "
        << std::setw(30) << problem::vel_S[1] << std::endl;
    com << std::setw(30) << problem::vel_P[2] << " "
        << std::setw(30) << problem::vel_S[2] << std::endl;
    com << std::setw(30) << problem::t_ff_P << " "
        << std::setw(30) << problem::t_ff_S << std::endl;

    com.close();



    std::ofstream rotation;
    rotation.open(checkpoint_dir + "/Rotation", std::ios::out);

    rotation << std::scientific;
    rotation.precision(19);

    rotation << std::setw(30) << castro::rotational_period << std::endl;

    rotation.close();



    std::ofstream extrema;
    extrema.open(checkpoint_dir + "/Extrema", std::ios::out);

    extrema << std::scientific;
    extrema.precision(19);

    extrema << std::setw(30) << problem::T_global_max << std::endl;
    extrema << std::setw(30) << problem::rho_global_max << std::endl;
    extrema << std::setw(30) << problem::ts_te_global_max << std::endl;

    extrema.close();



    std::ofstream energy;
    energy.open(checkpoint_dir + "/Energy", std::ios::out);

    energy << std::scientific;
    energy.precision(19);

    for (int i = 0; i < problem::num_previous_ener_timesteps; ++i) {
        energy << std::setw(30) << problem::total_ener_array[i] << std::endl;
    }

    energy.close();



    std::ofstream relaxation;
    relaxation.open(checkpoint_dir + "/Relaxation", std::ios::out);

    relaxation << std::setw(30) << problem::relaxation_is_done << std::endl;

    relaxation.close();



    // If the job is done, write a file in the checkpoint indicating this.

    if (problem::jobIsDone) {

        std::ofstream jobisdone;
        jobisdone.open(checkpoint_dir + "/jobIsDone", std::ios::out);

        jobisdone.close();

    }
}

#endif
