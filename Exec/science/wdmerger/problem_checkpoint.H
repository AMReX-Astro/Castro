#ifndef problem_checkpoint_H
#define problem_checkpoint_H

#include <iostream>
#include <fstream>

#include <prob_parameters.H>

AMREX_INLINE
void problem_checkpoint (const std::string& checkpoint_dir)
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



    std::ofstream relaxation;
    relaxation.open(checkpoint_dir + "/Relaxation", std::ios::out);

    relaxation << std::setw(30) << problem::relaxation_is_done << std::endl;

    relaxation.close();

    std::ofstream radial_damping;
    radial_damping.open(checkpoint_dir + "/RadialDamping", std::ios::out);

    radial_damping << std::setw(30) << problem::radial_damping_is_done << std::endl;

    radial_damping.close();
}

#endif
