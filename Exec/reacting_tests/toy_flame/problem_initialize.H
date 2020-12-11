#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>
#include <network.H>

AMREX_INLINE
void problem_initialize ()
{

    problem::ifuel = network_spec_index("fuel");
    problem::iash = network_spec_index("ash");

    // output flame speed and width estimates

    eos_t eos_state;
    eos_state.rho = rho_burn_ref;
    eos_state.T = T_burn_ref;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = 0.0_rt;
    }
    eos_state.xn[0] = 1.0_rt;

    eos(eos_input_rt, eos_state);

    Real lambda_f = std::sqrt(const_conductivity * T_burn_ref/
                              (rho_burn_ref * specific_q_burn * nu * rtilde));

    Real v_f = std::sqrt(const_conductivity * specific_q_burn * nu * rtilde/
                         (rho_burn_ref * eos_state.cp * eos_state.cp * T_burn_ref));

    std::cout << "flame width = " << lambda_f << std::endl;
    std::cout << "flame speed = " << v_f << std::endl;

    // Set up Castro data logs for this problem

    if (castro::sum_interval > 0 && amrex::ParallelDescriptor::IOProcessor()) {

        Castro::problem_data_logs.resize(1);

        Castro::problem_data_logs[0].reset(new std::fstream);
        Castro::problem_data_logs[0]->open("toy_flame.log", std::ios::out | std::ios::app);
        if (!Castro::problem_data_logs[0]->good()) {
            amrex::FileOpenFailed("toy_flame.log");
        }

    }

}

#endif
