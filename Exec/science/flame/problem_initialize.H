#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>
#include <network.H>

AMREX_INLINE
void problem_initialize ()
{

    int ifuel1 = network_spec_index(problem::fuel1_name);
    int iash1 = network_spec_index(problem::ash1_name);

    int ifuel2 = network_spec_index(problem::fuel2_name);
    int iash2 = network_spec_index(problem::ash2_name);

    int ifuel3 = network_spec_index(problem::fuel3_name);
    int iash3 = network_spec_index(problem::ash3_name);

    int ifuel4 = network_spec_index(problem::fuel4_name);
    int iash4 = network_spec_index(problem::ash4_name);

    if (iash1 < 0 && iash2 < 0 && iash3 < 0 && iash4 < 0) {
        amrex::Error("no valid ash state defined");
    }

    // fuel state

    for (int n = 0; n < NumSpec; n++) {
        problem::xn_fuel[n] = problem::smallx_init;
        problem::xn_ash[n] = problem::smallx_init;
    }

    if (ifuel1 >= 0) {
        problem::xn_fuel[ifuel1] = problem::X_fuel1;
    }

    if (ifuel2 >= 0) {
        problem::xn_fuel[ifuel2] = problem::X_fuel2;
    }

    if (ifuel3 >= 0) {
        problem::xn_fuel[ifuel3] = problem::X_fuel3;
    }

    if (ifuel4 >= 0) {
        problem::xn_fuel[ifuel4] = problem::X_fuel4;
    }

    if (iash1 >= 0) {
        problem::xn_ash[iash1] = problem::X_ash1;
    }

    if (iash2 >= 0) {
        problem::xn_ash[iash2] = problem::X_ash2;
    }

    if (iash3 >= 0) {
        problem::xn_ash[iash3] = problem::X_ash3;
    }

    if (iash4 >= 0) {
        problem::xn_ash[iash4] = problem::X_ash4;
    }

    // normalize

    Real sum_Xf = 0.0_rt;
    Real sum_Xa = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
        sum_Xf += problem::xn_fuel[n];
        sum_Xa += problem::xn_ash[n];
    }
    for (int n = 0; n < NumSpec; n++) {
        problem::xn_fuel[n] /= sum_Xf;
        problem::xn_ash[n] /= sum_Xa;
    }

    eos_t eos_state;
    eos_state.rho = problem::rho_fuel;
    eos_state.T = problem::T_fuel;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = problem::xn_fuel[n];
    }

    eos(eos_input_rt, eos_state);

    problem::e_fuel = eos_state.e;
    problem::p_fuel = eos_state.p;

    // compute the ash state -- this should be hot but in pressure
    // equilibrium

    eos_state.rho = problem::rho_fuel;  // initial guess
    eos_state.p = problem::p_fuel;      // pressure equilibrium
    eos_state.T = problem::T_ash;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = problem::xn_ash[n];
    }

    eos(eos_input_tp, eos_state);

    problem::rho_ash = eos_state.rho;
    problem::e_ash = eos_state.e;

    // mass flux will be constant across the flame

    problem::mass_flux = problem::rho_fuel * problem::v_inflow;

    // Set up Castro data logs for this problem

    if (castro::sum_interval > 0 && amrex::ParallelDescriptor::IOProcessor()) {

        Castro::problem_data_logs.resize(1);

        Castro::problem_data_logs[0] = std::make_unique<std::fstream>();
        Castro::problem_data_logs[0]->open("toy_flame.log", std::ios::out | std::ios::app);
        if (!Castro::problem_data_logs[0]->good()) {
            amrex::FileOpenFailed("toy_flame.log");
        }

    }

}
#endif
