#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>
#include <network.H>
#ifdef NSE_THERMO
#include <nse.H>
#endif

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

    int ifuel5 = network_spec_index(problem::fuel5_name);
    int iash5 = network_spec_index(problem::ash5_name);

    int ifuel6 = network_spec_index(problem::fuel6_name);
    int iash6 = network_spec_index(problem::ash6_name);

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

    if (ifuel5 >= 0) {
        problem::xn_fuel[ifuel5] = problem::X_fuel5;
    }

    if (ifuel6 >= 0) {
        problem::xn_fuel[ifuel6] = problem::X_fuel6;
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

    if (iash5 >= 0) {
        problem::xn_ash[iash5] = problem::X_ash5;
    }

    if (iash6 >= 0) {
        problem::xn_ash[iash6] = problem::X_ash6;
    }

    if (iash1 < 0 && iash2 < 0 && iash3 < 0 && iash4 < 0 && iash5 < 0 && iash6 < 0) {
        amrex::Error("no valid ash state defined");
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

    // get the species indices

    problem::ihe4 = network_spec_index("helium-4");
    problem::ic12 = network_spec_index("carbon-12");
    problem::io16 = network_spec_index("oxygen-16");

    if (problem::ihe4 < 0 || problem::ic12 < 0 || problem::io16 < 0) {
        amrex::Error("ERROR: species indices not found");
    }

    // make sure that the carbon fraction falls between 0 and 1

    if (problem::cfrac > 1.0_rt || problem::cfrac < 0.0_rt) {
        amrex::Error("ERROR: cfrac must fall between 0 and 1");
    }

    // make sure that the oxygen fraction falls between 0 and 1

    if (problem::ofrac > 1.0_rt || problem::ofrac < 0.0_rt) {
        amrex::Error("ERROR: ofrac must fall between 0 and 1");
    }

    // make sure that the C/O fraction sums to no more than 1

    if (problem::cfrac + problem::ofrac > 1.0_rt) {
        amrex::Error("ERROR: cfrac + ofrac cannot exceed 1.");
    }

    // set the default mass fractions

    for (int n = 0; n < NumSpec; n++) {
        problem::xn[n] = problem::smallx;
    }
    problem::xn[problem::ic12] = amrex::max(problem::cfrac, problem::smallx);
    problem::xn[problem::io16] = amrex::max(problem::ofrac, problem::smallx);
    problem::xn[problem::ihe4] = 1.0_rt - problem::cfrac - problem::ofrac - (NumSpec - 2) * problem::smallx;

    // Set the ambient material

    problem::ambient_dens = problem::dens_ash;
    for (int n = 0; n < NumSpec; n++) {
        problem::ambient_comp[n] = problem::xn[n];
    }

    eos_t eos_state;


    eos_state.rho = problem::ambient_dens;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n]  = problem::ambient_comp[n];
    }

    eos_state.T = problem::T_ash;

    eos(eos_input_rt, eos_state);

    problem::ambient_e_ash = eos_state.e;

    eos_state.T = problem::T_fuel;

    eos(eos_input_rt, eos_state);

    problem::ambient_e_fuel = eos_state.e;

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
    eos_state.p = problem::p_fuel;      // pressure equilibrum
    eos_state.T = problem::T_ash;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = problem::xn_ash[n];
    }

    eos(eos_input_tp, eos_state);

    problem::rho_ash = eos_state.rho;
    problem::e_ash = eos_state.e;

    // mass flux will be constant across the flame

    problem::mass_flux = problem::rho_fuel * problem::v_inflow;

}

#endif
