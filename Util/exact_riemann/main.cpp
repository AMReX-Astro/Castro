#include <iostream>
#include <iomanip>
#include <cstring>
#include <fstream>
#include <vector>

#include <AMReX_ParmParse.H>

#include <extern_parameters.H>
#include <prob_parameters.H>
#include <eos.H>
#include <network.H>

#include <exact_riemann.H>

#include <castro_params.H>
#include <params_type.H>

// a global struct to hold the params
#include <struct_params.H>

int main(int argc, char *argv[]) {

    amrex::Initialize(argc, argv);

    std::cout << "starting the exact Riemann solver..." << std::endl;

    // initialize the Castro runtime parameters

    amrex::ParmParse pp("castro");
#include <castro_queries.H>

    // initialize the external runtime parameters in C++

    init_prob_parameters();

    init_extern_parameters();

    // now initialize the C++ Microphysics
#ifdef REACTIONS
    network_init();
#endif

    amrex::Real small_temp = 1.e-200;
    amrex::Real small_dens = 1.e-200;
    eos_init(small_temp, small_dens);

    exact_riemann();

}
