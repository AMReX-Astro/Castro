#include <iostream>
#include <iomanip>
#include <cstring>
#include <fstream>
#include <vector>

#include <AMReX_ParmParse.H>

#include <extern_parameters.H>
#include <prob_parameters.H>
#ifdef MICROPHYSICS_FORT
#include <microphysics_F.H>
#endif
#include <eos.H>
#include <network.H>
#include <castro_params.H>
#include <exact_riemann.H>

int main(int argc, char *argv[]) {

  std::cout << "starting the exact Riemann solver..." << std::endl;
  std::cout << argv[1] << std::endl;
  std::cout << strlen(argv[1]) << std::endl;

  // initialize the external runtime parameters in C++

  init_extern_parameters();

  // now initialize the C++ Microphysics
#ifdef REACTIONS
  network_init();
#endif

  eos_init();

  exact_riemann();

}
