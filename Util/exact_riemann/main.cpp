#include <iostream>
#include <cstring>
#include <vector>

#include <AMReX_ParmParse.H>

#include <riemann_F.H>
#include <extern_parameters_F.H>
#include <extern_parameters.H>
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

  // initialize the external runtime parameters

  const int probin_file_length = strlen(argv[1]);
  std::vector<int> probin_file_name(probin_file_length);

  for (int i = 0; i < probin_file_length; i++)
    probin_file_name[i] = argv[1][i];

  // read them in in Fortran from the probin file

  runtime_init(probin_file_name.data(), &probin_file_length);

  // grab them from Fortran to C++; then read any C++ parameters directly
  // from inputs (via ParmParse)

  init_extern_parameters();

#ifdef MICROPHYSICS_FORT
  // finally, update the Fortran side via ParmParse to override the
  // values of any parameters that were set in inputs

  update_fortran_extern_after_cxx();
#endif


  // initialize microphysics
#ifdef MICROPHYSICS_FORT
#if !defined(NETWORK_HAS_CXX_IMPLEMENTATION)
  // Initialize the Fortran Microphysics
  microphysics_initialize(castro::small_temp, castro::small_dens);
#endif
#endif

  // now initialize the C++ Microphysics
#ifdef REACTIONS
  network_init();
#endif

  eos_init();

  riemann_exact();

}
