#ifdef SPONGE

#include "Castro.H"
#include "Castro_F.H"

void Castro::update_sponge_params(const amrex::Real time) {

  // This function is designed to permit you to update
  // the sponge parameters as a time-dependent process.
  // By default it does nothing, meaning that the probin
  // parameters are used.

}

#endif
