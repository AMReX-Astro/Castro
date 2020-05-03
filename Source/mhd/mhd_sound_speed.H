#ifndef _MHD_SOUND_SPEED_H_
#define _MHD_SOUND_SPEED_H_

#include "Castro_util.H"

using namespace amrex;


// sound speed calc for ideal MHD
void
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
eos_soundspeed_mhd(Real& c, Real as, Real ca, Real bd) {

  // inputs:
  // as = P_g * gam1/rho
  // ca = B^2/rho
  // bd = B_direction^2/rho

  // Fast Magneto-Sonic Wave
  c = 0.5_rt * ((as + ca) + std::sqrt((as + ca)*(as + ca) - 4.0_rt*as*bd));
  c = std::sqrt(c);
}


#endif
