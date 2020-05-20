#ifndef CASTRO_MHD_SOUND_SPEED_H
#define CASTRO_MHD_SOUND_SPEED_H

#include "Castro_util.H"

using namespace amrex;

constexpr int UMAGX = NUM_STATE;
constexpr int UMAGY = NUM_STATE+1;
constexpr int UMAGZ = NUM_STATE+2;

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

void
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Castro::qflux(Real* qflx, Real* flx, Real* q_zone) {

  // Calculate the C to P Jacobian applied to the fluxes

  // this is step 10 in the paper, just after Eq. 48

  // this implements dW/dU . qflux, where dW/dU is the Jacobian of
  // the primitive quantities (W) with respect to conserved quantities (U)

  qflx[QRHO] = flx[URHO];
  qflx[QU] = (flx[UMX] - flx[URHO] * q_zone[QU]) / q_zone[QRHO];
  qflx[QV] = (flx[UMY] - flx[URHO] * q_zone[QV]) / q_zone[QRHO];
  qflx[QW] = (flx[UMZ] - flx[URHO] * q_zone[QW]) / q_zone[QRHO];

  for (int n = 0; n < NumSpec; n++) {
    qflx[QFS+n] = (flx[UFS+n] - flx[URHO] * q_zone[QFS+n]) / q_zone[QRHO];
  }

  eos_t eos_state;
  eos_state.rho = q_zone[QRHO];
  eos_state.p = q_zone[QPRES];
  eos_state.T = 100.0_rt; // dummy initial guess
  for (int n = 0; n < NumSpec; n++) {
    eos_state.xn[n] = q_zone[QFS+n];
  }

  eos(eos_input_rp, eos_state);

  Real dedrho = eos_state.dedr - eos_state.dedT * eos_state.dpdr * 1.0_rt / eos_state.dpdT;
  Real dedp = eos_state.dedT * 1.0_rt / eos_state.dpdT;

  qflx[QPRES] = (-q_zone[QMAGX] * flx[UMAGX] - q_zone[QMAGY] * flx[UMAGY] -
                  q_zone[QMAGZ] * flx[UMAGZ] + flx[UEDEN] -
                 flx[UMX] * q_zone[QU] - flx[UMY] * q_zone[QV] - flx[UMZ] * q_zone[QW] +
                 flx[URHO] * (0.5_rt * (q_zone[QU] * q_zone[QU] +
                                        q_zone[QV] * q_zone[QV] +
                                        q_zone[QW] * q_zone[QW]) -
                              eos_state.e - q_zone[QRHO] * dedrho)) /
    (dedp * q_zone[QRHO]);

  qflx[QMAGX] = flx[UMAGX];
  qflx[QMAGY] = flx[UMAGY];
  qflx[QMAGZ] = flx[UMAGZ];

  qflx[QTEMP] = 0.0_rt;
}
#endif

