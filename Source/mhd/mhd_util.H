#ifndef CASTRO_MHD_SOUND_SPEED_H
#define CASTRO_MHD_SOUND_SPEED_H

#include "Castro_util.H"
#include "mhd_sizes.H"

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


void
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Castro::electric(Real* q_zone, Real& E_zone, const int comp) {

  // this takes the cell-center primitive state, q_zone, and computes the cell-center
  // electric field, E_zone, using Faraday's law:
  // E = -v X B

  if (comp == 0) {
    E_zone = -q_zone[QV] * q_zone[QMAGZ] + q_zone[QW] * q_zone[QMAGY];

  } else if (comp == 1) {
    E_zone = -q_zone[QW] * q_zone[QMAGX] + q_zone[QU] * q_zone[QMAGZ];

  } else {
    E_zone = -q_zone[QU] * q_zone[QMAGY] + q_zone[QV] * q_zone[QMAGX];
  }

}


void
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Castro::ConsToPrim(Real* q_zone, Real* U_zone) {

  // calculate the primitive variables from the conserved
  // U has NUM_STATE+3 components
  // q has NQ components


  q_zone[QRHO] = U_zone[URHO];
  q_zone[QU] = U_zone[UMX] / q_zone[QRHO];
  q_zone[QV] = U_zone[UMY] / q_zone[QRHO];
  q_zone[QW] = U_zone[UMZ] / q_zone[QRHO];
  q_zone[QREINT] = U_zone[UEDEN] -
    0.5_rt * q_zone[QRHO] * (q_zone[QU] * q_zone[QU] +
                             q_zone[QV] * q_zone[QV] +
                             q_zone[QW] * q_zone[QW]) -
    0.5_rt * (U_zone[UMAGX] * U_zone[UMAGX] +
              U_zone[UMAGY] * U_zone[UMAGY] +
              U_zone[UMAGZ] * U_zone[UMAGZ]);

  // species

  for (int n = 0; n < NumSpec; n++) {
    q_zone[QFS+n] = U_zone[UFS+n] / U_zone[URHO];
  }

  eos_t eos_state;

  eos_state.rho = q_zone[QRHO];
  eos_state.e = q_zone[QREINT] / eos_state.rho;
  for (int n = 0; n < NumSpec; n++) {
    eos_state.xn[n] = q_zone[QFS+n];
  }
  eos_state.T = 100.0_rt; // initial guess

  eos(eos_input_re, eos_state);

  q_zone[QTEMP] = eos_state.T;
  q_zone[QPRES] = eos_state.p;

  q_zone[QMAGX] = U_zone[UMAGX];
  q_zone[QMAGY] = U_zone[UMAGY];
  q_zone[QMAGZ] = U_zone[UMAGZ];

}



int
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
epsilon_ijk(const int i, const int j, const int k) {

  int val;

  if (i == j || j == k || i == k) {
    val = 0;
  } else if ((i == 0 && j == 1 && k == 2) ||
             (i == 1 && j == 2 && k == 0) ||
             (i == 2 && j == 0 && k == 1)) {
    val = 1;
  } else {
    val = -1;
  }

  return val;
}

#endif

