#include <cmath>

#ifndef _rad_util_H_
#define _rad_util_H_

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
compute_ptot_ctot(const Real* lam, const Real* q,
                  int is_comoving, int limiter, int closure,
                  Real cg,
                  Real& ptot, Real& ctot, Real& gamc_tot) {

  Real csrad2 = 0.0_rt;
  Real prad = 0.0_rt;

  for (int g = 0; g < NGROUPS; g++) {
    Real gamr = 0.0;
    if (is_comoving) {
      Real Eddf = Edd_factor(lam[g], limiter, closure);
      gamr = (3.0_rt - Eddf)/2.0_rt;
    } else {
      gamr = lam[g] + 1.0_rt;
    }

    prad += lam[g]*q[QRAD+g];
    csrad2 += gamr * (lam[g]*q[QRAD+g]) / q[QRHO];
  }

  ptot = q[QPRES] + prad;

  ctot = cg*cg + csrad2;
  gamc_tot = ctot * q[QRHO] / ptot;

  ctot = std::sqrt(ctot);
}

