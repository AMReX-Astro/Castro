#ifndef PROB_UTIL_H
#define PROB_UTIL_H

#include "prob_parameters.H"
#include "extern_parameters.H"
#include "gravity_params.H"

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
gamma_law_initial_model(Real* const pressure, Real* const density, Real* const temp, Real* const eint,
                        const int npts, const GpuArray<Real,AMREX_SPACEDIM>& dx) {

  pressure[0] = problem::pres_base;
  density[0]  = problem::dens_base;

  // only initialize the first species
  Real xn[NumSpec];
  xn[0] = 1.0_rt;

  // compute the pressure scale height (for an isothermal, ideal-gas
  // atmosphere)
  Real H = problem::pres_base / problem::dens_base / std::abs(gravity::const_grav);

  for (int j = 0; j < npts; j++) {

    // initial guess
     temp[j] = T_guess;

     if (problem::do_isentropic) {
       Real z = static_cast<Real>(j) * dx[AMREX_SPACEDIM-1];
       density[j] = problem::dens_base *
         std::pow((gravity::const_grav * problem::dens_base * (eos_rp::eos_gamma - 1.0_rt) * z/
                   (eos_rp::eos_gamma * problem::pres_base) + 1.0_rt), 1.0_rt/(eos_rp::eos_gamma - 1.0_rt));
     } else {
       Real z = (static_cast<Real>(j) + 0.5_rt) * dx[AMREX_SPACEDIM-1];
       density[j] = problem::dens_base * std::exp(-z/H);
     }

     if (j > 0) {
        pressure[j] = pressure[j-1] -
          dx[AMREX_SPACEDIM-1] * 0.5_rt *
          (density[j] + density[j-1]) * std::abs(gravity::const_grav);
     }

     eos_t eos_state;

     eos_state.rho = density[j];
     eos_state.T = temp[j];
     for (int n = 0; n < NumSpec; n++) {
       eos_state.xn[n] = xn[n];
     }
     eos_state.p = pressure[j];

     eos(eos_input_rp, eos_state);

     eint[j] = eos_state.e;
     temp[j] = eos_state.T;

  }

}

#endif
