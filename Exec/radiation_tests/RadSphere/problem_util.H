#ifndef PROBLEM_UTIL_H
#define PROBLEM_UTIL_H

#include "fundamental_constants.H"

// problem parameters
//
// we will specify the opacity as kappa = kapp_0 (nu_0 / nu)**3
// where nu_0 is a reference frequency
//
// Swesty and Myra (2009) only say that the "opacity of the material
// is proportional to 1/E**3", but they don't give the proportionality
// constant.  The output depends critically on the constant.  Eric
// says that they used kappa = 1.e13 * (1.5e-6 MeV / E)**3
//
// Converting 1.5e-6 MeV to a frequency gives nu_0 = 3.63e14 Hz

constexpr Real nu_0 = 3.6e14_rt; //  ref. freq (Hz)
constexpr Real kappa_0 = 1.e13_rt; // scattering opacity (1/cm)

// geometry parameters
constexpr Real R_sphere = 0.02_rt; // sphere radius (cm)
constexpr Real r_obs = 0.06_rt; // observer location (cm)

constexpr Real T_sphere = 1500.0_rt * C::ev2erg / C::k_B; // sphere temp (K)

// to guard against over/underflows, we define some 'safe' functions
// (they come from Eric)

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real
sferfc(const Real arg1) {
  return std::erfc(amrex::max(-26.0_rt, amrex::min(26.0_rt, arg1)));
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real
sfexp(const Real arg1) {
  return std::exp(amrex::max(-650.0_rt, amrex::min(650.0_rt, arg1)));
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real
planck(const Real nu, const Real T) {

  // the Planck function for a Blackbody (actually, energy density,
  // B = (4 pi / c) I, where I is the normal Planck function
  //
  // nu = frequency (Hz)
  // T  = temperature (K)
  //
  // see Swesty and Myra (2009), eq. 23, but note that we are working
  // in Hz here, so we have units of erg / cm^3 / Hz, where they
  // have units of erg / cm^3 / MeV.  As a result, we have one
  // less factor of hplanck.

  Real B = (8.0_rt * M_PI * C::hplanck * nu * nu * nu  / (C::c_light * C::c_light * C::c_light)) /
    (std::exp(C::hplanck * nu / (C::k_B * T)) - 1.0_rt);
  B = amrex::min(amrex::max(B, 1.e-50_rt), 1.e200_rt);

  return B;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real
F_radsphere(const Real r, const Real t, const Real nu) {

  // the function F(r,t) as defined in Swesty and Myra
  //
  // r      = position
  // t      = time (s)
  // nu     = frequency (Hz)


  Real kappa = kappa_0 * (nu_0 / nu) * (nu_0 / nu) * (nu_0 / nu);

  Real erfc_term1 = sferfc(std::sqrt(3.0_rt * kappa / (4.0_rt * C::c_light * amrex::max(t, 1.e-50_rt))) *
                           (r - R_sphere) - std::sqrt(C::c_light * t * kappa));

  Real erfc_term2 = sferfc(std::sqrt(3.0_rt * kappa / (4.0_rt * C::c_light * amrex::max(t, 1.e-50_rt))) *
                           (r - R_sphere) + std::sqrt(C::c_light * t * kappa));

  Real F = 0.5_rt * (sfexp(-std::sqrt(3.0_rt) * kappa * (r - R_sphere)) * erfc_term1 +
                     sfexp( std::sqrt(3.0_rt) * kappa * (r - R_sphere)) * erfc_term2);

  return amrex::min(amrex::max(F, 1.e-50_rt), 1.e200_rt);
}

#endif
