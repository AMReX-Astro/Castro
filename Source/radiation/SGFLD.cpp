#include "Radiation.H"

#include "RAD_F.H"

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace amrex;

void Radiation::SGFLD_compute_rosseland(MultiFab& kappa_r, const MultiFab& state)
{
  BL_PROFILE("Radiation::SGFLD_compute_rosseland (MultiFab)");

  if (use_opacity_table_module) {
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(kappa_r,true); mfi.isValid(); ++mfi) {
	  const Box& kbox = mfi.growntilebox();
#pragma gpu box(kbox)
	  ca_compute_rosseland(AMREX_INT_ANYD(kbox.loVect()), AMREX_INT_ANYD(kbox.hiVect()),
			       BL_TO_FORTRAN_ANYD(kappa_r[mfi]),
                               BL_TO_FORTRAN_ANYD(state[mfi]));
      }
  }
  else if (const_scattering > 0.0) {
      Real k_exp_p = 0.0;
      Real s_exp_p = 0.0;
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(kappa_r,true); mfi.isValid(); ++mfi) {
	  const Box& bx = mfi.growntilebox();
	  ca_compute_powerlaw_kappa_s(bx.loVect(), bx.hiVect(), 
				      BL_TO_FORTRAN(kappa_r[mfi]), BL_TO_FORTRAN(state[mfi]),
				      &const_kappa_r, &kappa_r_exp_m, &kappa_r_exp_n, &k_exp_p,
				      &const_scattering, &scattering_exp_m, &scattering_exp_n, &s_exp_p,
				      &prop_temp_floor, &kappa_r_floor);	 
      }
  }
  else {
      Real k_exp_p = 0.0;
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(kappa_r,true); mfi.isValid(); ++mfi) {
	  const Box& bx = mfi.growntilebox();
	  ca_compute_powerlaw_kappa(bx.loVect(), bx.hiVect(), 
				    BL_TO_FORTRAN(kappa_r[mfi]), BL_TO_FORTRAN(state[mfi]),
				    &const_kappa_r, &kappa_r_exp_m, &kappa_r_exp_n, &k_exp_p,
				    &prop_temp_floor, &kappa_r_floor);	 
      }
  }
}

void Radiation::SGFLD_compute_rosseland(FArrayBox& kappa_r, const FArrayBox& state)
{
  BL_PROFILE("Radiation::SGFLD_compute_rosseland (FArrayBox)");

  const Box& kbox = kappa_r.box();

  if (use_opacity_table_module) {
#pragma gpu box(kbox) sync
      ca_compute_rosseland(AMREX_INT_ANYD(kbox.loVect()), AMREX_INT_ANYD(kbox.hiVect()),
                           BL_TO_FORTRAN_ANYD(kappa_r),
                           BL_TO_FORTRAN_ANYD(state));
  }
  else if (const_scattering > 0.0) {
    Real k_exp_p = 0.0;
    Real s_exp_p = 0.0;
    ca_compute_powerlaw_kappa_s(kbox.loVect(), kbox.hiVect(), 
				BL_TO_FORTRAN(kappa_r), BL_TO_FORTRAN(state),
				&const_kappa_r, &kappa_r_exp_m, &kappa_r_exp_n, &k_exp_p,
				&const_scattering, &scattering_exp_m, &scattering_exp_n, &s_exp_p, 
				&prop_temp_floor, &kappa_r_floor);	 
  }
  else {
    Real k_exp_p = 0.0;
    ca_compute_powerlaw_kappa(kbox.loVect(), kbox.hiVect(), 
			      BL_TO_FORTRAN(kappa_r), BL_TO_FORTRAN(state),
			      &const_kappa_r, &kappa_r_exp_m, &kappa_r_exp_n, &k_exp_p,
			      &prop_temp_floor, &kappa_r_floor);
  }
}

