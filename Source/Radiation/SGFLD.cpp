#include "Radiation.H"

#include "RAD_F.H"

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

void Radiation::SGFLD_compute_rosseland(MultiFab& kappa_r, const MultiFab& state)
{
  BL_PROFILE("Radiation::SGFLD_compute_rosseland (MultiFab)");

  if (use_opacity_table_module) {
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(kappa_r,true); mfi.isValid(); ++mfi) {
	  const Box& kbox = mfi.growntilebox();
	  BL_FORT_PROC_CALL(CA_COMPUTE_ROSSELAND, ca_compute_rosseland)
	      (kbox.loVect(), kbox.hiVect(), 
	       BL_TO_FORTRAN(kappa_r[mfi]), BL_TO_FORTRAN(state[mfi]));
      }
  }
  else if (const_scattering[0] > 0.0) {
      Real k_exp_p = 0.0;
      Real s_exp_p = 0.0;
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(kappa_r,true); mfi.isValid(); ++mfi) {
	  const Box& bx = mfi.growntilebox();
	  BL_FORT_PROC_CALL(CA_COMPUTE_POWERLAW_KAPPA_S, ca_compute_powerlaw_kappa_s)
	      (bx.loVect(), bx.hiVect(), 
	       BL_TO_FORTRAN(kappa_r[mfi]), BL_TO_FORTRAN(state[mfi]),
	       &const_kappa_r[0], &kappa_r_exp_m[0], &kappa_r_exp_n[0], &k_exp_p,
	       &const_scattering[0], &scattering_exp_m[0], &scattering_exp_n[0], &s_exp_p,
	       &prop_temp_floor[0], &kappa_r_floor);	 
      }
  }
  else {
      Real k_exp_p = 0.0;
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(kappa_r,true); mfi.isValid(); ++mfi) {
	  const Box& bx = mfi.growntilebox();
	  BL_FORT_PROC_CALL(CA_COMPUTE_POWERLAW_KAPPA, ca_compute_powerlaw_kappa)
	      (bx.loVect(), bx.hiVect(), 
	       BL_TO_FORTRAN(kappa_r[mfi]), BL_TO_FORTRAN(state[mfi]),
	       &const_kappa_r[0], &kappa_r_exp_m[0], &kappa_r_exp_n[0], &k_exp_p,
	       &prop_temp_floor[0], &kappa_r_floor);	 
      }
  }
}

void Radiation::SGFLD_compute_rosseland(FArrayBox& kappa_r, const FArrayBox& state)
{
  BL_PROFILE("Radiation::SGFLD_compute_rosseland (FArrayBox)");

  const Box& kbox = kappa_r.box();

  if (use_opacity_table_module) {
    BL_FORT_PROC_CALL(CA_COMPUTE_ROSSELAND, ca_compute_rosseland)
	(kbox.loVect(), kbox.hiVect(), 
	 BL_TO_FORTRAN(kappa_r), BL_TO_FORTRAN(state));
  }
  else if (const_scattering[0] > 0.0) {
    Real k_exp_p = 0.0;
    Real s_exp_p = 0.0;
    BL_FORT_PROC_CALL(CA_COMPUTE_POWERLAW_KAPPA_S, ca_compute_powerlaw_kappa_s)
	(kbox.loVect(), kbox.hiVect(), 
	 BL_TO_FORTRAN(kappa_r), BL_TO_FORTRAN(state),
	 &const_kappa_r[0], &kappa_r_exp_m[0], &kappa_r_exp_n[0], &k_exp_p,
	 &const_scattering[0], &scattering_exp_m[0], &scattering_exp_n[0], &s_exp_p, 
	 &prop_temp_floor[0], &kappa_r_floor);	 
  }
  else {
    Real k_exp_p = 0.0;
    BL_FORT_PROC_CALL(CA_COMPUTE_POWERLAW_KAPPA, ca_compute_powerlaw_kappa)
	(kbox.loVect(), kbox.hiVect(), 
	 BL_TO_FORTRAN(kappa_r), BL_TO_FORTRAN(state),
	 &const_kappa_r[0], &kappa_r_exp_m[0], &kappa_r_exp_n[0], &k_exp_p,
	 &prop_temp_floor[0], &kappa_r_floor);	 
  }
}

