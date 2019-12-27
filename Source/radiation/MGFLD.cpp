#include "Radiation.H"
#include "Castro_F.H"
#include "RAD_F.H"

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace amrex;

void Radiation::check_convergence_er(Real& relative, Real& absolute, Real& err_er,
				     const MultiFab& Er_new, const MultiFab& Er_pi, 
				     const MultiFab& kappa_p,
				     const MultiFab& etaTz, const MultiFab& etaYz, 
				     const MultiFab& theTz, const MultiFab& theYz,
				     const MultiFab& temp_new, const MultiFab& Ye_new,
				     const BoxArray& grids, Real delta_t)
{
  BL_PROFILE("Radiation::check_convergence_er (MGFLD)");

  relative = 0.0;
  absolute = 0.0;
  err_er = 0.0;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
      Real priv_relative = 0.0;
      Real priv_absolute = 0.0;
      Real priv_err_er   = 0.0;

      for (MFIter mfi(Er_new,true); mfi.isValid(); ++mfi) {
	  const Box& bx = mfi.tilebox();
#ifdef NEUTRINO
	  BL_FORT_PROC_CALL(CA_CHECK_CONV_ER_NEUT, ca_check_conv_er_neut)
	      (bx.loVect(), bx.hiVect(),
	       BL_TO_FORTRAN(Er_new[mfi]),
	       BL_TO_FORTRAN(Er_pi[mfi]),
	       BL_TO_FORTRAN(kappa_p[mfi]),
	       BL_TO_FORTRAN(etaTz[mfi]),
	       BL_TO_FORTRAN(etaYz[mfi]),
	       BL_TO_FORTRAN(theTz[mfi]),
	       BL_TO_FORTRAN(theYz[mfi]),
	       BL_TO_FORTRAN(temp_new[mfi]),
	       BL_TO_FORTRAN(Ye_new[mfi]),
	       &priv_relative, &priv_absolute, &priv_err_er, &delta_t);
#else
	  BL_FORT_PROC_CALL(CA_CHECK_CONV_ER, ca_check_conv_er)
	      (bx.loVect(), bx.hiVect(),
	       BL_TO_FORTRAN(Er_new[mfi]),
	       BL_TO_FORTRAN(Er_pi[mfi]),
	       BL_TO_FORTRAN(kappa_p[mfi]),
	       BL_TO_FORTRAN(etaTz[mfi]),
	       BL_TO_FORTRAN(temp_new[mfi]),
	       &priv_relative, &priv_absolute, &priv_err_er, &delta_t);
#endif
      }

#ifdef _OPENMP
#pragma omp critical(rad_check_conv_er)
#endif
      {
	  relative = std::max(relative, priv_relative);
	  absolute = std::max(absolute, priv_absolute);
	  err_er   = std::max(err_er  , priv_err_er);
      }
  }

  Real data[3] = {relative, absolute, err_er};

  ParallelDescriptor::ReduceRealMax(data, 3);

  relative = data[0];
  absolute = data[1];
  err_er   = data[2];
}


void Radiation::check_convergence_matt(const MultiFab& rhoe_new, const MultiFab& rhoe_star, 
				       const MultiFab& rhoe_step, const MultiFab& Er_new, 
				       const MultiFab& temp_new, const MultiFab& temp_star, 
				       const MultiFab& rhoYe_new, const MultiFab& rhoYe_star, 
				       const MultiFab& rhoYe_step, const MultiFab& rho,
				       const MultiFab& kappa_p, const MultiFab& jg,
				       const MultiFab& dedT, const MultiFab& dedY,
				       Real& rel_rhoe, Real& abs_rhoe, 
				       Real& rel_FT,   Real& abs_FT, 
				       Real& rel_T,    Real& abs_T, 
				       Real& rel_FY,   Real& abs_FY, 
				       Real& rel_Ye,   Real& abs_Ye,
				       const BoxArray& grids, Real delta_t)
{
  rel_rhoe = abs_rhoe = 0.0;
  rel_FT   = abs_FT   = 0.0;
  rel_T    = abs_T    = 0.0;
  rel_FY   = abs_FY   = 0.0;
  rel_Ye   = abs_Ye   = 0.0;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
      Real priv_rel_rhoe = 0.0;
      Real priv_abs_rhoe = 0.0;
      Real priv_rel_FT   = 0.0;
      Real priv_abs_FT   = 0.0;
      Real priv_rel_T    = 0.0;
      Real priv_abs_T    = 0.0;
#ifdef NEUTRINO
      Real priv_rel_FY   = 0.0;
      Real priv_abs_FY   = 0.0;
      Real priv_rel_Ye   = 0.0;
      Real priv_abs_Ye   = 0.0;
#endif      

      for (MFIter mfi(rhoe_new,true); mfi.isValid(); ++mfi) {
	  const Box& bx = mfi.tilebox();
#ifdef NEUTRINO    
	  BL_FORT_PROC_CALL(CA_CHECK_CONV_NEUT, ca_check_conv_neut)
	      (bx.loVect(), bx.hiVect(),
	       BL_TO_FORTRAN(rhoe_new[mfi]),
	       BL_TO_FORTRAN(rhoe_star[mfi]),
	       BL_TO_FORTRAN(rhoe_step[mfi]),
	       BL_TO_FORTRAN(Er_new[mfi]),
	       BL_TO_FORTRAN(temp_new[mfi]),
	       BL_TO_FORTRAN(temp_star[mfi]),
	       BL_TO_FORTRAN(rhoYe_new[mfi]),
	       BL_TO_FORTRAN(rhoYe_star[mfi]),
	       BL_TO_FORTRAN(rhoYe_step[mfi]),
	       BL_TO_FORTRAN(rho[mfi]),
	       BL_TO_FORTRAN(kappa_p[mfi]),
	       BL_TO_FORTRAN(jg[mfi]),
	       BL_TO_FORTRAN(dedT[mfi]),
	       BL_TO_FORTRAN(dedY[mfi]),
	       &priv_rel_rhoe, &priv_abs_rhoe, 
	       &priv_rel_FT,   &priv_abs_FT, 
	       &priv_rel_T,    &priv_abs_T, 
	       &priv_rel_FY,   &priv_abs_FY, 
	       &priv_rel_Ye,   &priv_abs_Ye,
	       &delta_t);
#else
	  BL_FORT_PROC_CALL(CA_CHECK_CONV, ca_check_conv)
	      (bx.loVect(), bx.hiVect(),
	       BL_TO_FORTRAN(rhoe_new[mfi]),
	       BL_TO_FORTRAN(rhoe_star[mfi]),
	       BL_TO_FORTRAN(rhoe_step[mfi]),
	       BL_TO_FORTRAN(Er_new[mfi]),
	       BL_TO_FORTRAN(temp_new[mfi]),
	       BL_TO_FORTRAN(temp_star[mfi]),
	       BL_TO_FORTRAN(rho[mfi]),
	       BL_TO_FORTRAN(kappa_p[mfi]),
	       BL_TO_FORTRAN(jg[mfi]),
	       BL_TO_FORTRAN(dedT[mfi]),
	       &priv_rel_rhoe, &priv_abs_rhoe, 
	       &priv_rel_FT,   &priv_abs_FT, 
	       &priv_rel_T,    &priv_abs_T, 
	       &delta_t);
#endif
      }

#ifdef _OPENMP
#pragma omp critical(rad_check_conv)
#endif
      {
	  rel_rhoe = std::max(rel_rhoe, priv_rel_rhoe);  
	  abs_rhoe = std::max(abs_rhoe, priv_abs_rhoe); 
	  rel_FT   = std::max(rel_FT  , priv_rel_FT  ); 
	  abs_FT   = std::max(abs_FT  , priv_abs_FT  ); 
	  rel_T    = std::max(rel_T   , priv_rel_T   ); 
	  abs_T    = std::max(abs_T   , priv_abs_T   ); 
#ifdef NEUTRINO
	  rel_FY   = std::max(rel_FY  , priv_rel_FY  ); 
	  abs_FY   = std::max(abs_FY  , priv_abs_FY  ); 
	  rel_Ye   = std::max(rel_Ye  , priv_rel_Ye  ); 
	  abs_Ye   = std::max(abs_Ye  , priv_abs_Ye  ); 
#endif      
      }
  }

#ifdef NEUTRINO
  int ndata = 10;
Real data[10] = {rel_rhoe, abs_rhoe, rel_FT, abs_FT, rel_T, abs_T, rel_FY, abs_FY, rel_Ye, abs_Ye};
#else
  int ndata = 6;
  Real data[6 ] = {rel_rhoe, abs_rhoe, rel_FT, abs_FT, rel_T, abs_T};
#endif

  ParallelDescriptor::ReduceRealMax(data,ndata);

  rel_rhoe = data[0]; 
  abs_rhoe = data[1];
  rel_FT   = data[2];
  abs_FT   = data[3];
  rel_T    = data[4];
  abs_T    = data[5];
#ifdef NEUTRINO
  rel_FY   = data[6];
  abs_FY   = data[7];
  rel_Ye   = data[8];
  abs_Ye   = data[9];
#endif      
}


void Radiation::compute_coupling(MultiFab& coupT, MultiFab& coupY, 
				 const MultiFab& kpp, const MultiFab& Eg,
				 const MultiFab& jg)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(kpp,true); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
#ifdef NEUTRINO 
	BL_FORT_PROC_CALL(CA_COMPUTE_COUPTY, ca_compute_coupty)
	    (bx.loVect(), bx.hiVect(),
	     BL_TO_FORTRAN(coupT[mfi]),
	     BL_TO_FORTRAN(coupY[mfi]),
	     BL_TO_FORTRAN(kpp[mfi]),
	     BL_TO_FORTRAN(Eg[mfi]),    
	     BL_TO_FORTRAN(jg[mfi]));
#else
	BL_FORT_PROC_CALL(CA_COMPUTE_COUPT, ca_compute_coupt)
	    (bx.loVect(), bx.hiVect(),
	     BL_TO_FORTRAN(coupT[mfi]),
	     BL_TO_FORTRAN(kpp[mfi]),
	     BL_TO_FORTRAN(Eg[mfi]),    
	     BL_TO_FORTRAN(jg[mfi]));
#endif
    }
}


void Radiation::compute_eta_theta(MultiFab& etaT, MultiFab& etaTz, 
				  MultiFab& etaY, MultiFab& etaYz, 
				  MultiFab& eta1,
				  MultiFab& thetaT, MultiFab& thetaTz, 
				  MultiFab& thetaY, MultiFab& thetaYz, 
				  MultiFab& theta1,
				  MultiFab& djdT, MultiFab& djdY, 
				  const MultiFab& dkdT, const MultiFab& dkdY, 
				  const MultiFab& dedT, const MultiFab& dedY, 
				  const MultiFab& Er_star, const MultiFab& rho, 
				  const BoxArray& grids, Real delta_t, Real ptc_tau)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(rho,true); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

#ifdef NEUTRINO
	BL_FORT_PROC_CALL(CA_COMPUTE_ETA_THE, ca_compute_eta_the)
	    (bx.loVect(), bx.hiVect(),
	      BL_TO_FORTRAN(etaT[mfi]),
	      BL_TO_FORTRAN(etaTz[mfi]),
	      BL_TO_FORTRAN(etaY[mfi]),
	      BL_TO_FORTRAN(etaYz[mfi]),
	      BL_TO_FORTRAN(eta1[mfi]),
	      BL_TO_FORTRAN(thetaT[mfi]),
	      BL_TO_FORTRAN(thetaTz[mfi]),
	      BL_TO_FORTRAN(thetaY[mfi]),
	      BL_TO_FORTRAN(thetaYz[mfi]),
	      BL_TO_FORTRAN(theta1[mfi]),
	      BL_TO_FORTRAN(djdT[mfi]),
	      BL_TO_FORTRAN(djdY[mfi]),
	      BL_TO_FORTRAN(dkdT[mfi]),
	      BL_TO_FORTRAN(dkdY[mfi]),
	      BL_TO_FORTRAN(dedT[mfi]),
	      BL_TO_FORTRAN(dedY[mfi]),
	      BL_TO_FORTRAN(Er_star[mfi]),
	      BL_TO_FORTRAN(rho[mfi]),
	      &delta_t, &ptc_tau);
#else
	BL_FORT_PROC_CALL(CA_COMPUTE_ETAT, ca_compute_etat)
	    (bx.loVect(), bx.hiVect(),
	      BL_TO_FORTRAN(etaT[mfi]),
	      BL_TO_FORTRAN(etaTz[mfi]),
	      BL_TO_FORTRAN(eta1[mfi]),
	      BL_TO_FORTRAN(djdT[mfi]),
	      BL_TO_FORTRAN(dkdT[mfi]),
	      BL_TO_FORTRAN(dedT[mfi]),
	      BL_TO_FORTRAN(Er_star[mfi]),
	      BL_TO_FORTRAN(rho[mfi]),
	      &delta_t, &ptc_tau);
#endif
    }
}


void Radiation::eos_opacity_emissivity(const MultiFab& S_new, 
				       const MultiFab& temp_new, const MultiFab& Ye_new, 
				       const MultiFab& temp_star, const MultiFab& Ye_star,
				       MultiFab& kappa_p, MultiFab& kappa_r, MultiFab& jg, 
				       MultiFab& djdT, MultiFab& djdY, 
				       MultiFab& dkdT, MultiFab& dkdY, 
				       MultiFab& dedT, MultiFab& dedY, 
				       int level, const BoxArray& grids, int it, int ngrow)
{
  int star_is_valid = 1 - ngrow;

  int lag_opac;
  if (it == 1) {
    lag_opac = 0;
  }
  else if (it <= update_opacity) {
    lag_opac = 0;
  }
  else {
    lag_opac = 1;
  }

  const Geometry& geom = parent->Geom(level);

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi) {
      const Box& box = mfi.tilebox();
#ifdef NEUTRINO
      if (radiation_type == Neutrino) {
	  BL_FORT_PROC_CALL(CA_COMPUTE_DEDX, ca_compute_dedx)
	      (box.loVect(), box.hiVect(),
	       BL_TO_FORTRAN(S_new[mfi]),
	       BL_TO_FORTRAN(temp_new[mfi]),
	       BL_TO_FORTRAN(Ye_new[mfi]),
	       BL_TO_FORTRAN(temp_star[mfi]),
	       BL_TO_FORTRAN(Ye_star[mfi]),
	       BL_TO_FORTRAN(dedT[mfi]),
	       BL_TO_FORTRAN(dedY[mfi]),
	       &star_is_valid);
      }
      else {
	  dedY[mfi].setVal(0.0,box,0);
#endif
	  if (do_real_eos == 1) {
#pragma gpu box(box) sync
	    ca_compute_c_v
                (AMREX_INT_ANYD(box.loVect()), AMREX_INT_ANYD(box.hiVect()),
                 BL_TO_FORTRAN_ANYD(dedT[mfi]),
                 BL_TO_FORTRAN_ANYD(temp_new[mfi]),
                 BL_TO_FORTRAN_ANYD(S_new[mfi]));
	  }
	  else if (c_v_exp_m == 0.0 && c_v_exp_n == 0.0) {
	      dedT[mfi].setVal(const_c_v,box,0);
	  }
	  else if (const_c_v > 0.0) {
#pragma gpu box(box) sync
	      gcv(AMREX_INT_ANYD(box.loVect()), AMREX_INT_ANYD(box.hiVect()),
		  BL_TO_FORTRAN_ANYD(dedT[mfi]), 
		  BL_TO_FORTRAN_ANYD(temp_new[mfi]), 
		  const_c_v, c_v_exp_m, c_v_exp_n,
		  prop_temp_floor,
		  BL_TO_FORTRAN_ANYD(S_new[mfi]));
	  }
	  else {
	      amrex::Error("ERROR Radiation::eos_opacity_emissivity");
	  }
#ifdef NEUTRINO
      }
#endif
  }

  if (dedT_fac > 1.0) {
    dedT.mult(dedT_fac);
  }
#ifdef NEUTRINO
  if (dedY_fac > 1.0) {
    dedY.mult(dedY_fac);
  }
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.growntilebox(ngrow);
#ifdef NEUTRINO
      if (radiation_type == Neutrino) {
	  BL_FORT_PROC_CALL(CA_OPAC_EMIS_NEUT, ca_opac_emis_neut)
	      (bx.loVect(), bx.hiVect(),
	       BL_TO_FORTRAN(S_new[mfi]),
	       BL_TO_FORTRAN(temp_new[mfi]),
	       BL_TO_FORTRAN(Ye_new[mfi]),
	       BL_TO_FORTRAN(temp_star[mfi]),
	       BL_TO_FORTRAN(Ye_star[mfi]),
	       BL_TO_FORTRAN(kappa_p[mfi]),
	       BL_TO_FORTRAN(kappa_r[mfi]),
	       BL_TO_FORTRAN(jg[mfi]),
	       BL_TO_FORTRAN(djdT[mfi]),
	       BL_TO_FORTRAN(djdY[mfi]),
	       BL_TO_FORTRAN(dkdT[mfi]),
	       BL_TO_FORTRAN(dkdY[mfi]),
	       &use_dkdT, &star_is_valid, &lag_opac);
      }
      else {
	  djdY[mfi].setVal(0.0,bx,0,nGroups);
	  dkdY[mfi].setVal(0.0,bx,0,nGroups);
#endif

	  if (use_opacity_table_module) {

	      BL_FORT_PROC_CALL(CA_OPACS, ca_opacs)
		  (bx.loVect(), bx.hiVect(),
		   BL_TO_FORTRAN(S_new[mfi]),
		   BL_TO_FORTRAN(temp_new[mfi]),
		   BL_TO_FORTRAN(temp_star[mfi]),
		   BL_TO_FORTRAN(kappa_p[mfi]),
		   BL_TO_FORTRAN(kappa_r[mfi]),
		   BL_TO_FORTRAN(dkdT[mfi]),
		   &use_dkdT, &star_is_valid, &lag_opac);

	  }
	  else {  // use power-law 
	      
	      BL_FORT_PROC_CALL(CA_COMPUTE_KAPPAS, ca_compute_kappas)
		  (bx.loVect(), bx.hiVect(),
		   BL_TO_FORTRAN(S_new[mfi]),
		   BL_TO_FORTRAN(temp_new[mfi]),
		   BL_TO_FORTRAN(kappa_p[mfi]),
		   BL_TO_FORTRAN(kappa_r[mfi]),
		   BL_TO_FORTRAN(dkdT[mfi]), 
		   &do_kappa_stm_emission, &use_dkdT,
		   &const_kappa_p, &kappa_p_exp_m, 
		   &kappa_p_exp_n, &kappa_p_exp_p,
		   &const_kappa_r, &kappa_r_exp_m, 
		   &kappa_r_exp_n, &kappa_r_exp_p,
		   &const_scattering, &scattering_exp_m, 
		   &scattering_exp_n, &scattering_exp_p,
		   &prop_temp_floor);
	      
	  }

	  const Box& reg = mfi.tilebox();
      
	  Vector<Real> PFcoef(nGroups, -1.0); // picket-fence model coefficients
#ifdef MG_SU_OLSON
	  PFcoef[0] = 0.5;
	  PFcoef[1] = 0.5;
#endif
	  BL_FORT_PROC_CALL(CA_COMPUTE_EMISSIVITY, ca_compute_emissivity)
	      (reg.loVect(), reg.hiVect(),
	       BL_TO_FORTRAN(jg[mfi]),  
	       BL_TO_FORTRAN(djdT[mfi]),  
	       BL_TO_FORTRAN(temp_new[mfi]),
	       BL_TO_FORTRAN(kappa_p[mfi]),
	       BL_TO_FORTRAN(dkdT[mfi]),
	       PFcoef.dataPtr(), 
	       use_WiensLaw, integrate_Planck, Tf_Wien);

#ifdef NEUTRINO
      }
#endif
  }    

  if (ngrow == 0 && !lag_opac) {
      kappa_r.FillBoundary(geom.periodicity());
  }
}


void Radiation::gray_accel(MultiFab& Er_new, MultiFab& Er_pi, 
			   MultiFab& kappa_p, MultiFab& kappa_r,
			   MultiFab& etaT, MultiFab& etaY, MultiFab& eta1,
			   MultiFab& thetaT, MultiFab& thetaY, 
			   MultiFab& mugT, MultiFab& mugY, 
			   Array<MultiFab, BL_SPACEDIM>& lambda,
			   RadSolve& solver, MGRadBndry& mgbd, 
			   const BoxArray& grids, int level, Real time, 
			   Real delta_t, Real ptc_tau)
{
  const Geometry& geom = parent->Geom(level);
  const Real* dx = parent->Geom(level).CellSize();
  const Castro *castro = dynamic_cast<Castro*>(&parent->getLevel(level));
  const DistributionMapping& dmap = castro->DistributionMap();

  if (nGroups > 1) {
    solver.setHypreMulti(1.0);
  }
  else {
    solver.setHypreMulti(0.0);
  }
  mgbd.setCorrection();

  MultiFab Er_zero(grids, dmap, 1, 0);
  Er_zero.setVal(0.0);
  getBndryDataMG_ga(mgbd, Er_zero, level);

  MultiFab spec(grids, dmap, nGroups, 1);
#ifdef _OPENMP
#pragma omp parallel
#endif 
  for (MFIter mfi(spec,true); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
#ifdef NEUTRINO
    BL_FORT_PROC_CALL(CA_ACCEL_SPEC_NEUT, ca_accel_spec_neut) 
      (bx.loVect(), bx.hiVect(),
       BL_TO_FORTRAN(Er_new[mfi]),
       BL_TO_FORTRAN(Er_pi[mfi]),
       BL_TO_FORTRAN(kappa_p[mfi]),
       BL_TO_FORTRAN(etaT[mfi]),
       BL_TO_FORTRAN(etaY[mfi]),
       BL_TO_FORTRAN(thetaT[mfi]),
       BL_TO_FORTRAN(thetaY[mfi]),
       BL_TO_FORTRAN(mugT[mfi]),
       BL_TO_FORTRAN(mugY[mfi]),
       BL_TO_FORTRAN(spec[mfi]),
       &delta_t, &ptc_tau);
#else
    BL_FORT_PROC_CALL(CA_ACCEL_SPEC, ca_accel_spec) 
      (bx.loVect(), bx.hiVect(),
       BL_TO_FORTRAN(kappa_p[mfi]),
       BL_TO_FORTRAN(mugT[mfi]),
       BL_TO_FORTRAN(spec[mfi]),
       &delta_t, &ptc_tau);
#endif
  }

  // Extrapolate spectrum out one cell
  for (int indx = 0; indx < nGroups; indx++) {
    extrapolateBorders(spec, indx);
  }
  // Overwrite all extrapolated components with values from
  // neighboring fine grids where they exist:
  spec.FillBoundary(parent->Geom(level).periodicity());

  // set boundary condition
  solver.levelBndry(mgbd,0);

  // A coefficients
  MultiFab acoefs(grids, dmap, 1, 0);
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(acoefs,true); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
#ifdef NEUTRINO
      BL_FORT_PROC_CALL(CA_ACCEL_ACOE_NEUT, ca_accel_acoe_neut)
	  (bx.loVect(), bx.hiVect(),
	   BL_TO_FORTRAN(eta1[mfi]),
	   BL_TO_FORTRAN(thetaT[mfi]),
	   BL_TO_FORTRAN(thetaY[mfi]),
	   BL_TO_FORTRAN(spec[mfi]),
	   BL_TO_FORTRAN(kappa_p[mfi]),
	   BL_TO_FORTRAN(acoefs[mfi]),
	   &delta_t, &ptc_tau);    
#else
      BL_FORT_PROC_CALL(CA_ACCEL_ACOE, ca_accel_acoe)
	  (bx.loVect(), bx.hiVect(),
	   BL_TO_FORTRAN(eta1[mfi]),
	   BL_TO_FORTRAN(spec[mfi]),
	   BL_TO_FORTRAN(kappa_p[mfi]),
	   BL_TO_FORTRAN(acoefs[mfi]),
	   &delta_t, &ptc_tau);    
#endif
  }  

  solver.cellCenteredApplyMetrics(level, acoefs);
  solver.setLevelACoeffs(level, acoefs);

  const DistributionMapping& dm = castro->DistributionMap();

  // B & C coefficients
  Array<MultiFab, BL_SPACEDIM> bcoefs, ccoefs, bcgrp;
  for (int idim = 0; idim < BL_SPACEDIM; idim++) {
    const BoxArray& edge_boxes = castro->getEdgeBoxArray(idim);

    bcoefs[idim].define(edge_boxes, dm, 1, 0);
    bcoefs[idim].setVal(0.0);

    bcgrp [idim].define(edge_boxes, dm, 1, 0);

    if (nGroups > 1) {
	ccoefs[idim].define(edge_boxes, dm, 2, 0);
	ccoefs[idim].setVal(0.0);
    }
  }

  for (int igroup = 0; igroup < nGroups; igroup++) {
    for (int idim=0; idim<BL_SPACEDIM; idim++) {
      solver.computeBCoeffs(bcgrp[idim], idim, kappa_r, igroup,
			    lambda[idim], igroup, c, geom);
      // metrics is already in bcgrp

#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(spec,true); mfi.isValid(); ++mfi) {
	  const Box&  bx  = mfi.nodaltilebox(idim);
	  const Box& bbox = bcoefs[idim][mfi].box();

	  lbcoefna(bcoefs[idim][mfi].dataPtr(),
		   bcgrp[idim][mfi].dataPtr(),
		   ARLIM(bbox.loVect()), ARLIM(bbox.hiVect()),
		   ARLIM(bx.loVect()), ARLIM(bx.hiVect()),
		   BL_TO_FORTRAN_N(spec[mfi], igroup), 
		   idim);
	  
	  if (nGroups > 1) {
	      BL_FORT_PROC_CALL(CA_ACCEL_CCOE, ca_accel_ccoe)
		  (bx.loVect(), bx.hiVect(),
		   BL_TO_FORTRAN(bcgrp[idim][mfi]),
		   BL_TO_FORTRAN(spec[mfi]),
		   BL_TO_FORTRAN(ccoefs[idim][mfi]),
		   dx, &idim, &igroup);
	  }
      }
    }
  }

  for (int idim = 0; idim < BL_SPACEDIM; idim++) {
    solver.setLevelBCoeffs(level, bcoefs[idim], idim);

    if (nGroups > 1) {
      solver.setLevelCCoeffs(level, ccoefs[idim], idim);
    }
  }

  // rhs
  MultiFab rhs(grids,dmap,1,0);
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(rhs,true); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
#ifdef NEUTRINO
      BL_FORT_PROC_CALL(CA_ACCEL_RHS_NEUT, ca_accel_rhs_neut) 
	  (bx.loVect(), bx.hiVect(),
	   BL_TO_FORTRAN(Er_new[mfi]),
	   BL_TO_FORTRAN(Er_pi[mfi]),
	   BL_TO_FORTRAN(kappa_p[mfi]),
	   BL_TO_FORTRAN(etaT[mfi]),
	   BL_TO_FORTRAN(etaY[mfi]),
	   BL_TO_FORTRAN(thetaT[mfi]),
	   BL_TO_FORTRAN(thetaY[mfi]),
	   BL_TO_FORTRAN(rhs[mfi]),
	   &delta_t);
#else
      BL_FORT_PROC_CALL(CA_ACCEL_RHS, ca_accel_rhs) 
	  (bx.loVect(), bx.hiVect(),
	   BL_TO_FORTRAN(Er_new[mfi]),
	   BL_TO_FORTRAN(Er_pi[mfi]),
	   BL_TO_FORTRAN(kappa_p[mfi]),
	   BL_TO_FORTRAN(etaT[mfi]),
	   BL_TO_FORTRAN(rhs[mfi]),
	   &delta_t);
#endif
  }

  // must apply metrics to rhs here
  solver.cellCenteredApplyMetrics(level, rhs);

  // solve
  MultiFab accel(grids,dmap,1,0);
  accel.setVal(0.0);
  solver.levelSolve(level, accel, 0, rhs, 0.01);

  // update Er_new
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(spec,true); mfi.isValid(); ++mfi) {
      const Box& reg  = mfi.tilebox();

      ljupna(BL_TO_FORTRAN(Er_new[mfi]), 
	     ARLIM(reg.loVect()), ARLIM(reg.hiVect()),
	     BL_TO_FORTRAN(spec[mfi]),
	     BL_TO_FORTRAN(accel[mfi]), 
	     nGroups);
  }

  mgbd.unsetCorrection();
  getBndryDataMG(mgbd, Er_new, time, level);

  solver.restoreHypreMulti();
}


void Radiation::local_accel(MultiFab& Er_new, const MultiFab& Er_pi,
			    const MultiFab& kappa_p, 
			    const MultiFab& etaT, const MultiFab& etaY, 
			    const MultiFab& thetaT, const MultiFab& thetaY, 
			    const MultiFab& mugT, const MultiFab& mugY, 
			    Real delta_t, Real ptc_tau)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(Er_new,true); mfi.isValid(); ++mfi) {
	const Box& bx = mfi.tilebox();

#ifdef NEUTRINO
	BL_FORT_PROC_CALL(CA_LOCAL_ACCEL_NEUT, ca_local_accel_neut)
	    (bx.loVect(), bx.hiVect(),
	     BL_TO_FORTRAN(Er_new[mfi]),
	     BL_TO_FORTRAN(Er_pi[mfi]),
	     BL_TO_FORTRAN(kappa_p[mfi]),
	     BL_TO_FORTRAN(etaT[mfi]),
	     BL_TO_FORTRAN(etaY[mfi]),
	     BL_TO_FORTRAN(thetaT[mfi]),
	     BL_TO_FORTRAN(thetaY[mfi]),
	     BL_TO_FORTRAN(mugT[mfi]),
	     BL_TO_FORTRAN(mugY[mfi]),
	     &delta_t, &ptc_tau);
#else
	BL_FORT_PROC_CALL(CA_LOCAL_ACCEL, ca_local_accel)
	    (bx.loVect(), bx.hiVect(),
	     BL_TO_FORTRAN(Er_new[mfi]),
	     BL_TO_FORTRAN(Er_pi[mfi]),
	     BL_TO_FORTRAN(kappa_p[mfi]),
	     BL_TO_FORTRAN(etaT[mfi]),
	     BL_TO_FORTRAN(mugT[mfi]),
	     &delta_t, &ptc_tau);
#endif
    }
}


void Radiation::state_energy_update(MultiFab& state, const MultiFab& rhoe, 
				    const MultiFab& Ye, const MultiFab& temp, 
				    const BoxArray& grids,
				    Real& derat, Real& dT, Real& dye, int level)
{
  BL_PROFILE("Radiation::state_energy_update (MGFLD)");

  const DistributionMapping& dmap = state.DistributionMap();

  MultiFab msk(grids,dmap,1,0);

  {
      BoxArray baf;

      if (level < parent->finestLevel()) {
	  baf = parent->boxArray(level+1);
	  baf.coarsen(parent->refRatio(level));
      }

#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(msk); mfi.isValid(); ++mfi) 
      {
	  msk[mfi].setVal(1.0);

	  // Now mask off finer level also:
	  if (level < parent->finestLevel()) {
	      std::vector< std::pair<int,Box> > isects = baf.intersections(mfi.validbox());
	      for (int ii = 0; ii < isects.size(); ii++) {
		  msk[mfi].setVal(0.0, isects[ii].second, 0);
	      }
	  } 
      }
  }

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(state,true); mfi.isValid(); ++mfi) 
  {
      const Box& reg  = mfi.tilebox();
#ifdef NEUTRINO
      BL_FORT_PROC_CALL(CA_STATE_UPDATE_NEUT, ca_state_update_neut)
	  (reg.loVect(), reg.hiVect(),
	   BL_TO_FORTRAN(state[mfi]),
	   BL_TO_FORTRAN(rhoe[mfi]),
	   BL_TO_FORTRAN(Ye[mfi]),
	   BL_TO_FORTRAN(temp[mfi]),
	   BL_TO_FORTRAN(msk[mfi]),
	   &derat, &dT, &dye);
#else
      BL_FORT_PROC_CALL(CA_STATE_UPDATE, ca_state_update)
	  (reg.loVect(), reg.hiVect(),
	   BL_TO_FORTRAN(state[mfi]),
	   BL_TO_FORTRAN(rhoe[mfi]),
	   BL_TO_FORTRAN(temp[mfi]),
	   BL_TO_FORTRAN(msk[mfi]),
	   &derat, &dT);
#endif
  }

  ParallelDescriptor::ReduceRealMax(derat);
  ParallelDescriptor::ReduceRealMax(dT);
#ifdef NEUTRINO
  ParallelDescriptor::ReduceRealMax(dye);
#endif
}


void Radiation::update_matter(MultiFab& rhoe_new, MultiFab& temp_new, 
			      MultiFab& rhoYe_new, MultiFab& Ye_new, 
			      const MultiFab& Er_new, const MultiFab& Er_pi,
			      const MultiFab& rhoe_star, const MultiFab& rhoYe_star, 
			      const MultiFab& rhoe_step, const MultiFab& rhoYe_step,
			      const MultiFab& etaT, const MultiFab& etaTz,
			      const MultiFab& etaY, const MultiFab& etaYz,
			      const MultiFab& eta1,
			      const MultiFab& thetaT, const MultiFab& thetaTz, 
			      const MultiFab& thetaY, const MultiFab& thetaYz, 
			      const MultiFab& theta1,
			      const MultiFab& coupT, const MultiFab& coupY, 
			      const MultiFab& kappa_p, const MultiFab& jg, 
			      const MultiFab& mugT, const MultiFab& mugY, 
			      const MultiFab& S_new, 
			      int level, Real delta_t, Real ptc_tau,
			      int it, bool conservative_update)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(rhoe_new,true); mfi.isValid(); ++mfi) {
	const Box& bx = mfi.tilebox(); 

#ifdef NEUTRINO    
	if (conservative_update) {
	    BL_FORT_PROC_CALL(CA_UPDATE_MATTER_NEUT, ca_update_matter_neut)
		(bx.loVect(), bx.hiVect(),
		 BL_TO_FORTRAN(rhoe_new[mfi]),
		 BL_TO_FORTRAN(rhoYe_new[mfi]),
		 BL_TO_FORTRAN(Ye_new[mfi]),
		 BL_TO_FORTRAN(Er_new[mfi]),
		 BL_TO_FORTRAN(Er_pi[mfi]),
		 BL_TO_FORTRAN(rhoe_star[mfi]),
		 BL_TO_FORTRAN(rhoYe_star[mfi]),
		 BL_TO_FORTRAN(rhoe_step[mfi]),
		 BL_TO_FORTRAN(rhoYe_step[mfi]),
		 BL_TO_FORTRAN(etaT[mfi]),
		 BL_TO_FORTRAN(etaY[mfi]),
		 BL_TO_FORTRAN(eta1[mfi]),
		 BL_TO_FORTRAN(thetaT[mfi]),
		 BL_TO_FORTRAN(thetaY[mfi]),
		 BL_TO_FORTRAN(theta1[mfi]),
		 BL_TO_FORTRAN(coupT[mfi]),
		 BL_TO_FORTRAN(coupY[mfi]),
		 BL_TO_FORTRAN(kappa_p[mfi]),
		 BL_TO_FORTRAN(mugT[mfi]),
		 BL_TO_FORTRAN(mugY[mfi]),
		 BL_TO_FORTRAN(S_new[mfi]),
		 &delta_t, &ptc_tau);

	    if (radiation_type == Neutrino) {
#pragma gpu box(bx) sync
		ca_compute_temp_given_reye
		    (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
		     BL_TO_FORTRAN_ANYD(temp_new[mfi]), 
		     BL_TO_FORTRAN_ANYD(rhoe_new[mfi]),
		     BL_TO_FORTRAN_ANYD(Ye_new[mfi]),
		     BL_TO_FORTRAN_ANYD(S_new[mfi]));
	    }
	    else {
		temp_new[mfi].copy(rhoe_new[mfi],bx);
		
		if (do_real_eos > 0) {
#pragma gpu box(bx) sync
		  ca_compute_temp_given_rhoe
                      (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()), 
                       BL_TO_FORTRAN_ANYD(temp_new[mfi]), 
                       BL_TO_FORTRAN_ANYD(S_new[mfi]));
		}
		else if (do_real_eos == 0) {
#pragma gpu box(bx) sync
		  ca_compute_temp_given_cv
                      (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                       BL_TO_FORTRAN_ANYD(temp_new[mfi]), 
                       BL_TO_FORTRAN_ANYD(S_new[mfi]),
                       const_c_v, c_v_exp_m, c_v_exp_n);
		}
		else {
		    amrex::Error("ERROR Radiation::do_real_eos < 0");
		}
	    }
	}
	else {
	    BL_FORT_PROC_CALL(CA_NCUPDATE_MATTER_NEUT, ca_ncupdate_matter_neut)
		(bx.loVect(), bx.hiVect(),
		 BL_TO_FORTRAN(temp_new[mfi]),
		 BL_TO_FORTRAN(Ye_new[mfi]),
		 BL_TO_FORTRAN(Er_new[mfi]),
		 BL_TO_FORTRAN(rhoe_star[mfi]),
		 BL_TO_FORTRAN(rhoYe_star[mfi]),
		 BL_TO_FORTRAN(rhoe_step[mfi]),
		 BL_TO_FORTRAN(rhoYe_step[mfi]),
		 BL_TO_FORTRAN(etaTz[mfi]),
		 BL_TO_FORTRAN(etaYz[mfi]),
		 BL_TO_FORTRAN(thetaTz[mfi]),
		 BL_TO_FORTRAN(thetaYz[mfi]),
		 BL_TO_FORTRAN(kappa_p[mfi]),
		 BL_TO_FORTRAN(jg[mfi]),
		 &delta_t);
	    
	    BL_FORT_PROC_CALL(CA_COMPUTE_REYE_GIVEN_TY,ca_compute_reye_given_ty)
		(bx.loVect(), bx.hiVect(),
		 BL_TO_FORTRAN(rhoe_new[mfi]),
		 BL_TO_FORTRAN(rhoYe_new[mfi]),
		 BL_TO_FORTRAN(temp_new[mfi]), 
		 BL_TO_FORTRAN(Ye_new[mfi]),
		 BL_TO_FORTRAN(S_new[mfi]));
	}
#else
	if (conservative_update) {
	    BL_FORT_PROC_CALL(CA_UPDATE_MATTER, ca_update_matter)
		(bx.loVect(), bx.hiVect(),
		 BL_TO_FORTRAN(rhoe_new[mfi]),
		 BL_TO_FORTRAN(Er_new[mfi]),
		 BL_TO_FORTRAN(Er_pi[mfi]),
		 BL_TO_FORTRAN(rhoe_star[mfi]),
		 BL_TO_FORTRAN(rhoe_step[mfi]),
		 BL_TO_FORTRAN(eta1[mfi]),
		 BL_TO_FORTRAN(coupT[mfi]),
		 BL_TO_FORTRAN(kappa_p[mfi]),
		 BL_TO_FORTRAN(mugT[mfi]),
		 BL_TO_FORTRAN(S_new[mfi]),
		 &delta_t, &ptc_tau);
	    
	    temp_new[mfi].copy(rhoe_new[mfi],bx);

	    if (do_real_eos > 0) {
#pragma gpu box(bx) sync
	      ca_compute_temp_given_rhoe
                  (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()), 
                   BL_TO_FORTRAN_ANYD(temp_new[mfi]), 
                   BL_TO_FORTRAN_ANYD(S_new[mfi]));
	    }
	    else if (do_real_eos == 0) {
#pragma gpu box(bx) sync
	      ca_compute_temp_given_cv
                  (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                   BL_TO_FORTRAN_ANYD(temp_new[mfi]), 
                   BL_TO_FORTRAN_ANYD(S_new[mfi]),
                   const_c_v, c_v_exp_m, c_v_exp_n);
	    }
	    else {
		amrex::Error("ERROR Radiation::do_real_eos < 0");
	    }
	}
	else {
	    BL_FORT_PROC_CALL(CA_NCUPDATE_MATTER, ca_ncupdate_matter)
		(bx.loVect(), bx.hiVect(),
		 BL_TO_FORTRAN(temp_new[mfi]),
		 BL_TO_FORTRAN(Er_new[mfi]),
		 BL_TO_FORTRAN(rhoe_star[mfi]),
		 BL_TO_FORTRAN(rhoe_step[mfi]),
		 BL_TO_FORTRAN(etaTz[mfi]),
		 BL_TO_FORTRAN(kappa_p[mfi]),
		 BL_TO_FORTRAN(jg[mfi]),
		 &delta_t);

#pragma gpu box(bx) sync
	    ca_get_rhoe
                (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                 BL_TO_FORTRAN_ANYD(rhoe_new[mfi]),
                 BL_TO_FORTRAN_ANYD(temp_new[mfi]), 
                 BL_TO_FORTRAN_ANYD(S_new[mfi]));
	}
#endif
    }  
}

// ========================================================================
// for the hyperbolic solver

void Radiation::compute_limiter(int level, const BoxArray& grids,
				const MultiFab &Sborder, 
				const MultiFab &Erborder,
				MultiFab &lamborder)
{ // it works for both single- and multi-group
  int ngrow = lamborder.nGrow();
  const DistributionMapping& dmap = lamborder.DistributionMap();

  if (filter_lambda_T) {
    BL_ASSERT(ngrow == 4);
  }

  if (limiter == 0) {

    lamborder.setVal(1./3., ngrow);

  }
  else {

    MultiFab kpr(grids,dmap,Radiation::nGroups,ngrow);  

    if (do_multigroup) {
	MGFLD_compute_rosseland(kpr, Sborder); 
    }
    else {
      SGFLD_compute_rosseland(kpr, Sborder); 
    }

    MultiFab Er_wide(grids, dmap, nGroups, ngrow+1);
    Er_wide.setVal(-1.0);
    MultiFab::Copy(Er_wide, Erborder, 0, 0, nGroups, 0);
    
    Er_wide.FillBoundary(parent->Geom(level).periodicity());
    
    const Real* dx = parent->Geom(level).CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif    
    for (MFIter mfi(Er_wide,false); mfi.isValid(); ++mfi) {
      BL_FORT_PROC_CALL(CA_COMPUTE_LAMBORDER, ca_compute_lamborder)
	(BL_TO_FORTRAN(Er_wide[mfi]), 
	 BL_TO_FORTRAN(kpr[mfi]),
	 BL_TO_FORTRAN(lamborder[mfi]), 
	 dx, &ngrow, &limiter, &filter_lambda_T, &filter_lambda_S);
    }

    if (filter_lambda_T) {
	lamborder.FillBoundary(parent->Geom(level).periodicity());
    }
  }
}


void Radiation::estimate_gamrPr(const FArrayBox& state, const FArrayBox& Er, 
				FArrayBox& gPr, const Real*dx, const Box& box)
{
  if (limiter == 0) {
    BL_FORT_PROC_CALL(CA_EST_GPR0, ca_est_gpr0)
      (BL_TO_FORTRAN(Er),
       BL_TO_FORTRAN(gPr));
  }
  else {

    FArrayBox kappa_r(gPr.box(), nGroups);

    if (do_multigroup) {
      MGFLD_compute_rosseland(kappa_r, state);
    }
    else {
      SGFLD_compute_rosseland(kappa_r, state);
    }

    BL_FORT_PROC_CALL(CA_EST_GPR2, ca_est_gpr2)
      (BL_TO_FORTRAN(kappa_r),
       BL_TO_FORTRAN(Er),
       BL_TO_FORTRAN(gPr), 
       box.loVect(),box.hiVect(),
       dx, &limiter, &comoving);
  }
}


void Radiation::MGFLD_compute_rosseland(FArrayBox& kappa_r, const FArrayBox& state)
{
  BL_PROFILE("Radiation::MGFLD_compute_rosseland (FArrayBox)");

  const Box& kbox = kappa_r.box();

#ifdef NEUTRINO
  if (radiation_type == Neutrino) {
    BL_FORT_PROC_CALL(CA_COMPUTE_ROSSELAND_NEUT, ca_compute_rosseland_neut)
	(kbox.loVect(), kbox.hiVect(),
	 BL_TO_FORTRAN(kappa_r), BL_TO_FORTRAN(state));
  }
  else {
#endif
    
    if (use_opacity_table_module) {
#pragma gpu box(kbox) sync
        ca_compute_rosseland(AMREX_INT_ANYD(kbox.loVect()), AMREX_INT_ANYD(kbox.hiVect()),
                             BL_TO_FORTRAN_ANYD(kappa_r),
                             BL_TO_FORTRAN_ANYD(state));
    }
    else if (const_kappa_r < 0.0) {
      ca_compute_powerlaw_kappa_s(kbox.loVect(), kbox.hiVect(),
				  BL_TO_FORTRAN(kappa_r), BL_TO_FORTRAN(state),
				  &const_kappa_p, &kappa_p_exp_m, &kappa_p_exp_n, &kappa_p_exp_p, 
				  &const_scattering, &scattering_exp_m, &scattering_exp_n, &scattering_exp_p, 
				  &prop_temp_floor, &kappa_r_floor);	 
    }
    else {
      ca_compute_powerlaw_kappa(kbox.loVect(), kbox.hiVect(),
				BL_TO_FORTRAN(kappa_r), BL_TO_FORTRAN(state),
				&const_kappa_r, &kappa_r_exp_m, &kappa_r_exp_n, &kappa_r_exp_p, 
				&prop_temp_floor, &kappa_r_floor);
    }
#ifdef NEUTRINO
  }
#endif
}


void Radiation::MGFLD_compute_rosseland(MultiFab& kappa_r, const MultiFab& state)
{
    BL_PROFILE("Radiation::MGFLD_compute_rosseland (MultiFab)");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(kappa_r,true); mfi.isValid(); ++mfi) {
	const Box& bx = mfi.growntilebox();
#ifdef NEUTRINO
	if (radiation_type == Neutrino) {
	    BL_FORT_PROC_CALL(CA_COMPUTE_ROSSELAND_NEUT, ca_compute_rosseland_neut)
		(bx.loVect(), bx.hiVect(),
		 BL_TO_FORTRAN(kappa_r[mfi]), BL_TO_FORTRAN(state[mfi]));
	}
	else {
#endif
	    if (use_opacity_table_module) {
#pragma gpu box(bx) sync
                ca_compute_rosseland(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                                     BL_TO_FORTRAN_ANYD(kappa_r[mfi]),
                                     BL_TO_FORTRAN_ANYD(state[mfi]));
	    }
	    else if (const_kappa_r < 0.0) {
	      ca_compute_powerlaw_kappa_s(bx.loVect(), bx.hiVect(),
					  BL_TO_FORTRAN(kappa_r[mfi]), BL_TO_FORTRAN(state[mfi]),
					  &const_kappa_p, &kappa_p_exp_m, &kappa_p_exp_n, &kappa_p_exp_p, 
					  &const_scattering, &scattering_exp_m, &scattering_exp_n, &scattering_exp_p, 
					  &prop_temp_floor, &kappa_r_floor);	 
	    }
	    else {
	      ca_compute_powerlaw_kappa(bx.loVect(), bx.hiVect(),
					BL_TO_FORTRAN(kappa_r[mfi]), BL_TO_FORTRAN(state[mfi]),
					&const_kappa_r, &kappa_r_exp_m, &kappa_r_exp_n, &kappa_r_exp_p, 
					&prop_temp_floor, &kappa_r_floor);	 
	    }
#ifdef NEUTRINO
	}
#endif
    }
}

void Radiation::MGFLD_compute_scattering(FArrayBox& kappa_s, const FArrayBox& state)
{
    BL_PROFILE("Radiation::MGFLD_compute_scattering");

    const Box& kbox = kappa_s.box();

#ifdef NEUTRINO
    amrex::Abort("MGFLD_compute_scattering: not supposted to be here");
#else
    if (use_opacity_table_module) {
	BL_FORT_PROC_CALL(CA_COMPUTE_SCATTERING, ca_compute_scattering)
		(ARLIM_3D(kbox.loVect()), ARLIM_3D(kbox.hiVect()),
		 BL_TO_FORTRAN_ANYD(kappa_s), 
		 BL_TO_FORTRAN_ANYD(state));
    } else {
	BL_ASSERT(kappa_r_exp_p == 0.0 && kappa_p_exp_p == 0.0 && scattering_exp_p == 0.0);
	if (const_kappa_r < 0.0) {
	  ca_compute_powerlaw_kappa(kbox.loVect(), kbox.hiVect(),
				    BL_TO_FORTRAN(kappa_s), BL_TO_FORTRAN(state),
				    &const_scattering, &scattering_exp_m, &scattering_exp_n, &scattering_exp_p, 
				    &prop_temp_floor, &kappa_r_floor);	 
	    
	} else {
	    BL_FORT_PROC_CALL(CA_COMPUTE_SCATTERING_2, ca_compute_scattering_2)
		(ARLIM_3D(kbox.loVect()), ARLIM_3D(kbox.hiVect()),
		 BL_TO_FORTRAN_ANYD(kappa_s), BL_TO_FORTRAN_ANYD(state),
		 &const_kappa_p, &kappa_p_exp_m, &kappa_p_exp_n, 
		 &const_kappa_r, &kappa_r_exp_m, &kappa_r_exp_n,
		 &prop_temp_floor, &kappa_r_floor);	 
	}
    }
#endif
}

void Radiation::bisect_matter(MultiFab& rhoe_new, MultiFab& temp_new, 
			      MultiFab& rhoYe_new, MultiFab& Ye_new, 
			      const MultiFab& rhoe_star, const MultiFab& temp_star, 
			      const MultiFab& rhoYe_star, const MultiFab& Ye_star, 
			      const MultiFab& S_new, const BoxArray& grids, int level)
{
  temp_new.plus(temp_star, 0, 1, 0);
  temp_new.mult(0.5, 0);
#ifdef NEUTRINO
  Ye_new.plus(Ye_star, 0, 1, 0);
  Ye_new.mult(0.5, 0);
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(rhoe_new,true); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

#ifdef NEUTRINO
#pragma gpu box(bx) sync
      ca_compute_reye_given_ty
	  (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
	   BL_TO_FORTRAN_ANYD(rhoe_new[mfi]),
	   BL_TO_FORTRAN_ANYD(rhoYe_new[mfi]),
	   BL_TO_FORTRAN_ANYD(temp_new[mfi]), 
	   BL_TO_FORTRAN_ANYD(Ye_new[mfi]),
	   BL_TO_FORTRAN_ANYD(S_new[mfi]));
#else
      if (do_real_eos > 0) {
#pragma gpu box(bx) sync
	ca_get_rhoe
            (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
             BL_TO_FORTRAN_ANYD(rhoe_new[mfi]),
             BL_TO_FORTRAN_ANYD(temp_new[mfi]), 
             BL_TO_FORTRAN_ANYD(S_new[mfi]));
      }
      else {
	  amrex::Abort("do_real_eos == 0 not supported in bisect_matter");
      }
#endif
  }
}


void Radiation::rhstoEr(MultiFab& rhs, Real dt, int level)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    {    
	Vector<Real> r, s;

	for (MFIter ri(rhs,true); ri.isValid(); ++ri) 
	{
	    const Box &reg = ri.tilebox();	    

	    RadSolve::getCellCenterMetric(parent->Geom(level), reg, r, s);

	    BL_FORT_PROC_CALL(CA_RHSTOER, ca_rhstoer)
		(reg.loVect(), reg.hiVect(),
		 BL_TO_FORTRAN(rhs[ri]), r.dataPtr(), &dt);
	}
    }
}

void Radiation::inelastic_scattering(int level)
{
    if (do_inelastic_scattering) {
	Real dt = parent->dtLevel(level);
	Castro *castro = dynamic_cast<Castro*>(&parent->getLevel(level));
	MultiFab& S_new = castro->get_new_data(State_Type);
	MultiFab& Er_new = castro->get_new_data(Rad_Type);

#ifdef _OPENMP
#pragma omp parallel
#endif
	{
	    FArrayBox kps;
	    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
	    {
		const Box& bx = mfi.tilebox();

		kps.resize(bx,1); // we assume scattering is independent of nu
		MGFLD_compute_scattering(kps, S_new[mfi]);

		ca_inelastic_sct(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
				 BL_TO_FORTRAN_ANYD(S_new[mfi]),
				 BL_TO_FORTRAN_ANYD(Er_new[mfi]),
				 BL_TO_FORTRAN_ANYD(kps),
				 dt);		
	    }
	}

	if (do_real_eos > 0) {
	  castro->computeTemp(S_new, castro->state[State_Type].curTime(), S_new.nGrow());
	}
    }
}
