#include "Radiation.H"
#include "Castro_F.H"
#include "RAD_F.H"

#include <Using.H>

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

  for (MFIter mfi(Er_new); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box& bx = grids[i]; 

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
       &relative, &absolute, & err_er, &delta_t);
#else
    BL_FORT_PROC_CALL(CA_CHECK_CONV_ER, ca_check_conv_er)
      (bx.loVect(), bx.hiVect(),
       BL_TO_FORTRAN(Er_new[mfi]),
       BL_TO_FORTRAN(Er_pi[mfi]),
       BL_TO_FORTRAN(kappa_p[mfi]),
       BL_TO_FORTRAN(etaTz[mfi]),
       BL_TO_FORTRAN(temp_new[mfi]),
       &relative, &absolute, & err_er, &delta_t);
#endif
  }

  ParallelDescriptor::ReduceRealMax(relative);
  ParallelDescriptor::ReduceRealMax(absolute);
  ParallelDescriptor::ReduceRealMax(err_er);
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

  for (MFIter mfi(rhoe_new); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box& bx = grids[i]; 

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
       &rel_rhoe, &abs_rhoe, 
       &rel_FT,   &abs_FT, 
       &rel_T,    &abs_T, 
       &rel_FY,   &abs_FY, 
       &rel_Ye,   &abs_Ye,
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
       &rel_rhoe, &abs_rhoe, 
       &rel_FT,   &abs_FT, 
       &rel_T,    &abs_T, 
       &delta_t);
#endif
  }

  ParallelDescriptor::ReduceRealMax(rel_rhoe);
  ParallelDescriptor::ReduceRealMax(abs_rhoe);
  ParallelDescriptor::ReduceRealMax(rel_FT);
  ParallelDescriptor::ReduceRealMax(abs_FT);
  ParallelDescriptor::ReduceRealMax(rel_T);
  ParallelDescriptor::ReduceRealMax(abs_T);
#ifdef NEUTRINO
  ParallelDescriptor::ReduceRealMax(rel_FY);
  ParallelDescriptor::ReduceRealMax(abs_FY);
  ParallelDescriptor::ReduceRealMax(rel_Ye);
  ParallelDescriptor::ReduceRealMax(abs_Ye);
#endif
}


void Radiation::compute_coupling(MultiFab& coupT, MultiFab& coupY, 
				 const MultiFab& kpp, const MultiFab& Eg,
				 const MultiFab& jg)
{
  for (MFIter mfi(kpp); mfi.isValid(); ++mfi) {
#ifdef NEUTRINO 
    BL_FORT_PROC_CALL(CA_COMPUTE_COUPTY, ca_compute_coupty)
      (BL_TO_FORTRAN(coupT[mfi]),
       BL_TO_FORTRAN(coupY[mfi]),
       BL_TO_FORTRAN(kpp[mfi]),
       BL_TO_FORTRAN(Eg[mfi]),    
       BL_TO_FORTRAN(jg[mfi]));
#else
    BL_FORT_PROC_CALL(CA_COMPUTE_COUPT, ca_compute_coupt)
      (BL_TO_FORTRAN(coupT[mfi]),
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
  for (MFIter mfi(rho); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box& bx = grids[i];

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

  for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
    int i = mfi.index();
#ifdef NEUTRINO
    if (radiation_type == Neutrino) {
      BL_FORT_PROC_CALL(CA_COMPUTE_DEDX, ca_compute_dedx)
	(BL_TO_FORTRAN(S_new[mfi]),
	 BL_TO_FORTRAN(temp_new[mfi]),
	 BL_TO_FORTRAN(Ye_new[mfi]),
	 BL_TO_FORTRAN(temp_star[mfi]),
	 BL_TO_FORTRAN(Ye_star[mfi]),
	 BL_TO_FORTRAN(dedT[mfi]),
	 BL_TO_FORTRAN(dedY[mfi]),
	 &star_is_valid);
    }
    else {
      dedY[mfi].setVal(0.0);
#endif
      const Box& reg = grids[i];
      if (do_real_eos == 1) {
	BL_FORT_PROC_CALL(CA_COMPUTE_C_V,ca_compute_c_v)
	  (reg.loVect(), reg.hiVect(),
	   BL_TO_FORTRAN(dedT[mfi]), BL_TO_FORTRAN(temp_new[mfi]), BL_TO_FORTRAN(S_new[mfi]));
      }
      else if (c_v_exp_m[0] == 0.0 && c_v_exp_n[0] == 0.0) {
	dedT[mfi].setVal(const_c_v[0]);
      }
      else if (const_c_v[0] > 0.0) {
	  FORT_GCV(dedT[mfi].dataPtr(), dimlist(reg),
		   temp_new[mfi].dataPtr(), dimlist(temp_new[mfi].box()),
		   const_c_v.dataPtr(), c_v_exp_m.dataPtr(), c_v_exp_n.dataPtr(),
		   prop_temp_floor.dataPtr(),
		   S_new[mfi].dataPtr(), dimlist(S_new[mfi].box()));
      }
      else {
	BoxLib::Error("ERROR Radiation::eos_opacity_emissivity");
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

  for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box& bx = BoxLib::grow(grids[i], ngrow);
    
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
      djdY[mfi].setVal(0.0);
      dkdY[mfi].setVal(0.0);
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
	   &const_kappa_p[0], &kappa_p_exp_m[0], 
	   &kappa_p_exp_n[0], &kappa_p_exp_p[0],
	   &const_kappa_r[0], &kappa_r_exp_m[0], 
	   &kappa_r_exp_n[0], &kappa_r_exp_p[0],
	   &const_scattering[0], &scattering_exp_m[0], 
	   &scattering_exp_n[0], &scattering_exp_p[0],
	   &prop_temp_floor[0]);
	
      }

      const Box& reg = grids[i];
      
      Array<Real> PFcoef(nGroups, -1.0); // picket-fence model coefficients
#ifdef MG_SU_OLSON
      PFcoef.set(0, 0.5);
      PFcoef.set(1, 0.5);
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
    kappa_r.FillBoundary();
    if (geom.isAnyPeriodic()) {
      geom.FillPeriodicBoundary(kappa_r, false);
    }
  }
}


void Radiation::gray_accel(MultiFab& Er_new, MultiFab& Er_pi, 
			   MultiFab& kappa_p, MultiFab& kappa_r,
			   MultiFab& etaT, MultiFab& etaY, MultiFab& eta1,
			   MultiFab& thetaT, MultiFab& thetaY, 
			   MultiFab& mugT, MultiFab& mugY, 
			   Tuple<MultiFab, BL_SPACEDIM>& lambda,
			   RadSolve& solver, MGRadBndry& mgbd, 
			   const BoxArray& grids, int level, Real time, 
			   Real delta_t, Real ptc_tau)
{
  const Geometry& geom = parent->Geom(level);
  const Real* dx = parent->Geom(level).CellSize();

  if (nGroups > 1) {
    solver.setHypreMulti(1.0);
  }
  else {
    solver.setHypreMulti(0.0);
  }
  mgbd.setCorrection();

  MultiFab Er_zero(grids, 1, 0);
  Er_zero.setVal(0.0);
  getBndryDataMG_ga(mgbd, Er_zero, level);

  MultiFab spec(grids, nGroups, 1); 
  for (MFIter mfi(spec); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box& bx = grids[i];
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
  spec.FillBoundary();
  if (parent->Geom(level).isAnyPeriodic()) {
    parent->Geom(level).FillPeriodicBoundary(spec, false);
  }

  // set boundary condition
  solver.levelBndry(mgbd,0);

  // A coefficients
  MultiFab acoefs(grids, 1, 0);
  for (MFIter mfi(acoefs); mfi.isValid(); ++mfi) {
#ifdef NEUTRINO
    BL_FORT_PROC_CALL(CA_ACCEL_ACOE_NEUT, ca_accel_acoe_neut)
      (BL_TO_FORTRAN(eta1[mfi]),
       BL_TO_FORTRAN(thetaT[mfi]),
       BL_TO_FORTRAN(thetaY[mfi]),
       BL_TO_FORTRAN(spec[mfi]),
       BL_TO_FORTRAN(kappa_p[mfi]),
       BL_TO_FORTRAN(acoefs[mfi]),
       &delta_t, &ptc_tau);    
#else
    BL_FORT_PROC_CALL(CA_ACCEL_ACOE, ca_accel_acoe)
      (BL_TO_FORTRAN(eta1[mfi]),
       BL_TO_FORTRAN(spec[mfi]),
       BL_TO_FORTRAN(kappa_p[mfi]),
       BL_TO_FORTRAN(acoefs[mfi]),
       &delta_t, &ptc_tau);    
#endif
  }  

  solver.cellCenteredApplyMetrics(level, acoefs);
  solver.setLevelACoeffs(level, acoefs);

  // B & C coefficients
  Tuple<MultiFab, BL_SPACEDIM> bcoefs, ccoefs, bcgrp;
  for (int idim = 0; idim < BL_SPACEDIM; idim++) {
    BoxArray edge_boxes(grids);
    edge_boxes.surroundingNodes(idim);

    bcoefs[idim].define(edge_boxes, 1, 0, Fab_allocate);
    bcoefs[idim].setVal(0.0);

    bcgrp [idim].define(edge_boxes, 1, 0, Fab_allocate);

    if (nGroups > 1) {
      ccoefs[idim].define(edge_boxes, 2, 0, Fab_allocate);
      ccoefs[idim].setVal(0.0);
    }
  }

  for (int igroup = 0; igroup < nGroups; igroup++) {
    for (int idim=0; idim<BL_SPACEDIM; idim++) {
      solver.computeBCoeffs(bcgrp[idim], idim, kappa_r, igroup,
			    lambda[idim], igroup, c, geom);
      // metrics is already in bcgrp

      for (MFIter mfi(spec); mfi.isValid(); ++mfi) {
	int i = mfi.index();
	const Box& reg  = grids[i];
	const Box& bbox = bcoefs[idim][mfi].box();
	const Box& sbox = spec[mfi].box();
	FORT_LBCOEFNA(bcoefs[idim][mfi].dataPtr(),
		      bcgrp[idim][mfi].dataPtr(),
		      dimlist(bbox), dimlist(reg),
		      spec[mfi].dataPtr(igroup), dimlist(sbox),
		      idim);

	if (nGroups > 1) {
	  BL_FORT_PROC_CALL(CA_ACCEL_CCOE, ca_accel_ccoe)
	    (BL_TO_FORTRAN(bcgrp[idim][mfi]),
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
  MultiFab rhs(grids,1,0);
  for (MFIter mfi(spec); mfi.isValid(); ++mfi) {
#ifdef NEUTRINO
    BL_FORT_PROC_CALL(CA_ACCEL_RHS_NEUT, ca_accel_rhs_neut) 
      (BL_TO_FORTRAN(Er_new[mfi]),
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
      (BL_TO_FORTRAN(Er_new[mfi]),
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
  MultiFab accel(grids,1,0);
  accel.setVal(0.0);
  solver.levelSolve(level, accel, 0, rhs, 0.01);

  // update Er_new
  for (MFIter mfi(spec); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box& reg  = grids[i];
    const Box& jbox = Er_new[mfi].box();
    const Box& sbox = spec[mfi].box();
    FORT_LJUPNA(Er_new[mfi].dataPtr(), dimlist(jbox),
		dimlist(reg),
		spec[mfi].dataPtr(), dimlist(sbox),
		accel[mfi].dataPtr(), nGroups);
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
			    const BoxArray& grids, Real delta_t, Real ptc_tau)
{
  for (MFIter mfi(Er_new); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box& bx = grids[i];

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

  BoxArray baf;

  if (level < parent->finestLevel()) {
    baf = parent->boxArray(level+1);
    baf.coarsen(parent->refRatio(level));
  }

  for (MFIter mfi(state); mfi.isValid(); ++mfi) {
    int i = mfi.index();

    const Box& reg  = state.box(i); // interior region

    const Box& sbox = state[mfi].box();
    FArrayBox msk(sbox,1);
    {
      msk.setVal(0.0, 0);
      msk.setVal(1.0, reg, 0); // msk is 1 in interior, 0 in ghost cells

      // Now mask off finer level also:
      if (level < parent->finestLevel()) {
        std::vector< std::pair<int,Box> > isects = baf.intersections(reg);

        for (int ii = 0; ii < isects.size(); ii++) {
          msk.setVal(0.0, isects[ii].second, 0);
        }
      }

      // The mask here is for the derat and dye calculation,
      // to limit it to active cells on each level.
    }

#ifdef NEUTRINO
    BL_FORT_PROC_CALL(CA_STATE_UPDATE_NEUT, ca_state_update_neut)
      (reg.loVect(), reg.hiVect(),
       BL_TO_FORTRAN(state[mfi]),
       BL_TO_FORTRAN(rhoe[mfi]),
       BL_TO_FORTRAN(Ye[mfi]),
       BL_TO_FORTRAN(temp[mfi]),
       BL_TO_FORTRAN(msk),
       &derat, &dT, &dye);
#else
    BL_FORT_PROC_CALL(CA_STATE_UPDATE, ca_state_update)
      (reg.loVect(), reg.hiVect(),
       BL_TO_FORTRAN(state[mfi]),
       BL_TO_FORTRAN(rhoe[mfi]),
       BL_TO_FORTRAN(temp[mfi]),
       BL_TO_FORTRAN(msk),
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
			      int level, const BoxArray& grids, Real delta_t, Real ptc_tau,
			      int it, bool conservative_update)
{
  for (MFIter mfi(rhoe_new); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box& bx = grids[i]; 

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
	BL_FORT_PROC_CALL(CA_COMPUTE_TEMP_GIVEN_REYE,ca_compute_temp_given_reye)
	  (bx.loVect(), bx.hiVect(),
	   BL_TO_FORTRAN(temp_new[mfi]), 
	   BL_TO_FORTRAN(rhoe_new[mfi]),
	   BL_TO_FORTRAN(Ye_new[mfi]),
	   BL_TO_FORTRAN(S_new[mfi]));
      }
      else {
	temp_new[mfi].copy(rhoe_new[mfi]);

	if (do_real_eos > 0) {
	  BL_FORT_PROC_CALL(CA_COMPUTE_TEMP_GIVEN_RHOE, ca_compute_temp_given_rhoe)
	    (bx.loVect(), bx.hiVect(), 
	     BL_TO_FORTRAN(temp_new[mfi]), 
	     BL_TO_FORTRAN(S_new[mfi]));
	}
	else if (do_real_eos == 0) {
	  BL_FORT_PROC_CALL(CA_COMPUTE_TEMP_GIVEN_CV, ca_compute_temp_given_cv)
	    (bx.loVect(), bx.hiVect(), 
	     BL_TO_FORTRAN(temp_new[mfi]), 
	     BL_TO_FORTRAN(S_new[mfi]),
	     &const_c_v[0], &c_v_exp_m[0], &c_v_exp_n[0]);
	}
	else {
	  BoxLib::Error("ERROR Radiation::do_real_eos < 0");
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
      
      temp_new[mfi].copy(rhoe_new[mfi]);

      if (do_real_eos > 0) {
	BL_FORT_PROC_CALL(CA_COMPUTE_TEMP_GIVEN_RHOE, ca_compute_temp_given_rhoe)
	  (bx.loVect(), bx.hiVect(), 
	   BL_TO_FORTRAN(temp_new[mfi]), 
	   BL_TO_FORTRAN(S_new[mfi]));
      }
      else if (do_real_eos == 0) {
	BL_FORT_PROC_CALL(CA_COMPUTE_TEMP_GIVEN_CV, ca_compute_temp_given_cv)
	  (bx.loVect(), bx.hiVect(), 
	   BL_TO_FORTRAN(temp_new[mfi]), 
	   BL_TO_FORTRAN(S_new[mfi]),
	   &const_c_v[0], &c_v_exp_m[0], &c_v_exp_n[0]);
      }
      else {
	BoxLib::Error("ERROR Radiation::do_real_eos < 0");
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

      BL_FORT_PROC_CALL(CA_GET_RHOE,ca_get_rhoe)
      	(bx.loVect(), bx.hiVect(),
      	 BL_TO_FORTRAN(rhoe_new[mfi]),
      	 BL_TO_FORTRAN(temp_new[mfi]), 
      	 BL_TO_FORTRAN(S_new[mfi]));
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
  int ngrow = Sborder.nGrow();

  if (filter_lambda_T) {
    BL_ASSERT(ngrow == 4);
  }

  if (limiter == 0) {

    lamborder.setVal(1./3., ngrow);

  }
  else {

    MultiFab kpr(grids,Radiation::nGroups,ngrow);  

    if (do_multigroup) {
      MGFLD_compute_rosseland(kpr, Sborder); 
    }
    else {
      SGFLD_compute_rosseland(kpr, Sborder); 
    }

    MultiFab Er_wide(grids, nGroups, ngrow+1);
    for (MFIter mfi(Er_wide); mfi.isValid(); ++mfi) {
      int i = mfi.index();
      const Box &bx = grids[i];
      Er_wide[mfi].copy(Erborder[mfi], 0, 0, nGroups);
      Er_wide[mfi].setComplement(-1.0, bx, 0, nGroups);
    }
    
    Er_wide.FillBoundary();
    
    if (parent->Geom(level).isAnyPeriodic()) {
      parent->Geom(level).FillPeriodicBoundary(Er_wide, true);
    }
    
    const Real* dx = parent->Geom(level).CellSize();
    
    for (MFIter mfi(Er_wide); mfi.isValid(); ++mfi) {
      BL_FORT_PROC_CALL(CA_COMPUTE_LAMBORDER, ca_compute_lamborder)
	(BL_TO_FORTRAN(Er_wide[mfi]), 
	 BL_TO_FORTRAN(kpr[mfi]),
	 BL_TO_FORTRAN(lamborder[mfi]), 
	 dx, &ngrow, &limiter, &filter_lambda_T, &filter_lambda_S);
    }

    if (filter_lambda_T) {
      lamborder.FillBoundary();

      if (parent->Geom(level).isAnyPeriodic()) {
	parent->Geom(level).FillPeriodicBoundary(lamborder, true);
      }
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

    FArrayBox kappa_r(box, nGroups);

    if (do_multigroup) {
      MGFLD_compute_rosseland(kappa_r, state);
    }
    else {
      SGFLD_compute_rosseland(kappa_r, state);
    }

    BL_FORT_PROC_CALL(CA_EST_GPR2, ca_est_gpr2)
      (BL_TO_FORTRAN(kappa_r),
       BL_TO_FORTRAN(Er),
       BL_TO_FORTRAN(gPr), dx, &limiter, &comoving);
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
      BL_FORT_PROC_CALL(CA_COMPUTE_ROSSELAND, ca_compute_rosseland)
	  (kbox.loVect(), kbox.hiVect(),
	   BL_TO_FORTRAN(kappa_r), BL_TO_FORTRAN(state));
    }
    else if (const_kappa_r[0] < 0.0) {
      BL_FORT_PROC_CALL(CA_COMPUTE_POWERLAW_KAPPA_S, ca_compute_powerlaw_kappa_s)
	  (kbox.loVect(), kbox.hiVect(),
	   BL_TO_FORTRAN(kappa_r), BL_TO_FORTRAN(state),
	   &const_kappa_p[0], &kappa_p_exp_m[0], &kappa_p_exp_n[0], &kappa_p_exp_p[0], 
	   &const_scattering[0], &scattering_exp_m[0], &scattering_exp_n[0], &scattering_exp_p[0], 
	   &prop_temp_floor[0], &kappa_r_floor);	 
    }
    else {
      BL_FORT_PROC_CALL(CA_COMPUTE_POWERLAW_KAPPA, ca_compute_powerlaw_kappa)
	  (kbox.loVect(), kbox.hiVect(),
	   BL_TO_FORTRAN(kappa_r), BL_TO_FORTRAN(state),
	   &const_kappa_r[0], &kappa_r_exp_m[0], &kappa_r_exp_n[0], &kappa_r_exp_p[0], 
	   &prop_temp_floor[0], &kappa_r_floor);	       
    }
#ifdef NEUTRINO
  }
#endif
}


void Radiation::MGFLD_compute_rosseland(MultiFab& kappa_r, const MultiFab& state)
{
  BL_PROFILE("Radiation::MGFLD_compute_rosseland (MultiFab)");

  for (MFIter mfi(kappa_r); mfi.isValid(); ++mfi) {
    const Box& kbox = kappa_r[mfi].box();
#ifdef NEUTRINO
    if (radiation_type == Neutrino) {
      BL_FORT_PROC_CALL(CA_COMPUTE_ROSSELAND_NEUT, ca_compute_rosseland_neut)
	  (kbox.loVect(), kbox.hiVect(),
	   BL_TO_FORTRAN(kappa_r[mfi]), BL_TO_FORTRAN(state[mfi]));
    }
    else {
#endif
      if (use_opacity_table_module) {
	BL_FORT_PROC_CALL(CA_COMPUTE_ROSSELAND, ca_compute_rosseland)
	    (kbox.loVect(), kbox.hiVect(),
	     BL_TO_FORTRAN(kappa_r[mfi]), BL_TO_FORTRAN(state[mfi]));
      }
      else if (const_kappa_r[0] < 0.0) {
	BL_FORT_PROC_CALL(CA_COMPUTE_POWERLAW_KAPPA_S, ca_compute_powerlaw_kappa_s)
	    (kbox.loVect(), kbox.hiVect(),
	     BL_TO_FORTRAN(kappa_r[mfi]), BL_TO_FORTRAN(state[mfi]),
	     &const_kappa_p[0], &kappa_p_exp_m[0], &kappa_p_exp_n[0], &kappa_p_exp_p[0], 
	     &const_scattering[0], &scattering_exp_m[0], &scattering_exp_n[0], &scattering_exp_p[0], 
	     &prop_temp_floor[0], &kappa_r_floor);	 
      }
      else {
	BL_FORT_PROC_CALL(CA_COMPUTE_POWERLAW_KAPPA, ca_compute_powerlaw_kappa)
	    (kbox.loVect(), kbox.hiVect(),
	     BL_TO_FORTRAN(kappa_r[mfi]), BL_TO_FORTRAN(state[mfi]),
	     &const_kappa_r[0], &kappa_r_exp_m[0], &kappa_r_exp_n[0], &kappa_r_exp_p[0], 
	     &prop_temp_floor[0], &kappa_r_floor);	 
      }
#ifdef NEUTRINO
    }
#endif
  }
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

  for (MFIter mfi(rhoe_new); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box& bx = grids[i]; 

#ifdef NEUTRINO    
    BL_FORT_PROC_CALL(CA_COMPUTE_REYE_GIVEN_TY,ca_compute_reye_given_ty)
      (bx.loVect(), bx.hiVect(),
       BL_TO_FORTRAN(rhoe_new[mfi]),
       BL_TO_FORTRAN(rhoYe_new[mfi]),
       BL_TO_FORTRAN(temp_new[mfi]), 
       BL_TO_FORTRAN(Ye_new[mfi]),
       BL_TO_FORTRAN(S_new[mfi]));
#else
    if (do_real_eos > 0) {
      BL_FORT_PROC_CALL(CA_GET_RHOE,ca_get_rhoe)
      	(bx.loVect(), bx.hiVect(),
      	 BL_TO_FORTRAN(rhoe_new[mfi]),
      	 BL_TO_FORTRAN(temp_new[mfi]), 
      	 BL_TO_FORTRAN(S_new[mfi]));
    }
    else {
      BoxLib::Abort("do_real_eos == 0 not supported in bisect_matter");
    }
#endif
  }
}


void Radiation::rhstoEr(MultiFab& rhs, const BoxArray& grids, Real dt, int level)
{
    Array<Real> r, s;

    for (MFIter ri(rhs); ri.isValid(); ++ri) 
    {
	int i = ri.index();
	const Box &reg = grids[i];

	const int I = (BL_SPACEDIM >= 2) ? 1 : 0;
	if (CoordSys::IsCartesian()) 
	{
	    r.resize(reg.length(0), 1);
	    s.resize(reg.length(I), 1);
	}
	else if (CoordSys::IsRZ()) 
	{
	    parent->Geom(level).GetCellLoc(r, reg, 0);
	    s.resize(reg.length(I), 1);
	}
	else 
	{
	    parent->Geom(level).GetCellLoc(r, reg, 0);
	    parent->Geom(level).GetCellLoc(s, reg, I);
	    const Real *dx = parent->Geom(level).CellSize();
	    FORT_SPHC(r.dataPtr(), s.dataPtr(), dimlist(reg), dx);
	}

	BL_FORT_PROC_CALL(CA_RHSTOER, ca_rhstoer)
	    (BL_TO_FORTRAN(rhs[ri]), r.dataPtr(), &dt);
    }
}
