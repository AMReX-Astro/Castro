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
				     const MultiFab& etaTz,
				     const MultiFab& temp_new,
				     Real delta_t)
{
  BL_PROFILE("Radiation::check_convergence_er (MGFLD)");

  relative = 0.0;
  absolute = 0.0;
  err_er = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(max:relative, absolute, err_er)
#endif
  for (MFIter mfi(Er_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
      ca_check_conv_er
          (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
           BL_TO_FORTRAN_ANYD(Er_new[mfi]),
           BL_TO_FORTRAN_ANYD(Er_pi[mfi]),
           BL_TO_FORTRAN_ANYD(kappa_p[mfi]),
           BL_TO_FORTRAN_ANYD(etaTz[mfi]),
           BL_TO_FORTRAN_ANYD(temp_new[mfi]),
           AMREX_MFITER_REDUCE_MAX(&relative),
           AMREX_MFITER_REDUCE_MAX(&absolute),
           AMREX_MFITER_REDUCE_MAX(&err_er),
           delta_t);
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
				       const MultiFab& rho,
				       const MultiFab& kappa_p, const MultiFab& jg,
				       const MultiFab& dedT,
				       Real& rel_rhoe, Real& abs_rhoe, 
				       Real& rel_FT,   Real& abs_FT, 
				       Real& rel_T,    Real& abs_T, 
				       Real delta_t)
{
  rel_rhoe = 0.0;
  rel_FT   = 0.0;
  rel_T    = 0.0;
  abs_rhoe = 0.0;
  abs_FT   = 0.0;
  abs_T    = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(max:rel_rhoe, abs_rhoe, rel_FT, abs_FT, rel_T, abs_T)
#endif
  for (MFIter mfi(rhoe_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

      ca_check_conv
          (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
           BL_TO_FORTRAN_ANYD(rhoe_new[mfi]),
           BL_TO_FORTRAN_ANYD(rhoe_star[mfi]),
           BL_TO_FORTRAN_ANYD(rhoe_step[mfi]),
           BL_TO_FORTRAN_ANYD(Er_new[mfi]),
           BL_TO_FORTRAN_ANYD(temp_new[mfi]),
           BL_TO_FORTRAN_ANYD(temp_star[mfi]),
           BL_TO_FORTRAN_ANYD(rho[mfi]),
           BL_TO_FORTRAN_ANYD(kappa_p[mfi]),
           BL_TO_FORTRAN_ANYD(jg[mfi]),
           BL_TO_FORTRAN_ANYD(dedT[mfi]),
           AMREX_MFITER_REDUCE_MAX(&rel_rhoe),
           AMREX_MFITER_REDUCE_MAX(&abs_rhoe),
           AMREX_MFITER_REDUCE_MAX(&rel_FT),
           AMREX_MFITER_REDUCE_MAX(&abs_FT),
           AMREX_MFITER_REDUCE_MAX(&rel_T),
           AMREX_MFITER_REDUCE_MAX(&abs_T),
           delta_t);
  }

  int ndata = 6;
  Real data[6] = {rel_rhoe, abs_rhoe, rel_FT, abs_FT, rel_T, abs_T};

  ParallelDescriptor::ReduceRealMax(data, ndata);

  rel_rhoe = data[0]; 
  abs_rhoe = data[1];
  rel_FT   = data[2];
  abs_FT   = data[3];
  rel_T    = data[4];
  abs_T    = data[5];
}


void Radiation::compute_coupling(MultiFab& coupT,
				 const MultiFab& kpp, const MultiFab& Eg,
				 const MultiFab& jg)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(kpp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
	ca_compute_coupt
	    (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
	     BL_TO_FORTRAN_ANYD(coupT[mfi]),
	     BL_TO_FORTRAN_ANYD(kpp[mfi]),
	     BL_TO_FORTRAN_ANYD(Eg[mfi]),    
	     BL_TO_FORTRAN_ANYD(jg[mfi]));
    }
}


void Radiation::compute_etat(MultiFab& etaT, MultiFab& etaTz, 
                             MultiFab& eta1, MultiFab& djdT,
                             const MultiFab& dkdT, const MultiFab& dedT,
                             const MultiFab& Er_star, const MultiFab& rho,
                             Real delta_t, Real ptc_tau)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(rho, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
	ca_compute_etat
	    (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
             BL_TO_FORTRAN_ANYD(etaT[mfi]),
             BL_TO_FORTRAN_ANYD(etaTz[mfi]),
             BL_TO_FORTRAN_ANYD(eta1[mfi]),
             BL_TO_FORTRAN_ANYD(djdT[mfi]),
             BL_TO_FORTRAN_ANYD(dkdT[mfi]),
             BL_TO_FORTRAN_ANYD(dedT[mfi]),
             BL_TO_FORTRAN_ANYD(Er_star[mfi]),
             BL_TO_FORTRAN_ANYD(rho[mfi]),
             delta_t, ptc_tau);
    }
}


void Radiation::eos_opacity_emissivity(const MultiFab& S_new, 
				       const MultiFab& temp_new,
				       const MultiFab& temp_star,
				       MultiFab& kappa_p, MultiFab& kappa_r, MultiFab& jg, 
				       MultiFab& djdT, MultiFab& dkdT, MultiFab& dedT,
				       int level, int it, int ngrow)
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
  for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& box = mfi.tilebox();

#pragma gpu box(box)
      ca_compute_c_v
          (AMREX_INT_ANYD(box.loVect()), AMREX_INT_ANYD(box.hiVect()),
           BL_TO_FORTRAN_ANYD(dedT[mfi]),
           BL_TO_FORTRAN_ANYD(temp_new[mfi]),
           BL_TO_FORTRAN_ANYD(S_new[mfi]));
  }

  if (dedT_fac > 1.0) {
    dedT.mult(dedT_fac);
  }

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.growntilebox(ngrow);

#pragma gpu box(bx) sync
      ca_opacs
          (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
           BL_TO_FORTRAN_ANYD(S_new[mfi]),
           BL_TO_FORTRAN_ANYD(temp_new[mfi]),
           BL_TO_FORTRAN_ANYD(temp_star[mfi]),
           BL_TO_FORTRAN_ANYD(kappa_p[mfi]),
           BL_TO_FORTRAN_ANYD(kappa_r[mfi]),
           BL_TO_FORTRAN_ANYD(dkdT[mfi]),
           use_dkdT, star_is_valid, lag_opac);

      const Box& reg = mfi.tilebox();
      
      Vector<Real> PFcoef(nGroups, -1.0); // picket-fence model coefficients
#ifdef MG_SU_OLSON
      PFcoef[0] = 0.5;
      PFcoef[1] = 0.5;
#endif
      ca_compute_emissivity
          (AMREX_INT_ANYD(reg.loVect()), AMREX_INT_ANYD(reg.hiVect()),
           BL_TO_FORTRAN_ANYD(jg[mfi]),  
           BL_TO_FORTRAN_ANYD(djdT[mfi]),  
           BL_TO_FORTRAN_ANYD(temp_new[mfi]),
           BL_TO_FORTRAN_ANYD(kappa_p[mfi]),
           BL_TO_FORTRAN_ANYD(dkdT[mfi]),
           PFcoef.dataPtr(), 
           use_WiensLaw, integrate_Planck, Tf_Wien);
  }    

  if (ngrow == 0 && !lag_opac) {
      kappa_r.FillBoundary(geom.periodicity());
  }
}


void Radiation::gray_accel(MultiFab& Er_new, MultiFab& Er_pi, 
			   MultiFab& kappa_p, MultiFab& kappa_r,
			   MultiFab& etaT, MultiFab& eta1,
			   MultiFab& mugT,
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

    BL_FORT_PROC_CALL(CA_ACCEL_SPEC, ca_accel_spec) 
      (bx.loVect(), bx.hiVect(),
       BL_TO_FORTRAN(kappa_p[mfi]),
       BL_TO_FORTRAN(mugT[mfi]),
       BL_TO_FORTRAN(spec[mfi]),
       &delta_t, &ptc_tau);
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

      BL_FORT_PROC_CALL(CA_ACCEL_ACOE, ca_accel_acoe)
	  (bx.loVect(), bx.hiVect(),
	   BL_TO_FORTRAN(eta1[mfi]),
	   BL_TO_FORTRAN(spec[mfi]),
	   BL_TO_FORTRAN(kappa_p[mfi]),
	   BL_TO_FORTRAN(acoefs[mfi]),
	   &delta_t, &ptc_tau);    
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

      BL_FORT_PROC_CALL(CA_ACCEL_RHS, ca_accel_rhs) 
	  (bx.loVect(), bx.hiVect(),
	   BL_TO_FORTRAN(Er_new[mfi]),
	   BL_TO_FORTRAN(Er_pi[mfi]),
	   BL_TO_FORTRAN(kappa_p[mfi]),
	   BL_TO_FORTRAN(etaT[mfi]),
	   BL_TO_FORTRAN(rhs[mfi]),
	   &delta_t);
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
			    const MultiFab& etaT,
			    const MultiFab& mugT,
			    Real delta_t, Real ptc_tau)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(Er_new,true); mfi.isValid(); ++mfi) {
	const Box& bx = mfi.tilebox();

	BL_FORT_PROC_CALL(CA_LOCAL_ACCEL, ca_local_accel)
	    (bx.loVect(), bx.hiVect(),
	     BL_TO_FORTRAN(Er_new[mfi]),
	     BL_TO_FORTRAN(Er_pi[mfi]),
	     BL_TO_FORTRAN(kappa_p[mfi]),
	     BL_TO_FORTRAN(etaT[mfi]),
	     BL_TO_FORTRAN(mugT[mfi]),
	     &delta_t, &ptc_tau);
    }
}


void Radiation::state_energy_update(MultiFab& state, const MultiFab& rhoe, 
				    const MultiFab& temp, 
				    const BoxArray& grids,
				    Real& derat, Real& dT, int level)
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

      BL_FORT_PROC_CALL(CA_STATE_UPDATE, ca_state_update)
	  (reg.loVect(), reg.hiVect(),
	   BL_TO_FORTRAN(state[mfi]),
	   BL_TO_FORTRAN(rhoe[mfi]),
	   BL_TO_FORTRAN(temp[mfi]),
	   BL_TO_FORTRAN(msk[mfi]),
	   &derat, &dT);
  }

  ParallelDescriptor::ReduceRealMax(derat);
  ParallelDescriptor::ReduceRealMax(dT);
}


void Radiation::update_matter(MultiFab& rhoe_new, MultiFab& temp_new, 
			      const MultiFab& Er_new, const MultiFab& Er_pi,
			      const MultiFab& rhoe_star,
			      const MultiFab& rhoe_step,
			      const MultiFab& etaT, const MultiFab& etaTz,
			      const MultiFab& eta1,
			      const MultiFab& coupT,
			      const MultiFab& kappa_p, const MultiFab& jg, 
			      const MultiFab& mugT,
			      const MultiFab& S_new, 
			      int level, Real delta_t, Real ptc_tau,
			      int it, bool conservative_update)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(rhoe_new,true); mfi.isValid(); ++mfi) {
	const Box& bx = mfi.tilebox(); 

	if (conservative_update) {
#pragma gpu box(bx) sync
	    ca_update_matter
		(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
		 BL_TO_FORTRAN_ANYD(rhoe_new[mfi]),
		 BL_TO_FORTRAN_ANYD(Er_new[mfi]),
		 BL_TO_FORTRAN_ANYD(Er_pi[mfi]),
		 BL_TO_FORTRAN_ANYD(rhoe_star[mfi]),
		 BL_TO_FORTRAN_ANYD(rhoe_step[mfi]),
		 BL_TO_FORTRAN_ANYD(eta1[mfi]),
		 BL_TO_FORTRAN_ANYD(coupT[mfi]),
		 BL_TO_FORTRAN_ANYD(kappa_p[mfi]),
		 delta_t, ptc_tau);

	    temp_new[mfi].copy(rhoe_new[mfi],bx);

#pragma gpu box(bx) sync
            ca_compute_temp_given_rhoe
                (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()), 
                 BL_TO_FORTRAN_ANYD(temp_new[mfi]), 
                 BL_TO_FORTRAN_ANYD(S_new[mfi]),
                 0);
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

#pragma gpu box(kbox) sync
  ca_compute_rosseland(AMREX_INT_ANYD(kbox.loVect()), AMREX_INT_ANYD(kbox.hiVect()),
                       BL_TO_FORTRAN_ANYD(kappa_r),
                       BL_TO_FORTRAN_ANYD(state),
                       0, nGroups-1, nGroups);

}


void Radiation::MGFLD_compute_rosseland(MultiFab& kappa_r, const MultiFab& state)
{
    BL_PROFILE("Radiation::MGFLD_compute_rosseland (MultiFab)");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(kappa_r, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
	const Box& bx = mfi.growntilebox();

#pragma gpu box(bx) sync
        ca_compute_rosseland(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                             BL_TO_FORTRAN_ANYD(kappa_r[mfi]),
                             BL_TO_FORTRAN_ANYD(state[mfi]),
                             0, nGroups-1, nGroups);
    }
}

void Radiation::MGFLD_compute_scattering(FArrayBox& kappa_s, const FArrayBox& state)
{
    BL_PROFILE("Radiation::MGFLD_compute_scattering");

    const Box& kbox = kappa_s.box();

#pragma gpu box(kbox) sync
    ca_compute_scattering
        (AMREX_INT_ANYD(kbox.loVect()), AMREX_INT_ANYD(kbox.hiVect()),
         BL_TO_FORTRAN_ANYD(kappa_s), 
         BL_TO_FORTRAN_ANYD(state));
}

void Radiation::bisect_matter(MultiFab& rhoe_new, MultiFab& temp_new, 
			      const MultiFab& rhoe_star, const MultiFab& temp_star, 
			      const MultiFab& S_new, const BoxArray& grids, int level)
{
  temp_new.plus(temp_star, 0, 1, 0);
  temp_new.mult(0.5, 0);

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(rhoe_new,true); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

#pragma gpu box(bx) sync
      ca_get_rhoe
          (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
           BL_TO_FORTRAN_ANYD(rhoe_new[mfi]),
           BL_TO_FORTRAN_ANYD(temp_new[mfi]), 
           BL_TO_FORTRAN_ANYD(S_new[mfi]));
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

        castro->computeTemp(S_new, castro->state[State_Type].curTime(), S_new.nGrow());
    }
}
