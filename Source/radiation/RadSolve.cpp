
#include <AMReX_ParmParse.H>
#include <AMReX_AmrLevel.H>

#include <AMReX_LO_BCTYPES.H>
//#include <CompSolver.H>

#include "RadSolve.H"
#include "Radiation.H"  // for access to static physical constants only

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "RAD_F.H"

#include "HABEC_F.H"    // only for nonsymmetric flux; may be changed?

using namespace amrex;

#include <radsolve_defaults.H>

RadSolve::RadSolve (Amr* Parent, int level, const BoxArray& grids, const DistributionMapping& dmap)
    : parent(Parent)
{
    read_params();

    if (level_solver_flag < 100) {
        hd.reset(new HypreABec(grids, dmap, parent->Geom(level), level_solver_flag));
    }
    else {
        if (use_hypre_nonsymmetric_terms == 0) {
            hm.reset(new HypreMultiABec(level, level, level_solver_flag));
        }
        else {
            hm.reset(new HypreExtMultiABec(level, level, level_solver_flag));
            HypreExtMultiABec *hem = (HypreExtMultiABec*)hm.get();
            cMulti  = hem->cMultiplier();
            d1Multi = hem->d1Multiplier();
            d2Multi = hem->d2Multiplier();
        }
        hm->addLevel(level, parent->Geom(level), grids, dmap,
                     IntVect::TheUnitVector());
        hm->buildMatrixStructure();
    }
}

void
RadSolve::read_params ()
{
    ParmParse pp("radsolve");

    // Override some defaults manually.

    if (BL_SPACEDIM == 1) {
        // pfmg will not work in 1D
        level_solver_flag = 0;
    }

    if (Radiation::SolverType == Radiation::SGFLDSolver
        && Radiation::Er_Lorentz_term) { 
        use_hypre_nonsymmetric_terms = 1;
    }

    if (Radiation::SolverType == Radiation::MGFLDSolver && 
        Radiation::accelerate == 2 && Radiation::nGroups > 1) {
        use_hypre_nonsymmetric_terms = 1;
    }

    if (Radiation::SolverType == Radiation::SGFLDSolver ||
        Radiation::SolverType == Radiation::MGFLDSolver) {
        abstol = 0.0;
    }

#include "radsolve_queries.H"

    // Check for unsupported options.

    if (BL_SPACEDIM == 1) {
        if (level_solver_flag == 1) {
            amrex::Error("radsolve.level_solver_flag = 1 is not supported in 1D");
        }
    }

    if (Radiation::SolverType == Radiation::SGFLDSolver
        && Radiation::Er_Lorentz_term) { 

        if (level_solver_flag < 100) {
            amrex::Error("To do Lorentz term implicitly level_solver_flag must be >= 100.");
        }
    }

    if (Radiation::SolverType == Radiation::MGFLDSolver && 
        Radiation::accelerate == 2 && Radiation::nGroups > 1) {

        if (level_solver_flag < 100) {
            amrex::Error("When accelerate is 2, level_solver_flag must be >= 100.");
        }
    }

}

void RadSolve::levelInit(int level)
{
  BL_PROFILE("RadSolve::levelInit");
}

void RadSolve::levelBndry(RadBndry& bd)
{
  BL_PROFILE("RadSolve::levelBndry");

  if (hd) {
    hd->setBndry(bd);
  }
  else if (hm) {
    hm->setBndry(hm->crseLevel(), bd);
  }
}

// update multigroup version
void RadSolve::levelBndry(MGRadBndry& mgbd, const int comp)
{
  BL_PROFILE("RadSolve::levelBndryMG (updated)");

  if (hd) {
    hd->setBndry(mgbd, comp);
  }
  else if (hm) {
    hm->setBndry(hm->crseLevel(), mgbd, comp);
  }
}

void RadSolve::cellCenteredApplyMetrics(int level, MultiFab& cc)
{
  BL_PROFILE("RadSolve::cellCenteredApplyMetrics");
  BL_ASSERT(cc.nGrow() == 0);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
      Vector<Real> r, s;

      for (MFIter mfi(cc,true); mfi.isValid(); ++mfi) 
      {
	  const Box &reg = mfi.tilebox();

	  getCellCenterMetric(parent->Geom(level), reg, r, s);
	  
	  multrs(BL_TO_FORTRAN(cc[mfi]), 
		 ARLIM(reg.loVect()), ARLIM(reg.hiVect()), 
		 r.dataPtr(), s.dataPtr());
      }
  }
}

void RadSolve::setLevelACoeffs(int level, const MultiFab& acoefs)
{
    if (hd) {
	hd->aCoefficients(acoefs);
    }
    else if (hm) {
	hm->aCoefficients(level, acoefs);
    }
}

void RadSolve::setLevelBCoeffs(int level, const MultiFab& bcoefs, int dir)
{
    if (hd) {
	hd->bCoefficients(bcoefs, dir);
    }
    else if (hm) {
	hm->bCoefficients(level, bcoefs, dir);
    }
}

void RadSolve::setLevelCCoeffs(int level, const MultiFab& ccoefs, int dir)
{
  if (hm) {
    HypreExtMultiABec *hem = dynamic_cast<HypreExtMultiABec*>(hm.get());
    if (hem) {
      hem->cCoefficients(level, ccoefs, dir);
    }
  }
}

void RadSolve::levelACoeffs(int level,
                            MultiFab& fkp, MultiFab& eta, MultiFab& etainv,
                            Real c, Real delta_t, Real theta)
{
  BL_PROFILE("RadSolve::levelACoeffs");
  const BoxArray& grids = parent->boxArray(level);
  const DistributionMapping& dmap = parent->DistributionMap(level);
  const Geometry& geom = parent->Geom(level);
  const Real* dx       = geom.CellSize();

  // Allocate space for ABecLapacian acoeffs, fill with values

  int Ncomp  = 1;
  int Nghost = 0;

  MultiFab acoefs(grids, dmap, Ncomp, Nghost);

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(fkp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box &bx = mfi.tilebox();

#pragma gpu box(bx)
      lacoef(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
             BL_TO_FORTRAN_ANYD(acoefs[mfi]),
             BL_TO_FORTRAN_ANYD(fkp[mfi]),
             BL_TO_FORTRAN_ANYD(eta[mfi]),
             BL_TO_FORTRAN_ANYD(etainv[mfi]),
             AMREX_REAL_ANYD(dx),
             c, delta_t, theta);
  }

  if (hd) {
    hd->aCoefficients(acoefs);
  }
  else if (hm) {
    hm->aCoefficients(level, acoefs);
  }
}

void RadSolve::levelSPas(int level, Array<MultiFab, BL_SPACEDIM>& lambda, int igroup, 
			 int lo_bc[3], int hi_bc[3])
{
  const BoxArray& grids = parent->boxArray(level);
  const DistributionMapping& dmap = parent->DistributionMap(level);
  const Geometry& geom = parent->Geom(level);
  const Box& domainBox = geom.Domain();

  MultiFab spa(grids, dmap, 1, 0);
#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(spa,true); mfi.isValid(); ++mfi) {
      const Box& reg  = mfi.tilebox();
    
      spa[mfi].setVal(1.e210,reg,0);
    
      bool nexttoboundary=false;
      for (int idim=0; idim<BL_SPACEDIM; idim++) {
	  if (lo_bc[idim] == LO_SANCHEZ_POMRANING &&
	      reg.smallEnd(idim) == domainBox.smallEnd(idim)) {
	      nexttoboundary=true;
	      break;
	  }
	  if (hi_bc[idim] == LO_SANCHEZ_POMRANING &&
	      reg.bigEnd(idim) == domainBox.bigEnd(idim)) {
	      nexttoboundary=true;
	      break;
	  }
      }
    
      if (nexttoboundary) {
	  ca_spalpha(reg.loVect(), reg.hiVect(),
		     BL_TO_FORTRAN(spa[mfi]),
		     D_DECL(BL_TO_FORTRAN(lambda[0][mfi]),
			    BL_TO_FORTRAN(lambda[1][mfi]),
			    BL_TO_FORTRAN(lambda[2][mfi])),
		     &igroup);
      }
  }

  if (hm) {
    hm->SPalpha(level, spa);
  }
  else if (hd) {
    hd->SPalpha(spa);
  }
  else {
    amrex::Abort("Should not be in RadSolve::levelSPas");    
  }
}

void RadSolve::levelBCoeffs(int level,
                            Array<MultiFab, BL_SPACEDIM>& lambda,
                            MultiFab& kappa_r, int kcomp,
                            Real c, int lamcomp)
{
  BL_PROFILE("RadSolve::levelBCoeffs");
  BL_ASSERT(kappa_r.nGrow() == 1);

  const Geometry& geom = parent->Geom(level);
  const Real* dx       = geom.CellSize();

  for (int idim = 0; idim < BL_SPACEDIM; ++idim) {

    MultiFab bcoefs(lambda[idim].boxArray(), lambda[idim].DistributionMap(), 1, 0);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(lambda[idim], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
        bclim(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
              BL_TO_FORTRAN_ANYD(bcoefs[mfi]), 
              BL_TO_FORTRAN_N_ANYD(lambda[idim][mfi], lamcomp),
              idim, 
              BL_TO_FORTRAN_N_ANYD(kappa_r[mfi], kcomp), 
              c, AMREX_REAL_ANYD(dx));
    }

    if (hd) {
        hd->bCoefficients(bcoefs, idim);
    }
    else if (hm) {
      hm->bCoefficients(level, bcoefs, idim);
    }
  } // -->> over dimension
}

void RadSolve::levelDCoeffs(int level, Array<MultiFab, BL_SPACEDIM>& lambda,
			    MultiFab& vel, MultiFab& dcf)
{
    BL_PROFILE("RadSolve::levelDCoeffs");
    const Castro *castro = dynamic_cast<Castro*>(&parent->getLevel(level));
    const DistributionMapping& dm = castro->DistributionMap();

    for (int idim=0; idim<BL_SPACEDIM; idim++) {

	MultiFab dcoefs(castro->getEdgeBoxArray(idim), dm, 1, 0);

#ifdef _OPENMP
#pragma omp parallel
#endif
	{
	    Vector<Real> r, s;
	    
	    for (MFIter mfi(dcoefs,true); mfi.isValid(); ++mfi) {
		const Box& ndbx = mfi.tilebox();

		getEdgeMetric(idim, parent->Geom(level), ndbx, r, s);

		ca_compute_dcoefs(ndbx.loVect(), ndbx.hiVect(),
				  BL_TO_FORTRAN(dcoefs[mfi]), BL_TO_FORTRAN(lambda[idim][mfi]),
				  BL_TO_FORTRAN(vel[mfi]), BL_TO_FORTRAN(dcf[mfi]), 
				  r.dataPtr(), &idim);
	    }
	}

	HypreExtMultiABec *hem = (HypreExtMultiABec*)hm.get();
	hem->d2Coefficients(level, dcoefs, idim);
	hem->d2Multiplier() = 1.0;
    }
}

void RadSolve::levelRhs(int level, MultiFab& rhs,
                        MultiFab& temp,
                        MultiFab& fkp, MultiFab& eta, MultiFab& etainv,
                        MultiFab& rhoem, MultiFab& rhoes,
                        MultiFab& dflux_old, MultiFab& Er_old, MultiFab& Edot,
                        Real delta_t, Real sigma, Real c, Real theta,
                        FluxRegister* fine_corr, Real scale,
                        int igroup, Real nu, Real dnu)
{
  BL_PROFILE("RadSolve::levelRhs");
  BL_ASSERT(rhs.nGrow() == 0);

  const Geometry& geom = parent->Geom(level);
  const Real* dx       = geom.CellSize();

  rhs.setVal(0.0);
  if (fine_corr) {
    // This works trivially for a multilevel solve since the finer level is
    // present in the rhs to overwrite the junk produced under it.
    // In a single-level version we have to be sure that fine_corr
    // has been cleaned up using ClearInternalBorders.

    // Hack:  For the single group case igroup defaults to -1, which is
    // significant later in this routine.  So we have to construct the
    // correct component number here:

    int igrouptmp = (igroup < 0) ? 0 : igroup;

    fine_corr->Reflux(rhs, scale, igrouptmp, 0, 1, parent->Geom(level));
  }

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(rhs, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
      lrhs(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
           BL_TO_FORTRAN_ANYD(rhs[mfi]), 
           BL_TO_FORTRAN_ANYD(temp[mfi]),
           BL_TO_FORTRAN_ANYD(fkp[mfi]),
           BL_TO_FORTRAN_ANYD(eta[mfi]),
           BL_TO_FORTRAN_ANYD(etainv[mfi]),
           BL_TO_FORTRAN_ANYD(rhoem[mfi]),
           BL_TO_FORTRAN_ANYD(rhoes[mfi]),
           BL_TO_FORTRAN_ANYD(dflux_old[mfi]),
           BL_TO_FORTRAN_N_ANYD(Er_old[mfi], 0), 
           BL_TO_FORTRAN_ANYD(Edot[mfi]),
           delta_t, AMREX_REAL_ANYD(dx), sigma, c, theta);
  }
}


void RadSolve::levelSolve(int level,
                          MultiFab& Er, int igroup, MultiFab& rhs,
                          Real sync_absres_factor)
{
  BL_PROFILE("RadSolve::levelSolve");

  // Set coeffs, build solver, solve
  if (hd) {
    hd->setScalars(alpha, beta);
  }
  else if (hm) {
    hm->setScalars(alpha, beta);
  }

  if (hd) {
    hd->setupSolver(reltol, abstol, maxiter);
    hd->solve(Er, igroup, rhs, Inhomogeneous_BC);
    Real res = hd->getAbsoluteResidual();
    if (verbose >= 2 && ParallelDescriptor::IOProcessor()) {
      int oldprec = std::cout.precision(20);
      std::cout << "Absolute residual = " << res << std::endl;
      std::cout.precision(oldprec);
    }
    res *= sync_absres_factor;
    hd->clearSolver();
  }
  else if (hm) {
    hm->loadMatrix();
    hm->finalizeMatrix();
    hm->loadLevelVectors(level, Er, igroup, rhs, Inhomogeneous_BC);
    hm->finalizeVectors();
    hm->setupSolver(reltol, abstol, maxiter);
    hm->solve();
    hm->getSolution(level, Er, igroup);
    Real res = hm->getAbsoluteResidual();
    if (verbose >= 2 && ParallelDescriptor::IOProcessor()) {
      int oldprec = std::cout.precision(20);
      std::cout << "Absolute residual = " << res << std::endl;
      std::cout.precision(oldprec);
    }
    res *= sync_absres_factor;
    hm->clearSolver();
  }
}

void RadSolve::levelFluxFaceToCenter(int level, const Array<MultiFab, BL_SPACEDIM>& Flux,
				     MultiFab& flx, int iflx)
{
    int nflx = flx.nComp();
    
    const Geometry& geom = parent->Geom(level);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	Vector<Real> r, s;
    
	for (int idim = 0; idim < BL_SPACEDIM; idim++) {
	    for (MFIter mfi(flx,true); mfi.isValid(); ++mfi) 
	    {
		const Box &ccbx  = mfi.tilebox();
		const Box &ndbx = amrex::surroundingNodes(ccbx, idim);

		getEdgeMetric(idim, geom, ndbx, r, s);

		int rlo = ndbx.smallEnd(0);
		int rhi = rlo + r.size() - 1;

		ca_flux_face2center(ccbx.loVect(), ccbx.hiVect(),
				    BL_TO_FORTRAN(flx[mfi]),
				    BL_TO_FORTRAN(Flux[idim][mfi]),
				    r.dataPtr(), &rlo, &rhi, 
				    &nflx, &idim, &iflx);
	    }
	}
    }
}

void RadSolve::levelFlux(int level,
                         Array<MultiFab, BL_SPACEDIM>& Flux,
                         MultiFab& Er, int igroup)
{
  BL_PROFILE("RadSolve::levelFlux");
  const BoxArray& grids = parent->boxArray(level);
  const DistributionMapping& dmap = parent->DistributionMap(level);

  // grow a larger MultiFab to hold Er so we can difference across faces
  MultiFab Erborder(grids, dmap, 1, 1);
  Erborder.setVal(0.0);
  MultiFab::Copy(Erborder, Er, igroup, 0, 1, 0);

  Erborder.FillBoundary(parent->Geom(level).periodicity()); // zeroes left in off-level boundaries

  const Real* dx = parent->Geom(level).CellSize();

  for (int n = 0; n < BL_SPACEDIM; n++) {

      const MultiFab *bp;

      if (hd) {
          bp = &hd->bCoefficients(n);
      }
      else if (hm) {
          bp = &hm->bCoefficients(level, n);
      }

      MultiFab &bcoef = *(MultiFab*)bp;

#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(Flux[n], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
          set_abec_flux(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                        n,
                        BL_TO_FORTRAN_ANYD(Erborder[mfi]), 
                        BL_TO_FORTRAN_ANYD(bcoef[mfi]), 
                        beta, AMREX_REAL_ANYD(dx),
                        BL_TO_FORTRAN_ANYD(Flux[n][mfi]));
      }

  }

  // Correct fluxes at physical and coarse-fine boundaries.

  // Note: It would be good to move much of this section into the
  // linear solver drivers themselves (one layer down), to enable
  // an implementation that does not compute fluxes everywhere in
  // an EdgeVar.  We can't just move the nonsymmetric pieces easily
  // by themselves, though, because the current implementation
  // trashes the boundary fluxes before fixing them.

  if (hd) {
    hd->boundaryFlux(&Flux[0], Er, igroup, Inhomogeneous_BC);
  }
  else if (hm) {
    hm->boundaryFlux(level, &Flux[0], Er, igroup, Inhomogeneous_BC);
  }
}

void RadSolve::levelFluxReg(int level,
                            FluxRegister* flux_in, FluxRegister* flux_out,
                            const Array<MultiFab, BL_SPACEDIM>& Flux,
                            int igroup)
{
  BL_PROFILE("RadSolve::levelFluxReg");

  const Real* dx = parent->Geom(level).CellSize();

  const Real volume = D_TERM(dx[0], * dx[1], * dx[2]);

  if (flux_in) {
    for (int n = 0; n < BL_SPACEDIM; n++) {
      const Real scale = volume / dx[n];
      flux_in->CrseInit(Flux[n], n, 0, igroup, 1, scale);
    }
  }
  if (flux_out) {
    for (OrientationIter face; face; ++face) {
      Orientation ori = face();
      (*flux_out)[ori].setVal(0.0, igroup, 1);
    }
    for (int n = 0; n < BL_SPACEDIM; n++) {
      const Real scale = volume / dx[n];
      flux_out->FineAdd(Flux[n], n, 0, igroup, 1, scale);
    }
  }
}

void RadSolve::levelDterm(int level, MultiFab& Dterm, MultiFab& Er, int igroup)
{
  BL_PROFILE("RadSolve::levelDterm");
  const BoxArray& grids = parent->boxArray(level);
  const DistributionMapping& dmap = parent->DistributionMap(level);
  const Geometry& geom = parent->Geom(level);
  const Real* dx = parent->Geom(level).CellSize();
  const Castro *castro = dynamic_cast<Castro*>(&parent->getLevel(level));

  Array<MultiFab, BL_SPACEDIM> Dterm_face;
  for (int idim=0; idim<BL_SPACEDIM; idim++) {
      Dterm_face[idim].define(castro->getEdgeBoxArray(idim), dmap, 1, 0);
  }

  // grow a larger MultiFab to hold Er so we can difference across faces
  MultiFab Erborder(grids, dmap, 1, 1);
  Erborder.setVal(0.0);
  MultiFab::Copy(Erborder, Er, igroup, 0, 1, 0);

  Erborder.FillBoundary(parent->Geom(level).periodicity()); // zeroes left in off-level boundaries

  HypreExtMultiABec *hem = (HypreExtMultiABec*)hm.get();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (int n = 0; n < BL_SPACEDIM; n++) {
      const MultiFab *dp;

      dp = &hem->d2Coefficients(level, n);
      MultiFab &dcoef = *(MultiFab*)dp;
      
      for (MFIter fi(dcoef,true); fi.isValid(); ++fi) {
	  const Box& bx = fi.tilebox();
	  ca_set_dterm_face(bx.loVect(), bx.hiVect(),
			    BL_TO_FORTRAN(Erborder[fi]),
			    BL_TO_FORTRAN(dcoef[fi]), 
			    BL_TO_FORTRAN(Dterm_face[n][fi]), 
			    dx, &n);
      }
  }

  // Correct D terms at physical and coarse-fine boundaries.
  hem->boundaryDterm(level, &Dterm_face[0], Er, igroup);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
      Vector<Real> rc, re, s;
      
      if (geom.IsSPHERICAL()) {
	  for (MFIter fi(Dterm_face[0]); fi.isValid(); ++fi) {  // omp over boxes
	      int i = fi.index();
	      const Box &reg = grids[i];
	      parent->Geom(level).GetEdgeLoc(re, reg, 0);
	      parent->Geom(level).GetCellLoc(rc, reg, 0);
	      parent->Geom(level).GetCellLoc(s, reg, 0);
	      const Box &dbox = Dterm_face[0][fi].box();
	      sphe(re.dataPtr(), s.dataPtr(), 0,
		   ARLIM(dbox.loVect()), ARLIM(dbox.hiVect()), dx);
	      
	      ca_correct_dterm(D_DECL(BL_TO_FORTRAN(Dterm_face[0][fi]),
				      BL_TO_FORTRAN(Dterm_face[1][fi]),
				      BL_TO_FORTRAN(Dterm_face[2][fi])),
			       re.dataPtr(), rc.dataPtr());
	  }
#ifdef _OPENMP
#pragma omp barrier
#endif
      }
      else if (geom.IsRZ()) {
	  for (MFIter fi(Dterm_face[0]); fi.isValid(); ++fi) {  // omp over boxes
	      int i = fi.index();
	      const Box &reg = grids[i];
	      parent->Geom(level).GetEdgeLoc(re, reg, 0);
	      parent->Geom(level).GetCellLoc(rc, reg, 0);
	      
	      ca_correct_dterm(D_DECL(BL_TO_FORTRAN(Dterm_face[0][fi]),
				      BL_TO_FORTRAN(Dterm_face[1][fi]),
				      BL_TO_FORTRAN(Dterm_face[2][fi])),
			       re.dataPtr(), rc.dataPtr());
	  }
#ifdef _OPENMP
#pragma omp barrier
#endif
      }

      for (MFIter fi(Dterm,true); fi.isValid(); ++fi) {
	  const Box& bx = fi.tilebox();
	  int scomp = 0;
	  int dcomp = 0;
	  int ncomp = 1;
	  int nf = 1;
	  int nc = 1;
	  ca_face2center(bx.loVect(), bx.hiVect(), scomp, dcomp, ncomp, nf, nc,
			 D_DECL(BL_TO_FORTRAN(Dterm_face[0][fi]),
				BL_TO_FORTRAN(Dterm_face[1][fi]),
				BL_TO_FORTRAN(Dterm_face[2][fi])),
			 BL_TO_FORTRAN(Dterm[fi]));
      }
  }
}

// <MGFLD routines>
void RadSolve::computeBCoeffs(MultiFab& bcoefs, int idim,
                              MultiFab& kappa_r, int kcomp,
                              MultiFab& lambda, int lamcomp,
                              Real c, const Geometry& geom)
{
  BL_PROFILE("RadSolve::computeBCoeffs (MGFLD)");
  BL_ASSERT(kappa_r.nGrow() == 1);

  const Real* dx = geom.CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(lambda, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
      bclim(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
            BL_TO_FORTRAN_ANYD(bcoefs[mfi]),
            BL_TO_FORTRAN_N_ANYD(lambda[mfi], lamcomp),
            idim, 
            BL_TO_FORTRAN_N_ANYD(kappa_r[mfi], kcomp), 
            c, AMREX_REAL_ANYD(dx));
  }
}

void RadSolve::levelACoeffs(int level, MultiFab& kpp, 
			    Real delta_t, Real c, int igroup, Real ptc_tau)
{
  BL_PROFILE("RadSolve::levelACoeffs (MGFLD)");
  const BoxArray& grids = parent->boxArray(level);
  const DistributionMapping& dmap = parent->DistributionMap(level);

  // allocate space for ABecLaplacian acoeffs, fill with values

  int Ncomp = 1;
  int Nghost = 0;
  MultiFab acoefs(grids, dmap, Ncomp, Nghost);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
      Vector<Real> r, s;

      for (MFIter mfi(kpp,true); mfi.isValid(); ++mfi) {
	  const Box &reg = mfi.tilebox();
	  
	  getCellCenterMetric(parent->Geom(level), reg, r, s);
	  
	  Real dt_ptc = delta_t/(1.0+ptc_tau);
	  lacoefmgfld(BL_TO_FORTRAN(acoefs[mfi]), 
		      ARLIM(reg.loVect()), ARLIM(reg.hiVect()), 
		      BL_TO_FORTRAN_N(kpp[mfi], igroup), 
		      r.dataPtr(), s.dataPtr(), dt_ptc, c);
      }
  }

  // set a coefficients
  if (hd) {
    hd->aCoefficients(acoefs);
  }
  else if (hm) {
    hm->aCoefficients(level,acoefs);
  }
}


void RadSolve::levelRhs(int level, MultiFab& rhs, const MultiFab& jg, 
			const MultiFab& mugT, const MultiFab& mugY, 
			const MultiFab& coupT, const MultiFab& coupY, 
			const MultiFab& etaT, const MultiFab& etaY, 
			const MultiFab& thetaT, const MultiFab& thetaY, 
			const MultiFab& Er_step, const MultiFab& rhoe_step, const MultiFab& rhoYe_step, 
			const MultiFab& Er_star, const MultiFab& rhoe_star, const MultiFab& rhoYe_star,
			Real delta_t, int igroup, int it, Real ptc_tau)
{
  BL_PROFILE("RadSolve::levelRhs (MGFLD version)");
  Castro *castro = dynamic_cast<Castro*>(&parent->getLevel(level));
  Real time = castro->get_state_data(Rad_Type).curTime();

#ifdef _OPENMP
#pragma omp parallel
#endif
  {  
      Vector<Real> r, s;

      for (MFIter ri(rhs,true); ri.isValid(); ++ri) {

	  const Box &reg = ri.tilebox();

#ifdef MG_SU_OLSON

	  parent->Geom(level).GetCellLoc(r, reg, 0);

	  ca_compute_rhs_so(reg.loVect(), reg.hiVect(),
			    BL_TO_FORTRAN(rhs[ri]),
			    BL_TO_FORTRAN(jg[ri]),
			    BL_TO_FORTRAN(mugT[ri]),
			    BL_TO_FORTRAN(coupT[ri]),
			    BL_TO_FORTRAN(etaT[ri]),
			    BL_TO_FORTRAN(Er_step[ri]),
			    BL_TO_FORTRAN(rhoe_step[ri]),
			    BL_TO_FORTRAN(rhoe_star[ri]),
			    r.dataPtr(), 
			    &time, &delta_t, &igroup);

#else
	  getCellCenterMetric(parent->Geom(level), reg, r, s);

	  ca_compute_rhs(reg.loVect(), reg.hiVect(),
			 BL_TO_FORTRAN(rhs[ri]),
			 BL_TO_FORTRAN(jg[ri]),
			 BL_TO_FORTRAN(mugT[ri]),
			 BL_TO_FORTRAN(coupT[ri]),
			 BL_TO_FORTRAN(etaT[ri]),
			 BL_TO_FORTRAN(Er_step[ri]),
			 BL_TO_FORTRAN(rhoe_step[ri]),
			 BL_TO_FORTRAN(Er_star[ri]),
			 BL_TO_FORTRAN(rhoe_star[ri]),
			 r.dataPtr(), 
			 &delta_t, &igroup, &ptc_tau);
#endif
      }
  }
}

// </ MGFLD routines>

void RadSolve::setHypreMulti(Real cMul, Real d1Mul, Real d2Mul)
{
  HypreExtMultiABec *hem = dynamic_cast<HypreExtMultiABec*>(hm.get());
  if (hem) {
    hem-> cMultiplier() =  cMul;
    hem->d1Multiplier() = d1Mul;
    hem->d2Multiplier() = d2Mul;
  }
}

void RadSolve::restoreHypreMulti()
{
  HypreExtMultiABec *hem = dynamic_cast<HypreExtMultiABec*>(hm.get());
  if (hem) {
    hem-> cMultiplier() =  cMulti;
    hem->d1Multiplier() = d1Multi;
    hem->d2Multiplier() = d2Multi;  
  }
}

void RadSolve::getCellCenterMetric(const Geometry& geom, const Box& reg, Vector<Real>& r, Vector<Real>& s)
{
    const int I = (BL_SPACEDIM >= 2) ? 1 : 0;
    if (geom.IsCartesian()) {
	r.resize(reg.length(0), 1);
	s.resize(reg.length(I), 1);
    }
    else if (geom.IsRZ()) {
	geom.GetCellLoc(r, reg, 0);
	s.resize(reg.length(I), 1);
    }
    else {
	geom.GetCellLoc(r, reg, 0);
	geom.GetCellLoc(s, reg, I);
	const Real *dx = geom.CellSize();
	sphc(r.dataPtr(), s.dataPtr(),
	     ARLIM(reg.loVect()), ARLIM(reg.hiVect()), dx);
    }
}
	
void RadSolve::getEdgeMetric(int idim, const Geometry& geom, const Box& edgebox, 
			     Vector<Real>& r, Vector<Real>& s)
{
    const Box& reg = amrex::enclosedCells(edgebox);
    const int I = (BL_SPACEDIM >= 2) ? 1 : 0;
    if (geom.IsCartesian()) {
	r.resize(reg.length(0)+1, 1);
	s.resize(reg.length(I)+1, 1);
    }
    else if (geom.IsRZ()) {
	if (idim == 0) {
	    geom.GetEdgeLoc(r, reg, 0);
	}
	else {
	    geom.GetCellLoc(r, reg, 0);
	}
	s.resize(reg.length(I)+1, 1);
    }
    else {
      if (idim == 0) {
        geom.GetEdgeLoc(r, reg, 0);
        geom.GetCellLoc(s, reg, I);
      }
      else {
        geom.GetCellLoc(r, reg, 0);
        geom.GetEdgeLoc(s, reg, I);
      }
      const Real *dx = geom.CellSize();
      sphe(r.dataPtr(), s.dataPtr(), idim,
	   ARLIM(edgebox.loVect()), ARLIM(edgebox.hiVect()), dx);
    }
}

