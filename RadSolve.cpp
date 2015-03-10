
#include <ParmParse.H>
#include <AmrLevel.H>

#include <LO_BCTYPES.H>
//#include <CompSolver.H>

#include "RadSolve.H"
#include "Radiation.H"  // for access to static physical constants only

#include <Using.H>

#undef BL_USE_ARLIM

#include "RAD_F.H"

#include "HABEC_F.H"    // only for nonsymmetric flux; may be changed?

Array<Real> RadSolve::absres(0);

RadSolve::RadSolve(Amr* Parent) : parent(Parent),
  hd(NULL), hm(NULL)
{
  ParmParse pp("radsolve");

  if (BL_SPACEDIM == 1) {
    // pfmg will not work in 1D
    level_solver_flag            = 0;
  }
  else {
    level_solver_flag            = 1;
  }
  pp.query("level_solver_flag",            level_solver_flag);

  use_hypre_nonsymmetric_terms = 0;
  pp.query("use_hypre_nonsymmetric_terms", use_hypre_nonsymmetric_terms);

  if (Radiation::SolverType == Radiation::SGFLDSolver 
      && Radiation::Er_Lorentz_term) { 

    use_hypre_nonsymmetric_terms = 1;

    if (level_solver_flag < 100) {
      BoxLib::Error("To do Lorentz term implicitly level_solver_flag must be >= 100.");
    }
  }

  if (Radiation::SolverType == Radiation::MGFLDSolver && 
      Radiation::accelerate == 2 && Radiation::nGroups > 1) {
    use_hypre_nonsymmetric_terms = 1;

    if (level_solver_flag < 100) {
      BoxLib::Error("When accelerate is 2, level_solver_flag must be >= 100.");
    }
  }

  ParmParse ppr("radiation");

  reltol     = 1.0e-10;   pp.query("reltol",  reltol);
  if (Radiation::SolverType == Radiation::SGFLDSolver ||
      Radiation::SolverType == Radiation::MGFLDSolver) {
    abstol = 0.0;
  }
  else {
    abstol     = 1.0e-10;   
  }
  pp.query("abstol",  abstol);
  maxiter    = 40;        pp.query("maxiter", maxiter);

  // For the radiation problem these are always +1:
  alpha = 1.0; pp.query("alpha",alpha);
  beta  = 1.0; pp.query("beta",beta);

  verbose = 0; pp.query("v", verbose); pp.query("verbose", verbose);

  {
    // Putting this here is a kludge, but I make the factors static and
    // enter them here for both kinds of solvers so that any solver
    // objects created by, for example, the CompSolver, will get the
    // right values.  They are set each time this constructor is called
    // to allow for the fact that we might conceivably have two different
    // radiation-like sets of equations being solved with different
    // conventions about flux_factor (photons and neutrinos, for example).
    // The assumption then is that the RadSolve object will only persist
    // for the duration of one type of physics update.
    ParmParse pp1("radiation");
    Real c = Radiation::clight;
    pp1.query("c", c);
    HypreABec::fluxFactor() = c;
    HypreMultiABec::fluxFactor() = c;
  }

  static int first = 1;
  if (verbose >= 1 && first && ParallelDescriptor::IOProcessor()) {
    first = 0;
    cout << "radsolve.level_solver_flag      = " << level_solver_flag << endl;
    cout << "radsolve.maxiter                = " << maxiter << endl;
    cout << "radsolve.reltol                 = " << reltol << endl;
    cout << "radsolve.abstol                 = " << abstol << endl;
    cout << "radsolve.use_hypre_nonsymmetric_terms = "
         << use_hypre_nonsymmetric_terms << endl;
    cout << "radsolve.verbose                = " << verbose << endl;
  }

  // Static initialization:
  if (absres.size() == 0) {
    absres.resize(parent->maxLevel() + 1, 0.0);
  }
}

void RadSolve::levelInit(int level)
{
  BL_PROFILE("RadSolve::levelInit");
  const BoxArray& grids = parent->boxArray(level);
//  const Real *dx = parent->Geom(level).CellSize();

  if (level_solver_flag < 100) {
      hd = new HypreABec(grids, parent->Geom(level), level_solver_flag);
  }
  else {
      if (use_hypre_nonsymmetric_terms == 0) {
	  hm = new HypreMultiABec(level, level, level_solver_flag);
      }
      else {
	  hm = new HypreExtMultiABec(level, level, level_solver_flag);
	  HypreExtMultiABec *hem = (HypreExtMultiABec*)hm;
	  cMulti  = hem->cMultiplier();
	  d1Multi = hem->d1Multiplier();
	  d2Multi = hem->d2Multiplier();
      }
      hm->addLevel(level, parent->Geom(level), grids,
                   IntVect::TheUnitVector());
      hm->buildMatrixStructure();
  }
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

void RadSolve::levelClear()
{
  if (hd) {
    delete hd;
    hd = NULL;
  }
  else if (hm) {
    delete hm;
    hm = NULL;
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
      Array<Real> r, s;

      for (MFIter mfi(cc,true); mfi.isValid(); ++mfi) 
      {
	  const Box &reg = mfi.tilebox();

	  getCellCenterMetric(parent->Geom(level), reg, r, s);
	  
	  FORT_MULTRS(cc[mfi].dataPtr(), dimlist(cc[mfi].box()), dimlist(reg),
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
    HypreExtMultiABec *hem = dynamic_cast<HypreExtMultiABec*>(hm);
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

  // Allocate space for ABecLapacian acoeffs, fill with values

  int Ncomp  = 1;
  int Nghost = 0;

  MultiFab acoefs(grids, Ncomp, Nghost, Fab_allocate);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
      Array<Real> r, s;

      for (MFIter mfi(fkp,true); mfi.isValid(); ++mfi) {
	  const Box &abox = acoefs[mfi].box();
	  const Box &reg  = mfi.tilebox();

	  getCellCenterMetric(parent->Geom(level), reg, r, s);

	  FORT_LACOEF(acoefs[mfi].dataPtr(), dimlist(abox), dimlist(reg),
		      fkp[mfi].dataPtr(), eta[mfi].dataPtr(), etainv[mfi].dataPtr(),
		      r.dataPtr(), s.dataPtr(), c, delta_t, theta);
      }
  }

  if (hd) {
    hd->aCoefficients(acoefs);
  }
  else if (hm) {
    hm->aCoefficients(level, acoefs);
  }
}

void RadSolve::levelSPas(int level, Tuple<MultiFab, BL_SPACEDIM>& lambda, int igroup, 
			 int lo_bc[3], int hi_bc[3])
{
  const BoxArray& grids = parent->boxArray(level);
  const Geometry& geom = parent->Geom(level);
  const Box& domainBox = geom.Domain();

  MultiFab spa(grids, 1, 0);
#ifdef _OPENMP
#pragma omp
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
	  BL_FORT_PROC_CALL(CA_SPALPHA, ca_spalpha)
	      (reg.loVect(), reg.hiVect(),
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
    BoxLib::Abort("Should not be in RadSolve::levelSPas");    
  }
}

void RadSolve::levelBCoeffs(int level,
                            Tuple<MultiFab, BL_SPACEDIM>& lambda,
                            MultiFab& kappa_r, int kcomp,
                            Real c, int lamcomp)
{
  BL_PROFILE("RadSolve::levelBCoeffs");
  BL_ASSERT(kappa_r.nGrow() == 1);

  const Geometry& geom = parent->Geom(level);
  const Real* dx       = geom.CellSize();

  for (int idim = 0; idim < BL_SPACEDIM; idim++) {

    MultiFab bcoefs(lambda[idim].boxArray(), 1, 0, Fab_allocate);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	Array<Real> r, s;

	for (MFIter mfi(lambda[idim],true); mfi.isValid(); ++mfi) {
	    const Box &bbox = lambda[idim][mfi].box();
	    const Box &kbox = kappa_r[mfi].box();

	    const Box &ndbox  = mfi.tilebox();
	    getEdgeMetric(idim, geom, ndbox, r, s);

	    const Box& reg = BoxLib::enclosedCells(ndbox);
	    FORT_BCLIM(bcoefs[mfi].dataPtr(), lambda[idim][mfi].dataPtr(lamcomp),
		       dimlist(bbox), dimlist(reg),
		       idim, kappa_r[mfi].dataPtr(kcomp), dimlist(kbox),
		       r.dataPtr(), s.dataPtr(), c, dx);
	}
    }

    if (hd) {
	hd->bCoefficients(bcoefs, idim);
    }
    else if (hm) {
	hm->bCoefficients(level, bcoefs, idim);
    }
  } // -->> over dimension
}

void RadSolve::levelDCoeffs(int level, Tuple<MultiFab, BL_SPACEDIM>& lambda,
			    MultiFab& vel, MultiFab& dcf)
{
    BL_PROFILE("RadSolve::levelDCoeffs");
    const BoxArray& grids = parent->boxArray(level);

    for (int idim=0; idim<BL_SPACEDIM; idim++) {

	BoxArray edge_boxes(grids);
	edge_boxes.surroundingNodes(idim);
	MultiFab dcoefs(edge_boxes, 1, 0, Fab_allocate);

#ifdef _OPENMP
#pragma omp parallel
#endif
	{
	    Array<Real> r, s;
	    
	    for (MFIter mfi(dcoefs,true); mfi.isValid(); ++mfi) {
		const Box& ndbx = mfi.tilebox();

		getEdgeMetric(idim, parent->Geom(level), ndbx, r, s);

		BL_FORT_PROC_CALL(CA_COMPUTE_DCOEFS, ca_compute_dcoefs)
		    (ndbx.loVect(), ndbx.hiVect(),
		     BL_TO_FORTRAN(dcoefs[mfi]), BL_TO_FORTRAN(lambda[idim][mfi]),
		     BL_TO_FORTRAN(vel[mfi]), BL_TO_FORTRAN(dcf[mfi]), 
		     r.dataPtr(), &idim);
	    }
	}

	HypreExtMultiABec *hem = (HypreExtMultiABec*)hm;
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

  rhs.setVal(0.0);
  if (fine_corr) {
    // This works trivially for a multilevel solve since the finer level is
    // present in the rhs to overwrite the junk produced under it.
    // In a single-level version we have to be sure that fine_corr
    // has been cleaned up using clear_internal_borders.

    // Hack:  For the single group case igroup defaults to -1, which is
    // significant later in this routine.  So we have to construct the
    // correct component number here:

    int igrouptmp = (igroup < 0) ? 0 : igroup;

    fine_corr->Reflux(rhs, scale, igrouptmp, 0, 1, parent->Geom(level));
  }

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
      Array<Real> r, s;
      
      for (MFIter ri(rhs,true); ri.isValid(); ++ri) {
	  const Box &rbox = rhs[ri].box();
	  const Box &ebox = Er_old[ri].box();
	  const Box &reg  = ri.tilebox();
	  
	  getCellCenterMetric(parent->Geom(level), reg, r, s);
	  
	  FORT_LRHS(rhs[ri].dataPtr(), dimlist(rbox), dimlist(reg),
		    temp[ri].dataPtr(),
		    fkp[ri].dataPtr(), eta[ri].dataPtr(), etainv[ri].dataPtr(),
		    rhoem[ri].dataPtr(), rhoes[ri].dataPtr(),
		    dflux_old[ri].dataPtr(),
		    Er_old[ri].dataPtr(0), dimlist(ebox),
		    Edot[ri].dataPtr(),
		    r.dataPtr(), s.dataPtr(), delta_t, sigma, c, theta);
      }
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
      int oldprec = cout.precision(20);
      cout << "Absolute residual = " << res << endl;
      cout.precision(oldprec);
    }
    res *= sync_absres_factor;
    absres[level] = (absres[level] > res) ? absres[level] : res;
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
      int oldprec = cout.precision(20);
      cout << "Absolute residual = " << res << endl;
      cout.precision(oldprec);
    }
    res *= sync_absres_factor;
    absres[level] = (absres[level] > res) ? absres[level] : res;
    hm->clearSolver();
  }
}


void RadSolve::levelFlux(int level,
                         FluxRegister* flux_in, FluxRegister* flux_out,
                         MultiFab& Er, int igroup)
{
  BL_PROFILE("RadSolve::levelFlux");
  const BoxArray& grids = parent->boxArray(level);

  Tuple<MultiFab, BL_SPACEDIM> Flux;
  for (int n = 0; n < BL_SPACEDIM; n++) {
    BoxArray edge_boxes(grids);
    edge_boxes.surroundingNodes(n);
    Flux[n].define(edge_boxes, 1, 0, Fab_allocate);
  }

  levelFlux(level, Flux, Er, igroup);

  levelFluxReg(level, flux_in, flux_out, Flux, igroup);
}

void RadSolve::levelFlux(int level,
                         FluxRegister* flux_in, FluxRegister* flux_out,
                         MultiFab& Er, int igroup, MultiFab& flx)
{
  BL_PROFILE("RadSolve::levelFlux");
  const BoxArray& grids = parent->boxArray(level);

  Tuple<MultiFab, BL_SPACEDIM> Flux;
  for (int n = 0; n < BL_SPACEDIM; n++) {
    BoxArray edge_boxes(grids);
    edge_boxes.surroundingNodes(n);
    Flux[n].define(edge_boxes, 1, 0, Fab_allocate);
  }

  levelFlux(level, Flux, Er, igroup);

  levelFluxReg(level, flux_in, flux_out, Flux, igroup);

  levelFluxFaceToCenter(level, Flux, flx, igroup);
}

void RadSolve::levelFluxFaceToCenter(int level, Tuple<MultiFab, BL_SPACEDIM>& Flux,
				     MultiFab& flx, int igroup)
{
    int nflx = flx.nComp();
    
    const Geometry& geom = parent->Geom(level);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	Array<Real> r, s;
    
	for (int idim = 0; idim < BL_SPACEDIM; idim++) {
	    for (MFIter mfi(flx,true); mfi.isValid(); ++mfi) 
	    {
		const Box &ccbx  = mfi.tilebox();
		const Box &ndbx = BoxLib::surroundingNodes(ccbx, idim);

		getEdgeMetric(idim, geom, ndbx, r, s);

		int rlo = ndbx.smallEnd(0);
		int rhi = rlo + r.size() - 1;

		BL_FORT_PROC_CALL(CA_TEST_TYPE_FLUX, ca_test_type_flux)
		    (ccbx.loVect(), ccbx.hiVect(),
		     BL_TO_FORTRAN(flx[mfi]),
		     BL_TO_FORTRAN(Flux[idim][mfi]),
		     r.dataPtr(), &rlo, &rhi, 
		     &nflx, &idim, &igroup);
	    }
	}
    }
}

void RadSolve::levelFluxFaceToCenter(int level, MultiFab& state,
				     Tuple<MultiFab, BL_SPACEDIM>& lambda,
				     MultiFab& Er, Tuple<MultiFab, BL_SPACEDIM>& Flux,
				     MultiFab& flx)
{
  int nflx = flx.nComp();

  const Geometry& geom = parent->Geom(level);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	Array<Real> r, s;

	for (int idim = 0; idim < BL_SPACEDIM; idim++) {
	    for (MFIter mfi(flx,true); mfi.isValid(); ++mfi)
	    {
		const Box &ccbx  = mfi.tilebox();
		const Box &ndbx = BoxLib::surroundingNodes(ccbx, idim);

		getEdgeMetric(idim, geom, ndbx, r, s);

		int rlo = ndbx.smallEnd(0);
		int rhi = rlo + r.size() - 1;

		BL_FORT_PROC_CALL(CA_TEST_TYPE_FLUX_LAB, ca_test_type_flux_lab)
		    (ccbx.loVect(), ccbx.hiVect(),
		     BL_TO_FORTRAN(flx[mfi]),
		     BL_TO_FORTRAN(Flux[idim][mfi]),
		     D_DECL(BL_TO_FORTRAN(lambda[0][mfi]),
			    BL_TO_FORTRAN(lambda[1][mfi]),
			    BL_TO_FORTRAN(lambda[2][mfi])),
		     BL_TO_FORTRAN(Er[mfi]),
		     BL_TO_FORTRAN(state[mfi]),
		     r.dataPtr(), &rlo, &rhi, 
		     &nflx, &idim);
	    }
	}
  }
}

void RadSolve::levelFlux(int level,
                         Tuple<MultiFab, BL_SPACEDIM>& Flux,
                         MultiFab& Er, int igroup)
{
  BL_PROFILE("RadSolve::levelFlux");
  const BoxArray& grids = parent->boxArray(level);

  // grow a larger MultiFab to hold Er so we can difference across faces
  MultiFab Erborder(grids, 1, 1);
  Erborder.setVal(0.0);
  MultiFab::Copy(Erborder, Er, igroup, 0, 1, 0);

  Erborder.FillBoundary(); // zeroes left in off-level boundaries
  if (parent->Geom(level).isAnyPeriodic()) {
      parent->Geom(level).FillPeriodicBoundary(Erborder, true);
  }

  const Real* dx = parent->Geom(level).CellSize();

#ifdef _OPENMP
#pragam omp parallel
#endif
  for (int n = 0; n < BL_SPACEDIM; n++) {
    const MultiFab *bp; //, *cp;
    if (hd) {
      bp = &hd->bCoefficients(n);
    }
    else if (hm) {
      bp = &hm->bCoefficients(level, n);
    }
    // w.z. I commented this out because we may not always have ccoef 
    //      when use_hypre_nonsymmetric_terms == 1.
    //      And ccoef is not being used anyway.
    // if (use_hypre_nonsymmetric_terms == 1) {
    //   HypreExtMultiABec *hem = (HypreExtMultiABec*)hm;
    //   cp = &hem->cCoefficients(level, n);
    // }
    MultiFab &bcoef = *(MultiFab*)bp;
    //    MultiFab &ccoef = *(MultiFab*)cp;

    for (MFIter fi(Flux[n],true); fi.isValid(); ++fi) {
	const Box& reg = fi.tilebox();
	FORT_SET_ABEC_FLUX(dimlist(reg), &n,
			   Erborder[fi].dataPtr(), dimlist(Erborder[fi].box()),
			   bcoef[fi].dataPtr(),    dimlist(bcoef[fi].box()),
			   &beta,
			   dx,
			   Flux[n][fi].dataPtr(),  dimlist(Flux[n][fi].box()));
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
  if (use_hypre_nonsymmetric_terms == 1) {
    //HypreExtMultiABec *hem = (HypreExtMultiABec*)hm;
    //hem->boundaryFlux(level, &Flux[0], Er);
  }
}

void RadSolve::levelFluxReg(int level,
                            FluxRegister* flux_in, FluxRegister* flux_out,
                            Tuple<MultiFab, BL_SPACEDIM>& Flux,
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
  const Real* dx = parent->Geom(level).CellSize();

  Tuple<MultiFab, BL_SPACEDIM> Dterm_face;
  for (int idim=0; idim<BL_SPACEDIM; idim++) {
    BoxArray edge_boxes(grids);
    edge_boxes.surroundingNodes(idim);
    Dterm_face[idim].define(edge_boxes, 1, 0, Fab_allocate);
  }

  // grow a larger MultiFab to hold Er so we can difference across faces
  MultiFab Erborder(grids, 1, 1);
  Erborder.setVal(0.0);
  MultiFab::Copy(Erborder, Er, igroup, 0, 1, 0);

  Erborder.FillBoundary(); // zeroes left in off-level boundaries
  if (parent->Geom(level).isAnyPeriodic()) {
    parent->Geom(level).FillPeriodicBoundary(Erborder, true);
  }

  HypreExtMultiABec *hem = (HypreExtMultiABec*)hm;

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (int n = 0; n < BL_SPACEDIM; n++) {
      const MultiFab *dp;

      dp = &hem->d2Coefficients(level, n);
      MultiFab &dcoef = *(MultiFab*)dp;
      
      for (MFIter fi(dcoef,true); fi.isValid(); ++fi) {
	  const Box& bx = fi.tilebox();
	  BL_FORT_PROC_CALL(CA_SET_DTERM_FACE, ca_set_dterm_face)
	      (bx.loVect(), bx.hiVect(),
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
      Array<Real> rc, re, s;
      
      if (Geometry::IsSPHERICAL()) {
	  for (MFIter fi(Dterm_face[0]); fi.isValid(); ++fi) {  // omp over boxes
	      int i = fi.index();
	      const Box &reg = grids[i];
	      parent->Geom(level).GetEdgeLoc(re, reg, 0);
	      parent->Geom(level).GetCellLoc(rc, reg, 0);
	      parent->Geom(level).GetCellLoc(s, reg, 0);
	      const Box &dbox = Dterm_face[0][fi].box();
	      FORT_SPHE(re.dataPtr(), s.dataPtr(), 0, dimlist(dbox), dx);
	      
	      BL_FORT_PROC_CALL(CA_CORRECT_DTERM, ca_correct_dterm)
		  (D_DECL(BL_TO_FORTRAN(Dterm_face[0][fi]),
			  BL_TO_FORTRAN(Dterm_face[1][fi]),
			  BL_TO_FORTRAN(Dterm_face[2][fi])),
		   re.dataPtr(), rc.dataPtr());
	  }
#ifdef _OPENMP
#pragma omp barrier
#endif
      }
      else if (Geometry::IsRZ()) {
	  for (MFIter fi(Dterm_face[0]); fi.isValid(); ++fi) {  // omp over boxes
	      int i = fi.index();
	      const Box &reg = grids[i];
	      parent->Geom(level).GetEdgeLoc(re, reg, 0);
	      parent->Geom(level).GetCellLoc(rc, reg, 0);
	      
	      BL_FORT_PROC_CALL(CA_CORRECT_DTERM, ca_correct_dterm)
		  (D_DECL(BL_TO_FORTRAN(Dterm_face[0][fi]),
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
	  BL_FORT_PROC_CALL(CA_FACE2CENTER, ca_face2center)
	      (bx.loVect(), bx.hiVect(),
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
  {
      Array<Real> r, s;

      for (MFIter mfi(lambda,true); mfi.isValid(); ++mfi) {
	  const Box &bbox = lambda[mfi].box();
	  const Box &kbox = kappa_r[mfi].box();

	  const Box &ndbx  = mfi.tilebox();
	  const Box &reg   = BoxLib::enclosedCells(ndbx);

	  getEdgeMetric(idim, geom, ndbx, r, s);
    
	  FORT_BCLIM(bcoefs[mfi].dataPtr(), lambda[mfi].dataPtr(lamcomp),
		     dimlist(bbox), dimlist(reg),
		     idim, kappa_r[mfi].dataPtr(kcomp), dimlist(kbox),
		     r.dataPtr(), s.dataPtr(), c, dx);
      }
  }
}

void RadSolve::levelACoeffs(int level, MultiFab& kpp, 
			    Real delta_t, Real c, int igroup, Real ptc_tau)
{
  BL_PROFILE("RadSolve::levelACoeffs (MGFLD)");
  const BoxArray& grids = parent->boxArray(level);

  // allocate space for ABecLaplacian acoeffs, fill with values

  int Ncomp = 1;
  int Nghost = 0;
  MultiFab acoefs(grids, Ncomp, Nghost, Fab_allocate);

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
      Array<Real> r, s;

      for (MFIter mfi(kpp,true); mfi.isValid(); ++mfi) {
	  const Box &abox = acoefs[mfi].box();
	  const Box &reg = mfi.tilebox();
	  
	  getCellCenterMetric(parent->Geom(level), reg, r, s);
	  
	  const Box &kbox = kpp[mfi].box();
	  Real dt_ptc = delta_t/(1.0+ptc_tau);
	  FORT_LACOEFMGFLD(acoefs[mfi].dataPtr(), dimlist(abox), dimlist(reg),
			   kpp[mfi].dataPtr(igroup), dimlist(kbox),
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
      Array<Real> r, s;

      for (MFIter ri(rhs,true); ri.isValid(); ++ri) {

	  const Box &reg = ri.tilebox();

#ifdef MG_SU_OLSON

	  parent->Geom(level).GetCellLoc(r, reg, 0);

	  BL_FORT_PROC_CALL(CA_COMPUTE_RHS_SO, ca_compute_rhs_so)
	      (reg.loVect(), reg.hiVect(),
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

#ifdef NEUTRINO
	  BL_FORT_PROC_CALL(CA_COMPUTE_RHS_NEUT, ca_compute_rhs_neut)
	      (reg.loVect(), reg.hiVect(),
	       BL_TO_FORTRAN(rhs[ri]),
	       BL_TO_FORTRAN(jg[ri]),
	       BL_TO_FORTRAN(mugT[ri]),
	       BL_TO_FORTRAN(mugY[ri]),
	       BL_TO_FORTRAN(coupT[ri]),
	       BL_TO_FORTRAN(coupY[ri]),
	       BL_TO_FORTRAN(etaT[ri]),
	       BL_TO_FORTRAN(etaY[ri]),
	       BL_TO_FORTRAN(thetaT[ri]),
	       BL_TO_FORTRAN(thetaY[ri]),
	       BL_TO_FORTRAN(Er_step[ri]),
	       BL_TO_FORTRAN(rhoe_step[ri]),
	       BL_TO_FORTRAN(rhoYe_step[ri]),
	       BL_TO_FORTRAN(Er_star[ri]),
	       BL_TO_FORTRAN(rhoe_star[ri]),
	       BL_TO_FORTRAN(rhoYe_star[ri]),
	       r.dataPtr(), 
	       &delta_t, &igroup, &ptc_tau);
#else
	  BL_FORT_PROC_CALL(CA_COMPUTE_RHS, ca_compute_rhs)
	      (reg.loVect(), reg.hiVect(),
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
#endif
      }
  }
}

// </ MGFLD routines>

void RadSolve::setHypreMulti(Real cMul, Real d1Mul, Real d2Mul)
{
  HypreExtMultiABec *hem = dynamic_cast<HypreExtMultiABec*>(hm);
  if (hem) {
    hem-> cMultiplier() =  cMul;
    hem->d1Multiplier() = d1Mul;
    hem->d2Multiplier() = d2Mul;
  }
}

void RadSolve::restoreHypreMulti()
{
  HypreExtMultiABec *hem = dynamic_cast<HypreExtMultiABec*>(hm);
  if (hem) {
    hem-> cMultiplier() =  cMulti;
    hem->d1Multiplier() = d1Multi;
    hem->d2Multiplier() = d2Multi;  
  }
}

void RadSolve::getCellCenterMetric(const Geometry& geom, const Box& reg, Array<Real>& r, Array<Real>& s)
{
    const int I = (BL_SPACEDIM >= 2) ? 1 : 0;
    if (Geometry::IsCartesian()) {
	r.resize(reg.length(0), 1);
	s.resize(reg.length(I), 1);
    }
    else if (Geometry::IsRZ()) {
	geom.GetCellLoc(r, reg, 0);
	s.resize(reg.length(I), 1);
    }
    else {
	geom.GetCellLoc(r, reg, 0);
	geom.GetCellLoc(s, reg, I);
	const Real *dx = geom.CellSize();
	FORT_SPHC(r.dataPtr(), s.dataPtr(), dimlist(reg), dx);
    }
}
	
void RadSolve::getEdgeMetric(int idim, const Geometry& geom, const Box& edgebox, 
			     Array<Real>& r, Array<Real>& s)
{
    const Box& reg = BoxLib::enclosedCells(edgebox);
    const int I = (BL_SPACEDIM >= 2) ? 1 : 0;
    if (Geometry::IsCartesian()) {
	r.resize(reg.length(0)+1, 1);
	s.resize(reg.length(I)+1, 1);
    }
    else if (Geometry::IsRZ()) {
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
      FORT_SPHE(r.dataPtr(), s.dataPtr(), idim, dimlist(edgebox), dx);
    }
}

