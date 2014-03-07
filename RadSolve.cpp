
#include <ParmParse.H>
#include <AmrLevel.H>

#include <LO_BCTYPES.H>
#include <CompSolver.H>

#include "RadSolve.H"
#include "Radiation.H"  // for access to static physical constants only

#include <Using.H>

#undef BL_USE_ARLIM

#include "RAD_F.H"

#include "HABEC_F.H"    // only for nonsymmetric flux; may be changed?

Array<Real> RadSolve::absres(0);

RadSolve::RadSolve(Amr* Parent) : parent(Parent),
  hd(NULL), hm(NULL), solver(NULL)
{
  ParmParse pp("radsolve");

  if (BL_SPACEDIM == 1) {
    // pfmg will not work in 1D
    use_hypre_level              = 1;
    use_hypre_multilevel         = 0;
    use_hypre_nonsymmetric_terms = 0;
    level_solver_flag            = 0;
    multilevel_solver_flag       = 0;
    pp.query("use_hypre_level",              use_hypre_level);
    pp.query("use_hypre_multilevel",         use_hypre_multilevel);
    pp.query("use_hypre_nonsymmetric_terms", use_hypre_nonsymmetric_terms);
    pp.query("level_solver_flag",            level_solver_flag);
    pp.query("multilevel_solver_flag",       multilevel_solver_flag);
  }
  else {
    use_hypre_level              = 1;
    use_hypre_multilevel         = 0;
    use_hypre_nonsymmetric_terms = 0;
    level_solver_flag            = 1;
    multilevel_solver_flag       = 1;
    pp.query("use_hypre_level",              use_hypre_level);
    pp.query("use_hypre_multilevel",         use_hypre_multilevel);
    pp.query("use_hypre_nonsymmetric_terms", use_hypre_nonsymmetric_terms);
    pp.query("level_solver_flag",            level_solver_flag);
    pp.query("multilevel_solver_flag",       multilevel_solver_flag);
  }

  if (!use_hypre_level) {
    if ( ParallelDescriptor::IOProcessor() ) {
      cout << "radsolve.use_hypre_level = 1 is currently required" << endl;
    }
    exit(0);
  }

  if (Radiation::SolverType == Radiation::SGFLDSolver 
      && Radiation::Er_Lorentz_term) { 

    use_hypre_nonsymmetric_terms = 1;

    //    static int first = 1;
    
    if (level_solver_flag < 100) {
      BoxLib::Error("To do Lorentz term implicitly level_solver_flag must be >= 100.");
      // int old_flag = level_solver_flag;
      // level_solver_flag = 109;  // PFMG-preconditioned GMRES
      // use_hypre_nonsymmetric_terms = 1;
      // if ( ParallelDescriptor::IOProcessor() && first ) {
      // 	cout << "To do Lorentz term implicitly level_solver_flag must be >= 100." << endl;
      // 	cout << "level_solver_flag has been reset from "<<old_flag<<" to "<<level_solver_flag<<endl;
      // 	first = 0;
      // }
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

  // additional CompSolver parameters
  use_harmonic_avg   = 1; pp.query("use_harmonic_avg",     use_harmonic_avg);
  multilevel_version = 1; pp.query("multilevel_version",   multilevel_version);
  cs_reltol_mult     = 100.0;   pp.query("cs_reltol_mult", cs_reltol_mult);
  bottomnumiter      = 1;       pp.query("bottomnumiter",  bottomnumiter);
  bottomtol          = 1.0e-6;  pp.query("bottomtol",      bottomtol);
  secondtol          = 1.0e-10; pp.query("secondtol",      secondtol);
  presmooth          = -1;      pp.query("presmooth",      presmooth);
  postsmooth         = -1;      pp.query("postsmooth",     postsmooth);

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
    cout << "radsolve.use_hypre_level        = " << use_hypre_level << endl;
    cout << "radsolve.level_solver_flag      = " << level_solver_flag << endl;
    cout << "radsolve.use_hypre_multilevel   = "
         << use_hypre_multilevel << endl;
    cout << "radsolve.multilevel_solver_flag = "
         << multilevel_solver_flag << endl;
    cout << "radsolve.maxiter                = " << maxiter << endl;
    cout << "radsolve.reltol                 = " << reltol << endl;
    cout << "radsolve.abstol                 = " << abstol << endl;
    cout << "radsolve.use_hypre_nonsymmetric_terms = "
         << use_hypre_nonsymmetric_terms << endl;
    if (use_hypre_multilevel == 0) {
      cout << "radsolve.use_harmonic_avg       = " << use_harmonic_avg << endl;
      cout << "radsolve.multilevel_version     = " << multilevel_version<<endl;
      cout << "radsolve.cs_reltol_mult         = " << cs_reltol_mult << endl;
      cout << "radsolve.bottomnumiter          = " << bottomnumiter << endl;
      cout << "radsolve.bottomtol              = " << bottomtol << endl;
      cout << "radsolve.secondtol              = " << secondtol << endl;
      if (presmooth > 0 || postsmooth > 0) {
        cout << "radsolve.presmooth              = " << presmooth << endl;
        cout << "radsolve.postsmooth             = " << postsmooth << endl;
      }
    }
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
  const Real *dx = parent->Geom(level).CellSize();

  if (use_hypre_level) {
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
  const BoxArray& grids = parent->boxArray(level);

  Array<Real> r, s;

  BL_ASSERT(cc.nGrow() == 0);

  for (MFIter mfi(cc); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box &reg = grids[i];
    const int I = (BL_SPACEDIM >= 2) ? 1 : 0;
    if (Geometry::IsCartesian()) {
      r.resize(reg.length(0), 1);
      s.resize(reg.length(I), 1);
    }
    else if (Geometry::IsRZ()) {
      parent->Geom(level).GetCellLoc(r, reg, 0);
      s.resize(reg.length(I), 1);
    }
    else {
      parent->Geom(level).GetCellLoc(r, reg, 0);
      parent->Geom(level).GetCellLoc(s, reg, I);
      const Real *dx = parent->Geom(level).CellSize();
      FORT_SPHC(r.dataPtr(), s.dataPtr(), dimlist(reg), dx);
    }

    FORT_MULTRS(cc[i].dataPtr(), dimlist(reg), dimlist(reg),
                r.dataPtr(), s.dataPtr());
  }
}

void RadSolve::setLevelACoeffs(int level, const MultiFab& acoefs)
{
  if (hd || hm) {
    if (hd) {
      hd->aCoefficients(acoefs);
    }
    else if (hm) {
      hm->aCoefficients(level, acoefs);
    }
  }
  else if (solver) {
    solver->aCoefficients(level - base, acoefs);
  }
}

void RadSolve::setLevelBCoeffs(int level, const MultiFab& bcoefs, int dir)
{
  if (hd || hm) {
    if (hd) {
      hd->bCoefficients(bcoefs, dir);
    }
    else if (hm) {
      hm->bCoefficients(level, bcoefs, dir);
    }
  }
  else if (solver) {
    solver->bCoefficients(level - base, bcoefs, dir);
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

const MultiFab& RadSolve::getLevelACoeffs(int level)
{
  const MultiFab *ap;
  if (hd || hm) {
    if (hd) {
      ap = &hd->aCoefficients();
    }
    else if (hm) {
      ap = &hm->aCoefficients(level);
    }
  }
  else if (solver) {
    BoxLib::Error("CompSolver can not return coefficients");
    //ap = &solver->Acoef();
  }
  return *ap;
}

const MultiFab& RadSolve::getLevelBCoeffs(int level, int dir)
{
  const MultiFab *bp;
  if (hd || hm) {
    if (hd) {
      bp = &hd->bCoefficients(dir);
    }
    else if (hm) {
      bp = &hm->bCoefficients(level, dir);
    }
  }
  else if (solver) {
    BoxLib::Error("CompSolver can not return coefficients");
    //bp = &solver->Bcoef();
  }
  return *bp;
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

  Array<Real> r, s;

  MultiFab acoefs(grids, Ncomp, Nghost, Fab_allocate);

  for (MFIter mfi(fkp); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box &abox = acoefs[i].box();
    const Box &reg  = grids[i];
    const int I = (BL_SPACEDIM >= 2) ? 1 : 0;
    if (Geometry::IsCartesian()) {
      r.resize(reg.length(0), 1);
      s.resize(reg.length(I), 1);
    }
    else if (Geometry::IsRZ()) {
      parent->Geom(level).GetCellLoc(r, reg, 0);
      s.resize(reg.length(I), 1);
    }
    else {
      parent->Geom(level).GetCellLoc(r, reg, 0);
      parent->Geom(level).GetCellLoc(s, reg, I);
      const Real *dx = parent->Geom(level).CellSize();
      FORT_SPHC(r.dataPtr(), s.dataPtr(), dimlist(reg), dx);
    }
    FORT_LACOEF(acoefs[i].dataPtr(), dimlist(abox), dimlist(reg),
		fkp[i].dataPtr(), eta[i].dataPtr(), etainv[i].dataPtr(),
		r.dataPtr(), s.dataPtr(), c, delta_t, theta);
  }

  if (hd) {
    hd->aCoefficients(acoefs);
  }
  else if (hm) {
    hm->aCoefficients(level, acoefs);
  }
}

void RadSolve::computeBCoeffs(MultiFab& bcoefs, int idim,
                              MultiFab& kappa_r, int kcomp,
                              MultiFab& Erborder, int igroup,
                              Real c, int limiter,
                              const Geometry& geom)
{
  BL_PROFILE("RadSolve::computeBCoeffs");
  const BoxArray& grids = kappa_r.boxArray(); // valid region only

  BL_ASSERT(kappa_r.nGrow() == 1);

  Array<Real> q, r, s;

  const Real *dx = geom.CellSize();

  for (MFIter mfi(bcoefs); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box &bbox = bcoefs[i].box();
    const Box &reg  = grids[i];
    const Box &kbox = kappa_r[i].box();
    const int I = (BL_SPACEDIM >= 2) ? 1 : 0;
    if (Geometry::IsCartesian()) {
      q.resize(reg.length(0)+1, 1);
      r.resize(reg.length(0)+1, 1);
      s.resize(reg.length(I)+1, 1);
    }
    else if (Geometry::IsRZ()) {
      q.resize(reg.length(0)+1, 1);
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
        geom.GetEdgeLoc(q, reg, 0);
        geom.GetEdgeLoc(r, reg, 0);
        geom.GetCellLoc(s, reg, I);
      }
      else {
        geom.GetCellLoc(q, reg, 0);
        geom.GetCellLoc(r, reg, 0);
        geom.GetEdgeLoc(s, reg, I);
      }
      FORT_SPHE(r.dataPtr(), s.dataPtr(), idim, dimlist(bbox), dx);
    }
    if (limiter == 0) {
      FORT_BCLIM0(bcoefs[i].dataPtr(), dimlist(bbox), dimlist(reg),
		  idim, kappa_r[i].dataPtr(kcomp), dimlist(kbox),
		  q.dataPtr(), r.dataPtr(), s.dataPtr(), c, dx);
    }
    else if (limiter%10 == 1) {
      FORT_BCLIM1(bcoefs[i].dataPtr(), dimlist(bbox), dimlist(reg),
                  idim, kappa_r[i].dataPtr(kcomp), dimlist(kbox),
                  Erborder[i].dataPtr(igroup),
                  q.dataPtr(), r.dataPtr(), s.dataPtr(), c, dx, limiter);
    }
    else if (limiter%10 == 2) {
#if (BL_SPACEDIM >= 2)
      Fab dtmp(kbox, BL_SPACEDIM - 1);
#endif
      FORT_BCLIM2(bcoefs[i].dataPtr(), dimlist(bbox), dimlist(reg),
                  idim, kappa_r[i].dataPtr(kcomp), dimlist(kbox),
        D_DECL(Erborder[i].dataPtr(igroup), dtmp.dataPtr(0), dtmp.dataPtr(1)),
                  q.dataPtr(), r.dataPtr(), s.dataPtr(), c, dx, limiter);
    }
    else {
#if (BL_SPACEDIM >= 2)
      Fab dtmp(kbox, BL_SPACEDIM - 1);
#endif
      FORT_BCLIM3(bcoefs[i].dataPtr(), dimlist(bbox), dimlist(reg),
                  idim, kappa_r[i].dataPtr(kcomp), dimlist(kbox),
        D_DECL(Erborder[i].dataPtr(igroup), dtmp.dataPtr(0), dtmp.dataPtr(1)),
                  q.dataPtr(), r.dataPtr(), s.dataPtr(), c, dx, limiter);
    }
  }
}

void RadSolve::levelSPas(int level, Tuple<MultiFab, BL_SPACEDIM>& lambda, int igroup, 
			 int lo_bc[3], int hi_bc[3])
{
  const BoxArray& grids = parent->boxArray(level);
  const Geometry& geom = parent->Geom(level);
  const Box& domainBox = geom.Domain();

  MultiFab spa(grids, 1, 0);
  for (MFIter mfi(spa); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box& reg  = grids[i]; 
    
    spa[i].setVal(1.e210);
    
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
	(BL_TO_FORTRAN(spa[i]),
	 D_DECL(BL_TO_FORTRAN(lambda[0][i]),
		BL_TO_FORTRAN(lambda[1][i]),
		BL_TO_FORTRAN(lambda[2][i])),
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
  const BoxArray& grids = parent->boxArray(level);

  BL_ASSERT(kappa_r.nGrow() == 1);

  Array<Real> r, s;

  const Geometry& geom = parent->Geom(level);
  const Real* dx       = geom.CellSize();

  for (int idim = 0; idim < BL_SPACEDIM; idim++) {

    MultiFab bcoefs(lambda[idim].boxArray(), 1, 0, Fab_allocate);

    for (MFIter mfi(lambda[idim]); mfi.isValid(); ++mfi) {
      int i = mfi.index();
      const Box &bbox = lambda[idim][i].box();
      const Box &reg  = grids[i];
      const Box &kbox = kappa_r[i].box();

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
      else { // support only 1D spherical here
        geom.GetEdgeLoc(r, reg, 0);
        geom.GetCellLoc(s, reg, I);
        FORT_SPHE(r.dataPtr(), s.dataPtr(), idim, dimlist(bbox), dx);
      }

      FORT_BCLIM(bcoefs[i].dataPtr(), lambda[idim][i].dataPtr(lamcomp),
                 dimlist(bbox), dimlist(reg),
                 idim, kappa_r[i].dataPtr(kcomp), dimlist(kbox),
                 r.dataPtr(), s.dataPtr(), c, dx);
    }

    if (hd || hm) {
      if (hd) {
	hd->bCoefficients(bcoefs, idim);
      }
      else if (hm) {
	hm->bCoefficients(level, bcoefs, idim);
      }
    }
    else if (solver) {
      solver->bCoefficients(level - base, bcoefs, idim);
    }

  } // -->> over dimension
}

void RadSolve::levelBCoeffs(int level,
                            MultiFab& kappa_r, int kcomp,
                            MultiFab& Er, int igroup,
                            Real c, int limiter)
{
  BL_PROFILE("RadSolve::levelBCoeffs");
  const BoxArray& grids = parent->boxArray(level);

  BL_ASSERT(kappa_r.nGrow() == 1);

  int idim;

  MultiFab Erborder;
  if (limiter > 0) {
    Erborder.define(grids,1,1,Fab_allocate);
    for (MFIter mfi(Erborder); mfi.isValid(); ++mfi) {
      int i = mfi.index();
      Erborder[i].setVal(-1.0);
      Erborder[i].copy(Er[i], igroup, 0, 1);
    }
    // Values in ghost cells are set to -1, indicating that one-sided
    // differences should be used in computing the gradient term for
    // the flux limiter.  In order to make the solution independent
    // of the grid layout, we now go back and overwrite values in
    // those cells bordering grids at the same level:

    Erborder.FillBoundary();

    if (parent->Geom(level).isAnyPeriodic()) {
      parent->Geom(level).FillPeriodicBoundary(Erborder, true);
    }
  }

  // Allocate space for ABecLapacian coeffs, fill with values

  int Ncomp  = 1;
  int Nghost = 0;

  Array<Real> q, r, s;

  for (idim = 0; idim < BL_SPACEDIM; idim++) {

    BoxArray bsC(grids);
    MultiFab bcoefs(bsC.surroundingNodes(idim), Ncomp, Nghost, Fab_allocate);

    computeBCoeffs(bcoefs, idim, kappa_r, kcomp, Erborder, 0,
                   c, limiter, parent->Geom(level));

    // this routine is called for both level and multilevel solves,
    // so we must install coefs in the right place:

    if (hd || hm) {
      if (hd) {
	hd->bCoefficients(bcoefs, idim);
      }
      else if (hm) {
	hm->bCoefficients(level, bcoefs, idim);
      }
    }
    else if (solver) {
      solver->bCoefficients(level - base, bcoefs, idim);
    }

  } // -->> over dimension
}

void RadSolve::levelDCoeffs(int level, Tuple<MultiFab, BL_SPACEDIM>& lambda,
			    MultiFab& vel, MultiFab& dcf)
{
  BL_PROFILE("RadSolve::levelDCoeffs");
  const BoxArray& grids = parent->boxArray(level);

  Array<Real> r, s;

  for (int idim=0; idim<BL_SPACEDIM; idim++) {

    BoxArray edge_boxes(grids);
    edge_boxes.surroundingNodes(idim);
    MultiFab dcoefs(edge_boxes, 1, 0, Fab_allocate);

    for (MFIter mfi(dcoefs); mfi.isValid(); ++mfi) {
      int i = mfi.index();
      const Box &reg = grids[i];

      // metric terms
      const int I = (BL_SPACEDIM >= 2) ? 1 : 0;
      if (Geometry::IsCartesian()) {
	r.resize(reg.length(0)+1, 1);
      }
      else if (Geometry::IsRZ()) {
	if (idim == 0) {
	  parent->Geom(level).GetEdgeLoc(r, reg, 0);
	}
	else {
	  parent->Geom(level).GetCellLoc(r, reg, 0);
	}
      }
      else {
	parent->Geom(level).GetEdgeLoc(r, reg, 0);
	parent->Geom(level).GetCellLoc(s, reg, I);
	const Real *dx = parent->Geom(level).CellSize();
	const Box &dbox = dcoefs[i].box();
	FORT_SPHE(r.dataPtr(), s.dataPtr(), idim, dimlist(dbox), dx);
      }

      BL_FORT_PROC_CALL(CA_COMPUTE_DCOEFS, ca_compute_dcoefs)
    	(BL_TO_FORTRAN(dcoefs[i]), BL_TO_FORTRAN(lambda[idim][i]),
    	 BL_TO_FORTRAN(vel[i]), BL_TO_FORTRAN(dcf[i]), 
	 r.dataPtr(), &idim);
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
  const BoxArray& grids = parent->boxArray(level);

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

  Array<Real> r, s, rf;

  for (MFIter ri(rhs); ri.isValid(); ++ri) {
    int i = ri.index();
    const Box &rbox = rhs[i].box();
    const Box &ebox = Er_old[i].box();
    const Box &reg  = grids[i];
    const Real *dx  = parent->Geom(level).CellSize();

    const int I = (BL_SPACEDIM >= 2) ? 1 : 0;
    if (Geometry::IsCartesian()) {
      r.resize(reg.length(0), 1);
      s.resize(reg.length(I), 1);
    }
    else if (Geometry::IsRZ()) {
      parent->Geom(level).GetCellLoc(r, reg, 0);
      s.resize(reg.length(I), 1);
    }
    else {
      parent->Geom(level).GetCellLoc(r, reg, 0);
      parent->Geom(level).GetCellLoc(s, reg, I);
      FORT_SPHC(r.dataPtr(), s.dataPtr(), dimlist(reg), dx);
    }

    FORT_LRHS(rhs[i].dataPtr(), dimlist(rbox), dimlist(reg),
	      temp[i].dataPtr(),
	      fkp[i].dataPtr(), eta[i].dataPtr(), etainv[i].dataPtr(),
	      rhoem[i].dataPtr(), rhoes[i].dataPtr(),
	      dflux_old[i].dataPtr(),
	      Er_old[i].dataPtr(0), dimlist(ebox),
	      Edot[i].dataPtr(),
	      r.dataPtr(), s.dataPtr(), delta_t, sigma, c, theta);
  }
}


void RadSolve::levelSolve(int level,
                          MultiFab& Er, int igroup, MultiFab& rhs,
                          Real sync_absres_factor)
{
  BL_PROFILE("RadSolve::levelSolve");
  const BoxArray& grids = parent->boxArray(level);

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
  const BoxArray& grids = parent->boxArray(level);

  int nflx = flx.nComp();

  Array<Real> r, s;
  const Geometry& geom = parent->Geom(level);
  const Real *dx = geom.CellSize();

  for (int idim = 0; idim < BL_SPACEDIM; idim++) {
    for (MFIter mfi(flx); mfi.isValid(); ++mfi) {
      int i = mfi.index();

      const Box &Fbox = Flux[idim][i].box();
      const int* Flo = Fbox.loVect();
      const int* Fhi = Fbox.hiVect();

      const Box &reg  = grids[i];
      const int* reglo = reg.loVect();
      const int* reghi = reg.hiVect();

      int rlo, rhi;

      const int I = (BL_SPACEDIM >= 2) ? 1 : 0;
      if (Geometry::IsCartesian()) {
	r.resize(reg.length(0)+1, 1);
	rlo = Flo[0];
	rhi = Fhi[0];
      }
      else if (Geometry::IsRZ()) {
	if (idim == 0) {
	  geom.GetEdgeLoc(r, reg, 0);
	  rlo = Flo[0];
	  rhi = Fhi[0];
	}
	else {
	  geom.GetCellLoc(r, reg, 0);
	  rlo = reglo[0];
	  rhi = reghi[0];
	}
      }
      else {
	geom.GetEdgeLoc(r, reg, 0);
	geom.GetCellLoc(s, reg, I);
	rlo = Flo[0];
	rhi = Fhi[0];

	FORT_SPHE(r.dataPtr(), s.dataPtr(), idim, dimlist(Fbox), dx);
      }

      BL_FORT_PROC_CALL(CA_TEST_TYPE_FLUX, ca_test_type_flux)
	(BL_TO_FORTRAN(flx[i]),
	 BL_TO_FORTRAN(Flux[idim][i]),
	 r.dataPtr(), &rlo, &rhi, 
	 &nflx, &idim, &igroup);
    }
  }
}

void RadSolve::levelFluxFaceToCenter(int level, MultiFab& state,
				     Tuple<MultiFab, BL_SPACEDIM>& lambda,
				     MultiFab& Er, Tuple<MultiFab, BL_SPACEDIM>& Flux,
				     MultiFab& flx)
{
  const BoxArray& grids = parent->boxArray(level);

  int nflx = flx.nComp();

  Array<Real> r, s;
  const Geometry& geom = parent->Geom(level);
  const Real *dx = geom.CellSize();

  for (int idim = 0; idim < BL_SPACEDIM; idim++) {
    for (MFIter mfi(flx); mfi.isValid(); ++mfi) {
      int i = mfi.index();

      const Box &Fbox = Flux[idim][i].box();
      const int* Flo = Fbox.loVect();
      const int* Fhi = Fbox.hiVect();

      const Box &reg  = grids[i];
      const int* reglo = reg.loVect();
      const int* reghi = reg.hiVect();

      int rlo, rhi;

      const int I = (BL_SPACEDIM >= 2) ? 1 : 0;
      if (Geometry::IsCartesian()) {
	r.resize(reg.length(0)+1, 1);
	rlo = Flo[0];
	rhi = Fhi[0];
      }
      else if (Geometry::IsRZ()) {
	if (idim == 0) {
	  geom.GetEdgeLoc(r, reg, 0);
	  rlo = Flo[0];
	  rhi = Fhi[0];
	}
	else {
	  geom.GetCellLoc(r, reg, 0);
	  rlo = reglo[0];
	  rhi = reghi[0];
	}
      }
      else {
	geom.GetEdgeLoc(r, reg, 0);
	geom.GetCellLoc(s, reg, I);
	rlo = Flo[0];
	rhi = Fhi[0];

	FORT_SPHE(r.dataPtr(), s.dataPtr(), idim, dimlist(Fbox), dx);
      }

      BL_FORT_PROC_CALL(CA_TEST_TYPE_FLUX_LAB, ca_test_type_flux_lab)
	(BL_TO_FORTRAN(flx[i]),
	 BL_TO_FORTRAN(Flux[idim][i]),
	 D_DECL(BL_TO_FORTRAN(lambda[0][i]),
		BL_TO_FORTRAN(lambda[1][i]),
		BL_TO_FORTRAN(lambda[2][i])),
	 BL_TO_FORTRAN(Er[i]),
	 BL_TO_FORTRAN(state[i]),
	 r.dataPtr(), &rlo, &rhi, 
	 &nflx, &idim);
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
  for (MFIter ei(Er); ei.isValid(); ++ei) {
    int i = ei.index();
    Erborder[i].copy(Er[i], igroup, 0, 1);
  }

  if (hd || hm) {
    Erborder.FillBoundary(); // zeroes left in off-level boundaries
    if (parent->Geom(level).isAnyPeriodic()) {
      parent->Geom(level).FillPeriodicBoundary(Erborder, true);
    }
  }

  const Real* dx = parent->Geom(level).CellSize();

  for (int n = 0; n < BL_SPACEDIM; n++) {
    const MultiFab *bp, *cp;
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
    for (MFIter fi(Flux[n]); fi.isValid(); ++fi) {
      int i = fi.index();
      FORT_SET_ABEC_FLUX(&n,
                         Erborder[i].dataPtr(), dimlist(Erborder[i].box()),
                         bcoef[i].dataPtr(),    dimlist(bcoef[i].box()),
                         &beta,
                         dx,
                         Flux[n][i].dataPtr(),  dimlist(Flux[n][i].box()));
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
      for (MFIter fi(Flux[n]); fi.isValid(); ++fi) {
        int i = fi.index();
        flux_out->FineAdd(Flux[n][i], n, i, 0, igroup, 1, scale);
      }
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
  for (MFIter ei(Er); ei.isValid(); ++ei) {
    int i = ei.index();
    Erborder[i].copy(Er[i], igroup, 0, 1);
  }

  Erborder.FillBoundary(); // zeroes left in off-level boundaries
  if (parent->Geom(level).isAnyPeriodic()) {
    parent->Geom(level).FillPeriodicBoundary(Erborder, true);
  }

  HypreExtMultiABec *hem = (HypreExtMultiABec*)hm;

  for (int n = 0; n < BL_SPACEDIM; n++) {
    const MultiFab *dp;

    dp = &hem->d2Coefficients(level, n);
    MultiFab &dcoef = *(MultiFab*)dp;

    for (MFIter fi(dcoef); fi.isValid(); ++fi) {
      int i = fi.index();

      BL_FORT_PROC_CALL(CA_SET_DTERM_FACE, ca_set_dterm_face)
	(BL_TO_FORTRAN(Erborder[i]),
	 BL_TO_FORTRAN(dcoef[i]), 
	 BL_TO_FORTRAN(Dterm_face[n][i]), 
	 dx, &n);
    }
  }

  // Correct D terms at physical and coarse-fine boundaries.
  hem->boundaryDterm(level, &Dterm_face[0], Er, igroup);

  Array<Real> rc, re, s;
  if (Geometry::IsSPHERICAL()) {
    for (MFIter fi(Dterm_face[0]); fi.isValid(); ++fi) {
      int i = fi.index();
      const Box &reg = grids[i];
      parent->Geom(level).GetEdgeLoc(re, reg, 0);
      parent->Geom(level).GetCellLoc(rc, reg, 0);
      parent->Geom(level).GetCellLoc(s, reg, 0);
      const Box &dbox = Dterm_face[0][i].box();
      FORT_SPHE(re.dataPtr(), s.dataPtr(), 0, dimlist(dbox), dx);

      BL_FORT_PROC_CALL(CA_CORRECT_DTERM, ca_correct_dterm)
      (D_DECL(BL_TO_FORTRAN(Dterm_face[0][i]),
	      BL_TO_FORTRAN(Dterm_face[1][i]),
	      BL_TO_FORTRAN(Dterm_face[2][i])),
       re.dataPtr(), rc.dataPtr());
    }
  }
  else if (Geometry::IsRZ()) {
    for (MFIter fi(Dterm_face[0]); fi.isValid(); ++fi) {
      int i = fi.index();
      const Box &reg = grids[i];
      parent->Geom(level).GetEdgeLoc(re, reg, 0);
      parent->Geom(level).GetCellLoc(rc, reg, 0);

      BL_FORT_PROC_CALL(CA_CORRECT_DTERM, ca_correct_dterm)
      (D_DECL(BL_TO_FORTRAN(Dterm_face[0][i]),
	      BL_TO_FORTRAN(Dterm_face[1][i]),
	      BL_TO_FORTRAN(Dterm_face[2][i])),
       re.dataPtr(), rc.dataPtr());
    }
  }

  for (MFIter fi(Dterm); fi.isValid(); ++fi) {
    int i = fi.index();
    BL_FORT_PROC_CALL(CA_FACE2CENTER, ca_face2center)
      (D_DECL(BL_TO_FORTRAN(Dterm_face[0][i]),
	      BL_TO_FORTRAN(Dterm_face[1][i]),
	      BL_TO_FORTRAN(Dterm_face[2][i])),
       BL_TO_FORTRAN(Dterm[i]));
  }
}

void RadSolve::multilevelInit(int crse_level, int fine_level,
                              const BCRec& rad_bc, Real time)
{
  BL_PROFILE("RadSolve::multilevelInit");

  if (use_hypre_multilevel == 1) {
    RadBndry::setTime(time);

    if (use_hypre_nonsymmetric_terms == 0) {
      hm = new HypreMultiABec(crse_level, fine_level,
                              multilevel_solver_flag);
    }
    else {
      hm = new HypreExtMultiABec(crse_level, fine_level,
                                 multilevel_solver_flag);
    }

    // add grids in reverse order in case we are masking covered part of coarse
    for (int level = fine_level; level >= crse_level; level--) {
      IntVect ratio = ((level < fine_level) ? parent->refRatio(level)
                       : IntVect::TheUnitVector());
      hm->addLevel(level, parent->Geom(level),
                   parent->boxArray(level), ratio);
    }
    for (int level = crse_level; level <= fine_level; level++) {
      RadBndry *bdp = new RadBndry(parent->boxArray(level),
                                   parent->Geom(level));
      if (level > crse_level) {
        IntVect rat =  parent->refRatio(level-1);
        bdp->setBndryConds(rad_bc, parent->Geom(level), rat);
        bdp->setBndryFluxConds(rad_bc);
      }
      hm->setBndry(level, *bdp);
    }
    hm->buildMatrixStructure();
    hm->setScalars(alpha, beta);

    return;
  }

  int use_hypre = 1; // this is only option currently supported
  solver = new CompSolver(use_hypre, multilevel_solver_flag,
			  use_harmonic_avg, multilevel_version);

  bbld = new RadBndryBld;

  solver->SetBndryConds(bbld, rad_bc);
  RadBndry::setTime(time);

  base   = crse_level;
  finest = fine_level;

  for (int level = crse_level; level <= fine_level; level++) {
    IntVect ratio = (level == 0) ? IntVect::TheUnitVector()
                                 : parent->refRatio(level-1);
    solver->AddLevel(level - base, parent->Geom(level),
		     parent->boxArray(level), ratio);
  }

  solver->SetScalars(alpha, beta);
}

void RadSolve::multilevelClear()
{
  if (use_hypre_multilevel == 1) {
    for (int level = hm->crseLevel(); level <= hm->fineLevel(); level++) {
      hm->bndryClear(level);
    }
    delete hm;
    hm = NULL;
  }
  else {
    delete solver;
    solver = NULL;
    delete bbld;
  }
}

RadBndry& RadSolve::multilevelCrseBndryData()
{
  if (use_hypre_multilevel == 1) {
    return (RadBndry&)hm->bndryData(hm->crseLevel());
  }
  else {
    return *(RadBndry*)solver->GetCrseBndryData();
  }
}

void RadSolve::multilevelACoeffs(int level, MultiFab& fkp, Real c)
{
  BL_PROFILE("RadSolve::multilevelACoeffs");
  const BoxArray& grids = parent->boxArray(level);

  // Allocate space for ABecLapacian acoeffs, fill with values

  int Ncomp  = 1;
  int Nghost = 0;

  Array<Real> r, s;

  MultiFab acoefs(grids, Ncomp, Nghost, Fab_allocate);

  for (MFIter mfi(fkp); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box& abox = acoefs[i].box();
    const Box& reg  = grids[i];
    const int I = (BL_SPACEDIM >= 2) ? 1 : 0;
    if (Geometry::IsCartesian()) {
      r.resize(reg.length(0), 1);
      s.resize(reg.length(I), 1);
    }
    else if (Geometry::IsRZ()) {
      parent->Geom(level).GetCellLoc(r, reg, 0);
      s.resize(reg.length(I), 1);
    }
    else {
      parent->Geom(level).GetCellLoc(r, reg, 0);
      parent->Geom(level).GetCellLoc(s, reg, I);
      const Real *dx = parent->Geom(level).CellSize();
      FORT_SPHC(r.dataPtr(), s.dataPtr(), dimlist(reg), dx);
    }
    FORT_MACOEF(acoefs[i].dataPtr(), dimlist(abox), dimlist(reg),
		fkp[i].dataPtr(), r.dataPtr(), s.dataPtr(), c);
  }

  if (use_hypre_multilevel == 1) {
    hm->aCoefficients(level, acoefs);
  }
  else {
    solver->aCoefficients(level - base, acoefs);
  }
}

void RadSolve::multilevelRhs(int level, MultiFab& temp,
                             MultiFab& fkp, MultiFab& Ert, Real sigma)
{
  BL_PROFILE("RadSolve::multilevelRhs");
  int i, nf = fkp.nComp();
  const BoxArray& grids = parent->boxArray(level);

  int Ncomp  = 1;
  int Nghost = 0;

  Array<Real> r, s;

  MultiFab rhs(grids, Ncomp, Nghost, Fab_allocate);

  for (MFIter mfi(fkp); mfi.isValid(); ++mfi) {
    i = mfi.index();
    const Box &rbox = rhs[i].box();
    const Box &reg  = grids[i];
    const int I = (BL_SPACEDIM >= 2) ? 1 : 0;
    if (Geometry::IsCartesian()) {
      r.resize(reg.length(0), 1);
      s.resize(reg.length(I), 1);
    }
    else if (Geometry::IsRZ()) {
      parent->Geom(level).GetCellLoc(r, reg, 0);
      s.resize(reg.length(I), 1);
    }
    else {
      parent->Geom(level).GetCellLoc(r, reg, 0);
      parent->Geom(level).GetCellLoc(s, reg, I);
      const Real *dx = parent->Geom(level).CellSize();
      FORT_SPHC(r.dataPtr(), s.dataPtr(), dimlist(reg), dx);
    }
    FORT_MRHS(rhs[i].dataPtr(), dimlist(rbox), dimlist(reg),
	      temp[i].dataPtr(), fkp[i].dataPtr(),
	      Ert[i].dataPtr(), r.dataPtr(), s.dataPtr(),
	      sigma);
  }

  if (use_hypre_multilevel == 1) {
    hm->loadLevelVectorB(level, rhs, Inhomogeneous_BC);
  }
  else {
    solver->SetRhs(level - base, rhs);
  }
}

void RadSolve::multilevelGuess(int level, MultiFab& guess)
{
  if (use_hypre_multilevel == 1) {
    hm->loadLevelVectorX(level, guess, 0);
  }
  else {
    solver->SetInitialGuess(level - base, guess);
  }
}

void RadSolve::multilevelSolve(int crse_level,
                               int fine_active_level,
                               int is_sync)
{
  BL_PROFILE("RadSolve::multilevelSolve");

  if (use_hypre_multilevel == 1) {

    // Hypre multilevel solver in normal operation:

    hm->loadMatrix();
    hm->finalizeMatrix();
    hm->finalizeVectors();
    if (is_sync) {
      Real abstoltmp = 0.0;
      for (int lev = crse_level; lev <= fine_active_level; lev++) {
        abstoltmp = (abstoltmp > absres[lev]) ? abstoltmp : absres[lev];
        if (crse_level == 0)
          absres[lev] = 0.0;
      }
      if (verbose > 0 && ParallelDescriptor::IOProcessor()) {
        int oldprec = cout.precision(20);
        cout << "Using absolute tolerance = " << abstoltmp << endl;
        cout.precision(oldprec);
      }

      BL_ASSERT(crse_level == hm->crseLevel());
      // We assert this because it is all that is currently supported:
      BL_ASSERT(fine_active_level == hm->fineLevel());

      hm->setupSolver(reltol, abstoltmp, maxiter);
    }
    else {

      // The intent is that with hypre the "secondary solver"
      // business will not be necessary, so this should just work.

      hm->setupSolver(reltol, abstol, maxiter);
    }

    hm->setVerbose(verbose);
    hm->solve();
    hm->clearSolver();
    return;
  }

  // We're not using hypre multilevel, so use CompSolver:

  Real BottomTol    = bottomtol;
  int BottomMaxIter = maxiter;
  int BottomNumIter = bottomnumiter;

  if (finest == crse_level || fine_active_level == crse_level)
    BottomTol = reltol;

  solver->SetBottomParams(BottomTol, BottomMaxIter, BottomNumIter);

  int cycle_type = 1;    // 1 = V-cycle, 2 = W-cycle
  // these two now set in inputs file:
  //int presmooth = 1;     // smoothing sweeps before coarse correction
  //int postsmooth = 1;    // smoothing sweeps after coarse correction
  int interp_order = 0;  // 0 = pw-constant, 1 = pw-BL_SPACEDIM-linear
  int restrict_order = 0;// 0 = pw-constant, 1 = pw-BL_SPACEDIM-linear
  solver->SetParms(cycle_type, presmooth, postsmooth,
		   interp_order, restrict_order);

  BL_ASSERT(crse_level >= base);

  if (is_sync) {
    Real abstoltmp = 0.0;
    for (int lev = crse_level; lev <= fine_active_level; lev++) {
      abstoltmp = (abstoltmp > absres[lev]) ? abstoltmp : absres[lev];
      if (crse_level == 0)
	absres[lev] = 0.0;
    }

    if (verbose > 0 && ParallelDescriptor::IOProcessor()) {
      int oldprec = cout.precision(20);
      cout << "Using absolute tolerance = " << abstoltmp << endl;
      cout.precision(oldprec);
    }

    // This is the original tolerance used in the sync for years:
    //solver->Solve(reltol, 0.0, maxiter, verbose,
    //              crse_level - base, fine_active_level - base);

    // Pure abstol version, usually faster, in some cases too restrictive
    //solver->Solve(0.0, abstoltmp, maxiter, verbose,
    //		  crse_level - base, fine_active_level - base);

    // Compromise version, less restrictive (faster) than either of the above
    solver->Solve(cs_reltol_mult * reltol, abstoltmp, maxiter, verbose,
                  crse_level - base, fine_active_level - base);
  }
  else {
    solver->BuildSecondarySolvers(crse_level - base,
				  fine_active_level - base);
    solver->Solve(cs_reltol_mult * reltol, abstol, maxiter, verbose,
		  crse_level - base, fine_active_level - base);
    solver->SecondarySolve(secondtol, maxiter, verbose,
			   crse_level - base, fine_active_level - base);
  }
}

void RadSolve::multilevelSolution(int level, MultiFab& solution)
{
  BL_PROFILE("RadSolve::multilevelSolution");
  if (use_hypre_multilevel == 1) {
    hm->getSolution(level, solution, 0);
  }
  else {
    solver->GetSolution(level - base, solution);
  }
}

void RadSolve::multilevelFineCorrection(int level, MultiFab& diff)
{
  BL_PROFILE("RadSolve::multilevelCorrection");
  if (use_hypre_multilevel == 1) {
    cout << "Warning: multilevelFineCorrection not implemented" << endl;
    cout << "         for radsolve.use_hypre_multilevel = 1" << endl;
  }
  else {
    solver->GetFineCorrection(level - base, diff);
  }
}

void RadSolve::multilevelFlux(int level, MultiFab& dflux,
                              FluxRegister* flux_in, FluxRegister* flux_out,
                              int igroup, const BCRec& rad_bc)
{
  BL_PROFILE("RadSolve::multilevelFlux");
  const BoxArray& grids = parent->boxArray(level);

  if (use_hypre_multilevel == 1) {

    // This still doesn't do dflux, but seems to handle both flux_in
    // and flux_out now.  The part that looks at the coarse level
    // solution is ugly and not tested.
    // Values in this case would only matter if theta < 1.

    MultiFab Er(grids, 1, 0);
    multilevelSolution(level, Er);

    if (level > hm->crseLevel()) {
      // This section seems to work but has not been debugged.
      MultiFab Er_crse(parent->boxArray(level-1), 1, 0);
      multilevelSolution(level-1, Er_crse);
      BoxArray cgrids(grids);
      IntVect crse_ratio = parent->refRatio(level-1);
      cgrids.coarsen(crse_ratio);
      BndryRegister crse_br(cgrids, 0, 1, 1, 1);
      crse_br.setVal(1.0e30);
      BL_ASSERT(!parent->Geom(level).isAnyPeriodic());
      crse_br.copyFrom(Er_crse, 0, 0, 0, 1);

      RadBndry *bdp = (RadBndry*) &hm->bndryData(level);
      bdp->setBndryValues(crse_br, 0, Er, 0, 0, 1, crse_ratio, rad_bc);
      bdp->setBndryFluxConds(rad_bc);
    }

    levelFlux(level, flux_in, flux_out, Er, igroup);
  }
  else {

    // dflux gets only the flux portion of the operator, so turn off absorption
    solver->SetScalars(0.0, 1.0);
    solver->GetFlux(level - base, dflux, flux_in, flux_out);
    solver->SetScalars(alpha, beta);

    if (Geometry::IsRZ()) {
      Array<Real> r;
      for (MFIter di(dflux); di.isValid(); ++di) {
        int i = di.index();
        const Box &dbox = dflux[i].box();
        const Box &reg  = grids[i];
        parent->Geom(level).GetCellLoc(r, dbox, 0);
        FORT_DIVR(dflux[i].dataPtr(), dimlist(dbox), dimlist(reg),
                  r.dataPtr());
      }
    }
    else if (Geometry::IsSPHERICAL()) {
      const int I = (BL_SPACEDIM >= 2) ? 1 : 0;
      Array<Real> r, s;
      for (MFIter di(dflux); di.isValid(); ++di) {
        int i = di.index();
        const Box &dbox = dflux[i].box();
        const Box &reg  = grids[i];
        parent->Geom(level).GetCellLoc(r, dbox, 0);
        parent->Geom(level).GetCellLoc(s, dbox, I);
        const Real *dx = parent->Geom(level).CellSize();
        FORT_SPHC(r.dataPtr(), s.dataPtr(), dimlist(dbox), dx);
        FORT_DIVRS(dflux[i].dataPtr(), dimlist(dbox), dimlist(reg),
                   r.dataPtr(), s.dataPtr());
      }
    }
  }
}

void RadSolve::syncACoeffs(int level,
                           MultiFab& fkp, MultiFab& etainv,
                           Real c, Real delta_t)
{
  BL_PROFILE("RadSolve::syncACoeffs");
  const BoxArray& grids = parent->boxArray(level);

  // Allocate space for ABecLapacian acoeffs, fill with values

  int Ncomp  = 1;
  int Nghost = 0;

  Array<Real> r, s;

  MultiFab acoefs(grids, Ncomp, Nghost, Fab_allocate);

  for (MFIter mfi(fkp); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box &abox = acoefs[i].box();
    const Box &reg  = grids[i];
    const int I = (BL_SPACEDIM >= 2) ? 1 : 0;
    if (Geometry::IsCartesian()) {
      r.resize(reg.length(0), 1);
      s.resize(reg.length(I), 1);
    }
    else if (Geometry::IsRZ()) {
      parent->Geom(level).GetCellLoc(r, reg, 0);
      s.resize(reg.length(I), 1);
    }
    else {
      parent->Geom(level).GetCellLoc(r, reg, 0);
      parent->Geom(level).GetCellLoc(s, reg, I);
      const Real *dx = parent->Geom(level).CellSize();
      FORT_SPHC(r.dataPtr(), s.dataPtr(), dimlist(reg), dx);
    }
    FORT_SACOEF(acoefs[i].dataPtr(), dimlist(abox), dimlist(reg),
		fkp[i].dataPtr(), etainv[i].dataPtr(),
		r.dataPtr(), s.dataPtr(), c, delta_t);
  }

  if (use_hypre_multilevel == 1) {
    hm->aCoefficients(level, acoefs);
  }
  else {
    solver->aCoefficients(level - base, acoefs);
  }
}


void RadSolve::syncRhs(int level, FluxRegister* sync_flux,
                       Real scale, int igroup)
{
  BL_PROFILE("RadSolve::syncRhs");
  const BoxArray& grids = parent->boxArray(level);

  int Ncomp  = 1;
  int Nghost = 0;

  MultiFab rhs(grids, Ncomp, Nghost, Fab_allocate);
  rhs.setVal(0.0);

  if (sync_flux) {
    // This works directly for a multilevel solve since the finer level
    // is present in the rhs to overwrite the junk produced under it.
    // In the single-level version we have to be sure to clean up
    // fine-fine interfaces manually.  This is currently done by a call
    // to clear_internal_borders in Radiation::sync_solve.

    sync_flux->Reflux(rhs, scale, igroup, 0, 1, parent->Geom(level));
  }

  if (use_hypre_multilevel == 1) {
    hm->loadLevelVectorB(level, rhs, Homogeneous_BC);
  }
  else {
    solver->SetRhs(level - base, rhs);
  }
}


// <MGFLD routines>
void RadSolve::computeBCoeffs(MultiFab& bcoefs, int idim,
                              MultiFab& kappa_r, int kcomp,
                              MultiFab& lambda, int lamcomp,
                              Real c, const Geometry& geom)
{
  BL_PROFILE("RadSolve::computeBCoeffs (MGFLD)");
  const BoxArray& grids = kappa_r.boxArray(); // valid region only

  BL_ASSERT(kappa_r.nGrow() == 1);

  Array<Real> r, s;

  const Real* dx       = geom.CellSize();

  for (MFIter mfi(lambda); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box &bbox = lambda[i].box();
    const Box &reg  = grids[i];
    const Box &kbox = kappa_r[i].box();

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
    else { // support only 1D spherical here
      geom.GetEdgeLoc(r, reg, 0);
      geom.GetCellLoc(s, reg, I);
      FORT_SPHE(r.dataPtr(), s.dataPtr(), idim, dimlist(bbox), dx);
    }
    
    FORT_BCLIM(bcoefs[i].dataPtr(), lambda[i].dataPtr(lamcomp),
	       dimlist(bbox), dimlist(reg),
	       idim, kappa_r[i].dataPtr(kcomp), dimlist(kbox),
	       r.dataPtr(), s.dataPtr(), c, dx);
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

  Array<Real> r, s;

  MultiFab acoefs(grids, Ncomp, Nghost, Fab_allocate);

  for (MFIter mfi(kpp); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box &abox = acoefs[i].box();
    const Box &reg = grids[i];
    const int I = (BL_SPACEDIM >= 2) ? 1 : 0;
    if (Geometry::IsCartesian()) {
      r.resize(reg.length(0), 1);
      s.resize(reg.length(I), 1);
    }
    else if (Geometry::IsRZ()) {
      parent->Geom(level).GetCellLoc(r, reg, 0);
      s.resize(reg.length(I), 1);
    }
    else {
      parent->Geom(level).GetCellLoc(r, reg, 0);
      parent->Geom(level).GetCellLoc(s, reg, I);
      const Real *dx = parent->Geom(level).CellSize();
      FORT_SPHC(r.dataPtr(), s.dataPtr(), dimlist(reg), dx);
    }

    const Box &kbox = kpp[i].box();
    Real dt_ptc = delta_t/(1.0+ptc_tau);
    FORT_LACOEFMGFLD(acoefs[i].dataPtr(), dimlist(abox), dimlist(reg),
		     kpp[i].dataPtr(igroup), dimlist(kbox),
		     r.dataPtr(), s.dataPtr(), dt_ptc, c);
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
  const BoxArray& grids = parent->boxArray(level);
  Castro *castro = dynamic_cast<Castro*>(&parent->getLevel(level));
  Real time = castro->get_state_data(Rad_Type).curTime();
  
  Array<Real> r, s;

  for (MFIter ri(rhs); ri.isValid(); ++ri) {
    int i = ri.index();
    const Box &reg = grids[i];

#ifdef MG_SU_OLSON

    parent->Geom(level).GetCellLoc(r, reg, 0);

    BL_FORT_PROC_CALL(CA_COMPUTE_RHS_SO, ca_compute_rhs_so)
      (BL_TO_FORTRAN(rhs[i]),
       BL_TO_FORTRAN(jg[i]),
       BL_TO_FORTRAN(mugT[i]),
       BL_TO_FORTRAN(coupT[i]),
       BL_TO_FORTRAN(etaT[i]),
       BL_TO_FORTRAN(Er_step[i]),
       BL_TO_FORTRAN(rhoe_step[i]),
       BL_TO_FORTRAN(rhoe_star[i]),
       r.dataPtr(), 
       &time, &delta_t, &igroup);

#else

    const int I = (BL_SPACEDIM >= 2) ? 1 : 0;
    if (CoordSys::IsCartesian()) {
      r.resize(reg.length(0), 1);
      s.resize(reg.length(I), 1);
    }
    else if (CoordSys::IsRZ()) {
      parent->Geom(level).GetCellLoc(r, reg, 0);
      s.resize(reg.length(I), 1);
    }
    else {
      parent->Geom(level).GetCellLoc(r, reg, 0);
      parent->Geom(level).GetCellLoc(s, reg, I);
      const Real *dx = parent->Geom(level).CellSize();
      FORT_SPHC(r.dataPtr(), s.dataPtr(), dimlist(reg), dx);
    }

#ifdef NEUTRINO
    BL_FORT_PROC_CALL(CA_COMPUTE_RHS_NEUT, ca_compute_rhs_neut)
      (BL_TO_FORTRAN(rhs[i]),
       BL_TO_FORTRAN(jg[i]),
       BL_TO_FORTRAN(mugT[i]),
       BL_TO_FORTRAN(mugY[i]),
       BL_TO_FORTRAN(coupT[i]),
       BL_TO_FORTRAN(coupY[i]),
       BL_TO_FORTRAN(etaT[i]),
       BL_TO_FORTRAN(etaY[i]),
       BL_TO_FORTRAN(thetaT[i]),
       BL_TO_FORTRAN(thetaY[i]),
       BL_TO_FORTRAN(Er_step[i]),
       BL_TO_FORTRAN(rhoe_step[i]),
       BL_TO_FORTRAN(rhoYe_step[i]),
       BL_TO_FORTRAN(Er_star[i]),
       BL_TO_FORTRAN(rhoe_star[i]),
       BL_TO_FORTRAN(rhoYe_star[i]),
       r.dataPtr(), 
       &delta_t, &igroup, &ptc_tau);
#else
    BL_FORT_PROC_CALL(CA_COMPUTE_RHS, ca_compute_rhs)
      (BL_TO_FORTRAN(rhs[i]),
       BL_TO_FORTRAN(jg[i]),
       BL_TO_FORTRAN(mugT[i]),
       BL_TO_FORTRAN(coupT[i]),
       BL_TO_FORTRAN(etaT[i]),
       BL_TO_FORTRAN(Er_step[i]),
       BL_TO_FORTRAN(rhoe_step[i]),
       BL_TO_FORTRAN(Er_star[i]),
       BL_TO_FORTRAN(rhoe_star[i]),
       r.dataPtr(), 
       &delta_t, &igroup, &ptc_tau);
#endif
#endif
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
