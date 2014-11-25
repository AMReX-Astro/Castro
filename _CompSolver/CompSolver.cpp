#include <EdgeVar.H>
#include <CompSolver.H>

#include <Using.H>

CompSolver::CompSolver(int _use_hypre, int _solverflag,
		       int _use_harmonic_avg,
		       int _version)
  :CompLevel(PArrayManage),
   use_hypre(_use_hypre),
   solverflag(_solverflag),
   use_harmonic_avg(_use_harmonic_avg),
   version(_version),
   alpha(0.),
   beta(0.),
   ValidCoefAmrLevel(4096)
{
   if (solverflag > 200) {
      solverflag0= 103; // pcg with amg on coarsest
      solverflagi= 3;   // pcg with pfmg finer
   }
   else {
      // input solver for all levels
      solverflag0= solverflag; 
      solverflagi= solverflag; 
   }
}

CompSolver::~CompSolver()
{
}

void
CompSolver::SetBndryConds( const NGBndryBld * _BndryBld,
			   const BCRec & _PhysBcr )
{
  for( int n=0; n<BL_SPACEDIM; n++ ) {
    PhysBcr.setLo( n, _PhysBcr.lo(n) );
    PhysBcr.setHi( n, _PhysBcr.hi(n) );
  }
  BndryBld = _BndryBld;
}

void
CompSolver::SetParms( int cycle_type, int presmooth, int postsmooth,
		      int interp_order, int restrict_order )
{
  if( interp_order + restrict_order > 1 ) {
    BoxLib::Error("interp_order + restrict_order > 1");
  }

  for( int level=0; level<CompLevel.size(); level++ ) {
    CompLevel[level].SetParms( cycle_type, presmooth, postsmooth,
			       interp_order, restrict_order );
  }
}

static int log2local(int i)
{
  if (i <= 1)
    return 0;
  else
    return 1 + log2local(i/2);
}

void
CompSolver::AddLevel(int                   AmrLevel,
		     const Geometry &       AmrGeom,
		     const BoxArray &      AmrGrids,
		     IntVect           AmrCrseRatio)
{

  int OldLength = CompLevel.size();

  // If an existing level is being added again, first throw away
  // the existing level and all finer levels
  if( AmrLevel < AmrLevelIndex.size() ) {
    AmrLevelIndex.resize( AmrLevel );
    OldLength = (AmrLevel>0)? AmrLevelIndex[AmrLevel-1]+1: 0;
    CompLevel.resize( OldLength );
    // finer interpolated levels may exist without coefficients:
    //ValidCoefAmrLevel = Max( ValidCoefAmrLevel, AmrLevel+1 );
  }

  ValidCoefAmrLevel = min( ValidCoefAmrLevel, AmrLevel );

  int NumNewLevels;
  CompSolverLevel * CoarserAmrLevel;

  if( AmrLevel > 0 ) {
    int MaxRatio = 0;
    for( int n=0; n<BL_SPACEDIM; n++ ) {
      if( AmrCrseRatio[n] > MaxRatio ) MaxRatio = AmrCrseRatio[n];
    }
    NumNewLevels = log2local(MaxRatio);
    CoarserAmrLevel = &CompLevel[AmrLevelIndex[AmrLevel-1]];
  }
  else {
    NumNewLevels = 1;
    AmrCrseRatio = IntVect::TheUnitVector();
    CoarserAmrLevel = 0;
  }

  int NewLength = OldLength + NumNewLevels;

  AmrLevelIndex.resize(AmrLevelIndex.size() + 1);
  AmrLevelIndex[AmrLevel] = NewLength-1;

  CompLevel.resize(NewLength);

  BoxArray ba(AmrGrids);

  IntVect TwoVect(D_DECL(2,2,2));

  IntVect IntermediateCrseRatio = AmrCrseRatio;

  IntermediateCrseRatio.min(TwoVect);

  // if different solver on the coarsest grid, input a new solverflag
  int solver_type;
  solver_type= solverflagi;
  if (AmrLevel == 0)
  {
     solver_type= solverflag0;
  }

  // Add new AMR level
     CompLevel.set( NewLength-1, new CompSolverLevel(ba, IntermediateCrseRatio,
  	                                             CoarserAmrLevel,
						     AmrCrseRatio,
						     AmrGeom,
						     *BndryBld, PhysBcr,
						     true,
						     use_hypre,
						     solver_type,
						     use_harmonic_avg,
						     version) );

  if( AmrLevel > 0 ) {
    CoarserAmrLevel->SetFinerAmrLevel(&CompLevel[NewLength-1]);
  }

  // Create intermediate levels
  if( NumNewLevels > 1 && AmrLevel > 0 ) {

    for( int i=NumNewLevels-2; i>=0; i-- ) {
      int level = OldLength + i;

      // Make geometry for this level
      Box AmrDomain = AmrGeom.Domain();
      AmrDomain.coarsen(IntermediateCrseRatio);
      Geometry IntermediateGeom(AmrDomain);

      ba.coarsen(IntermediateCrseRatio);

      AmrCrseRatio /= 2;
      AmrCrseRatio.max(IntVect::TheUnitVector());

      // Compute ratio to the next coarser level
      IntermediateCrseRatio = AmrCrseRatio;
      IntermediateCrseRatio.min(TwoVect);

      CompLevel.set( level, new CompSolverLevel(ba, IntermediateCrseRatio,
						CoarserAmrLevel,
						AmrCrseRatio,
						IntermediateGeom,
						*BndryBld, PhysBcr,
						false,
						use_hypre,
						solver_type,  // note using solver_type
						use_harmonic_avg,
						version) );
      CompLevel[level].SetFinerAmrLevel(&CompLevel[NewLength-1]);
    }
  }

  int Start = max(OldLength,1);
  for( int level = Start; level<NewLength; level++ ) {
    CompLevel[level].SetNextCoarserLevel( &CompLevel[level-1] );
  }
}

void
CompSolver::Clear(void)
{
  CompLevel.clear();
  AmrLevelIndex.clear();
  alpha = 0.;
  beta = 0.;
  MLPrecondTol = -1.;
  MLPrecondMaxIter = -1;
}

void
CompSolver::SetScalars(Real _alpha,
		       Real  _beta)
{
  alpha = _alpha;
  beta  = _beta;

  for( int level=0; level<CompLevel.size(); level++ ) {
    CompLevel[level].SetScalars( alpha, beta );
  }
}

void
CompSolver::SetCoefficients(int           AmrLevel,
			    const MultiFab & acoef,
			    const MultiFab * bcoef)
{
  CompLevel[AmrLevelIndex[AmrLevel]].SetCoefficients(acoef, bcoef );
  ValidCoefAmrLevel = max( ValidCoefAmrLevel, AmrLevel );
}

void
CompSolver::aCoefficients(int           AmrLevel,
                          const MultiFab & acoef)
{
  CompLevel[AmrLevelIndex[AmrLevel]].aCoefficients(acoef);
  ValidCoefAmrLevel = max( ValidCoefAmrLevel, AmrLevel );
}

void
CompSolver::bCoefficients(int           AmrLevel,
                          const MultiFab & bcoef,
			  int                dir)
{
  CompLevel[AmrLevelIndex[AmrLevel]].bCoefficients( bcoef, dir );
  ValidCoefAmrLevel = max( ValidCoefAmrLevel, AmrLevel );
}

void
CompSolver::SetRhs(int             AmrLevel,
		   const MultiFab & Amr_rhs)
{
  MultiFab & rhs = CompLevel[AmrLevelIndex[AmrLevel]].Rhs();

  for(MFIter mfi(rhs); mfi.isValid(); ++mfi) {
    rhs[mfi].copy(Amr_rhs[mfi]);
  }
}

void
CompSolver::SetZeroInitialGuess( )
{
  for (int level = 0; level < CompLevel.size(); level++) {
    CompLevel[level].Solution().setVal(0.0);
  }
}

void
CompSolver::SetInitialGuess(int              AmrLevel,
			    const MultiFab & AmrGuess)
{
  MultiFab & solution = CompLevel[AmrLevelIndex[AmrLevel]].Solution();

  for(MFIter mfi(solution); mfi.isValid(); ++mfi) {
    if (solution.nGrow() > AmrGuess.nGrow()) {
      solution[mfi].setVal(0.0);
    }
    solution[mfi].copy(AmrGuess[mfi]);
  }

  // Average down initial guess to intermediate levels
  if( AmrLevel > 0 ) {
    int CurrentAmrIndex = AmrLevelIndex[AmrLevel];
    int NextAmrIndex = AmrLevelIndex[AmrLevel-1];
    for( int level = CurrentAmrIndex-1; level>NextAmrIndex; level-- ) {
      CompSolverLevel & Current = CompLevel[level];
      CompSolverLevel & Fine = CompLevel[level+1];
      IntVect Ratio = Fine.getCrseRatio();

      (Current.Solution()).setVal(0.);
      CellAvgDown( Fine.Solution(), Current.Solution(), Ratio );
    }
  }

}

void
CompSolver::GetSolution(int           AmrLevel,
			MultiFab & AmrSolution)
{
  int level = AmrLevelIndex[AmrLevel];
  MultiFab & solution = CompLevel[level].Solution();

  for(MFIter mfi(solution); mfi.isValid(); ++mfi) {
    AmrSolution[mfi].copy(solution[mfi]);
  }
}

void
CompSolver::GetFineCorrection(int       AmrLevel,
			      MultiFab & AmrDiff)
{
  int CurrentAmrIndex = AmrLevelIndex[AmrLevel];
  MultiFab & solution = CompLevel[CurrentAmrIndex].Solution();

  for(MFIter mfi(solution); mfi.isValid(); ++mfi) {
    AmrDiff[mfi].copy(solution[mfi]);
  }

  if (CurrentAmrIndex < CompLevel.size() - 1) {
    int FineAmrIndex = AmrLevelIndex[AmrLevel+1];
    CompSolverLevel & Fine = CompLevel[FineAmrIndex];
    IntVect Ratio = Fine.getAmrCrseRatio();
    CellAvgDown(Fine.Solution(), AmrDiff, Ratio);
  }
  else {
    BoxLib::Error("GetFineCorrection called at the finest level");
  }

  for(MFIter mfi(solution); mfi.isValid(); ++mfi) {
    AmrDiff[mfi].minus(solution[mfi]);
  }
}

void
CompSolver::GetFlux(int AmrLevel, MultiFab& dflux,
		    FluxRegister* flux_in, FluxRegister* flux_out)
{
  int level = AmrLevelIndex[AmrLevel];
  CompLevel[level].GetFlux(dflux, flux_in, flux_out);
}

void
CompSolver::AvgDownCoefs(int AmrBaseLevel)
{
  //int finest_level = CompLevel.size() - 1;
  int finest_level       = AmrLevelIndex[ValidCoefAmrLevel];
  int CompositeBaseLevel = AmrLevelIndex[AmrBaseLevel];

  for( int level=finest_level-1; level>=CompositeBaseLevel; level-- ) {

    CompSolverLevel & Current = CompLevel[level];

    const BoxArray & grids = Current.Grids();

    MultiFab acoef(grids,1,0);

    MultiFab bcoef[BL_SPACEDIM];
    for( int n=0; n<BL_SPACEDIM; n++ ) {
      BoxArray ba(grids);
      bcoef[n].define( ba.surroundingNodes(n), 1, 0, Fab_allocate );
    }

    if( Current.AmrLevel() ) {

      const MultiFab & Acoef = Current.Acoef();
      for(MFIter mfi(acoef); mfi.isValid(); ++mfi) {
	acoef[mfi].copy(Acoef[mfi]);
      }

      for( int n=0; n<BL_SPACEDIM; n++ ) {
	const MultiFab & Bcoef = Current.Bcoef(n);
	for(MFIter mfi(bcoef[n]); mfi.isValid(); ++mfi) {
	  bcoef[n][mfi].copy(Bcoef[mfi]);
	}
      }
    }
    else {

      acoef.setVal(0.);

      for( int n=0; n<BL_SPACEDIM; n++ ) {
	bcoef[n].setVal(0.);
      }
    }

    CompSolverLevel & Fine = CompLevel[level+1];
    IntVect Ratio = Fine.getCrseRatio();

    // Average down acoef
    const MultiFab & FineAcoef = Fine.Acoef();
    CellAvgDown( FineAcoef, acoef, Ratio );

    // Average down bcoefs
    for( int n=0; n<BL_SPACEDIM; n++ ) {
      const MultiFab & FineBcoef = Fine.Bcoef(n);
      if(use_harmonic_avg == 1) {
        EdgeHarmAvg( n, FineBcoef, bcoef[n], Ratio );
      } else {
        EdgeAvgDown( n, FineBcoef, bcoef[n], Ratio );
      }
    }

    Current.SetCoefficients( acoef, bcoef );
  }

  ValidCoefAmrLevel = AmrBaseLevel;
}

void
CompSolver::EdgeAvgDown(int               dir,
			const MultiFab & Fine,
			MultiFab &       Crse,
			IntVect &        nref)
{
  BoxArray crse_box( Fine.boxArray() ) ;
  crse_box.enclosedCells(dir);
  crse_box.coarsen(nref) ;
  crse_box.surroundingNodes(dir);
  MultiFab Coarse(crse_box,1,0) ;

  for (MFIter mfi(Coarse); mfi.isValid(); ++mfi) {
    FORT_EDGE_AVG_DOWN(&dir,
		       Fine[mfi].dataPtr(),   dimlist(Fine[mfi].box()),
		       nref.getVect(),
		       Coarse[mfi].dataPtr(), dimlist(Coarse[mfi].box()) );
  }
  Crse.copy(Coarse) ;
}

void
CompSolver::EdgeHarmAvg(int               dir,
			const MultiFab & Fine,
			MultiFab &       Crse,
			IntVect &        nref)
{
  const BoxArray & fine_ba = Fine.boxArray();

  BoxArray crse_ba( fine_ba ) ;
  crse_ba.enclosedCells(dir);
  crse_ba.coarsen(nref) ;
  crse_ba.surroundingNodes(dir);
  MultiFab Coarse(crse_ba,1,0) ;

  // Add a layer of ghost edges to the current fine MultiFab in
  // direction "dir" so that FORT_EDGE_HARM_AVG called below
  // does not have to handle special boundary cases
  BoxArray grown_fine_ba( fine_ba ) ;
  grown_fine_ba.grow(dir, 1);
  MultiFab GrownFine(grown_fine_ba,1,0);

  for (MFIter mfi(Fine); mfi.isValid(); ++mfi) {
    int i = mfi.index();

    // Fill ghost edges with valid data reflected from the
    // first interior edge
    FArrayBox Temp(fine_ba[i]);
    Temp.copy(Fine[mfi]);

    Temp.shift( dir, -2 );
    GrownFine[mfi].copy(Temp);
    Temp.shift( dir, 4 );
    GrownFine[mfi].copy(Temp);
  }

  // Fill ghost edges with valid data from other fine grids
  // and restore current grid data

  // copyToAll not needed because GrownFine uses extended BoxArray
  // instead of true ghost cells.

  //copyToAll(GrownFine, Fine, 0, 0, 1);
  GrownFine.copy(Fine);

  for (MFIter mfi(Fine); mfi.isValid(); ++mfi) {
    FORT_EDGE_HARM_AVG(&dir,
		       GrownFine[mfi].dataPtr(), dimlist(GrownFine[mfi].box()),
		       nref.getVect(),
		       Coarse[mfi].dataPtr(),    dimlist(Coarse[mfi].box()) );
  }

  Crse.copy(Coarse);
#if 0
  for ( int i=0; i<nfine; i++ ) {

    // Add a layer of ghost edges to the current fine Fab in
    // direction "dir" so that FORT_EDGE_HARM_AVG called below
    // does not have to handle special boundary cases
    Box grown_bx(fine_ba[i]);
    grown_bx.grow( dir, 1 );
    FArrayBox GrownFine(grown_bx);

    // Fill ghost edges with valid data reflected from the
    // first interior edge
    FArrayBox Temp(fine_ba[i]);
    Temp.copy(Fine[i]);

    Temp.shift( dir, -2 );
    GrownFine.copy(Temp);
    Temp.shift( dir, 4 );
    GrownFine.copy(Temp);

    // Fill ghost edges with valid data from other fine grids
    // and restore current grid data
    for ( int i1=0; i1<nfine; i1++ ) {
      GrownFine.copy(Fine[i1]);
    }

    FORT_EDGE_HARM_AVG(&dir,
		       GrownFine.dataPtr(), dimlist(GrownFine.box()),
		       nref.getVect(),
		       Coarse[i].dataPtr(), dimlist(Coarse[i].box()) );
    for( int j=0; j<ncoarse; j++ ) {
      Crse[j].copy(Coarse[i]) ;
    }
  }
#endif
}

void
CompSolver::CellAvgDown(const MultiFab & Fine,
			MultiFab &       Crse,
			IntVect          nref)
{
  int nscal = 1;

  BoxArray coarseGrids(Fine.boxArray());
  coarseGrids.coarsen(nref);

  int ncomp = 1, nghost = 0;
  MultiFab CoarseMF( coarseGrids, ncomp, nghost );

  for (MFIter mfi(CoarseMF); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box& ovlp = coarseGrids[i];
    const int* ovlo = ovlp.loVect();
    const int* ovhi = ovlp.hiVect();
    FORT_AVG_DOWN (CoarseMF[mfi].dataPtr(), dimlist(CoarseMF[mfi].box()),
		   &nscal,
		   Fine[mfi].dataPtr(),     dimlist(Fine[mfi].box()),
		   ovlo,ovhi,
		   nref.getVect());
  }
  Crse.copy(CoarseMF);
#if 0
  // average down cell centered data
  int num_crse = Crse.size();
  int num_fine = Fine.size();
  for (int crse = 0; crse < num_crse; crse++) {
    const Box & cbox = (Crse.boxArray())[crse];
    for (int fine = 0; fine < num_fine; fine++) {
      const BOX& fbox = (Fine.boxArray())[fine];
      Box ovlp(coarsen(fbox,nref));
      ovlp &= cbox;
      if (ovlp.ok()) {
	const int* ovlo = ovlp.loVect();
	const int* ovhi = ovlp.hiVect();
	FORT_AVG_DOWN (Crse[crse].dataPtr(), dimlist(Crse[crse].box()),
		       &nscal,
		       Fine[fine].dataPtr(), dimlist(Fine[fine].box()),
		       ovlo,ovhi,
		       nref.getVect());
      }
    }
  }
#endif
}

void
CompSolver::ZeroCoveredMF( int level, MultiFab & MF )
{
  cout << "CompSolver::ZeroCoveredMF not implemented in parallel" << endl;
  exit(1);

  if( level+1 < CompLevel.size() ) {

    CompSolverLevel & FineLevel = CompLevel[level+1];
    const BoxArray & FineGrids = FineLevel.Grids();

    BoxArray coarsenedFineGrids(FineGrids);
    coarsenedFineGrids.coarsen(FineLevel.getCrseRatio());

    for( int i=0; i<FineGrids.size(); i++ ) {

      FArrayBox Temp(coarsenedFineGrids[i]);
      Temp.setVal(0.);

      for( int j=0; j<MF.size(); j++ ) {
	MF[j].copy(Temp);
      }
    }
  }
}

void
CompSolver::SetBottomParams(Real BottomTol,
			    int BottomMaxIter, int BottomNumIter)
{
  for( int level=0; level<CompLevel.size(); level++ ) {
    CompSolverLevel & Current = CompLevel[level];

    Current.BottomSolveTol( BottomTol );
    Current.BottomSolveMaxIter( BottomMaxIter );
    Current.BottomSolveNumIter( BottomNumIter );
  }
}

void
CompSolver::BuildSecondarySolvers(int AmrSolveBase,
				  int AmrSolveTop)
{
  for (int lev = AmrSolveBase; lev < AmrSolveTop; lev++) {
    CompLevel[AmrLevelIndex[lev]].BuildSecondarySolver();
  }
}

void
CompSolver::SecondarySolve(Real reltol,
			   int  maxiter,
			   int  verbose,
			   int  AmrSolveBase,
			   int  AmrSolveTop)
{
  for (int lev = AmrSolveBase; lev < AmrSolveTop; lev++) {
    CompLevel[AmrLevelIndex[lev]].SecondarySolve(reltol,
						 maxiter, verbose);
  }

#if 1
  IntVect UnitRatio = IntVect::TheUnitVector();
  for (int lev = AmrSolveTop; lev > AmrSolveBase; lev--) {
    CompSolverLevel & Current = CompLevel[AmrLevelIndex[lev]];

    CompSolverLevel & NextAmrLevel = *Current.GetNextCoarserAmrLevel();
    NextAmrLevel.SetRegs( Current.getFluxReg(), NextAmrLevel.Solution(),
			  UpSolve, UnitRatio );
  }

  // Compute the residual norm, and finish setting the correct
  // coarse-fine differential in the flux registers:

  CompSolverLevel & finest = CompLevel[AmrLevelIndex[AmrSolveTop]];
  Real norm = finest.CompositeResidualNorm2();

  if( ParallelDescriptor::IOProcessor() ) {
    int oldprec = cout.precision(20);
    cout << "Multilevel Final Absolute Norm = " << norm << endl;
    cout.precision(oldprec);
  }
#endif
}

void
CompSolver::Solve(Real        reltol,
		  Real        abstol,
		  int       MaxCycle,
		  int        verbose,
		  int   AmrSolveBase,
		  int   AmrSolveTop)
{
  int SolveBase = AmrLevelIndex[AmrSolveBase];
  CompSolverLevel & SolveBaseLevel = CompLevel[SolveBase];

  int finest_level        = CompLevel.size() - 1;
  int finest_active_level =
    (AmrSolveTop >= 0) ? AmrLevelIndex[AmrSolveTop] : finest_level;

  // Make sure coefficients are valid for solve hierarchy
  if( AmrSolveBase < ValidCoefAmrLevel ) {
    AvgDownCoefs( AmrSolveBase );
  }

  // Average down and compute norm of the right hand side
  Real RHSnorm = 0.;
  int lev;
  for( lev=finest_active_level; lev>=SolveBase; lev-- ) {
    CompSolverLevel & Current = CompLevel[lev];
    MultiFab & Rhs = Current.Rhs();

    if( lev < finest_active_level ) {
      CompSolverLevel & Fine = CompLevel[lev+1];

      CellAvgDown( Fine.Rhs(), Rhs, Fine.getCrseRatio() );
    }

    if( Current.AmrLevel() ) {
      for(MFIter mfi(Rhs); mfi.isValid(); ++mfi) {
	Real tnorm = Rhs[mfi].norm(0,0,1);
	if( tnorm > RHSnorm ) RHSnorm = tnorm;
      }
    }
  }
  ParallelDescriptor::ReduceRealMax(RHSnorm);

  if( verbose && ParallelDescriptor::IOProcessor() ) {
    int oldprec = cout.precision(20);
    cout << "Multilevel RHS Norm = " << RHSnorm << endl;
    cout.precision(oldprec);
  }

  // Set the composite level numbers (relative to SolveBase)
  for( lev=SolveBase; lev<=finest_active_level; lev++ ) {
    CompSolverLevel & Current = CompLevel[lev];
    Current.SetCompLevel( lev - SolveBase );
  }

  // The multilevel (bottom) solvers used to be made here.  The hypre
  // version requires a fully-initialized boundary condition object,
  // so we now delay this step until after the bc has been finished.

#if 0
  // Make the bottom solver
  SolveBaseLevel.MakeBottomSolver(verbose ? verbose - 1 : 0);
  // Make multigrid objects at other levels if needed
  for (lev = SolveBase + 1; lev <= finest_active_level; lev++) {
    CompSolverLevel & Current = CompLevel[lev];
    Current.MakeBottomSolver(verbose ? verbose - 1 : 0);
  }
#endif

  // Perform composite solve

  CompSolverLevel & finest = CompLevel[finest_active_level];
  //int CompFinest = finest_active_level - SolveBase;

  int iter = 0;
  bool MoreIterations = true;

  while( MoreIterations ) {

    // Initialize the flux registers with the negative coarse fluxes
    // in preparation for computing the composite residual and its norm
    IntVect UnitRatio = IntVect::TheUnitVector();
    for( lev=finest_active_level; lev>SolveBase; lev-- ) {
      CompSolverLevel & Current = CompLevel[lev];
      if( !Current.AmrLevel() ) continue;

      CompSolverLevel & NextAmrLevel = *Current.GetNextCoarserAmrLevel();
      NextAmrLevel.SetRegs( Current.getFluxReg(), NextAmrLevel.Solution(),
			    UpSolve, UnitRatio );
    }

    // Compute the residual norm, and finish setting the correct
    // coarse-fine differential in the flux registers:

    Real norm, base_norm;
    if ( version == 1 ) {
      norm = finest.CompositeResidualNorm();
    }
    else {
      norm = finest.CompositeResidualNorm2();
    }

    if (iter == 0) {
      // We make the multilevel solvers here.  The hypre version has
      // to be done after the first residual calculation, since before
      // that the boundary condition object has not yet been fully
      // initialized.

      // Make the bottom solver
      SolveBaseLevel.MakeBottomSolver(verbose ? verbose - 1 : 0);
      // Make multigrid objects at other levels if needed
      for (lev = SolveBase + 1; lev <= finest_active_level; lev++) {
	CompSolverLevel & Current = CompLevel[lev];
	if ( version == 1 ) {
	  Current.MakeBottomSolver(verbose ? verbose - 1 : 0);
	}
	else {
	  if( Current.AmrLevel() ) {
	    Current.MakeBottomSolver(verbose ? verbose - 1 : 0);
	  }
	}
      }

      if( verbose && ParallelDescriptor::IOProcessor() ) {
        int oldprec = cout.precision(20);
	cout << "Multilevel Initial Absolute Norm = " << norm << endl;
        cout.precision(oldprec);
      }

      // This definition of base_norm is an improvement over RHSnorm
      // because it allows for the case of a problem where RHSnorm is
      // zero and the total RHS comes from the boundary conditions.
      // Even better would be the initial norm computed with the
      // initial guess set to zero, but we don't have that.

      base_norm = (norm > RHSnorm) ? norm : RHSnorm;
    }

    if( base_norm > 0.0 && verbose && ParallelDescriptor::IOProcessor() ) {
      int oldprec = cout.precision(20);
      cout << "Multilevel iter " << iter << " Relative Norm = "
           << norm/base_norm << endl;
      cout.precision(oldprec);
    }

    if( norm > reltol*base_norm && norm > abstol && iter < MaxCycle ) {

      // Zero the flux registers in preparation for relaxation of the
      // composite hierarchy
      for( lev=finest_active_level; lev>SolveBase; lev-- ) {
	CompSolverLevel & Current = CompLevel[lev];
	if( !Current.AmrLevel() ) continue;

	FluxRegister * Register = Current.getFluxReg();
	Register->setVal(0.);
      }

      // Do a v- or w-cycle:
      if ( version == 1 ) {
	finest.Relax();
      }
      else {
	finest.Relax2();
      }
      iter++;
    }
    else {
      MoreIterations = false;
    }
  }

  if (AmrSolveTop >= 0) {
    for (lev = AmrSolveTop + 1; lev < AmrLevelIndex.size(); lev++) {
      CompSolverLevel & Current = CompLevel[AmrLevelIndex[lev]];
      Current.InterpolateSolution();
    }
  }

  // Delete the bottom solver
  SolveBaseLevel.DeleteBottomSolver();

  for (lev = SolveBase + 1; lev <= finest_active_level; lev++) {
    CompSolverLevel & Current = CompLevel[lev];
    if( version == 1 || Current.AmrLevel() ) {
      Current.DeleteBottomSolver();
    }
  }
}

#if 0
#include <CG_F.H>

void
CompSolver::CGSolve(Real        Tolerance,
		    int           MaxIter,
		    int     PrecondMethod,
		    int           verbose,
		    int      AmrSolveBase)
{
  int finest_level = CompLevel.size() - 1;
  int SolveBase = AmrLevelIndex[AmrSolveBase];
  CompSolverLevel & SolveBaseLevel = CompLevel[SolveBase];

  const int ncomp = 1;

  // Make sure coefficients are valid for solve hierarchy
  if( AmrSolveBase < ValidCoefAmrLevel ) {
    AvgDownCoefs( AmrSolveBase );
  }

  // Set x=0, r=b and compute the L2 norm of b
  Real bnorm = 0.;
  for( int level=SolveBase; level<=finest_level; level++ ) {
    CompSolverLevel & Current = CompLevel[level];
    if( !Current.AmrLevel() ) continue;

    // Allocate CG temporaries
    const BoxArray & grids = Current.Grids();
    MultiFab * x = new MultiFab(grids,1,Fab_allocate);
    Current.SetCGx( x );
    MultiFab * p = new MultiFab(grids,1,Fab_allocate);
    Current.SetCGp( p );

    // Set initial guess to 0 and the initial residual to the
    // system right-hand side.
    x->setVal(0.);
    MultiFab & r = Current.Rhs();

    ZeroCoveredMF( level, r );

    const Real * dx = Current.Geom().CellSize();
    Real scale = 1.;
    for( int j=0; j<BL_SPACEDIM; j++ ) scale *= dx[j];
    for( int i=0; i<r.size(); i++ ) {
      Real trho;
      FORT_CGXDOTY( &trho, r[i].dataPtr(), ARLIM(r[i].loVect()),
		    ARLIM(r[i].hiVect()), r[i].dataPtr(),
		    ARLIM(r[i].loVect()), ARLIM(r[i].hiVect()),
		    BOXARG(grids[i]), &ncomp );
      bnorm += trho*scale;
    }
  }
  bnorm = sqrt(bnorm);

  Real rnorm = bnorm;
  int iter = 0;
  Real bknum, bkden;

  while( rnorm > Tolerance*bnorm && iter < MaxIter ) {

    // Solve Mz = r

    for( int level=SolveBase; level<=finest_level; level++ ) {
      CompSolverLevel & Current = CompLevel[level];
      if( !Current.AmrLevel() ) continue;

      MultiFab & r = Current.Rhs();
      MultiFab & z = Current.Solution();

      z.setVal(0.);

      if( PrecondMethod == 0 || PrecondMethod == 1 ) { // no or diagonal prec.
	for( int i=0; i<z.size(); i++ ) {
	  z[i].copy(r[i]);
	}
      }

      if( PrecondMethod == 1 ) { // diagonal preconditioner
	Current.CGDiagPrecond(z);
      }
      Current.ZeroPhysBndry( z );
    }

    if( PrecondMethod == 2 ) {

      Solve( MLPrecondTol, MLPrecondMaxIter, verbose, AmrSolveBase );

      for( int level=SolveBase; level<=finest_level; level++ ) {
	CompSolverLevel & Current = CompLevel[level];
	if( !Current.AmrLevel() ) continue;

	MultiFab & z = Current.Solution();

	ZeroCoveredMF( level, z );
	Current.ZeroPhysBndry( z );
      }
    }

    // Calculate bknum = (z,r) and p
    bknum = 0.;
    for( int level=SolveBase; level<=finest_level; level++ ) {
      CompSolverLevel & Current = CompLevel[level];
      if( !Current.AmrLevel() ) continue;

      MultiFab & r = Current.Rhs();
      MultiFab & z = Current.Solution();

      const BoxArray & grids = Current.Grids();
      const Real * dx = Current.Geom().CellSize();
      Real scale = 1.;
      for( int j=0; j<BL_SPACEDIM; j++ ) scale *= dx[j];
      for( int i=0; i<z.size(); i++ ) {
	Real trho;
	FORT_CGXDOTY( &trho, z[i].dataPtr(), ARLIM(z[i].loVect()),
		      ARLIM(z[i].hiVect()), r[i].dataPtr(),
		      ARLIM(r[i].loVect()), ARLIM(r[i].hiVect()),
		      BOXARG(grids[i]), &ncomp );
	bknum += trho*scale;
      }
    }

    if( iter == 0 ) {
      // p = z
      for( int level=SolveBase; level<=finest_level; level++ ) {
	CompSolverLevel & Current = CompLevel[level];
	if( !Current.AmrLevel() ) continue;

	MultiFab & p = Current.CGp();
	MultiFab & z = Current.Solution();

	for( int i=0; i<p.size(); i++ ) {
	  p[i].copy(z[i]);
	}
      }
    }
    else {
      // p = z + bk*p
      Real beta = bknum/bkden;
      for( int level=SolveBase; level<=finest_level; level++ ) {
	CompSolverLevel & Current = CompLevel[level];
	if( !Current.AmrLevel() ) continue;

	MultiFab & p = Current.CGp();
	MultiFab & z = Current.Solution();

	const BoxArray & grids = Current.Grids();
	for( int i=0; i<p.size(); i++ ) {
	  FORT_CGADVCP( p[i].dataPtr(), ARLIM(p[i].loVect()),
			ARLIM(p[i].hiVect()), z[i].dataPtr(),
			ARLIM(z[i].loVect()), ARLIM(z[i].hiVect()),
			&beta, BOXARG(grids[i]), &ncomp );
	}
      }
    }
    bkden = bknum;

    // Calculate z = Ap
    SolveBaseLevel.CGAtimes();

    // Compute akden = (p, z) and ak
    Real akden = 0.;
    for( int level=SolveBase; level<=finest_level; level++ ) {
      CompSolverLevel & Current = CompLevel[level];
      if( !Current.AmrLevel() ) continue;

      MultiFab & z = Current.Solution();
      MultiFab & p = Current.CGp();

      ZeroCoveredMF( level, z );

      const BoxArray & grids = Current.Grids();
      const Real * dx = Current.Geom().CellSize();
      Real scale = 1.;
      for( int j=0; j<BL_SPACEDIM; j++ ) scale *= dx[j];
      for( int i=0; i<z.size(); i++ ) {
	Real trho;
	FORT_CGXDOTY( &trho, z[i].dataPtr(), ARLIM(z[i].loVect()),
		      ARLIM(z[i].hiVect()), p[i].dataPtr(),
		      ARLIM(p[i].loVect()), ARLIM(p[i].hiVect()),
		      BOXARG(grids[i]), &ncomp );
	akden += trho*scale;
      }
    }

    Real alpha = bknum/akden;

    // Compute x += alpha p  and  r -= alpha z
    for( int level=SolveBase; level<=finest_level; level++ ) {
      CompSolverLevel & Current = CompLevel[level];
      if( !Current.AmrLevel() ) continue;

      MultiFab & p = Current.CGp();
      MultiFab & r = Current.Rhs();
      MultiFab & z = Current.Solution();
      MultiFab & x = Current.CGx();

      const BoxArray & grids = Current.Grids();
      for( int i=0; i<x.size(); i++ ) {
	FORT_CGUPDATE( x[i].dataPtr(), ARLIM(x[i].loVect()),
		       ARLIM(x[i].hiVect()), r[i].dataPtr(),
		       ARLIM(r[i].loVect()), ARLIM(r[i].hiVect()), &alpha,
		       z[i].dataPtr(), ARLIM(z[i].loVect()),
		       ARLIM(z[i].hiVect()), p[i].dataPtr(),
		       ARLIM(p[i].loVect()), ARLIM(p[i].hiVect()),
		       BOXARG(grids[i]), &ncomp );
      }
    }

    // Compute L2 norm of r
    rnorm = 0.;
    for( int level=SolveBase; level<=finest_level; level++ ) {
      CompSolverLevel & Current = CompLevel[level];
      if( !Current.AmrLevel() ) continue;

      MultiFab & r = Current.Rhs();

      const BoxArray & grids = Current.Grids();
      const Real * dx = Current.Geom().CellSize();
      Real scale = 1.;
      for( int j=0; j<BL_SPACEDIM; j++ ) scale *= dx[j];
      for( int i=0; i<r.size(); i++ ) {
	Real trho;
	FORT_CGXDOTY( &trho, r[i].dataPtr(), ARLIM(r[i].loVect()),
		      ARLIM(r[i].hiVect()), r[i].dataPtr(),
		      ARLIM(r[i].loVect()), ARLIM(r[i].hiVect()),
		      BOXARG(grids[i]), &ncomp );
	rnorm += trho*scale;
      }
    }
    rnorm = sqrt(rnorm);

    if( verbose ) {
      cout << "  CG iteration: " << iter << "  relative residual norm = "
	   << rnorm/bnorm << endl;
    }

    iter++;
  }

  if( iter >= MaxIter ) {
    BoxLib::Error( "Composite CG iteration failed to converge");
  }

  // Copy solution to the CompMLsolution array and deallocate temps
  for( int level=SolveBase; level<=finest_level; level++ ) {
    CompSolverLevel & Current = CompLevel[level];
    if( !Current.AmrLevel() ) continue;

    MultiFab & x = Current.CGx();
    MultiFab & z = Current.Solution();

    for( int i=0; i<z.size(); i++ ) {
      z[i].copy(x[i]);
    }

    delete &Current.CGx();
    delete &Current.CGp();
  }
}
#endif
