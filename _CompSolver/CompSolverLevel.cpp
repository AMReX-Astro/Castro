#include <climits>

#include <BCRec.H>
#include <CompSolverLevel.H>
#include <EdgeVar.H>
//#include <CGSolver.H>
//#include <MultiGrid.H>

#include <Using.H>
#include <ccse-mpi.H>

#define MG_BOTTOM

#include <COMPSOLVER_F.H>

// file scope declarations
static void ConservInterp( MultiFab &, MultiFab &,
			   const IntVect &, const Geometry& );
static void InterpolateAdd( MultiFab &, const MultiFab &, IntVect &,
			    const Geometry &, int );
static void Restrict( const MultiFab &, MultiFab &, const IntVect &, int );

CompSolverLevel::CompSolverLevel(const BoxArray &         _grids,
				 IntVect              _CrseRatio,
				 CompSolverLevel *  NextAmrLevel,
				 IntVect           _AmrCrseRatio,
				 const Geometry &          _geom,
				 const NGBndryBld &    _BndryBld,
				 const BCRec &          _PhysBcr,
				 bool              _AmrLevelFlag,
				 int               _use_hypre,
				 int               _solverflag,
			         int		   _use_harmonic_avg,
				 int               _version)
  : grids(_grids),
    CrseRatio(_CrseRatio),
    NextCoarserAmrLevel(NextAmrLevel),
    NextFinerAmrLevel(0),
    AmrCrseRatio(_AmrCrseRatio),
    geom(_geom),
    BndryBld(_BndryBld),
    AmrLevelFlag(_AmrLevelFlag),
    cycle_type(1),
    presmooth(1),
    postsmooth(1),
    interp_order(0),
    restrict_order(0),
    BottomSolver(0),
    BottomMaxIter(20),
    BottomTol(1.e-6),
    use_hypre(_use_hypre),
    solverflag(_solverflag),
    use_harmonic_avg(_use_harmonic_avg),
    version(_version),
    /* Op(NULL), */ Hd(NULL), Hd2(NULL), Hm(NULL), Hm2(NULL), Bd2(NULL)
{
  BuildBcr( _PhysBcr );

  Bd = BndryBld(grids, 1, geom);

  IntVect UnitRatio = IntVect::TheUnitVector();
  Bd->setBndryConds( PhysBcr, geom, UnitRatio );

  // Initialize boundary values here so that 3rd won't complain
  for (OrientationIter fi; fi; ++fi) {
    Orientation face = fi();
    for (FabSetIter bi((*Bd)[face]); bi.isValid(); ++bi) {
      (*Bd)[face][bi].setVal(0.0);
    }
  }

  if (use_hypre) {
    if (solverflag < 100) {
       Hd = new HypreABec( grids, geom, solverflag );
       Hd->setBndry( *(NGBndry*)Bd );
    }
    else
    {
       // patch solve only 1 level
       Hm = new HypreMultiABec(0, 0, solverflag);
       Hm->setBndry(0, *(NGBndry*)Bd );
       Hm->addLevel(0, geom, grids, IntVect::TheUnitVector());
       Hm->buildMatrixStructure();
    } 
  }

  solution.define(grids,1,1,Fab_allocate);
  correction.define(grids,1,1,Fab_allocate);
  rhs.define(grids,1,0,Fab_allocate);
  residual.define(grids,1,0,Fab_allocate);

  if( AmrLevelFlag ) {
    int AmrLevel = 0;    // This argument is not used

    FluxReg = new FluxRegister(grids,AmrCrseRatio,AmrLevel,1);
  }
  else {
    FluxReg = 0;
  }
}

CompSolverLevel::~CompSolverLevel()
{
  if (use_hypre) {
    if (solverflag < 100) {
       delete Hd;
       Hd = NULL;
    }
    else {
       delete Hm;
       Hm = NULL;
    }
  }
  delete Bd;
  if (Bd2 != NULL) {
    delete Bd2;
    if (Hd2 != NULL) {
      delete Hd2;
    }
    if (Hm2 != NULL) {
      delete Hm2;
    }
    delete rhs2;
  }
  if( AmrLevelFlag ) delete FluxReg;
}

void
CompSolverLevel::BuildSecondarySolver()
{
  // This routine should be called after all coefficients have been
  // initialized, but before CompSolver::Solve is called.  Solve calls
  // AvgDownCoefs, and we want the secondary solver to be based on the
  // original, non-averaged coefficient values.  The right hand side
  // should already be initialized, too.

  if (!use_hypre) {
    cout << "CompSolverLevel::BuildSecondarySolver requires use_hypre" << endl;
    exit(1);
  }

  if (NextFinerAmrLevel) {
    IntVect Ratio = NextFinerAmrLevel->getAmrCrseRatio();

    grids2 = NextFinerAmrLevel->Grids();
    grids2.coarsen(Ratio);

    Bd2 = BndryBld(grids2, 1, geom);

    IntVect UnitRatio = IntVect::TheUnitVector();
    Bd2->setBndryConds( PhysBcr, geom, UnitRatio );

    // Initialize boundary values here so that 3rd won't complain
    for (OrientationIter fi; fi; ++fi) {
      Orientation face = fi();
      for (FabSetIter bi((*Bd2)[face]); bi.isValid(); ++bi) {
	(*Bd2)[face][bi].setVal(0.0);
      }
    }

    if (use_hypre) {
      if (solverflag < 100) {
        Hd2 = new HypreABec( grids2, geom, solverflag );
        Hd2->setBndry( *(NGBndry*)Bd2 );
        Hd2->setScalars(Hd->getAlpha(), Hd->getBeta());
      }
      else {
        Hm2 = new HypreMultiABec( 0, 0, solverflag );
        Hm2->setBndry(0, *(NGBndry*)Bd2 );
        Hm2->addLevel(0, geom, grids2, IntVect::TheUnitVector());
        Hm2->buildMatrixStructure();
        Hm2->setScalars(Hm->getAlpha(), Hm->getBeta());
      }
    }

    // hypre makes internal copies of the coefficient arrays, so we don't
    // need to keep the acoefs and bcoefs around after giving them to Hd2.

    int ncomp=1;
    int ngrow=0;

    {
      MultiFab acoefs2(grids2, ncomp, ngrow);
      acoefs2.copy(Acoef());
      if (use_hypre) {
        if (solverflag < 100) {
          Hd2->aCoefficients(acoefs2);
        }
        else {
          Hm2->aCoefficients(0, acoefs2);
        }
      }
    }
 
    for (int idim = 0; idim < BL_SPACEDIM; idim++) {
      BoxArray edge_boxes(grids2);
      edge_boxes.surroundingNodes(idim);
      MultiFab bcoefs2(edge_boxes, ncomp, ngrow);
      bcoefs2.copy(Bcoef(idim));
      if (use_hypre) {
        if (solverflag < 100) {
          Hd2->bCoefficients(bcoefs2, idim);
        }
        else {
          Hm2->bCoefficients(0, bcoefs2, idim);
        }
      }
    }

    rhs2 = new MultiFab(grids2, ncomp, ngrow);
    rhs2->copy(rhs);
  }
}

void
CompSolverLevel::SecondarySolve(Real reltol, int maxiter, int verbose)
{
  // This routine should be called after CompSolver::Solve is called,
  // but before GetFlux is called.

  if (!use_hypre) {
    cout << "CompSolverLevel::SecondarySolve requires use_hypre" << endl;
    exit(1);
  }

  int ncomp = 1, ngrow = 0;

  MultiFab solution2(grids2, ncomp, ngrow);

  solution2.copy(solution);

  // Create "coarse" boundary register from valid region of solution

  int bndry_InRad = 0, bndry_OutRad = 1, bndry_Extent = 1;
  int nghost = 0, src_comp = 0, dst_comp = 0;

  BndryRegister cbr(grids2, bndry_InRad, bndry_OutRad, bndry_Extent, ncomp);
  cbr.setVal(0.0);

  for (OrientationIter face; face; ++face) {
    Orientation ori = face();
    FabSet& bnd_fs(cbr[ori]);
    bnd_fs.copyFrom(solution, nghost, src_comp, dst_comp, ncomp);
  }

  IntVect Ratio = IntVect::TheUnitVector(); // boundary data is from same level

  // "Interpolate" "coarse" data to "fine" boundary, where applicable.
  // (This is overkill, since the data are at the same resolution a
  // simple copy would be sufficient.  I want the code to be similar to
  // that used for interpolation just for consistency, though.)

  int cbr_Nstart = 0; int fine_Nstart = 0; int bndry_Nstart = 0;

  Bd2->setBndryValues( cbr, cbr_Nstart, solution2, fine_Nstart, bndry_Nstart,
		       ncomp, Ratio, PhysBcr );

  ((NGBndry*)Bd2)->setBndryFluxConds(PhysBcr, Inhomogeneous_BC);

  if (solverflag < 100) {
    Hd2->setVerbose(verbose);
    Hd2->setupSolver(reltol, -1., maxiter);
    Hd2->solve(solution2, 0, *rhs2, Inhomogeneous_BC);
    Hd2->clearSolver();
  }
  else {
    Hm2->setVerbose(verbose);
    Hm2->loadMatrix();
    Hm2->finalizeMatrix();
    Hm2->setupSolver(reltol, -1., maxiter);
    Hm2->loadLevelVectors(0, solution2, 0, *rhs2, Inhomogeneous_BC);
    Hm2->finalizeVectors();
    Hm2->solve();
    Hm2->getSolution(0, solution2, 0);
    Hm2->clearSolver();
  }

  solution.copy(solution2);

#if 1
  // restore coefficients and rhs for composite residual test

  {
    MultiFab acoefs(grids, ncomp, ngrow);
    if (use_hypre) {
      if (solverflag < 100) {
        acoefs.copy(Hd->aCoefficients());
        acoefs.copy(Hd2->aCoefficients());
        Hd->aCoefficients(acoefs);
      }
      else {
        acoefs.copy(Hm->aCoefficients(0));
        acoefs.copy(Hm2->aCoefficients(0));
        Hm->aCoefficients(0, acoefs);
      }
    }
  }

  for (int idim = 0; idim < BL_SPACEDIM; idim++) {
    BoxArray edge_boxes(grids);
    edge_boxes.surroundingNodes(idim);
    MultiFab bcoefs(edge_boxes, ncomp, ngrow);
    if (use_hypre) {
      if (solverflag < 100) {
        bcoefs.copy(Hd->bCoefficients(idim));
        bcoefs.copy(Hd2->bCoefficients(idim));
        Hd->bCoefficients(bcoefs, idim);
      }
      else {
        bcoefs.copy(Hm->bCoefficients(0, idim));
        bcoefs.copy(Hm2->bCoefficients(0, idim));
        Hm->bCoefficients(0, bcoefs, idim);
      }
    }
  }

  rhs.copy(*rhs2);

#endif
}

void
CompSolverLevel::SetScalars(Real alpha,
			    Real  beta)
{
  if (use_hypre) {
    if (solverflag < 100) {
       Hd->setScalars( alpha, beta );
    }
    else {
       Hm->setScalars( alpha, beta );
    }
  }
}

void
CompSolverLevel::SetCoefficients(
                                 const MultiFab & acoef,
				 const MultiFab * bcoef)
{
  if (use_hypre) {
    if (solverflag < 100) {
       Hd->aCoefficients( acoef );
       for (int idim = 0; idim < BL_SPACEDIM; idim++) {
          Hd->bCoefficients(bcoef[idim], idim);
       }
    }
    else {
       Hm->aCoefficients( 0, acoef ); // using only one amr level- level=0
       for (int idim = 0; idim < BL_SPACEDIM; idim++) {
          Hm->bCoefficients(0, bcoef[idim], idim); // using only one amr level- level=0
       }
    }
  }
}

void
CompSolverLevel::aCoefficients(const MultiFab & acoef )
{
  if (use_hypre) {
    if (solverflag < 100) {
       Hd->aCoefficients( acoef );
    }
    else {
       Hm->aCoefficients( 0, acoef ); // using only one amr level- level=0
    }
  }
}

void
CompSolverLevel::bCoefficients(const MultiFab & bcoef, int dir )
{
  if (use_hypre) {
    if (solverflag < 100) {
       Hd->bCoefficients( bcoef, dir );
    }
    else {
       Hm->bCoefficients( 0, bcoef, dir ); // one amr level, level= 0
    }
  }
}

void
CompSolverLevel::GetFlux(
                         MultiFab& dflux,
			 FluxRegister* flux_in,
			 FluxRegister* flux_out)
{
  int i, n;

  if (use_hypre) {
    if (solverflag < 100) {
      Hd->apply(dflux, solution, 0, Inhomogeneous_BC);
    }
    else {
      // one amr level, level = 0
      Hm->initializeApplyLevel(0, dflux, solution, 0, Inhomogeneous_BC);
      Hm->apply();
      Hm->getProduct(0, dflux); // one amr level, level = 0
    }
  }

  // set flux registers

  if (flux_in || flux_out) {
    EdgeVar Flux(grids, 1);

    Real beta;
    if (use_hypre) {
      if (solverflag < 100) {
        beta = Hd->getBeta();
      }
      else {
        beta = Hm->getBeta();
      }
    }

    const Real* dx = geom.CellSize();

    if (use_hypre) {
      solution.FillBoundary();
      if (geom.isAnyPeriodic()) {
	geom.FillPeriodicBoundary(solution, false);
      }
      // junk may be left in off-level boundaries (must be floats)
    }

    for (n = 0; n < BL_SPACEDIM; n++) {
      const MultiFab *bp;
      if (use_hypre) {
         if (solverflag < 100) {
	    bp = &Hd->bCoefficients(n);
         }
         else {
	    bp = &Hm->bCoefficients(0, n); // one amr level, level= 0 
         }
      }
      MultiFab &bcoef = *(MultiFab*)bp;
      for (MFIter fi(Flux[n]); fi.isValid(); ++fi) {
	FORT_SET_ABEC_FLUX(&n,
			   solution[fi].dataPtr(), dimlist(solution[fi].box()),
			   bcoef[fi].dataPtr(),    dimlist(bcoef[fi].box()),
			   &beta,
			   dx,
			   Flux[n][fi].dataPtr(),  dimlist(Flux[n][fi].box()));
      }
    }

    // correct fluxes at physical and coarse-fine boundaries:
    if (use_hypre && flux_out) {
      if (solverflag < 100) {
        Hd->boundaryFlux(Flux.FabPtr(), solution, 0, Inhomogeneous_BC);
      }
      else {
        // one amr level, level= 0
        Hm->boundaryFlux(0, Flux.FabPtr(), solution, 0, Inhomogeneous_BC);
      }
    }

    Flux.Extensive(dx);

    if (flux_in) {
      flux_in->setVal(0.0);
      // New (MultiFab) version, supposed to be more efficient:
      for (n = 0; n < BL_SPACEDIM; n++) {
        flux_in->CrseInit(Flux[n], n, 0, 0, 1, 1.0);
      }
    }
    if (flux_out) {
      flux_out->setVal(0.0);
      for (n = 0; n < BL_SPACEDIM; n++) {
	for (MFIter fi(Flux[n]); fi.isValid(); ++fi) {
	  i = fi.index();
	  flux_out->FineAdd(Flux[n][fi], n, i, 0, 0, 1, 1.0);
	}
      }
    }
  }
}

void
CompSolverLevel::SetRegs( 
                          FluxRegister *     Register,
			  MultiFab &              Phi,
			  int          SolveDirection,
			  const IntVect &       Ratio,
			  const BC_Mode  phys_bc_mode)
{
  BL_ASSERT(Phi.nGrow() > 0);
  BL_ASSERT(Phi.boxArray() == grids);

  EdgeVar Flux(grids, 1);

  Real beta;
  if (use_hypre) {
     if (solverflag < 100) {
        beta = Hd->getBeta();
     }
     else {
        beta = Hm->getBeta();
     }
  }
  const Real* dx = geom.CellSize();

  if (use_hypre) {
    Phi.FillBoundary(); // junk left in off-level boundaries (must be floats)
    if (geom.isAnyPeriodic()) {
      geom.FillPeriodicBoundary(Phi, false);
    }
  }

  for( int n=0; n<BL_SPACEDIM; n++ ) {
    const MultiFab *bp;
    if (use_hypre) {
       if (solverflag < 100) {
          bp = &Hd->bCoefficients(n);
       }
       else {
          bp = &Hm->bCoefficients(0, n); // one amr level, level= 0
       }
    }
    MultiFab &bcoef = *(MultiFab*)bp;
    for (MFIter fi(Flux[n]); fi.isValid(); ++fi) {
      FORT_SET_ABEC_FLUX(&n,
			 Phi[fi].dataPtr(),     dimlist(Phi[fi].box()),
			 bcoef[fi].dataPtr(),   dimlist(bcoef[fi].box()),
			 &beta,
			 dx,
			 Flux[n][fi].dataPtr(), dimlist(Flux[n][fi].box()));
    }
  }

  if( SolveDirection == UpSolve ) {

    // LHH: boundaryFlux not needed on UpSolve since not at boundary

    Flux.Extensive(dx);

    //  Initialize fine flux register
    Register->setVal(0.0);

    // New (MultiFab) version, supposed to be more efficient:
    for( int n=0; n<BL_SPACEDIM; n++ ) {
      Register->CrseInit(Flux[n], n, 0, 0, 1, -1.0 );
    }
  }

  if( SolveDirection == DownSolve ) {

    // LHH: this routine is always called with UnitRatio, so can all
    // this ratio stuff be dispensed with?

    //cout << "Wasteful implementation of SetRegs" << endl;

    // correct fluxes at physical and coarse-fine boundaries:

    if (use_hypre) {
       if (solverflag < 100) {
         Hd->boundaryFlux(Flux.FabPtr(), Phi, 0, phys_bc_mode);
       }
       else {
         // one amr level, level = 0
         Hm->boundaryFlux(0, Flux.FabPtr(), Phi, 0, phys_bc_mode); 
       }
    }

    BoxArray fine_grids(grids);
    fine_grids.refine(Ratio);
    EdgeVar FineFlux(fine_grids, 1);

    FineFlux.Interp( Flux, Ratio );

    Real finedx[BL_SPACEDIM];
    for( int n=0; n<BL_SPACEDIM; n++ ) {
      finedx[n] = dx[n] / Ratio[n];
    }

    FineFlux.Extensive(finedx);

    // Increment the current flux register
    for( int n=0; n<BL_SPACEDIM; n++ ) {
      for (MFIter fi(FineFlux[n]); fi.isValid(); ++fi) {
	int i = fi.index();
	Register->FineAdd( FineFlux[n][fi], n, i, 0, 0, 1, 1.0 );
      }
    }
  }
}

void
CompSolverLevel::SetBndryComposite(const MultiFab & CurrentPhi,
				   const MultiFab &  CoarsePhi,
				   IntVect               Ratio,
				   RadInterpBndryData &        bd,
				   const BC_Mode  phys_bc_mode)
{
  int ncomp = 1;
  int src_comp = 0;
  int dst_comp = 0;

  if( CompLevel > 0 ) {

    BL_ASSERT(Ratio == AmrCrseRatio);

    const Geometry& CoarseGeom = NextCoarserAmrLevel->Geom();

    BoxArray coarseGrids(grids); coarseGrids.coarsen(Ratio);

    int nghost = 1;
    MultiFab CoarseMF( coarseGrids, ncomp, nghost );
    CoarseMF.setVal(0.);

    //CoarseMF.copy(CoarsePhi); // not right---need ghost cells filled too
    copyToAll(CoarseMF, CoarsePhi, CoarseGeom, 0, 0, 1);
#if 0
    for( int i=0; i<grids.size(); i++ ) {
      for( int j=0; j<CoarsePhi.size(); j++ ) {
	CoarseMF[i].copy(CoarsePhi[j]);
      }
    }
#endif

    // Create coarse boundary register
    int bndry_InRad = 0; int bndry_OutRad = 1; int bndry_Extent = 1;
    BndryRegister cbr( coarseGrids, bndry_InRad, bndry_OutRad, bndry_Extent,
		       ncomp );
    for (OrientationIter face; face; ++face) {
        Orientation f = face();
        FabSet& bnd_fs(cbr[f]);
        bnd_fs.copyFrom( CoarseMF, nghost, src_comp, dst_comp, ncomp );
    }

    // Interpolate crse data to fine boundary, where applicable
    int cbr_Nstart = 0; int fine_Nstart = 0; int bndry_Nstart = 0;

    bd.setBndryValues( cbr, cbr_Nstart, CurrentPhi, fine_Nstart, bndry_Nstart,
		       ncomp, Ratio, PhysBcr );

  }
  else {   // composite solve base level

#if 0
    // LHH: 8/8/02
    // Cutting out this section, which used data fillpatched into the
    // ghost cells of the base level solution in cases where the base
    // level was not the coarsest AMR level.  This is not correct
    // because it is different from the way interpolated boundary
    // conditions are computed in the rest of the code.  In particular,
    // doing a 0-1-2 multilevel solve followed by a 1-2 multilevel
    // solve would change the solution, because the bc at the 0-1
    // interface was handled differently.

    // Instead, boundary conditions will be set into bd by the code
    // that calls the CompSolver before Solve is called.  For this to
    // work, we rely on the fact that only inhomogeneous boundary
    // conditions are used on the base level, so there is never a need
    // to reset the coarse bd based on phys_bc_mode.

    // Comments below are older, kept in case I want to revisit this:

    int mf_start = 0; int bnd_start = 0;
    bd.setBndryValues( CurrentPhi, mf_start, bnd_start, ncomp, PhysBcr );

    // According to Milo, this is used for coarser-level boundary data
    // in the case where the base level is not the coarsest AMR level.
    // This presupposes that we have loaded this boundary information
    // into the ghost cells of the base level already, and that there
    // isn't some better way to get the data in.  (Have it already in bd?)

    // Redo copy, since prior call does not copy ghost cells.
    int nghost = 1;
    for (OrientationIter fi; fi; ++fi) {
	bd[fi()].copyFrom( CurrentPhi, nghost, src_comp, dst_comp, ncomp);
    }
#endif
  }

  // We do this last, in case Er has ghost cells which get written into
  // the boundary values:

  ((NGBndry*)&bd)->setBndryFluxConds(PhysBcr, phys_bc_mode);
}

void
CompSolverLevel::ZeroPhysBndry(MultiFab & phi)
{
  cout << "CompSolverLevel::ZeroPhysBndry()---obsolete function, do not call"
       << endl;
  exit(1);

/*
  //const BndryData & bd = Op->bndryData();
  const BndryData bd;

  const int outside_flag = bd.outside_domain;

  Orientation orient_xlo(0,Orientation::low);
  const PArray<Mask>& mask_xlo = bd.bndryMasks(orient_xlo);

  Orientation orient_xhi(0,Orientation::high);
  const PArray<Mask>& mask_xhi = bd.bndryMasks(orient_xhi);

  Orientation orient_ylo(1,Orientation::low);
  const PArray<Mask>& mask_ylo = bd.bndryMasks(orient_ylo);

  Orientation orient_yhi(1,Orientation::high);
  const PArray<Mask>& mask_yhi = bd.bndryMasks(orient_yhi);

  for( int i=0; i<grids.size(); i++ ) {

    // The effect of this routine is only relevant if an inhomogeneous
    // Dirichlet boundary condition will be applied to phi.  Other
    // boundary conditions ignore the values previously in the ghost
    // cells, so it doesn't matter to them what we load here.  Thus
    // there is no need for this routine to use the Bcr (it can just
    // zero everything).  Conversely, if we are not using Dirichlet
    // boundaries there is no need to call this routine at all.

    FORT_ZERO_PHYS_BND(BOXARG(grids[i]),
		       BL_FAA(phi[i]),
		       FAA(mask_xlo[i]),
		       FAA(mask_xhi[i]),
		       FAA(mask_ylo[i]),
		       FAA(mask_yhi[i]),
		       Bcr[i].vect(),
		       &outside_flag);
  }
*/
}

//#define RECREATE 1

void
CompSolverLevel::MakeBottomSolver(
                                  int verbose)
{
  if (CompLevel > 0 && presmooth >= 0 && postsmooth >= 0)
    return;

#ifdef MG_BOTTOM

  if (use_hypre) {
     if (solverflag < 100) {
        Hd->setVerbose(verbose);
     }
     else {
         Hm->setVerbose(verbose);
     }
#ifndef RECREATE
    if (CompLevel == 0) {
       if (solverflag < 100) {
          Hd->setupSolver(BottomTol, -1., BottomNumIter);
       }
       else {
          Hm->loadMatrix();
          Hm->finalizeMatrix();
          Hm->setupSolver(BottomTol, -1., BottomNumIter);
       }
    }
    else {
      if ( version == 1 ) {
        if (solverflag < 100) {
	   Hd->setupSolver(BottomTol, -1., 1);
        }
        else {
           Hm->loadMatrix();
           Hm->finalizeMatrix();
	   Hm->setupSolver(BottomTol, -1., 1);
        }
      }
      else {
        if (solverflag < 100) {
	   Hd->setupSolver(BottomTol, -1., BottomNumIter);
        }
        else {
           Hm->loadMatrix();
           Hm->finalizeMatrix();
	   Hm->setupSolver(BottomTol, -1., BottomNumIter);
        }
      }
    }
#endif
  }

#else

  bool use_mg_precond = true;
  int solver_level = 0;
  CGSolver * cg = new CGSolver(*Op, use_mg_precond, solver_level);
  cg->setMaxIter( BottomMaxIter );

  BottomSolver = (void *)cg;

#endif
}

void
CompSolverLevel::DeleteBottomSolver(void)
{
#ifdef MG_BOTTOM
#ifndef RECREATE
  if (use_hypre) {
    if (solverflag < 100) {
       Hd->clearSolver();
    }
    else {
       Hm->clearSolver();
    }
  }
#endif
#endif

  if (BottomSolver == NULL)
    return;

#ifdef MG_BOTTOM

  //MultiGrid * mg = (MultiGrid *)BottomSolver;

  //delete mg;

#else

  CGSolver * cg = (CGSolver *)BottomSolver;

  delete cg;

#endif
}

void
CompSolverLevel::BuildBcr( const BCRec & _PhysBcr )
{
  for( int n=0; n<BL_SPACEDIM; n++ ) {
    PhysBcr.setLo( n, _PhysBcr.lo(n) );
    PhysBcr.setHi( n, _PhysBcr.hi(n) );
  }

  int ngrids = grids.size();
  Box domain = geom.Domain();

  Bcr.resize(ngrids);

  for( int i=0; i<ngrids; i++ ) {
    BCRec bcr;

    const int* bxlo = grids[i].loVect();
    const int* bxhi = grids[i].hiVect();
    const int* dlo = domain.loVect();
    const int* dhi = domain.hiVect();
    for ( int dir = 0; dir < BL_SPACEDIM; dir++) {
	bcr.setLo(dir, ( bxlo[dir]<=dlo[dir] ? PhysBcr.lo(dir) : INT_DIR ));
	bcr.setHi(dir, ( bxhi[dir]>=dhi[dir] ? PhysBcr.hi(dir) : INT_DIR ));
    }

    Bcr.set(i,bcr);
  }
}

static void clear_internal_borders(FluxRegister& fr)
{
  for (int dir = 0; dir < BL_SPACEDIM; dir++) {
    Orientation lo(dir, Orientation::low);
    Orientation hi(dir, Orientation::high);
    const BoxArray& grids = fr.boxes();
    for (int j = 0; j < grids.size(); j++) {
      Box jbox = BoxLib::bdryHi(grids[j], dir);
      for (FabSetIter fsi(fr[lo]); fsi.isValid(); ++fsi) {
	Box ibox = fr[lo][fsi].box();
	ibox &= jbox;
	if (ibox.ok()) {
	  fr[lo][fsi].setVal(0.0, ibox, 0);
	}
      }
      jbox = BoxLib::bdryLo(grids[j], dir);
      for (FabSetIter fsi(fr[hi]); fsi.isValid(); ++fsi) {
	Box ibox = fr[hi][fsi].box();
	ibox &= jbox;
	if (ibox.ok()) {
	  fr[hi][fsi].setVal(0.0, ibox, 0);
	}
      }
    }
  }
}

Real
CompSolverLevel::CompositeResidualNorm(
                                       FluxRegister * FineRegister,
				       MultiFab *     FineResidual,
				       IntVect *         FineRatio)
{
  // Notes:
  // (1) This function must be called recursively beginning with
  //     the finest level.
  // (2) The composite flux registers
  //     are *assumed* to have been previously initialized with the
  //     negative coarse fluxes.  These registers will be modified
  //     during the execution of this function, and on exit should
  //     contain the correct coarse-fine flux differential.

  MultiFab & phi = solution;
  MultiFab & R = residual;
  MultiFab & rho = rhs;
  RadInterpBndryData & bd = *Bd;

  IntVect UnitRatio = IntVect::TheUnitVector();

  // Compute the residual on the current level using boundary data
  // interpolated from the next coarser level (or data contained
  // in the ghost cells of the solution MultiFab at the composite
  // base level) ignoring any finer levels
  if( AmrLevelFlag ) {

    if( CompLevel > 0 ) {

      CompSolverLevel & CoarseAmr = *NextCoarserAmrLevel;
      MultiFab & Coarsephi = CoarseAmr.Solution();

      SetBndryComposite( phi, Coarsephi, AmrCrseRatio, bd );
    }
    else {

      SetBndryComposite( phi, phi, UnitRatio, bd );
    }

    // Copy current solution to a temporary to pass to the following
    // residual calculation, which will overwrite the boundary data
    // in the physical ghost cells.
    MultiFab phiTemp(grids,1,1);
    for (MFIter mfi(phiTemp); mfi.isValid(); ++mfi) {
      phiTemp[mfi].copy(phi[mfi]);
    }

    if (use_hypre) {
      if (solverflag < 100) {
        Hd->apply(R, phiTemp, 0, Inhomogeneous_BC);
      }
      else {
        // one amr level, level = 0
        Hm->initializeApplyLevel(0, R, phiTemp, 0, Inhomogeneous_BC);
        Hm->apply();
        Hm->getProduct(0, R); // one amr level, level = 0
      }
      for (MFIter mfi(R); mfi.isValid(); ++mfi) {
	R[mfi].minus(rho[mfi]);
	R[mfi].negate();
      }
    }

    // Put solution back
    Box dbox(geom.Domain());
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
      int i = mfi.index();
      Box CopyBox(grids[i]);
      CopyBox.grow(1);
      CopyBox &= dbox;
      phi[mfi].copy(phiTemp[mfi],CopyBox);
    }

    // Reflux the residual from the next finer AMR level
    if( FineRegister ) {

      Real scale = 1.0;
      FineRegister->Reflux( R, scale, 0, 0, 1, geom ) ;
    }
  }

  // Restrict the residual from the next finer level
  if( FineResidual ) {

    if( restrict_order == 1 ) {
      // need a layer of ghost cells for linear interpolation
      MultiFab MF( FineResidual->boxArray(), FineResidual->nComp(),
		   1, Fab_allocate );
      for (MFIter mfi(MF); mfi.isValid(); ++mfi) {
	MF[mfi].setVal(0.);
	MF[mfi].copy((*FineResidual)[mfi]);
      }
      Restrict( MF, R, *FineRatio, restrict_order );
    }
    else {
      Restrict( *FineResidual, R, *FineRatio, restrict_order );
    }
  }

  Real ResidualNorm = 0.;

  // Get residual norm from the coarser levels
  if( CompLevel > 0 ) {
    if( AmrLevelFlag ) {
      SetRegs( FluxReg, phi, DownSolve, UnitRatio );
      ResidualNorm = NextCoarserLevel->CompositeResidualNorm( 
                                                              FluxReg,
							      &R,
							      &CrseRatio);
    }
    else {
      ResidualNorm = NextCoarserLevel->CompositeResidualNorm( 
                                                              FineRegister,
							      &R,
							      &CrseRatio);
    }
  }

#if 0
  static int first = 1;
  if (first) {
    gopen(5);
    black();
    first = 0;
  }
  else if( CompLevel == 0 ) {
    fit(Geom().Domain());
    //Array<IntVect> a(2);
    //a.set(0, IntVect(2,2));
    //PArray<MultiFab> PR(2);
    //PR.set(0, &residual);
    //PR.set(1, &NextFinerAmrLevel->residual);
    //PR.set(0, &solution);
    //PR.set(1, &NextFinerAmrLevel->solution);
    contour(correction, unitvect, 11, 0);
    //contour(PR, a, 101, 0);
    //cout << CompLevel << " " << mfnorm(R) << endl;
    cin.get();
  }
#endif

  // Max the residuals on the current and coarser levels
  for (MFIter mfi(R); mfi.isValid(); ++mfi) {
    Real tnorm = R[mfi].norm(0,0,1);
    if( tnorm > ResidualNorm ) ResidualNorm = tnorm;
  }

  ParallelDescriptor::ReduceRealMax(ResidualNorm);

  return ResidualNorm;
}

Real
CompSolverLevel::CompositeResidualNorm2(
                                        FluxRegister * FineRegister,
				        MultiFab *     FineResidual,
					IntVect *         FineRatio)
{
  // Notes:
  // (1) This function must be called recursively beginning with
  //     the finest level.
  // (2) The composite flux registers
  //     are *assumed* to have been previously initialized with the
  //     negative coarse fluxes.  These registers will be modified
  //     during the execution of this function, and on exit should
  //     contain the correct coarse-fine flux differential.

  MultiFab & phi = solution;
  MultiFab & R = residual;
  MultiFab & rho = rhs;
  RadInterpBndryData & bd = *Bd;

  IntVect UnitRatio = IntVect::TheUnitVector();

  // Compute the residual on the current level using boundary data
  // interpolated from the next coarser level (or data contained
  // in the ghost cells of the solution MultiFab at the composite
  // base level) ignoring any finer levels
  if( AmrLevelFlag ) {

    if( CompLevel > 0 ) {

      CompSolverLevel & CoarseAmr = *NextCoarserAmrLevel;
      MultiFab & Coarsephi = CoarseAmr.Solution();

      SetBndryComposite( phi, Coarsephi, AmrCrseRatio, bd );
    }
    else {

      SetBndryComposite( phi, phi, UnitRatio, bd );
    }

    // Copy current solution to a temporary to pass to the following
    // residual calculation, which will overwrite the boundary data
    // in the physical ghost cells.
    MultiFab phiTemp(grids,1,1);
    for (MFIter mfi(phiTemp); mfi.isValid(); ++mfi) {
      phiTemp[mfi].copy(phi[mfi]);
    }

    if (use_hypre) {
      if (solverflag < 100) {
        Hd->apply(R, phiTemp, 0, Inhomogeneous_BC);
      }
      else {
        // one amr level, level = 0
        Hm->initializeApplyLevel(0, R, phiTemp, 0, Inhomogeneous_BC);
        Hm->apply();
        Hm->getProduct(0, R); // one amr level, level = 0
      }
      for (MFIter mfi(R); mfi.isValid(); ++mfi) {
	R[mfi].minus(rho[mfi]);
	R[mfi].negate();
      }
    }

    // Put solution back
    Box dbox(geom.Domain());
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
      int i = mfi.index();
      Box CopyBox(grids[i]);
      CopyBox.grow(1);
      CopyBox &= dbox;
      phi[mfi].copy(phiTemp[mfi],CopyBox);
    }

    // Reflux the residual from the next finer AMR level
    if( FineRegister ) {

#if 0
  // begin temporary
  if ( AmrLevelFlag ) {
    Real TmpResidualNorm = 0.;
    for (MFIter mfi(R); mfi.isValid(); ++mfi) {
      Real tnorm = R[mfi].norm(0,0,1);
      if( tnorm > TmpResidualNorm ) TmpResidualNorm = tnorm;
    }

    ParallelDescriptor::ReduceRealMax(TmpResidualNorm);

    int oldprec = cout.precision(20);
    cout << "TmpResidualNorm2 on level " << CompLevel
	 << " is " << TmpResidualNorm << endl;
    cout.precision(oldprec);
    //if (CompLevel == 0)
    //printRange(R[0],3,3,20,45,0);
  }

  // end temporary
#endif

      clear_internal_borders( *FineRegister );

      Real scale = 1.0;
      FineRegister->Reflux( R, scale, 0, 0, 1, geom ) ;
    }
  }
  else {
    R.setVal(0.0);
  }

  Real ResidualNorm = 0.;

  // Get residual norm from the coarser levels
  if( CompLevel > 0 ) {
    if( AmrLevelFlag ) {
      SetRegs( FluxReg, phi, DownSolve, UnitRatio );
      ResidualNorm = NextCoarserLevel->CompositeResidualNorm2( 
                                                               FluxReg,
							       &R,
							       &CrseRatio);
    }
    else {
      ResidualNorm = NextCoarserLevel->CompositeResidualNorm2( 
                                                               FineRegister,
							       &R,
							       &CrseRatio);
    }
  }

  if( AmrLevelFlag ) {

    // Max the residuals on the current and coarser AMR levels
    for (MFIter mfi(R); mfi.isValid(); ++mfi) {
      Real tnorm = R[mfi].norm(0,0,1);
      if( tnorm > ResidualNorm ) ResidualNorm = tnorm;
    }

    ParallelDescriptor::ReduceRealMax(ResidualNorm);
  }

  return ResidualNorm;
}

void
CompSolverLevel::InterpolateSolution()
{
  CompSolverLevel & Coarse = *NextCoarserAmrLevel;

  ConservInterp( solution, Coarse.Solution(), AmrCrseRatio,
		 Coarse.Geom() );
}

void
CompSolverLevel::Relax2(
                        FluxRegister * FineRegister)
{
#ifndef MG_BOTTOM
  cout << "MG_BOTTOM must be defined in CompSolverLevel::Relax2" << endl;
  exit(1);
#endif

  MultiFab & e = correction;
  MultiFab & R = residual;

  IntVect UnitRatio = IntVect::TheUnitVector();

  if( FineRegister == 0 ) {   // Finest composite level
    e.setVal(0.);
  }

  if( CompLevel > 0 ) {

    MultiFab eTemp(grids,1,1);
    MultiFab eTempSave(grids,1,1);
    MultiFab TempR(grids,1,restrict_order);
    TempR.setVal(0.);

    CompSolverLevel & Coarse = *NextCoarserLevel;
    MultiFab & CoarseR = Coarse.Residual();
    MultiFab & CoarseCorrection = Coarse.Correction();

    int ncycle = (AmrLevelFlag)? cycle_type: 1;

    for( int cycle=0; cycle<ncycle; cycle++ ) {

      // Initialize the error on the next coarser level
      CoarseCorrection.setVal(0.);

      eTemp.setVal(0.);

      if( AmrLevelFlag ) {

	if ( version == 1 ) {

	  // Presmoothing helps convergence even for version 2, but since
	  // version 2 does high-accuracy solves on fine levels the increased
	  // expense per iteration overwhelms the decrease in number of
          // iterations.

	  if (presmooth < 0) {
	    if (use_hypre) {
#ifdef RECREATE
              if (solverflag < 100) {
                Hd->setupSolver(BottomTol, -1., 1);
              }
              else {
                Hm->loadMatrix();
                Hm->finalizeMatrix();
                Hm->setupSolver(BottomTol, -1., 1);
              }
#endif
              if (solverflag < 100) {
                Hd->solve(eTemp, 0, R, Homogeneous_BC);
              }
              else {
                // one amr level, level = 0
                Hm->loadLevelVectors(0, eTemp, 0, R, Homogeneous_BC);
                Hm->finalizeVectors();
                Hm->solve();
                Hm->getSolution(0, eTemp, 0);
              }
#ifdef RECREATE
              if (solverflag < 100) {
                Hd->clearSolver();
              }
              else {
                Hm->clearSolver();
              }
#endif
	    }
	  }
	  else {
	    cout << "Use CompSolverLevel::Relax2 only with presmooth < 0"
		 << endl;
	    exit(1);
	  }
	}

	SetRegs( FluxReg, eTemp, DownSolve, UnitRatio, Homogeneous_BC );

	if (use_hypre) {
           if (solverflag < 100) {
             Hd->apply(TempR, eTemp, 0, Homogeneous_BC);
           }
           else {
             // one amr level, level = 0
             Hm->initializeApplyLevel(0, TempR, eTemp, 0, Homogeneous_BC);
             Hm->apply();
             Hm->getProduct(0, TempR);
           }

	  for (MFIter mfi(TempR); mfi.isValid(); ++mfi) {
	    TempR[mfi].minus(R[mfi]);
	    TempR[mfi].negate();
	  }
	}

	for (MFIter mfi(eTempSave); mfi.isValid(); ++mfi) {
	  eTempSave[mfi].copy(eTemp[mfi]);
	}
      }
      else {

	for (MFIter mfi(TempR); mfi.isValid(); ++mfi) {
	  TempR[mfi].copy(R[mfi]);
	}
      }

      if( Coarse.AmrLevel() ) {
	// Reflux the coarse residual
	Real scale = 1.0;
	const Geometry & Coarsegeom = Coarse.Geom();

	if( AmrLevelFlag ) {
	  clear_internal_borders( *FluxReg );
	  FluxReg->Reflux( CoarseR, scale, 0, 0, 1, Coarsegeom ) ;
	}
	else {
	  clear_internal_borders( *FineRegister );
	  FineRegister->Reflux( CoarseR, scale, 0, 0, 1, Coarsegeom ) ;
	}
      }

      // Restrict the coarse residual
      //Restrict( TempR, CoarseR, CrseRatio, restrict_order );

      // Relax the next coarser level
      if( AmrLevelFlag ) {
	Coarse.Relax2( FluxReg );
      }
      else {
	Coarse.Relax2( FineRegister );
      }

      //InterpolateAdd( eTemp, CoarseCorrection, CrseRatio,
      //                geom, interp_order );

      SetBndryComposite( eTemp, NextCoarserAmrLevel->Correction(),
			 AmrCrseRatio, *Bd, Homogeneous_BC);

      if( AmrLevelFlag ) {

      if (postsmooth < 0) {
	if (use_hypre) {
#ifdef RECREATE
           if (solverflag < 100) {
             Hd->setupSolver(BottomTol, -1., 1);
           }
           else {
             Hm->loadMatrix();
             Hm->finalizeMatrix();
             Hm->setupSolver(BottomTol, -1., 1);
           }
#endif
           if (solverflag < 100) {
             Hd->solve(eTemp, 0, R, Inhomogeneous_BC);
           }
           else {
             // one amr level, level = 0
             Hm->loadLevelVectors(0, eTemp, 0, R, Inhomogeneous_BC);
             Hm->finalizeVectors();
             Hm->solve();
             Hm->getSolution(0, eTemp, 0);
           }
#ifdef RECREATE
           if (solverflag < 100) {
             Hd->clearSolver();
           }
           else {
             Hm->clearSolver();
           }
#endif
	}
      }
      else {
	cout << "Use CompSolverLevel::Relax2 only with postsmooth < 0"
	     << endl;
	exit(1);
      }

      }

      if( AmrLevelFlag ) {

	if (use_hypre) {
          if (solverflag < 100) {
            Hd->apply(TempR, eTemp, 0, Inhomogeneous_BC);
          }
          else {
            // one amr level, level = 0
            Hm->initializeApplyLevel(0, TempR, eTemp, 0, Inhomogeneous_BC);
            Hm->apply();
            Hm->getProduct(0, TempR);
          }
	}

	for (MFIter mfi(R); mfi.isValid(); ++mfi) {
	  solution[mfi] += eTemp[mfi];
	  R[mfi] -= TempR[mfi];
	  e[mfi] += eTemp[mfi];
	  eTemp[mfi] -= eTempSave[mfi];
	}

	SetRegs( FluxReg, eTemp, DownSolve, UnitRatio );
      }
      else {

	for (MFIter mfi(e); mfi.isValid(); ++mfi) {
	  e[mfi] += eTemp[mfi];
	}
      }
    }
  }
  else {

    MultiFab DPTemp(grids,1,0);
    DPTemp.setVal(0.);

    if (use_hypre) {
#ifdef RECREATE
      if (solverflag < 100) {
        Hd->setupSolver(BottomTol, -1., BottomNumIter);
      }
      else {
        Hm->loadMatrix();
        Hm->finalizeMatrix();
        Hm->setupSolver(BottomTol, -1., BottomNumIter);
      }
#endif
      if (solverflag < 100) {
        Hd->solve(DPTemp, 0, R, Homogeneous_BC);
      }
      else {
        // one amr level, level = 0
        Hm->loadLevelVectors(0, DPTemp, 0, R, Homogeneous_BC);
        Hm->finalizeVectors();
        Hm->solve();
        Hm->getSolution(0, DPTemp, 0);
      }
#ifdef RECREATE
      if (solverflag < 100) {
        Hd->clearSolver();
      }
      else {
        Hm->clearSolver();
      }
#endif
    }

    for (MFIter mfi(e); mfi.isValid(); ++mfi) {
      e[mfi].copy(DPTemp[mfi]);
    }

    if (use_hypre) {
      if (solverflag < 100) {
        Hd->apply( DPTemp, e, 0, Homogeneous_BC );
      }
      else {
        // one amr level, level = 0
        Hm->initializeApplyLevel(0, DPTemp, e, 0, Homogeneous_BC);
        Hm->apply();
        Hm->getProduct(0, DPTemp);
      }
    }

    for (MFIter mfi(R); mfi.isValid(); ++mfi) {
      int i = mfi.index();
      solution[mfi].plus(e[mfi],grids[i],0,0,1);
      R[mfi] -= DPTemp[mfi];
    }

    // a FillBoundary has been done on this MultiFab before, when
    // SetRegs was called in CompSolver::Solve.  The FillBoundary
    // copy descriptor bug therefore kicks in at this point.
    solution.FillBoundary();
    if (geom.isAnyPeriodic()) {
      geom.FillPeriodicBoundary(solution, false);
    }
  }

  if( FineRegister && AmrLevelFlag ) {
    SetRegs( FineRegister, e, UpSolve, UnitRatio );
  }
}

void
CompSolverLevel::Relax(
                       FluxRegister * FineRegister)
{
  MultiFab & e = correction;
  MultiFab & R = residual;

  IntVect UnitRatio = IntVect::TheUnitVector();

  if( FineRegister == 0 ) {   // Finest composite level
    e.setVal(0.);
  }

  if( CompLevel > 0 ) {

    MultiFab eTemp(grids,1,1);
    MultiFab eTempSave(grids,1,1);
    MultiFab TempR(grids,1,restrict_order);
    TempR.setVal(0.);

    CompSolverLevel & Coarse = *NextCoarserLevel;
    MultiFab & CoarseR = Coarse.Residual();
    MultiFab & CoarseCorrection = Coarse.Correction();

    int ncycle = (AmrLevelFlag)? cycle_type: 1;

    for( int cycle=0; cycle<ncycle; cycle++ ) {

      // Initialize the error on the next coarser level
      CoarseCorrection.setVal(0.);

      eTemp.setVal(0.);

      if( AmrLevelFlag ) {

	for( int k=0; k<presmooth; k++ ) {
	  if (use_hypre) {
	    cout << "HypreABec has no smoothing capability" << endl;
	    cout << "Use Hypre only with presmooth < 0" << endl;
	    exit(1);
	  }
	}

#ifdef MG_BOTTOM
	if (presmooth < 0) {
	  if (use_hypre) {
#ifdef RECREATE
            if (solverflag < 100) {
              Hd->setupSolver(BottomTol, -1., 1);
            }
            else {
              Hm->loadMatrix();
              Hm->finalizeMatrix();
              Hm->setupSolver(BottomTol, -1., 1);
            }
#endif
            if (solverflag < 100) {
              Hd->solve(eTemp, 0, R, Homogeneous_BC);
            }
            else {
              // one amr level, level = 0
              Hm->loadLevelVectors(0, eTemp, 0, R, Homogeneous_BC);
              Hm->finalizeVectors();
              Hm->solve();
              Hm->getSolution(0, eTemp, 0);
            }
#ifdef RECREATE
            if (solverflag < 100) {
              Hd->clearSolver();
            }
            else {
              Hm->clearSolver();
            }
#endif
	  }
	}
#endif

	//	ZeroPhysBndry( eTemp );  // to reset bndry ghost cells

	SetRegs( FluxReg, eTemp, DownSolve, UnitRatio, Homogeneous_BC );

	if (use_hypre) {
          if (solverflag < 100) {
            Hd->apply(TempR, eTemp, 0, Homogeneous_BC);
          }
          else {
            // one amr level, level= 0
            Hm->initializeApplyLevel(0, TempR, eTemp, 0, Homogeneous_BC);
            Hm->apply();
            Hm->getProduct(0, TempR);
          }
	  for (MFIter mfi(TempR); mfi.isValid(); ++mfi) {
	    TempR[mfi].minus(R[mfi]);
	    TempR[mfi].negate();
	  }
	}

	for (MFIter mfi(eTempSave); mfi.isValid(); ++mfi) {
	  eTempSave[mfi].copy(eTemp[mfi]);
	}
      }
      else {

	for (MFIter mfi(TempR); mfi.isValid(); ++mfi) {
	  TempR[mfi].copy(R[mfi]);
	}
      }

      if( Coarse.AmrLevel() ) {
	// Reflux the coarse residual
	Real scale = 1.0;
	const Geometry & Coarsegeom = Coarse.Geom();

	if( AmrLevelFlag ) {
	  FluxReg->Reflux( CoarseR, scale, 0, 0, 1, Coarsegeom ) ;
	}
	else {
	  FineRegister->Reflux( CoarseR, scale, 0, 0, 1, Coarsegeom ) ;
	}
      }

      // Restrict the coarse residual
      Restrict( TempR, CoarseR, CrseRatio, restrict_order );

      // Relax the next coarser level
      if( AmrLevelFlag ) {
	Coarse.Relax( FluxReg );
      }
      else {
	Coarse.Relax( FineRegister );
      }

      InterpolateAdd( eTemp, CoarseCorrection, CrseRatio,
		      geom, interp_order );

      SetBndryComposite( eTemp, NextCoarserAmrLevel->Correction(),
			 AmrCrseRatio, *Bd, Homogeneous_BC);

      for( int j=0; j<postsmooth; j++ ) {
	if (use_hypre) {
	  cout << "HypreABec has no smoothing capability" << endl;
	  cout << "Use Hypre only with postsmooth < 0" << endl;
	  exit(1);
	}
      }

#ifdef MG_BOTTOM
      if (postsmooth < 0) {
	if (use_hypre) {
#ifdef RECREATE
          if (solverflag < 100) {
            Hd->setup_solver(BottomTol, -1., 1);
          }
          else {
            Hm->loadMatrix();
            Hm->finalizeMatrix();
            Hm->setupSolver(BottomTol, -1., 1);
          }
#endif
          if (solverflag < 100) {
            Hd->solve(eTemp, 0, R, Inhomogeneous_BC);
          }
          else {
            // one amr level, level = 0
            Hm->loadLevelVectors(0, eTemp, 0, R, Inhomogeneous_BC);
            Hm->finalizeVectors();
            Hm->solve();
            Hm->getSolution(0, eTemp, 0);
          }
#ifdef RECREATE
          if (solverflag < 100) {
            Hd->clearSolver();
          }
          else {
            Hm->clearSolver();
          }
#endif
	}
      }
#endif

      if( AmrLevelFlag ) {

	if (use_hypre) {
          if (solverflag < 100) {
            Hd->apply(TempR, eTemp, 0, Inhomogeneous_BC);
          }
          else {
            // one amr level, level= 0
            Hm->initializeApplyLevel(0, TempR, eTemp, 0, Inhomogeneous_BC);
            Hm->apply();
            Hm->getProduct(0, TempR);
          }
	}
	for (MFIter mfi(R); mfi.isValid(); ++mfi) {
	  solution[mfi] += eTemp[mfi];
	  R[mfi] -= TempR[mfi];
	  e[mfi] += eTemp[mfi];
	  eTemp[mfi] -= eTempSave[mfi];
	}

	SetRegs( FluxReg, eTemp, DownSolve, UnitRatio );
      }
      else {

	//ZeroPhysBndry( eTemp ); // LHH: unnecessary?

	for (MFIter mfi(e); mfi.isValid(); ++mfi) {
	  e[mfi] += eTemp[mfi];
	}
      }
    }
  }
  else {

    MultiFab DPTemp(grids,1,0);
    DPTemp.setVal(0.);

#ifdef MG_BOTTOM
    if (use_hypre) {
#ifdef RECREATE
      if (solverflag < 100) {
        Hd->setupSolver(BottomTol, -1., BottomNumIter);
      }
      else {
        Hm->loadMatrix();
        Hm->finalizeMatrix();
        Hm->setupSolver(BottomTol, -1., BottomNumIter);
      }
#endif
      if (solverflag < 100) {
        Hd->solve(DPTemp, 0, R, Homogeneous_BC);
      }
      else {
        // one amr level, level = 0
        Hm->loadLevelVectors(0, DPTemp, 0, R, Homogeneous_BC);
        Hm->finalizeVectors();
        Hm->solve();
        Hm->getSolution(0, DPTemp, 0);
      }
#ifdef RECREATE
      if (solverflag < 100) {
        Hd->clearSolver();
      }
      else {
        Hm->clearSolver();
      }
#endif
    }
#else

    CGSolver * cg = (CGSolver *)BottomSolver;
    cg->solve( DPTemp, R, BottomTol, -1., Homogeneous_BC );

#endif

    for (MFIter mfi(e); mfi.isValid(); ++mfi) {
      e[mfi].copy(DPTemp[mfi]);
    }
    // LHH: unnecessary since apply will do it?
    //e.FillBoundary();

    if (use_hypre) {
       if (solverflag < 100) {
         Hd->apply(DPTemp, e, 0, Homogeneous_BC);
       }
       else {
         // one amr level, level= 0
         Hm->initializeApplyLevel(0, DPTemp, e, 0, Homogeneous_BC);
         Hm->apply();
         Hm->getProduct(0, DPTemp);
       }
    }

    for (MFIter mfi(R); mfi.isValid(); ++mfi) {
      int i = mfi.index();
      solution[mfi].plus(e[mfi],grids[i],0,0,1);
      R[mfi] -= DPTemp[mfi];
    }

    // a FillBoundary has been done on this MultiFab before, when
    // SetRegs was called in CompSolver::Solve.  The FillBoundary
    // copy descriptor bug therefore kicks in at this point.
    solution.FillBoundary();
    if (geom.isAnyPeriodic()) {
      geom.FillPeriodicBoundary(solution, false);
    }
  }

  if( FineRegister && AmrLevelFlag ) {
    SetRegs( FineRegister, e, UpSolve, UnitRatio );
  }
}

void
CompSolverLevel::CGAtimes( 
                           MultiFab * Coarsep )
{
  BL_ASSERT(AmrLevelFlag);

  MultiFab & p = *cgp;
  RadInterpBndryData & bd = *Bd;

  IntVect UnitRatio = IntVect::TheUnitVector();

  if( Coarsep ) {
    SetBndryComposite( p, *Coarsep, AmrCrseRatio, bd );
  }
  else {
    SetBndryComposite( p, p, UnitRatio, bd );
  }

  //Op->bndryData( bd );

  MultiFab & Ap = solution;

  //Op->apply( Ap, p );

  ZeroPhysBndry(p);

  CompSolverLevel * FinerLevel = GetNextFinerAmrLevel();

  if( FinerLevel ) {

    FluxRegister * FineReg = FinerLevel->getFluxReg();

    SetRegs( FineReg, p, UpSolve, UnitRatio );

    FinerLevel->CGAtimes( &p );

    FineReg->Reflux( Ap, -1.0, 0, 0, 1, geom ) ;

    //    ZeroCoveredMF( Ap );
  }

  if( CompLevel > 0 ) {
    SetRegs( FluxReg, p, DownSolve, UnitRatio );
  }
}

void
CompSolverLevel::CGDiagPrecond(MultiFab & MF)
{
//  const Real * dx = geom.CellSize();

//  for( int i=0; i<grids.size(); i++ ) {
    //FORT_INVERT_DIAG( BOXARG(grids[i]), FAA((Op->aCoefficients())[i]),
    //		      FAA((Op->bCoefficients(0))[i]),
    //		      FAA((Op->bCoefficients(1))[i]), dx, FAA(MF[i]) );
//  }
}

// Interpolates conservatively from c to f using centered slopes.
// Ghost cells of c are modified by calling c.FillBoundary() and by
// linearly extrapolating at physical boundaries.  It is
// assumed that f is properly nested within c.

static void
ConservInterp(MultiFab &f, MultiFab &c,
	      const IntVect & ratio, const Geometry& cgeom)
{
  BoxArray coarseGrids(f.boxArray());
  coarseGrids.coarsen(ratio);

  int ncomp = 1, nghost = 1;
  MultiFab CoarseMF( coarseGrids, ncomp, nghost );
  //CoarseMF.setVal(0.0); // not necessary: proper nesting, cdom handled

  //CoarseMF.copy(c); // need copyToAll to fill ghost cells too
  copyToAll(CoarseMF, c, cgeom, 0, 0, 1);

  Box cdom(cgeom.Domain());
  for (int idim = 0; idim < BL_SPACEDIM; idim++) {
    if (cgeom.isPeriodic(idim))
      cdom.grow(idim, 1); // periodic bdys filled, don't extrapolate them
  }

  for (MFIter mfi(f); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box& ovlp = coarseGrids[i];
    FORT_CONS_INTERP(f[mfi].dataPtr(),        dimlist(f[mfi].box()),
		     CoarseMF[mfi].dataPtr(), dimlist(CoarseMF[mfi].box()),
		     ovlp.loVect(), ovlp.hiVect(),
		     ratio.getVect(),
		     cdom.loVect(), cdom.hiVect());
  }
}

static void
InterpolateAdd(MultiFab &f, const MultiFab &c, IntVect & ratio,
	       const Geometry& fgeom, int order )
{
  BL_ASSERT(f.nGrow() >= order);

  BoxArray coarseGrids(f.boxArray());
  coarseGrids.coarsen(ratio);

  int ncomp = 1, nghost = 0;
  MultiFab CoarseMF( coarseGrids, ncomp, nghost );
  //CoarseMF.setVal(0.0); // not necessary: copy fills all cells

  CoarseMF.copy(c);

  for (MFIter mfi(f); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box& ovlp = coarseGrids[i];
    switch( order )
      {
      case 0:  // Piecewise constant interpolation

	FORT_PC_INTERP(f[mfi].dataPtr(),        dimlist(f[mfi].box()),
		       CoarseMF[mfi].dataPtr(), dimlist(CoarseMF[mfi].box()),
		       ovlp.loVect(), ovlp.hiVect(),
		       ratio.getVect() );
	break;
      case 1:  // Piecewise linear interpolation

	FORT_PL_INTERP(f[mfi].dataPtr(),        dimlist(f[mfi].box()),
		       CoarseMF[mfi].dataPtr(), dimlist(CoarseMF[mfi].box()),
		       ovlp.loVect(), ovlp.hiVect(),
		       ratio.getVect() );
	break;
      default:
	BoxLib::Error("Illegal order passed to InterpolateAdd()");
      }
  }
  f.FillBoundary();
  if (fgeom.isAnyPeriodic()) {
    fgeom.FillPeriodicBoundary(f, false);
  }
}

static void
Restrict( const MultiFab &f, MultiFab &c, const IntVect & ratio, int order )
{
  BL_ASSERT(f.nGrow() >= order);

  BoxArray coarseGrids(f.boxArray());
  coarseGrids.coarsen(ratio);

  int ncomp = 1, nghost = 0;
  MultiFab CoarseMF( coarseGrids, ncomp, nghost );
  //CoarseMF.setVal(0.0); // not necessary: all cells will be filled

  for (MFIter mfi(CoarseMF); mfi.isValid(); ++mfi) {
    int i = mfi.index();
    const Box& ovlp = coarseGrids[i];
    switch( order )
      {
      case 0:  // Piecewise constant restriction

	FORT_PC_RESTRICT(f[mfi].dataPtr(),        dimlist(f[mfi].box()),
			 CoarseMF[mfi].dataPtr(), dimlist(CoarseMF[mfi].box()),
			 ovlp.loVect(), ovlp.hiVect(),
			 ratio.getVect() );
	break;
      case 1:  // Piecewise linear restriction

	FORT_PL_RESTRICT(f[mfi].dataPtr(),        dimlist(f[mfi].box()),
			 CoarseMF[mfi].dataPtr(), dimlist(CoarseMF[mfi].box()),
			 ovlp.loVect(), ovlp.hiVect(),
			 ratio.getVect() );
	break;
      default:
	BoxLib::Error("Illegal order passed to Restrict()");
      }
  }
  c.copy(CoarseMF);
}

// This copyToAll uses only the valid region of src but fills all
// possible cells of dest, including ghost cells, which intersect the
// valid src region or its periodic images.  We write this in terms
// of MultiFab rather than FabArray because FillBoundary operations
// are not defined for FabArray.

void copyToAll(MultiFab&       dest,
	       const MultiFab& src,
	       const Geometry& geom,
	       int             src_comp,
	       int             dest_comp,
	       int             num_comp)
{
  if (!geom.isAnyPeriodic() || dest.nGrow() == 0) {
    copyToAll(dest, src, src_comp, dest_comp, num_comp);
  }
  else {
    // This implementation is ugly, using two additional MultiFab's of
    // storage.  The idea is that some of the ghost cells of dest may
    // extend out of the domain into a periodic image of src.  Rather
    // than add periodic shifts into the already complex logic of
    // copyToAll, we create a new extended MultiFab gsrc in which the
    // potentially troublesome region has been filled with shifted
    // data already.  We do this using FillPeriodicBoundary, but this
    // requires yet another MultiFab since FillPeriodicBoundary operates
    // on ghost cells, and in gsrc these cells are considered part of
    // the valid region.

    BoxArray ggrids(src.boxArray());
    ggrids.grow(dest.nGrow());
    MultiFab gsrc(ggrids, num_comp, 0);

    {
      MultiFab psrc(src.boxArray(), num_comp, dest.nGrow());
      psrc.setVal(0.0); // avoid NaN's, even if some are in ghost cells of src

      for (MFIter mfi(src); mfi.isValid(); ++mfi) {
	int i = mfi.index();
	psrc[mfi].copy(src[mfi], src.box(i), src_comp, src.box(i), 0, num_comp);
      }
      psrc.FillBoundary();
      geom.FillPeriodicBoundary(psrc, false);

      for (MFIter mfi(src); mfi.isValid(); ++mfi) {
	gsrc[mfi].copy(psrc[mfi]);
      }
    }

    copyToAll(dest, gsrc, 0, dest_comp, num_comp);
  }
}

// copyToAll uses only the valid region of src but fills all possible cells
// of dest, including ghost cells, which intersect the valid src region.

template <class FAB>
void copyToAll(FabArray<FAB>&       dest,
	       const FabArray<FAB>& src,
	       int                  src_comp,
	       int                  dest_comp,
	       int                  num_comp)
{
  BL_PROFILE("CompSolverLevel::copyToAll");
  if (dest.nGrow() == 0) {
    dest.copy(src, src_comp, dest_comp, num_comp);
    return;
  }

  if (0 && src.boxArray() == dest.boxArray()) {
    // This would work if FabArray had the FillBoundary function
    // commented out below.
    for (MFIter fai(dest); fai.isValid(); ++fai) {
      int i = fai.index();

      Box intersect = dest.fabbox(i) & src.box(i);

      if (intersect.ok()) {
	dest[fai].copy(src[fai],
		     intersect,
		     src_comp,
		     intersect,
		     dest_comp,
		     num_comp);
      }
    }

    //dest.FillBoundary();
    return;
  }

  typedef typename FAB::value_type value_type;

  const int MyProc = ParallelDescriptor::MyProc();

  FabArrayBase::FabComTag tag;

  //Only did something for BSP, was a no-op for MPI:
  //ParallelDescriptor::SetMessageHeaderSize(sizeof(FabComTag));

  FAB fab;

  vector<FabArrayBase::FabComTag> sTags;
  Array<int>        msgs(ParallelDescriptor::NProcs(), 0);
  Array<int>        nrcv(ParallelDescriptor::NProcs(), 0);

  for (int i = 0; i < dest.size(); ++i) {
    if (dest.DistributionMap()[i] == MyProc) {
      for (int ii = 0; ii < src.size(); ++ii) {
	if (src.DistributionMap()[ii] == MyProc) {
	  //
	  // Both fabs are local
	  //
	  Box intersect = src.boxArray()[ii] & dest.fabbox(i);
	  if (intersect.ok()) {

	    dest[i].copy(src[ii],
			 intersect,
			 src_comp,
			 intersect,
			 dest_comp,
			 num_comp);
	  }
	}
      }
    }
    else {
      for (int ii = 0; ii < src.size(); ++ii) {
	if (src.DistributionMap()[ii] == MyProc) {
	  Box intersect = src.boxArray()[ii] & dest.fabbox(i);
	  if (intersect.ok()) {
	    tag.fromProc = MyProc;
	    tag.toProc   = dest.DistributionMap()[i];
	    tag.fabIndex = i;
	    tag.box      = intersect;
	    //
	    // Use tag.fineIndex to store index into `src'.
	    //
	    tag.fineIndex = ii;
	    sTags.push_back(tag);
	    msgs[dest.DistributionMap()[i]]++;
	  }
	}
      }
    }
  }

  //
  // Pass each processor # of IRecv()s it'll need to post.
  //
  int rc;

#if defined(BL_USE_MPI) || !(defined(BL_BGL) || defined(chaos_3_x86_64_ib) || defined(chaos_3_x86_64))
//#if defined(BL_USE_MPI)
  for (int i = 0; i < msgs.size(); i++) {
    if ((rc = MPI_Reduce(&msgs[i],
			 &nrcv[i],
			 1,
			 MPI_INT,
			 MPI_SUM,
			 i,
			 MPI_COMM_WORLD)) != MPI_SUCCESS)
      ParallelDescriptor::Abort(rc);
  }
  const int NumRecv = nrcv[MyProc];

  Array<MPI_Request> reqs(NumRecv);
  Array<MPI_Status>  stat(NumRecv);
  Array<ParallelDescriptor::CommData>    recv(NumRecv);
  PArray<FAB>        fabs(NumRecv,PArrayManage);
  //
  // First send/receive the box information.
  // I'll receive the NumRecv boxes in any order.
  //
  for (int i = 0; i < NumRecv; i++) {
    if ((rc = MPI_Irecv(recv[i].dataPtr(),
			recv[i].length(),
			MPI_INT,
			MPI_ANY_SOURCE,
			139,
			MPI_COMM_WORLD,
			&reqs[i])) != MPI_SUCCESS)
      ParallelDescriptor::Abort(rc);
  }

  for (int i = 0; i < sTags.size(); i++) {
    ParallelDescriptor::CommData senddata(0, // Not Used.
		      sTags[i].fabIndex,
		      MyProc,
		      //
		      // We use the index into loop over sTags as the ID.
		      // The combination of the loop index and the
		      // processor from which the message was sent forms
		      // a unique identifier.  We'll later use the
		      // combination of fromproc() and id() to match up
		      // the box()s being sent now with the FAB data on
		      // those box()s to be sent next.
		      //
		      i,
		      0, // Not Used.
		      0, // Not Used.
		      0, // Not Used.
		      sTags[i].box);

    if ((rc = MPI_Ssend(senddata.dataPtr(),
			senddata.length(),
			MPI_INT,
			sTags[i].toProc,
			139,
			MPI_COMM_WORLD)) != MPI_SUCCESS)
      ParallelDescriptor::Abort(rc);
  }

  if (NumRecv > 0) {
    if ((rc = MPI_Waitall(NumRecv,
			  reqs.dataPtr(),
			  stat.dataPtr())) != MPI_SUCCESS)
      ParallelDescriptor::Abort(rc);
  }
  //
  // Now the FAB data itself.
  //
  for (int i = 0; i < NumRecv; i++) {
    fabs.set(i, new FAB(recv[i].box(), num_comp));

    if ((rc = MPI_Irecv(fabs[i].dataPtr(),
			fabs[i].box().numPts() * num_comp,
			ParallelDescriptor::Mpi_typemap<value_type>::type(),
			recv[i].fromproc(),
			recv[i].id(),
			MPI_COMM_WORLD,
			&reqs[i])) != MPI_SUCCESS)
      ParallelDescriptor::Abort(rc);
  }
#endif

  for (int i = 0; i < sTags.size(); i++) {
    fab.resize(sTags[i].box, num_comp);

    fab.copy(src[sTags[i].fineIndex],
	     sTags[i].box,
	     src_comp,
	     sTags[i].box,
	     0,
	     num_comp);

    long count = sTags[i].box.numPts() * num_comp;

    BL_ASSERT(count < INT_MAX);
#if defined(BL_USE_MPI) || !(defined(BL_BGL) || defined(chaos_3_x86_64_ib) || defined(chaos_3_x86_64))
//#if defined(BL_USE_MPI)
    //
    // Use MPI_Ssend() to try and force the system not to buffer.
    //
    if ((rc = MPI_Ssend(fab.dataPtr(),
			int(count),
			ParallelDescriptor::Mpi_typemap<value_type>::type(),
			sTags[i].toProc,
			//
			// We use the index into loop over sTags as the ID.
			// The combination of the loop index and the
			// processor from which the message was sent forms
			// a unique identifier.
			//
			// Note that the form of this MPI_Ssend() MUST
			// match the MPI_Send() of the box()es
			// corresponding to this FAB above.
			//
			i,
			MPI_COMM_WORLD)) != MPI_SUCCESS)
      ParallelDescriptor::Abort(rc);
#endif
  }

#if defined(BL_USE_MPI) || !(defined(BL_BGL) || defined(chaos_3_x86_64_ib) || defined(chaos_3_x86_64))
//#if defined(BL_USE_MPI)
  if (NumRecv > 0) {
    if ((rc = MPI_Waitall(NumRecv,
			  reqs.dataPtr(),
			  stat.dataPtr())) != MPI_SUCCESS)
      ParallelDescriptor::Abort(rc);
  }

  for (int i = 0; i < NumRecv; i++) {
    BL_ASSERT(dest.DistributionMap()[recv[i].fabindex()] == MyProc);

    dest[recv[i].fabindex()].copy(fabs[i],
				  fabs[i].box(),
				  0,
				  fabs[i].box(),
				  dest_comp,
				  num_comp);
  }
#endif
}
