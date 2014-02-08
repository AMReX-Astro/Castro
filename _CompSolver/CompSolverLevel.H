#ifndef _COMP_SOLVER_LEVEL_H_
#define _COMP_SOLVER_LEVEL_H_

#include <BoxArray.H>
//#include <ABecLaplacian.H>
#include <BCRec.H>
#include <FluxRegister.H>
//#include <NGBndry.H>  // in HypreABec.H
#include <HypreABec.H>
#include <HypreMultiABec.H>
#include <RadBndry.H>

#include <COMPSOLVER_F.H>

void copyToAll(MultiFab&       dest,
	       const MultiFab& src,
	       const Geometry& geom,
	       int             src_comp,
	       int             dest_comp,
	       int             num_comp);

template <class FAB>
void copyToAll(FabArray<FAB>&       dest,
	       const FabArray<FAB>& src,
	       int                  src_comp,
	       int                  dest_comp,
	       int                  num_comp);

enum SolveType {DownSolve, UpSolve};

class CompSolverLevel
{

public:

  // constructors , destructor
  CompSolverLevel(const BoxArray & ba, IntVect CrseRatio,
		  CompSolverLevel * CoarserAmrLevel, IntVect AmrCrseRatio,
		  const Geometry & geom,
		  const NGBndryBld & _BndryBld, const BCRec & phys_bc,
		  bool AmrLevelFlag,
		  int use_hypre, int solverflag,
		  int use_harmonic_avg,
		  int version);
  ~CompSolverLevel();

  // operators

  // access functions
  void SetNextCoarserLevel( CompSolverLevel * Next ){NextCoarserLevel = Next;}

  void SetNextCoarserAmrLevel( CompSolverLevel * Next )
  {NextCoarserAmrLevel = Next;}

  void SetFinerAmrLevel( CompSolverLevel * Next )
  {NextFinerAmrLevel = Next;}

  CompSolverLevel * GetNextCoarserLevel() {return NextCoarserLevel;}

  CompSolverLevel * GetNextCoarserAmrLevel() {return NextCoarserAmrLevel;}

  CompSolverLevel * GetNextFinerAmrLevel() {return NextFinerAmrLevel;}

  void SetCompLevel( int _CompLevel ) {CompLevel = _CompLevel;}

  int GetCompLevel() {return CompLevel;}

  void SetScalars( Real alpha, Real beta );

  void SetCoefficients( const MultiFab & acoef, const MultiFab * bcoef );

  void aCoefficients( const MultiFab & acoef );

  void bCoefficients( const MultiFab & bcoef, int dir );

  void SetParms( int _cycle_type,
		 int _presmooth, int _postsmooth,
		 int _interp_order, int _restrict_order )
               { cycle_type = _cycle_type;
	         presmooth = _presmooth;
		 postsmooth = _postsmooth;
		 interp_order = _interp_order;
		 restrict_order = _restrict_order;
	       }

  MultiFab & Rhs() {return rhs;}

  MultiFab & CGx() {return *cgx;}

  void SetCGx( MultiFab * x ) {cgx = x;}

  MultiFab & CGp() {return *cgp;}

  void SetCGp( MultiFab * p ) {cgp = p;}

  MultiFab & Solution() {return solution;}

  MultiFab & Correction() {return correction;}

  MultiFab & Residual() {return residual;}

  RadInterpBndryData & getBd() {return *Bd;}

  Geometry & Geom() {return geom;}

  bool AmrLevel() {return AmrLevelFlag;}

  const BoxArray & Grids() {return grids;}

  const IntVect & getCrseRatio() {return CrseRatio;}

  const IntVect & getAmrCrseRatio() {return AmrCrseRatio;}

  FluxRegister * getFluxReg() {return FluxReg;}

  RadInterpBndryData * getBndryData() {return Bd;}

  // operator access functions

#if 0
  const MultiFab & Acoef() {
    return (use_hypre ? Hd->aCoefficients() : Op->aCoefficients());
  }

  const MultiFab & Bcoef( int dir ) {
    return (use_hypre ? Hd->bCoefficients(dir) : Op->bCoefficients(dir));
  }
#else
  const MultiFab & Acoef() {
    // HypreABec solver (hypre struct solver)
    if (solverflag < 100){ 
      return(Hd->aCoefficients());
    }
    // HypreMultiABec solver (hypre sstruct solver)
    else { 
      return(Hm->aCoefficients(0)); // one amr level, level= 0
    }
  }

  const MultiFab & Bcoef( int dir ) {
    // HypreABec solver (hypre struct solver)
    if (solverflag < 100){
      return(Hd->bCoefficients(dir));
    }
    // HypreMultiABec solver (hypre sstruct solver)
    else {
      return(Hm->bCoefficients(0, dir)); // one amr level, level= 0
    }
  }
#endif

  // work functions

  void BuildSecondarySolver();
  void SecondarySolve(Real reltol, int maxiter, int verbose);

  void GetFlux(MultiFab& dflux,
	       FluxRegister* flux_in, FluxRegister* flux_out);

  void SetRegs(FluxRegister * Reg, MultiFab & DeltaPhi,
	       int SolveDirection, const IntVect & Ratio,
	       const BC_Mode phys_bc_mode = Inhomogeneous_BC);

  void SetBndryComposite(const MultiFab & CurrentPhi,
			 const MultiFab & CoarsePhi, IntVect Ratio,
			 RadInterpBndryData & bd,
			 const BC_Mode phys_bc_mode = Inhomogeneous_BC);

  Real CompositeResidualNorm( 
                              FluxRegister * FineRegister = 0,
			      MultiFab * FineResidual = 0,
			      IntVect * FineRatio = 0 );

  Real CompositeResidualNorm2( 
                               FluxRegister * FineRegister = 0,
			       MultiFab * FineResidual = 0,
			       IntVect * FineRatio = 0 );

  void InterpolateSolution();

  void Relax2( 
               FluxRegister * FineRegister = 0 );

  void Relax( 
              FluxRegister * FineRegister = 0 );

  void ZeroPhysBndry(MultiFab & phi);

  void MakeBottomSolver(int verbose);

  void DeleteBottomSolver(void);

  void CGAtimes( MultiFab * Coarsep = 0 );

  void CGDiagPrecond(MultiFab & MF);

  void BottomSolveMaxIter( int MaxIter ) {BottomMaxIter = MaxIter;}
  int  BottomSolveMaxIter( void ) {return BottomMaxIter;}

  void BottomSolveNumIter( int NumIter ) {BottomNumIter = NumIter;}
  int  BottomSolveNumIter( void ) {return BottomNumIter;}

  void BottomSolveTol( Real Tol ) {BottomTol = Tol;}
  Real BottomSolveTol( void ) {return BottomTol;}


protected :

  int CompLevel;

  BoxArray grids;

  int use_hypre, solverflag, use_harmonic_avg;
  int version;

  //ABecLaplacian * Op;
  HypreABec     * Hd;
  HypreMultiABec * Hm;

  const NGBndryBld & BndryBld;
  RadInterpBndryData  * Bd;
  
  bool AmrLevelFlag;

  IntVect CrseRatio;

  IntVect AmrCrseRatio;

  Geometry geom;

  MultiFab rhs;
  MultiFab solution;
  MultiFab residual;
  MultiFab correction;

  MultiFab *cgx;
  MultiFab *cgp;

  Array<BCRec> Bcr;

  BCRec PhysBcr;

  FluxRegister * FluxReg;

  CompSolverLevel * NextCoarserLevel;

  CompSolverLevel * NextCoarserAmrLevel;

  CompSolverLevel * NextFinerAmrLevel;

  void *BottomSolver;

  // Secondary solver components
  BoxArray          grids2;
  HypreABec        *Hd2;
  HypreMultiABec   *Hm2;
  RadInterpBndryData  *Bd2;
  MultiFab         *rhs2;

  void BuildBcr( const BCRec & PhysBcr );

  int interp_order;

  int restrict_order;

  int presmooth;

  int postsmooth;

  int cycle_type;

  // Bottom solver parameters
  int BottomMaxIter;
  int BottomNumIter;
  Real BottomTol;

};

#endif
