#ifndef _COMP_SOLVER_H_
#define _COMP_SOLVER_H_

#include <PArray.H>
#include <BoxArray.H>

#include <CompSolverLevel.H>

class CompSolver
{

public:

  // constructors , destructor
  CompSolver( int _use_hypre = 1,
	      int _solverflag = 0,
	      int _use_harmonic_avg = 1,
	      int _version = 1);

  ~CompSolver();

  // access functions

  void SetBndryConds( const NGBndryBld * _BndryBld,
		      const BCRec & _PhysBcr );

  void AddLevel( int AmrLevel, const Geometry & AmrGeom,
		 const BoxArray & AmrGrids, IntVect AmrCrseRatio );

  RadInterpBndryData * GetCrseBndryData() {return CompLevel[0].getBndryData();}

  void Clear(void);

  void SetParms( int cycle_type, int presmooth, int postsmooth,
		 int interp_order, int restrict_order);

  void SetScalars( Real alpha, Real beta );

  void SetCoefficients( int Amr_level, const MultiFab & acoef,
			const MultiFab * bcoef );

  void aCoefficients( int Amr_level, const MultiFab & acoef );

  void bCoefficients( int Amr_level, const MultiFab & bcoef, int dir );

  void SetRhs( int Amr_level, const MultiFab & Rhs );

  void SetInitialGuess( int Amr_level, const MultiFab & Guess );

  void SetZeroInitialGuess( );

  void GetSolution( int Amr_level, MultiFab & Solution );

  void GetFineCorrection( int Amr_level, MultiFab & AmrDiff );

  void GetFlux(int AmrLevel, MultiFab& dflux,
	       FluxRegister* flux_in, FluxRegister* flux_out);

  void BuildSecondarySolvers( int AmrSolveBase, int AmrSolveTop );

  void Solve( Real reltol, Real abstol, int MaxCycle, int verbose,
	      int AmrSolveBase, int AmrSolveTop = -1);

  void SecondarySolve( Real reltol, int maxiter, int verbose,
		       int AmrSolveBase, int AmrSolveTop );

  void CGSolve( Real tol, int MaxCycle, int PrecondMethod, int verbose,
		int AmrSolveBase );

  void SetBottomParams( Real BottomTol,
			int BottomMaxIter, int BottomNumIter );

  void SetMLPrecondParams( Real Tol, int MaxIter ) {MLPrecondTol = Tol;
                                                  MLPrecondMaxIter = MaxIter;}

protected :

  void AvgDownCoefs(int AmrBaseLevel);

  void EdgeAvgDown( int dir, const MultiFab & Fine, MultiFab & Crse,
		    IntVect & nref );

  void EdgeHarmAvg( int dir, const MultiFab & Fine, MultiFab & Crse,
		    IntVect & nref );

  void CellAvgDown( const MultiFab & Fine, MultiFab & Crse, IntVect nref );

  void ZeroCoveredMF( int level, MultiFab & MF );

  PArray<CompSolverLevel> CompLevel;

  Array<int> AmrLevelIndex;

  int use_hypre, solverflag, use_harmonic_avg;
  int solverflag0, solverflagi;

  int version;

  Real alpha;
  Real beta;

  const NGBndryBld *BndryBld;
  BCRec PhysBcr;

  Real MLPrecondTol;
  int MLPrecondMaxIter;

  int ValidCoefAmrLevel;

};

#endif
