#ifndef _HypreABec_H_
#define _HypreABec_H_

#include <Tuple.H>
#include <MultiFab.H>

#include "NGBndry.H"

#include "_hypre_utilities.h"
#include "HYPRE_struct_ls.h"

class HypreABec {

 public:

  // solver_flag = 0 for SMG
  // solver_flag = 1 for PFMG

  HypreABec(const BoxArray& grids, const Geometry& geom,
	    int solver_flag = 0);
  ~HypreABec();

  void setVerbose(int v) {
    verbose = v;
  }

  void setScalars(Real alpha, Real beta);

  Real getAlpha() const {
    return alpha;
  }
  Real getBeta() const {
    return beta;
  }

  void aCoefficients(const MultiFab &a);
  void bCoefficients(const MultiFab &b, int dir);

  void SPalpha(const MultiFab &Spa);

  const MultiFab& aCoefficients() {
    return *acoefs;
  }
  const MultiFab& bCoefficients(int dir) {
    return *bcoefs[dir];
  }

  void setBndry(const NGBndry& bd, int _comp = 0) {
    bdp = &bd;
    bdcomp = _comp;
  }
  const NGBndry& getBndry() {
    return *bdp;
  }
  static Real& fluxFactor() {
    return flux_factor;
  }

  static void getFaceMetric(Array<Real>& r,
                            const Box& reg,
                            const Orientation& ori,
                            const Geometry& geom);

  // The argument inhom in the following methods formerly defaulted
  // to 1.  For greater type safety (to avoid confusion with icomp) it
  // is now an enum with no default.  The argument icomp is always a
  // component number for the independent variable, whether it is
  // called Er, vector, or dest.

  void boundaryFlux(MultiFab* Flux, MultiFab& Er, int icomp, BC_Mode inhom);

  void apply(MultiFab& product, MultiFab& vector, int icomp, BC_Mode inhom);

  // Three steps separated so that multiple calls to solve can be made
  void setupSolver(Real _reltol, Real _abstol, int maxiter);

  void solve(MultiFab& dest, int icomp, MultiFab& rhs, BC_Mode inhom);

  // This is the 2-norm of the complete rhs, including b.c. contributions
  Real getAbsoluteResidual();

  void clearSolver();

 protected:

  const Geometry& geom;

  std::unique_ptr<MultiFab> acoefs;
  std::unique_ptr<MultiFab> bcoefs[BL_SPACEDIM];
  Real alpha, beta;
  Real dx[BL_SPACEDIM];
  Real reltol, abstol;

  std::unique_ptr<MultiFab> SPa; // LO_SANCHEZ_POMRANING alpha

  const NGBndry *bdp;
  int bdcomp; // component number used for bdp

  int solver_flag, verbose, verbose_threshold, pfmg_relax_type, bho;

  HYPRE_StructGrid    grid;
  //HYPRE_StructStencil stencil;

  HYPRE_StructMatrix  A, A0;
  HYPRE_StructVector  b;
  HYPRE_StructVector  x;

  HYPRE_StructSolver  solver;
  HYPRE_StructSolver  precond;

  static Real flux_factor;
};

#endif
