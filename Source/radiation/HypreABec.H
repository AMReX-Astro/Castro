#ifndef _HypreABec_H_
#define _HypreABec_H_

#include <AMReX_Tuple.H>
#include <AMReX_MultiFab.H>

#include "NGBndry.H"

#include "_hypre_utilities.h"
#include "HYPRE_struct_ls.h"

class HypreABec {

 public:

  // solver_flag = 0 for SMG
  // solver_flag = 1 for PFMG

  HypreABec(const amrex::BoxArray& grids,
	    const amrex::DistributionMapping& dmap,
	    const amrex::Geometry& geom,
	    int solver_flag = 0);
  ~HypreABec();

  void setVerbose(int v) {
    verbose = v;
  }

  void setScalars(amrex::Real alpha, amrex::Real beta);

  amrex::Real getAlpha() const {
    return alpha;
  }
  amrex::Real getBeta() const {
    return beta;
  }

  void aCoefficients(const amrex::MultiFab &a);
  void bCoefficients(const amrex::MultiFab &b, int dir);

  void SPalpha(const amrex::MultiFab &Spa);

  const amrex::MultiFab& aCoefficients() {
    return *acoefs;
  }
  const amrex::MultiFab& bCoefficients(int dir) {
    return *bcoefs[dir];
  }

  void setBndry(const NGBndry& bd, int _comp = 0) {
    bdp = &bd;
    bdcomp = _comp;
  }
  const NGBndry& getBndry() {
    return *bdp;
  }
  static amrex::Real& fluxFactor() {
    return flux_factor;
  }

  static void getFaceMetric(amrex::Vector<amrex::Real>& r,
                            const amrex::Box& reg,
                            const amrex::Orientation& ori,
                            const amrex::Geometry& geom);

  // The argument inhom in the following methods formerly defaulted
  // to 1.  For greater type safety (to avoid confusion with icomp) it
  // is now an enum with no default.  The argument icomp is always a
  // component number for the independent variable, whether it is
  // called Er, vector, or dest.

  void boundaryFlux(amrex::MultiFab* Flux, amrex::MultiFab& Er, int icomp, BC_Mode inhom);

  void apply(amrex::MultiFab& product, amrex::MultiFab& vector, int icomp, BC_Mode inhom);

  // Three steps separated so that multiple calls to solve can be made
  void setupSolver(amrex::Real _reltol, amrex::Real _abstol, int maxiter);

  void solve(amrex::MultiFab& dest, int icomp, amrex::MultiFab& rhs, BC_Mode inhom);

  // This is the 2-norm of the complete rhs, including b.c. contributions
  amrex::Real getAbsoluteResidual();

  void clearSolver();

 protected:

  const amrex::Geometry& geom;

  std::unique_ptr<amrex::MultiFab> acoefs;
  std::unique_ptr<amrex::MultiFab> bcoefs[BL_SPACEDIM];
  amrex::Real alpha, beta;
  amrex::Real dx[BL_SPACEDIM];
  amrex::Real reltol, abstol;

  std::unique_ptr<amrex::MultiFab> SPa; // LO_SANCHEZ_POMRANING alpha

  const NGBndry *bdp;
  int bdcomp; // component number used for bdp

  int solver_flag, verbose, verbose_threshold, pfmg_relax_type, bho;

  HYPRE_StructGrid    hgrid;
  //HYPRE_StructStencil stencil;

  HYPRE_StructMatrix  A, A0;
  HYPRE_StructVector  b;
  HYPRE_StructVector  x;

  HYPRE_StructSolver  solver;
  HYPRE_StructSolver  precond;

  static amrex::Real flux_factor;
};

#endif
