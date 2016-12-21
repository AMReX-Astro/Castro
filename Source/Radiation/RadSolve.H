#ifndef _RadSolve_H_
#define _RadSolve_H_

#include <AMReX_Amr.H>

#include <AMReX_FluxRegister.H>

#include "LHH.H"
#include "RadBndry.H"
#include "MGRadBndry.H"

#include "HypreABec.H"
#include "HypreMultiABec.H"
#include "HypreExtMultiABec.H"

class RadSolve {

 public:

  RadSolve(Amr* Parent);
  virtual ~RadSolve() { }

  void levelInit(int level);
  void levelBndry(RadBndry& bd);
  // update multigroup version
  void levelBndry(MGRadBndry& mgbd, const int comp);

  void cellCenteredApplyMetrics(int level, MultiFab& acoefs);
  static void getCellCenterMetric(const Geometry& geom, const Box& reg, Array<Real>& r, Array<Real>& s);
  static void getEdgeMetric(int idim, const Geometry& geom, const Box& edgebox, 
			    Array<Real>& r, Array<Real>& s);

  void setLevelACoeffs(int level, const MultiFab &a);
  void setLevelBCoeffs(int level, const MultiFab &b, int dir);
  void setLevelCCoeffs(int level, const MultiFab &c, int dir);

  void levelACoeffs(int level,
		    MultiFab& fkp, MultiFab& eta, MultiFab& etainv,
		    Real c, Real delta_t, Real theta);

  void levelBCoeffs(int level,
                    Tuple<MultiFab, BL_SPACEDIM>& lambda,
		    MultiFab& kappa_r, int kcomp, Real c, int lamcomp=0);

  void levelRhs(int level, MultiFab& rhs,
		MultiFab& temp,
		MultiFab& fkp, MultiFab& eta, MultiFab& etainv,
		MultiFab& rhoem, MultiFab& rhoes,
		MultiFab& dflux_old, MultiFab& Er_old, MultiFab& Edot,
		Real delta_t, Real sigma, Real c, Real theta,
		FluxRegister* fine_corr, Real scale = 1.0,
                int igroup = -1, Real nu = -1.0, Real dnu = -1.0);

  void levelSolve(int level, MultiFab& Er, int igroup, MultiFab& rhs,
		  Real sync_absres_factor);

  void levelFlux(int level,
                 Tuple<MultiFab, BL_SPACEDIM>& Flux,
                 MultiFab& Er, int igroup);
  void levelFluxReg(int level,
                    FluxRegister* flux_in, FluxRegister* flux_out,
                    const Tuple<MultiFab, BL_SPACEDIM>& Flux,
                    int igroup);
  void levelFluxFaceToCenter(int level, const Tuple<MultiFab, BL_SPACEDIM>& Flux,
			     MultiFab& flx, int iflx);
  void levelDterm(int level, MultiFab& Dterm, MultiFab& Er, int igroup);
  void levelClear();

  // <MGFLD routines>
  void computeBCoeffs(MultiFab& bcoefs, int idim,
                      MultiFab& kappa_r, int kcomp,
                      MultiFab& lambda, int lamcomp,
                      Real c, const Geometry& geom);
  void levelACoeffs(int level, MultiFab& kappa_p, 
		    Real delta_t, Real c, int igroup, Real ptc_tau);  
  void levelRhs(int level, MultiFab& rhs, const MultiFab& jg,
		const MultiFab& muTg, const MultiFab& muYg, 
		const MultiFab& coupT, const MultiFab& coupY, 
		const MultiFab& etaT, const MultiFab& etaY, 
		const MultiFab& thetaT, const MultiFab& thetaY, 
		const MultiFab& Er_step, const MultiFab& rhoe_step, const MultiFab& rhoYe_step, 
		const MultiFab& Er_star, const MultiFab& rhoe_star, const MultiFab& rhoYe_star,
		Real delta_t, int igroup, int it, Real ptc_tau);
  void levelSPas(int level, Tuple<MultiFab, BL_SPACEDIM>& lambda, int igroup,
		 int lo_bc[], int hi_bc[]);
  // </ MGFLD routines>

  void levelDCoeffs(int level, Tuple<MultiFab, BL_SPACEDIM>& lambda,
		    MultiFab& vel, MultiFab& dcf);

  // temporarily change multipliers for C and D coefficients
  Real cMulti, d1Multi, d2Multi;
  void setHypreMulti(Real cMul, Real d1Mul=0.0, Real d2Mul=0.0);
  void restoreHypreMulti();

protected:

  int use_hypre_nonsymmetric_terms;
  int level_solver_flag;

  Real reltol, abstol;
  int maxiter;

  Real alpha, beta;
  Amr* parent;

  int verbose;

  HypreABec      *hd;
  HypreMultiABec *hm;

  // static storage for sync tolerance information
  static Array<Real> absres;
};

#endif
