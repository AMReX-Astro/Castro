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

  RadSolve(amrex::Amr* Parent);
  virtual ~RadSolve() { }

  void levelInit(int level);
  void levelBndry(RadBndry& bd);
  // update multigroup version
  void levelBndry(MGRadBndry& mgbd, const int comp);

  void cellCenteredApplyMetrics(int level, amrex::MultiFab& acoefs);
  static void getCellCenterMetric(const amrex::Geometry& geom, const amrex::Box& reg, amrex::Array<amrex::Real>& r, amrex::Array<amrex::Real>& s);
  static void getEdgeMetric(int idim, const amrex::Geometry& geom, const amrex::Box& edgebox, 
			    amrex::Array<amrex::Real>& r, amrex::Array<amrex::Real>& s);

  void setLevelACoeffs(int level, const amrex::MultiFab &a);
  void setLevelBCoeffs(int level, const amrex::MultiFab &b, int dir);
  void setLevelCCoeffs(int level, const amrex::MultiFab &c, int dir);

  void levelACoeffs(int level,
		    amrex::MultiFab& fkp, amrex::MultiFab& eta, amrex::MultiFab& etainv,
		    amrex::Real c, amrex::Real delta_t, amrex::Real theta);

  void levelBCoeffs(int level,
                    amrex::Tuple<amrex::MultiFab, BL_SPACEDIM>& lambda,
		    amrex::MultiFab& kappa_r, int kcomp, amrex::Real c, int lamcomp=0);

  void levelRhs(int level, amrex::MultiFab& rhs,
		amrex::MultiFab& temp,
		amrex::MultiFab& fkp, amrex::MultiFab& eta, amrex::MultiFab& etainv,
		amrex::MultiFab& rhoem, amrex::MultiFab& rhoes,
		amrex::MultiFab& dflux_old, amrex::MultiFab& Er_old, amrex::MultiFab& Edot,
		amrex::Real delta_t, amrex::Real sigma, amrex::Real c, amrex::Real theta,
		amrex::FluxRegister* fine_corr, amrex::Real scale = 1.0,
                int igroup = -1, amrex::Real nu = -1.0, amrex::Real dnu = -1.0);

  void levelSolve(int level, amrex::MultiFab& Er, int igroup, amrex::MultiFab& rhs,
		  amrex::Real sync_absres_factor);

  void levelFlux(int level,
                 amrex::Tuple<amrex::MultiFab, BL_SPACEDIM>& Flux,
                 amrex::MultiFab& Er, int igroup);
  void levelFluxReg(int level,
                    amrex::FluxRegister* flux_in, amrex::FluxRegister* flux_out,
                    const amrex::Tuple<amrex::MultiFab, BL_SPACEDIM>& Flux,
                    int igroup);
  void levelFluxFaceToCenter(int level, const amrex::Tuple<amrex::MultiFab, BL_SPACEDIM>& Flux,
			     amrex::MultiFab& flx, int iflx);
  void levelDterm(int level, amrex::MultiFab& Dterm, amrex::MultiFab& Er, int igroup);
  void levelClear();

  // <MGFLD routines>
  void computeBCoeffs(amrex::MultiFab& bcoefs, int idim,
                      amrex::MultiFab& kappa_r, int kcomp,
                      amrex::MultiFab& lambda, int lamcomp,
                      amrex::Real c, const amrex::Geometry& geom);
  void levelACoeffs(int level, amrex::MultiFab& kappa_p, 
		    amrex::Real delta_t, amrex::Real c, int igroup, amrex::Real ptc_tau);  
  void levelRhs(int level, amrex::MultiFab& rhs, const amrex::MultiFab& jg,
		const amrex::MultiFab& muTg, const amrex::MultiFab& muYg, 
		const amrex::MultiFab& coupT, const amrex::MultiFab& coupY, 
		const amrex::MultiFab& etaT, const amrex::MultiFab& etaY, 
		const amrex::MultiFab& thetaT, const amrex::MultiFab& thetaY, 
		const amrex::MultiFab& Er_step, const amrex::MultiFab& rhoe_step, const amrex::MultiFab& rhoYe_step, 
		const amrex::MultiFab& Er_star, const amrex::MultiFab& rhoe_star, const amrex::MultiFab& rhoYe_star,
		amrex::Real delta_t, int igroup, int it, amrex::Real ptc_tau);
  void levelSPas(int level, amrex::Tuple<amrex::MultiFab, BL_SPACEDIM>& lambda, int igroup,
		 int lo_bc[], int hi_bc[]);
  // </ MGFLD routines>

  void levelDCoeffs(int level, amrex::Tuple<amrex::MultiFab, BL_SPACEDIM>& lambda,
		    amrex::MultiFab& vel, amrex::MultiFab& dcf);

  // temporarily change multipliers for C and D coefficients
  amrex::Real cMulti, d1Multi, d2Multi;
  void setHypreMulti(amrex::Real cMul, amrex::Real d1Mul=0.0, amrex::Real d2Mul=0.0);
  void restoreHypreMulti();

protected:

  int use_hypre_nonsymmetric_terms;
  int level_solver_flag;

  amrex::Real reltol, abstol;
  int maxiter;

  amrex::Real alpha, beta;
  amrex::Amr* parent;

  int verbose;

  HypreABec      *hd;
  HypreMultiABec *hm;

  // static storage for sync tolerance information
  static amrex::Array<amrex::Real> absres;
};

#endif
