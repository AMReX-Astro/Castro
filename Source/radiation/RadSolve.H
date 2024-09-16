#ifndef CASTRO_RADSOLVE_H
#define CASTRO_RADSOLVE_H

#include <AMReX_Amr.H>

#include <AMReX_FluxRegister.H>

#include <RadBndry.H>
#include <MGRadBndry.H>

#include <HypreABec.H>
#include <HypreMultiABec.H>
#include <HypreExtMultiABec.H>

#include <radsolve_params.H>

///
/// @class RadSolve
/// @brief
///
class RadSolve {

 public:


///
/// @param Parent
///
  RadSolve (amrex::Amr* Parent, int level,
            const amrex::BoxArray& grids,
            const amrex::DistributionMapping& dmap);

///
/// query runtime parameters
///
  static void read_params ();

///
/// @param level
///
  void levelInit(int level);

///
/// @param bd
///
  void levelBndry(RadBndry& bd);

///
/// update multigroup version
///
/// @param mgbd
/// @param comp
///
  void levelBndry(MGRadBndry& mgbd, const int comp);


///
/// @param level
/// @param acoefs
///
  void cellCenteredApplyMetrics(int level, amrex::MultiFab& acoefs);


///
/// @param idim
/// @param geom
/// @param edgebox
/// @param r
/// @param s
///
  static void getEdgeMetric(int idim, const amrex::Geometry& geom, const amrex::Box& edgebox,
                            amrex::Vector<amrex::Real>& r, amrex::Vector<amrex::Real>& s);


///
/// @param level
/// @param &a
///
  void setLevelACoeffs(int level, const amrex::MultiFab &a);

///
/// @param level
/// @param &b
/// @param dir
///
  void setLevelBCoeffs(int level, const amrex::MultiFab &b, int dir);

///
/// @param level
/// @param &c
/// @param dir
///
  void setLevelCCoeffs(int level, const amrex::MultiFab &c, int dir);


///
/// @param level
/// @param fkp
/// @param eta
/// @param etainv
/// @param c
/// @param delta_t
/// @param theta
///
  void levelACoeffs(int level,
                    amrex::MultiFab& fkp, amrex::MultiFab& eta, amrex::MultiFab& etainv,
                    amrex::Real c, amrex::Real delta_t, amrex::Real theta);


///
/// @param level
/// @param amrex::Array<amrex::MultiFab
/// @param lambda
/// @param kappa_r
/// @param kcomp
/// @param c
/// @param lamcomp
///
  void levelBCoeffs(int level,
                    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& lambda,
                    amrex::MultiFab& kappa_r, int kcomp, amrex::Real c, int lamcomp=0);


///
/// @param level
/// @param rhs
/// @param temp
/// @param fkp
/// @param eta
/// @param etainv
/// @param rhoem
/// @param rhoes
/// @param dflux_old
/// @param Er_old
/// @param Edot
/// @param delta_t
/// @param sigma
/// @param c
/// @param theta
/// @param fine_corr
/// @param scale
/// @param igroup
/// @param nu
/// @param dnu
///
  void levelRhs(int level, amrex::MultiFab& rhs,
                amrex::MultiFab& temp,
                amrex::MultiFab& fkp, amrex::MultiFab& eta, amrex::MultiFab& etainv,
                amrex::MultiFab& rhoem, amrex::MultiFab& rhoes,
                amrex::MultiFab& dflux_old, amrex::MultiFab& Er_old, amrex::MultiFab& Edot,
                amrex::Real delta_t, amrex::Real sigma, amrex::Real c, amrex::Real theta,
                amrex::FluxRegister* fine_corr, amrex::Real scale = 1.0,
                int igroup = -1, amrex::Real nu = -1.0, amrex::Real dnu = -1.0);


///
/// @param level
/// @param Er
/// @param igroup
/// @param rhs
/// @param sync_absres_factor
///
  void levelSolve(int level, amrex::MultiFab& Er, int igroup, amrex::MultiFab& rhs,
                  amrex::Real sync_absres_factor);


///
/// @param level
/// @param amrex::Array<amrex::MultiFab
/// @param Flux
/// @param Er
/// @param igroup
///
  void levelFlux(int level,
                 amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& Flux,
                 amrex::MultiFab& Er, int igroup);

///
/// @param level
/// @param flux_in
/// @param flux_out
/// @param amrex::Array<amrex::MultiFab
/// @param Flux
/// @param igroup
///
  void levelFluxReg(int level,
                    amrex::FluxRegister* flux_in, amrex::FluxRegister* flux_out,
                    const amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& Flux,
                    int igroup);

///
/// @param level
/// @param amrex::Array<amrex::MultiFab
/// @param Flux
/// @param flx
/// @param iflx
///
  void levelFluxFaceToCenter(int level, const amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& Flux,
                             amrex::MultiFab& flx, int iflx);

///
/// @param level
/// @param Dterm
/// @param Er
/// @param igroup
///
  void levelDterm(int level, amrex::MultiFab& Dterm, amrex::MultiFab& Er, int igroup);
  void levelClear();


///
/// <MGFLD routines>
///
/// @param bcoefs
/// @param idim
/// @param kappa_r
/// @param kcomp
/// @param lambda
/// @param lamcomp
/// @param c
/// @param geom
///
  void computeBCoeffs(amrex::MultiFab& bcoefs, int idim,
                      amrex::MultiFab& kappa_r, int kcomp,
                      amrex::MultiFab& lambda, int lamcomp,
                      amrex::Real c, const amrex::Geometry& geom);

///
/// @param level
/// @param kappa_p
/// @param delta_t
/// @param c
/// @param igroup
/// @param ptc_tau
///
  void levelACoeffs(int level, amrex::MultiFab& kappa_p,
                    amrex::Real delta_t, amrex::Real c, int igroup, amrex::Real ptc_tau);

///
/// @param level
/// @param rhs
/// @param jg
/// @param muTg
/// @param muYg
/// @param coupT
/// @param etaT
/// @param Er_step
/// @param rhoe_step
/// @param Er_star
/// @param rhoe_star
/// @param delta_t
/// @param igroup
/// @param it
/// @param ptc_tau
///
  void levelRhs(int level, amrex::MultiFab& rhs, const amrex::MultiFab& jg,
                const amrex::MultiFab& muTg,
                const amrex::MultiFab& coupT,
                const amrex::MultiFab& etaT,
                const amrex::MultiFab& Er_step, const amrex::MultiFab& rhoe_step,
                const amrex::MultiFab& Er_star, const amrex::MultiFab& rhoe_star,
                amrex::Real delta_t, int igroup, int it, amrex::Real ptc_tau);

  void levelSPas(int level, amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& lambda, int igroup,
                 int lo_bc[], int hi_bc[]);

///
/// </ MGFLD routines>
///
/// @param level
/// @param amrex::Array<amrex::MultiFab
/// @param lambda
/// @param vel
/// @param dcf
///
  void levelDCoeffs(int level, amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& lambda,
                    amrex::MultiFab& vel, amrex::MultiFab& dcf);

///
/// temporarily change multipliers for C and D coefficients
///
  amrex::Real cMulti, d1Multi, d2Multi;

///
/// @param cMul
/// @param d1Mul
/// @param d2Mul
///
  void setHypreMulti(amrex::Real cMul, amrex::Real d1Mul=0.0, amrex::Real d2Mul=0.0);
  void restoreHypreMulti();

protected:

    amrex::Amr* parent;

    std::unique_ptr<HypreABec> hd;
    std::unique_ptr<HypreMultiABec> hm;
    std::unique_ptr<HypreExtMultiABec> hem;


};

#endif
