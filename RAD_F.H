#ifndef _RAD_F_H_
#define _RAD_F_H_

#include <ArrayLim.H>
#include <LHH.H>

#ifdef BL_FORT_USE_UNDERSCORE
#  define   FORT_SETLAYOUT setlayout_
#  define   FORT_MULTRS  multrs_
#  define   FORT_SPHC    sphc_
#  define   FORT_SPHE    sphe_
#  define   FORT_LACOEF  lacoef_
#  define   FORT_KAVG    kavg_
#  define   FORT_BCLIM   bclim_
#  define   FORT_FLXLIM  flxlim_
#  define   FORT_EDDFAC  eddfac_
#  define   FORT_SCGRD1  scgrd1_
#  define   FORT_SCGRD2  scgrd2_
#  define   FORT_SCGRD3  scgrd3_
#  define   FORT_LRHS    lrhs_
#  define   FORT_ANATW2  anatw2_
#  define   FORT_CFRHOE  cfrhoe_
#  define   FORT_GTEMP   gtemp_
#  define   FORT_GCV     gcv_
#  define   FORT_CEXCH   cexch_
#  define   FORT_CETA2   ceta2_
#  define   FORT_CEUP    ceup_
#  define   FORT_CEUPDTERM    ceupdterm_
#  define   FORT_NCEUP   nceup_
#  define   FORT_CETOT   cetot_
#  define   FORT_FKPN    fkpn_
#  define   FORT_ROSSE1  rosse1_
#  define   FORT_ROSSE1S  rosse1s_
#  define   FORT_NFLOOR  nfloor_
#  define   FORT_LACOEFMGFLD  lacoefmgfld_
#  define   FORT_RFFACE  rfface_
#  define   FORT_RADBNDRY  rbndry_
#  define   FORT_RADBNDRY2 rbndr2_
#  define   FORT_INIT_OPACITY_TABLE init_opacity_table_
#  define   FORT_BEXTRP bextrp_
#  define   FORT_LBCOEFNA lbcoefna_
#  define   FORT_LJUPNA ljupna_

#elif defined(BL_FORT_USE_UPPERCASE)

#  define   FORT_SETLAYOUT SETLAYOUT
#  define   FORT_MULTRS  MULTRS
#  define   FORT_SPHC    SPHC
#  define   FORT_SPHE    SPHE
#  define   FORT_LACOEF  LACOEF
#  define   FORT_KAVG    KAVG
#  define   FORT_BCLIM   BCLIM
#  define   FORT_FLXLIM  FLXLIM
#  define   FORT_EDDFAC  EDDFAC
#  define   FORT_SCGRD1  SCGRD1
#  define   FORT_SCGRD2  SCGRD2
#  define   FORT_SCGRD3  SCGRD3
#  define   FORT_LRHS    LRHS
#  define   FORT_ANATW2  ANATW2
#  define   FORT_CFRHOE  CFRHOE
#  define   FORT_GTEMP   GTEMP
#  define   FORT_GCV     GCV
#  define   FORT_CEXCH   CEXCH
#  define   FORT_CETA2   CETA2
#  define   FORT_CEUP    CEUP
#  define   FORT_CEUPDTERM    CEUPDTERM
#  define   FORT_NCEUP   NCEUP
#  define   FORT_CETOT   CETOT
#  define   FORT_FKPN    FKPN
#  define   FORT_ROSSE1  ROSSE1
#  define   FORT_ROSSE1S  ROSSE1S
#  define   FORT_NFLOOR  NFLOOR
#  define   FORT_LACOEFMGFLD  LACOEFMGFLD
#  define   FORT_RFFACE  RFFACE
#  define   FORT_RADBNDRY  RBNDRY
#  define   FORT_RADBNDRY2 RBNDR2
#  define   FORT_INIT_OPACITY_TABLE INIT_OPACITY_TABLE
#  define   FORT_BEXTRP BEXTRP
#  define   FORT_LBCOEFNA LBCOEFNA
#  define   FORT_LJUPNA LJUPNA

#endif

#include <BLFort.H>

BL_FORT_PROC_DECL(CA_INIT_RADHYDRO_PARS,ca_init_radhydro_pars)
    (const int & fspace_advection_type, const int & do_inelastic_scattering,
     const int & comoving,
     const int & first_order_hydro, const Real & flatten_pp_threshold);

BL_FORT_PROC_DECL(CA_INITRADCONSTANTS,ca_initradconstants)
  (const Real& pi, Real& clight, Real& hplanck,
   Real& kboltz, Real& stefbol, Real& avogadro,
   Real& mev2erg, int& J_is_used);

BL_FORT_PROC_DECL(CA_INITSINGLEGROUP,ca_initsinglegroup)
(const int& ngroups);

BL_FORT_PROC_DECL(CA_INITGROUPS,ca_initgroups)
  (const Real* nugroup, const Real* dnugroup,
   const int& ngroups, const int& ng0, const int& ng1);

BL_FORT_PROC_DECL(CA_INITGROUPS2,ca_initgroups2)
  (const Real* nugroup, const Real* dnugroup,
   const Real* xnu, const int& ngroups);

BL_FORT_PROC_DECL(CA_INITGROUPS3,ca_initgroups3)
  (const Real* nugroup, const Real* dnugroup, const Real* dlognugroup,
   const Real* xnu, const int& ngroups, const int& ng0, const int& ng1);

#ifdef NEUTRINO
BL_FORT_PROC_DECL(CA_DERNEUT,ca_derneut)
    (BL_FORT_FAB_ARG(der),const int* nvar,
     const BL_FORT_FAB_ARG(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const Real* delta, const Real* xlo,
     const Real* time, const Real* dt, const int* bcrec,
     const int* level, const int* grid_no);

BL_FORT_PROC_DECL(CA_DERRHOYL,ca_derrhoyl)
    (BL_FORT_FAB_ARG(der),const int* nvar,
     const BL_FORT_FAB_ARG(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const Real* delta, const Real* xlo,
     const Real* time, const Real* dt, const int* bcrec,
     const int* level, const int* grid_no);

BL_FORT_PROC_DECL(CA_DERYL,ca_deryl)
    (BL_FORT_FAB_ARG(der),const int* nvar,
     const BL_FORT_FAB_ARG(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const Real* delta, const Real* xlo,
     const Real* time, const Real* dt, const int* bcrec,
     const int* level, const int* grid_no);

BL_FORT_PROC_DECL(CA_DERYNUE,ca_derynue)
    (BL_FORT_FAB_ARG(der),const int* nvar,
     const BL_FORT_FAB_ARG(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const Real* delta, const Real* xlo,
     const Real* time, const Real* dt, const int* bcrec,
     const int* level, const int* grid_no);

BL_FORT_PROC_DECL(CA_DERYNUAE,ca_derynuae)
    (BL_FORT_FAB_ARG(der),const int* nvar,
     const BL_FORT_FAB_ARG(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const Real* delta, const Real* xlo,
     const Real* time, const Real* dt, const int* bcrec,
     const int* level, const int* grid_no);
#endif

BL_FORT_PROC_DECL(CA_COMPUTE_POWERLAW_KAPPA, ca_compute_powerlaw_kappa)
     (const int lo[], const int hi[],
      BL_FORT_FAB_ARG(kappa), const BL_FORT_FAB_ARG(u), 
      const Real* kappa0, const Real* m, const Real* n,  const Real* p, 
      const Real* Tfloor, const Real* kfloor);

BL_FORT_PROC_DECL(CA_COMPUTE_POWERLAW_KAPPA_S, ca_compute_powerlaw_kappa_s)
     (const int lo[], const int hi[],
      BL_FORT_FAB_ARG(kappa), const BL_FORT_FAB_ARG(u), 
      const Real* kappa0, const Real* m, const Real* n, const Real* p, 
      const Real* s0, const Real* sm, const Real* sn, const Real* sp, 
      const Real* Tfloor, const Real* kfloor);
 
BL_FORT_PROC_DECL(CA_COMPUTE_KAPKAP, ca_compute_kapkap)
     (BL_FORT_FAB_ARG(kapkap), const BL_FORT_FAB_ARG(kap_r)); 

BL_FORT_PROC_DECL(CA_COMPUTE_LAMBORDER, ca_compute_lamborder)
   (const BL_FORT_FAB_ARG(Er), const BL_FORT_FAB_ARG(kap), 
    BL_FORT_FAB_ARG(lam),  const Real* dx, const int* ngrow, const int* limiter,
    const int* filter_lambda_T, const int* filter_lambda_S);

BL_FORT_PROC_DECL(CA_GET_V_DCF, ca_get_v_dcf)
   (const int lo[], const int hi[],
    const BL_FORT_FAB_ARG(er), 
    const BL_FORT_FAB_ARG(state), 
    const BL_FORT_FAB_ARG(temp),
    const BL_FORT_FAB_ARG(c_v),
    const BL_FORT_FAB_ARG(kr),
    const BL_FORT_FAB_ARG(kp),
    const BL_FORT_FAB_ARG(kp2),
    const Real* dtemp, const Real* dtime, const Real* sigma, const Real* c,
    BL_FORT_FAB_ARG(v),
    BL_FORT_FAB_ARG(dcf));

BL_FORT_PROC_DECL(CA_COMPUTE_DCOEFS, ca_compute_dcoefs)
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(d), 
    const BL_FORT_FAB_ARG(lam), 
    const BL_FORT_FAB_ARG(v), const BL_FORT_FAB_ARG(dcf),
    const Real* r, const int* idir);
 
BL_FORT_PROC_DECL(CA_UPDATE_DCF, ca_update_dcf)
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(dcf),
    const BL_FORT_FAB_ARG(etainv),
    const BL_FORT_FAB_ARG(kp),
    const BL_FORT_FAB_ARG(kr));

BL_FORT_PROC_DECL(CA_EST_GPR0, ca_est_gpr0)
   (const BL_FORT_FAB_ARG(Er),
    BL_FORT_FAB_ARG(gPr));

BL_FORT_PROC_DECL(CA_EST_GPR2, ca_est_gpr2)
   (const BL_FORT_FAB_ARG(kap),
    const BL_FORT_FAB_ARG(Er),
    BL_FORT_FAB_ARG(gPr), 
    const int lo[], const int hi[],
    const Real* dx, const int* limiter, const int* comoving); 

BL_FORT_PROC_DECL(CA_SET_DTERM_FACE, ca_set_dterm_face)
   (const int lo[], const int hi[],
    const BL_FORT_FAB_ARG(Er),
    const BL_FORT_FAB_ARG(dcoef), 
    BL_FORT_FAB_ARG(dterm_face), 
    const Real* dx, const int* idir);

BL_FORT_PROC_DECL(CA_FACE2CENTER, ca_face2center)
    (const int lo[], const int hi[], 
     const int& scomp, const int& dcomp, const int& ncomp, const int& nf, const int& nc,
     D_DECL(const BL_FORT_FAB_ARG(foox),
	    const BL_FORT_FAB_ARG(fooy),
	    const BL_FORT_FAB_ARG(fooz)),
     BL_FORT_FAB_ARG(fooc));

BL_FORT_PROC_DECL(CA_CORRECT_DTERM, ca_correct_dterm)
    (D_DECL(BL_FORT_FAB_ARG(dfx),
	    BL_FORT_FAB_ARG(dfy),
	    BL_FORT_FAB_ARG(dfz)),
     const Real* re, const Real* rc);

// <MGFLD>
BL_FORT_PROC_DECL(CA_INITFLUXLIMITER, ca_initfluxlimiter)
   (const int* limiter, const int* closure);

#ifdef NEUTRINO
BL_FORT_PROC_DECL(CA_ACCEL_ACOE_NEUT, ca_accel_acoe_neut)
   (const int lo[], const int hi[],
    const BL_FORT_FAB_ARG(eta1),
    const BL_FORT_FAB_ARG(thetaT),
    const BL_FORT_FAB_ARG(thetaY),
    const BL_FORT_FAB_ARG(spc),
    const BL_FORT_FAB_ARG(kap),
    BL_FORT_FAB_ARG(aco),
    const Real* dt, const Real* ptc_tau);
BL_FORT_PROC_DECL(CA_ACCEL_RHS_NEUT, ca_accel_rhs_neut)
   (const int lo[], const int hi[],
    const BL_FORT_FAB_ARG(Ern),
    const BL_FORT_FAB_ARG(Erl),
    const BL_FORT_FAB_ARG(kap),
    const BL_FORT_FAB_ARG(etaT),
    const BL_FORT_FAB_ARG(etaY),
    const BL_FORT_FAB_ARG(thetaT),
    const BL_FORT_FAB_ARG(thetaY),
    BL_FORT_FAB_ARG(rhs),
    const Real* dt);
BL_FORT_PROC_DECL(CA_ACCEL_SPEC_NEUT, ca_accel_spec_neut)
   (const int lo[], const int hi[],
    const BL_FORT_FAB_ARG(Er_new),
    const BL_FORT_FAB_ARG(Er_pi),
    const BL_FORT_FAB_ARG(kap),
    const BL_FORT_FAB_ARG(etaT),
    const BL_FORT_FAB_ARG(etaY),
    const BL_FORT_FAB_ARG(thetaT),
    const BL_FORT_FAB_ARG(thetaY),
    const BL_FORT_FAB_ARG(mugT),
    const BL_FORT_FAB_ARG(mugY),
    BL_FORT_FAB_ARG(spec), 
    const Real* dt, const Real* ptc_tau);
BL_FORT_PROC_DECL(CA_CHECK_CONV_ER_NEUT, ca_check_conv_er_neut)
   (const int lo[], const int hi[],
    const BL_FORT_FAB_ARG(Ern),
    const BL_FORT_FAB_ARG(Erl),
    const BL_FORT_FAB_ARG(kap),
    const BL_FORT_FAB_ARG(etTz),
    const BL_FORT_FAB_ARG(etYz),
    const BL_FORT_FAB_ARG(thTz),
    const BL_FORT_FAB_ARG(thYz),
    const BL_FORT_FAB_ARG(temp),
    const BL_FORT_FAB_ARG(Ye),
    Real* rel, Real* abs, Real *errr, const Real* dt);
BL_FORT_PROC_DECL(CA_CHECK_CONV_NEUT, ca_check_conv_neut)
   (const int lo[], const int hi[],
    const BL_FORT_FAB_ARG(rhoe_new),
    const BL_FORT_FAB_ARG(rhoe_star),
    const BL_FORT_FAB_ARG(rhoe_step),
    const BL_FORT_FAB_ARG(Er_new),
    const BL_FORT_FAB_ARG(temp_new),
    const BL_FORT_FAB_ARG(temp_star),
    const BL_FORT_FAB_ARG(rhoYe_new),
    const BL_FORT_FAB_ARG(rhoYe_star),
    const BL_FORT_FAB_ARG(rhoYe_step),
    const BL_FORT_FAB_ARG(rho),
    const BL_FORT_FAB_ARG(kappa_p),
    const BL_FORT_FAB_ARG(jg),
    const BL_FORT_FAB_ARG(dedT),
    const BL_FORT_FAB_ARG(dedY),
    Real* rel_rhoe, Real* abs_rhoe, 
    Real* rel_FT  , Real* abs_FT, 
    Real* rel_T   , Real* abs_T, 
    Real* rel_FY  , Real* abs_FY, 
    Real* rel_Ye  , Real* abs_Ye,
    const Real* delta_t);
BL_FORT_PROC_DECL(CA_COMPUTE_COUPTY, ca_compute_coupty)
    (const int lo[], const int hi[],
     BL_FORT_FAB_ARG(coupt), 
     BL_FORT_FAB_ARG(coupy), 
     const BL_FORT_FAB_ARG(kpp),
     const BL_FORT_FAB_ARG(Eg),
     const BL_FORT_FAB_ARG(jg));
BL_FORT_PROC_DECL(CA_COMPUTE_DEDX, ca_compute_dedx)
    (const int* lo, const int *hi, 
     const BL_FORT_FAB_ARG(state),
     const BL_FORT_FAB_ARG(T),
     const BL_FORT_FAB_ARG(Ye),
     const BL_FORT_FAB_ARG(Ts),
     const BL_FORT_FAB_ARG(Yes),
     BL_FORT_FAB_ARG(dedT),
     BL_FORT_FAB_ARG(dedY),
     const int* validStar);
BL_FORT_PROC_DECL(CA_COMPUTE_ETA_THE, ca_compute_eta_the)
    (const int* lo, const int *hi, 
     BL_FORT_FAB_ARG(etaT), 
     BL_FORT_FAB_ARG(etaTz), 
     BL_FORT_FAB_ARG(etaY), 
     BL_FORT_FAB_ARG(etaYz), 
     BL_FORT_FAB_ARG(eta1), 
     BL_FORT_FAB_ARG(thetaT), 
     BL_FORT_FAB_ARG(thetaTz), 
     BL_FORT_FAB_ARG(thetaY), 
     BL_FORT_FAB_ARG(thetaYz), 
     BL_FORT_FAB_ARG(theta1), 
     BL_FORT_FAB_ARG(djdT),
     BL_FORT_FAB_ARG(djdY),
     const BL_FORT_FAB_ARG(dkdT),
     const BL_FORT_FAB_ARG(dkdY),
     const BL_FORT_FAB_ARG(dedT),
     const BL_FORT_FAB_ARG(dedY),
     const BL_FORT_FAB_ARG(Er_star),
     const BL_FORT_FAB_ARG(rho),
     const Real* delta_t, const Real* tau);
BL_FORT_PROC_DECL(CA_COMPUTE_RHS_NEUT, ca_compute_rhs_neut)
   (const int* lo, const int *hi, 
    BL_FORT_FAB_ARG(rhs),
    const BL_FORT_FAB_ARG(jg),
    const BL_FORT_FAB_ARG(mugT),
    const BL_FORT_FAB_ARG(mugY),
    const BL_FORT_FAB_ARG(coupT),
    const BL_FORT_FAB_ARG(coupY),
    const BL_FORT_FAB_ARG(etaT),
    const BL_FORT_FAB_ARG(etaY),
    const BL_FORT_FAB_ARG(thetaT),
    const BL_FORT_FAB_ARG(thetaY),
    const BL_FORT_FAB_ARG(Er_step),
    const BL_FORT_FAB_ARG(rhoe_step),
    const BL_FORT_FAB_ARG(rhoYe_step),
    const BL_FORT_FAB_ARG(Er_star),
    const BL_FORT_FAB_ARG(rhoe_star),    
    const BL_FORT_FAB_ARG(rhoYe_star),    
    const Real* r,
    const Real* dt, const int* igroup, const Real* ptc_tau);
BL_FORT_PROC_DECL(CA_COMPUTE_ROSSELAND_NEUT, ca_compute_rosseland_neut)
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(kpr),
    const BL_FORT_FAB_ARG(state));
BL_FORT_PROC_DECL(CA_COMPUTE_PLANCK_NEUT, ca_compute_planck_neut)
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(kpp),
    const BL_FORT_FAB_ARG(state));
BL_FORT_PROC_DECL(CA_LOCAL_ACCEL_NEUT, ca_local_accel_neut)
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(Ern),
    const BL_FORT_FAB_ARG(Erl),
    const BL_FORT_FAB_ARG(kap),
    const BL_FORT_FAB_ARG(etaT),
    const BL_FORT_FAB_ARG(etaY),
    const BL_FORT_FAB_ARG(theT),
    const BL_FORT_FAB_ARG(theY),
    const BL_FORT_FAB_ARG(mugT),
    const BL_FORT_FAB_ARG(mugY),
    const Real* dt, const Real* ptc_tau);
BL_FORT_PROC_DECL(CA_OPAC_EMIS_NEUT, ca_opac_emis_neut)
   (const int lo[], const int hi[],
    const BL_FORT_FAB_ARG(Snew),
    const BL_FORT_FAB_ARG(T),
    const BL_FORT_FAB_ARG(Ye),
    const BL_FORT_FAB_ARG(Ts),
    const BL_FORT_FAB_ARG(Yes),
    BL_FORT_FAB_ARG(kpp),
    BL_FORT_FAB_ARG(kpr),
    BL_FORT_FAB_ARG(jg),
    BL_FORT_FAB_ARG(djdT),
    BL_FORT_FAB_ARG(djdY),
    BL_FORT_FAB_ARG(dkdT),
    BL_FORT_FAB_ARG(dkdY),
    const int* use_dkdT, const int* validStar, const int* lag_opac); 
BL_FORT_PROC_DECL(CA_STATE_UPDATE_NEUT, ca_state_update_neut)
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(state),
    const BL_FORT_FAB_ARG(rhoe),
    const BL_FORT_FAB_ARG(Ye),
    const BL_FORT_FAB_ARG(temp),
    const BL_FORT_FAB_ARG(msk),
    Real* derat, Real* dTrat, Real* dye);
BL_FORT_PROC_DECL(CA_UPDATE_MATTER_NEUT, ca_update_matter_neut)
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(rhoe_new),
    BL_FORT_FAB_ARG(rhoYe_new),
    BL_FORT_FAB_ARG(Ye_new),
    const BL_FORT_FAB_ARG(Er_new),
    const BL_FORT_FAB_ARG(Er_pi),
    const BL_FORT_FAB_ARG(rhoe_star),
    const BL_FORT_FAB_ARG(rhoYe_star),
    const BL_FORT_FAB_ARG(rhoe_step),
    const BL_FORT_FAB_ARG(rhoYe_step),
    const BL_FORT_FAB_ARG(etaT),
    const BL_FORT_FAB_ARG(etaY),
    const BL_FORT_FAB_ARG(eta1),
    const BL_FORT_FAB_ARG(thetaT),
    const BL_FORT_FAB_ARG(thetaY),
    const BL_FORT_FAB_ARG(theta1),
    const BL_FORT_FAB_ARG(coupT),
    const BL_FORT_FAB_ARG(coupY),
    const BL_FORT_FAB_ARG(kappa_p),
    const BL_FORT_FAB_ARG(mugT),
    const BL_FORT_FAB_ARG(mugY),
    const BL_FORT_FAB_ARG(Snew),
    const Real* dt, const Real* ptc_tau);
BL_FORT_PROC_DECL(CA_NCUPDATE_MATTER_NEUT, ca_ncupdate_matter_neut)
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(temp_new),
    BL_FORT_FAB_ARG(Ye_new),
    const BL_FORT_FAB_ARG(Er_new),
    const BL_FORT_FAB_ARG(rhoe_star),
    const BL_FORT_FAB_ARG(rhoYe_star),
    const BL_FORT_FAB_ARG(rhoe_step),
    const BL_FORT_FAB_ARG(rhoYe_step),
    const BL_FORT_FAB_ARG(etaTz),
    const BL_FORT_FAB_ARG(etaYz),
    const BL_FORT_FAB_ARG(thetaTz),
    const BL_FORT_FAB_ARG(thetaYz),
    const BL_FORT_FAB_ARG(kappa_p),
    const BL_FORT_FAB_ARG(jg),
    const Real* dt);

#else   // photons

BL_FORT_PROC_DECL(CA_ACCEL_ACOE, ca_accel_acoe)
   (const int lo[], const int hi[],
    const BL_FORT_FAB_ARG(eta1),
    const BL_FORT_FAB_ARG(spc),
    const BL_FORT_FAB_ARG(kap),
    BL_FORT_FAB_ARG(aco),
    const Real* dt, const Real* ptc_tau);
BL_FORT_PROC_DECL(CA_ACCEL_RHS, ca_accel_rhs)
   (const int lo[], const int hi[],
    const BL_FORT_FAB_ARG(Ern),
    const BL_FORT_FAB_ARG(Erl),
    const BL_FORT_FAB_ARG(kap),
    const BL_FORT_FAB_ARG(etaT),
    BL_FORT_FAB_ARG(rhs),
    const Real* dt);
BL_FORT_PROC_DECL(CA_ACCEL_SPEC, ca_accel_spec)
   (const int lo[], const int hi[],
    const BL_FORT_FAB_ARG(kap),
    const BL_FORT_FAB_ARG(mugT),
    BL_FORT_FAB_ARG(spec), 
    const Real* dt, const Real* ptc_tau);
BL_FORT_PROC_DECL(CA_CHECK_CONV, ca_check_conv)
   (const int lo[], const int hi[],
    const BL_FORT_FAB_ARG(rhoe_new),
    const BL_FORT_FAB_ARG(rhoe_star),
    const BL_FORT_FAB_ARG(rhoe_step),
    const BL_FORT_FAB_ARG(Er_new),
    const BL_FORT_FAB_ARG(temp_new),
    const BL_FORT_FAB_ARG(temp_star),
    const BL_FORT_FAB_ARG(rho),
    const BL_FORT_FAB_ARG(kappa_p),
    const BL_FORT_FAB_ARG(jg),
    const BL_FORT_FAB_ARG(dedT),
    Real* rel_rhoe, Real* abs_rhoe, 
    Real* rel_FT  , Real* abs_FT, 
    Real* rel_T   , Real* abs_T, 
    const Real* delta_t);
BL_FORT_PROC_DECL(CA_CHECK_CONV_ER, ca_check_conv_er)
   (const int lo[], const int hi[],
    const BL_FORT_FAB_ARG(Ern),
    const BL_FORT_FAB_ARG(Erl),
    const BL_FORT_FAB_ARG(kap),
    const BL_FORT_FAB_ARG(etTz),
    const BL_FORT_FAB_ARG(temp),
    Real* rel, Real* abs, Real *errr, const Real* dt);
BL_FORT_PROC_DECL(CA_COMPUTE_COUPT, ca_compute_coupt)
    (const int lo[], const int hi[],
     BL_FORT_FAB_ARG(coupt), 
     const BL_FORT_FAB_ARG(kpp),
     const BL_FORT_FAB_ARG(Eg),
     const BL_FORT_FAB_ARG(jg));
BL_FORT_PROC_DECL(CA_COMPUTE_ETAT, ca_compute_etat)
    (const int* lo, const int *hi, 
     BL_FORT_FAB_ARG(etat), 
     BL_FORT_FAB_ARG(etaTz), 
     BL_FORT_FAB_ARG(eta1), 
     BL_FORT_FAB_ARG(djdT),
     const BL_FORT_FAB_ARG(dkdT),
     const BL_FORT_FAB_ARG(dedT),
     const BL_FORT_FAB_ARG(Er_star),
     const BL_FORT_FAB_ARG(rho),
     const Real* delta_t, const Real* tau);
BL_FORT_PROC_DECL(CA_COMPUTE_RHS, ca_compute_rhs)
   (const int* lo, const int *hi, 
    BL_FORT_FAB_ARG(rhs),
    const BL_FORT_FAB_ARG(jg),
    const BL_FORT_FAB_ARG(mugT),
    const BL_FORT_FAB_ARG(coupT),
    const BL_FORT_FAB_ARG(etaT),
    const BL_FORT_FAB_ARG(Er_step),
    const BL_FORT_FAB_ARG(rhoe_step),
    const BL_FORT_FAB_ARG(Er_star),
    const BL_FORT_FAB_ARG(rhoe_star),    
    const Real* r,
    const Real* dt, const int* igroup, const Real* ptc_tau);
BL_FORT_PROC_DECL(CA_LOCAL_ACCEL, ca_local_accel)
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(Ern),
    const BL_FORT_FAB_ARG(Erl),
    const BL_FORT_FAB_ARG(kap),
    const BL_FORT_FAB_ARG(etaT),
    const BL_FORT_FAB_ARG(mugT),
    const Real* dt, const Real* ptc_tau);
BL_FORT_PROC_DECL(CA_STATE_UPDATE, ca_state_update)
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(state),
    const BL_FORT_FAB_ARG(rhoe),
    const BL_FORT_FAB_ARG(temp),
    const BL_FORT_FAB_ARG(msk),
    Real* derat, Real* dTrat);
BL_FORT_PROC_DECL(CA_UPDATE_MATTER, ca_update_matter)
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(rhoe_new),
    const BL_FORT_FAB_ARG(Er_new),
    const BL_FORT_FAB_ARG(Er_pi),
    const BL_FORT_FAB_ARG(rhoe_star),
    const BL_FORT_FAB_ARG(rhoe_step),
    const BL_FORT_FAB_ARG(eta1),
    const BL_FORT_FAB_ARG(coupT),
    const BL_FORT_FAB_ARG(kappa_p),
    const BL_FORT_FAB_ARG(mugT),
    const BL_FORT_FAB_ARG(Snew),
    const Real* dt, const Real* ptc_tau);
BL_FORT_PROC_DECL(CA_NCUPDATE_MATTER, ca_ncupdate_matter)
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(temp_new),
    const BL_FORT_FAB_ARG(Er_new),
    const BL_FORT_FAB_ARG(rhoe_star),
    const BL_FORT_FAB_ARG(rhoe_step),
    const BL_FORT_FAB_ARG(etaTz),
    const BL_FORT_FAB_ARG(kappa_p),
    const BL_FORT_FAB_ARG(jg),
    const Real* dt);

#endif

BL_FORT_PROC_DECL(CA_COMPUTE_ROSSELAND, ca_compute_rosseland)
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(kpr),
    const BL_FORT_FAB_ARG(state));
BL_FORT_PROC_DECL(CA_COMPUTE_ROS_SCT, ca_compute_ros_sct)
   (const int * lo, const int* hi,
    BL_FORT_FAB_ARG_3D(kpr),
    BL_FORT_FAB_ARG_3D(kps),
    const BL_FORT_FAB_ARG_3D(state));
BL_FORT_PROC_DECL(CA_COMPUTE_PLANCK, ca_compute_planck)
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(kpp),
    const BL_FORT_FAB_ARG(state));

BL_FORT_PROC_DECL(CA_FILT_PRIM, ca_filt_prim)
    (const int lo[], const int hi[],
     BL_FORT_FAB_ARG(Stmp),
     BL_FORT_FAB_ARG(Snew),
     BL_FORT_FAB_ARG(mask),
     const int* filter_prim_T, const int* filter_prim_S,
     const int* domlo, const int* domhi,
     const Real* delta, const Real* xlo, const Real* problo,
     const Real* time, const int* level);

// photon problems solved by the neutrino solver need these
BL_FORT_PROC_DECL(CA_COMPUTE_EMISSIVITY, ca_compute_emissivity)
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(jg),
    BL_FORT_FAB_ARG(djdT),
    const BL_FORT_FAB_ARG(temp),
    const BL_FORT_FAB_ARG(kappa_p),
    const BL_FORT_FAB_ARG(dkdT),
    const Real* PFcoef, 
    const int& use_WiensLaw, const int& integrate_Planck,
    const Real& Tf);
BL_FORT_PROC_DECL(CA_OPACS, ca_opacs)
   (const int lo[], const int hi[],
    const BL_FORT_FAB_ARG(Snew),
    const BL_FORT_FAB_ARG(T),
    const BL_FORT_FAB_ARG(Ts),
    BL_FORT_FAB_ARG(kpp),
    BL_FORT_FAB_ARG(kpr),
    BL_FORT_FAB_ARG(dkdT),
    const int* use_dkdT, const int* validStar, const int* lag_opac); 
BL_FORT_PROC_DECL(CA_COMPUTE_KAPPAS, ca_compute_kappas)
   (const int lo[], const int hi[],
    const BL_FORT_FAB_ARG(state),
    const BL_FORT_FAB_ARG(temp),
    BL_FORT_FAB_ARG(kappa_p),
    BL_FORT_FAB_ARG(kappa_r),
    BL_FORT_FAB_ARG(dkdT),
    const int* do_kappa_stm_emission, const int* use_dkdT, 
    const Real* const_kappa_p,
    const Real* kappa_p_exp_m,
    const Real* kappa_p_exp_n,
    const Real* kappa_p_exp_p,
    const Real* const_kappa_r,
    const Real* kappa_r_exp_m,
    const Real* kappa_r_exp_n,
    const Real* kappa_r_exp_p,
    const Real* const_scattering,
    const Real* scattering_exp_m,
    const Real* scattering_exp_n,
    const Real* scattering_exp_p,
    const Real* prop_temp_floor);
#ifdef MG_SU_OLSON
BL_FORT_PROC_DECL(CA_COMPUTE_RHS_SO, ca_compute_rhs_so)
   (const int* lo, const int *hi, 
    BL_FORT_FAB_ARG(rhs),
    const BL_FORT_FAB_ARG(jg),
    const BL_FORT_FAB_ARG(mugT),
    const BL_FORT_FAB_ARG(coupT),
    const BL_FORT_FAB_ARG(etaT),
    const BL_FORT_FAB_ARG(Er_step),
    const BL_FORT_FAB_ARG(rhoe_step),
    const BL_FORT_FAB_ARG(rhoe_star),    
    const Real* x,
    const Real* t, const Real* dt, const int* igroup);
#endif

// other functions that work for both neutrino and photons
BL_FORT_PROC_DECL(CA_ACCEL_CCOE, ca_accel_ccoe)
   (const int* lo, const int *hi, 
    const BL_FORT_FAB_ARG(bcgr),
    const BL_FORT_FAB_ARG(spec),
    BL_FORT_FAB_ARG(ccoe),
    const Real* dx, const int* idim, const int* igroup);
BL_FORT_PROC_DECL(CA_FLUX_FACE2CENTER, ca_flux_face2center)
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(test), 
    const BL_FORT_FAB_ARG(flx),
    const Real* x, const int* xlo, const int* xhi, 
    const int* ntest, const int* idim, const int* iflx);
BL_FORT_PROC_DECL(RESET_EINT_COMPUTE_TEMP, reset_eint_compute_temp)
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(State),
    int* resetEint, Real* relative, Real* absolute);
BL_FORT_PROC_DECL(CA_SPALPHA, ca_spalpha)
   (const int* lo, const int* hi,
    BL_FORT_FAB_ARG(SPa),
    D_DECL(BL_FORT_FAB_ARG(lamx),
	   BL_FORT_FAB_ARG(lamy),
	   BL_FORT_FAB_ARG(lamz)),
    const int* igroup);
BL_FORT_PROC_DECL(CA_RHSTOER, ca_rhstoer)
   (const int* lo, const int* hi,
    BL_FORT_FAB_ARG(rhs),
    const Real* r, const Real* dt);

// in RadDerive_?d.f90
BL_FORT_PROC_DECL(CA_DERERTOT,ca_derertot)
    (BL_FORT_FAB_ARG(der),const int* nvar,
     const BL_FORT_FAB_ARG(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const Real* delta, const Real* xlo,
     const Real* time, const Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);
BL_FORT_PROC_DECL(CA_DERENUE,ca_derenue)
    (BL_FORT_FAB_ARG(der),const int* nvar,
     const BL_FORT_FAB_ARG(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const Real* delta, const Real* xlo,
     const Real* time, const Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);
BL_FORT_PROC_DECL(CA_DERENUAE,ca_derenuae)
    (BL_FORT_FAB_ARG(der),const int* nvar,
     const BL_FORT_FAB_ARG(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const Real* delta, const Real* xlo,
     const Real* time, const Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);
// </ MGFLD>

// RadEOS
BL_FORT_PROC_DECL(CA_COMPUTE_C_V,ca_compute_c_v)
  (const int lo[], const int hi[],
   const BL_FORT_FAB_ARG(cv),
   const BL_FORT_FAB_ARG(temp),
   const BL_FORT_FAB_ARG(state));

BL_FORT_PROC_DECL(CA_GET_RHOE,ca_get_rhoe)
  (const int lo[], const int hi[],
   const BL_FORT_FAB_ARG(rhoe),
   const BL_FORT_FAB_ARG(temp),
   const BL_FORT_FAB_ARG(state));

BL_FORT_PROC_DECL(CA_COMPUTE_TEMP_GIVEN_RHOE,ca_compute_temp_given_rhoe)
     (const int lo[], const int hi[],
     BL_FORT_FAB_ARG(temp),
     const BL_FORT_FAB_ARG(state));

BL_FORT_PROC_DECL(CA_COMPUTE_TEMP_GIVEN_CV,ca_compute_temp_given_cv)
     (const int lo[], const int hi[],
      BL_FORT_FAB_ARG(temp), 
      const BL_FORT_FAB_ARG(state),
      const Real* const_c_v,
      const Real* c_v_exp_m,
      const Real* c_v_exp_n);

#ifdef NEUTRINO
BL_FORT_PROC_DECL(CA_COMPUTE_TEMP_GIVEN_REYE,ca_compute_temp_given_reye)
  (const int lo[], const int hi[],
   BL_FORT_FAB_ARG(temp),
   const BL_FORT_FAB_ARG(rhoe),
   const BL_FORT_FAB_ARG(ye),
   const BL_FORT_FAB_ARG(state));

BL_FORT_PROC_DECL(CA_COMPUTE_REYE_GIVEN_TY,ca_compute_reye_given_ty)
  (const int lo[], const int hi[],
   BL_FORT_FAB_ARG(rhoe),
   BL_FORT_FAB_ARG(rhoYe),
   const BL_FORT_FAB_ARG(temp),
   const BL_FORT_FAB_ARG(ye),
   const BL_FORT_FAB_ARG(state));
#endif

BL_FORT_PROC_DECL(CA_ER_COM2LAB, ca_er_com2lab)
  (const int lo[], const int hi[],
   const BL_FORT_FAB_ARG(Snew),
   const BL_FORT_FAB_ARG(Ecom), 
   const BL_FORT_FAB_ARG(F), int& iflx, int& nflx, 
   BL_FORT_FAB_ARG(plotvar), int& ier, int& npv);

BL_FORT_PROC_DECL(CA_COMPUTE_FCC, ca_compute_fcc)
    (const int lo[], const int hi[], 
     D_DECL(const BL_FORT_FAB_ARG(lamx),
	    const BL_FORT_FAB_ARG(lamy),
	    const BL_FORT_FAB_ARG(lamz)), int& nlam,
     BL_FORT_FAB_ARG(fcc));

BL_FORT_PROC_DECL(CA_TRANSFORM_FLUX, ca_transform_flux)
    (const int lo[], const int hi[], const Real& flag,
     const BL_FORT_FAB_ARG(Snew),
     const BL_FORT_FAB_ARG(f),
     const BL_FORT_FAB_ARG(Er), 
     const BL_FORT_FAB_ARG(Fi), int& ifi, int& nfi, 
     BL_FORT_FAB_ARG(Fo), int& ifo, int& nfo);

extern "C" {
  void FORT_SETLAYOUT(int& Density, int& Xmom, int& Eden, int& Eint, int& Temp, 
		      int& FirstSpec, int& FirstAux, int& NUM_STATE);

  void FORT_MULTRS(Real* d, intS, intS, Real* r, Real* s);

  void FORT_SPHC(Real* r, Real* s, intS, const Real* dx);

  void FORT_SPHE(Real* r, Real* s, const int& idim, intS, const Real* dx);

  void FORT_LACOEF(Real* acoefs, intS, intS,
		   Real* fkp, Real* eta, Real* etainv,
		   Real* r, Real* s,
		   const Real& c, const Real& dt, const Real& theta);

  Real FORT_KAVG(const Real& a, const Real& b, const Real& d, const int& iopt);

  void FORT_BCLIM(Real* bcoefs, Real* lambda, intS, intS, const int& idim,
                  Real* kappa_r, intS,
                  Real* r, Real* s,
                  const Real& c, const Real* dx);

  void FORT_FLXLIM(Real* lambda, intS, intS, const int& limiter);

  void FORT_EDDFAC(Real* efact, intS, intS, const int& limiter, const int& idim);

  void FORT_SCGRD1(Real* R, intS, intS, const int& idim,
		   Real* kappa_r, intS, Real* Er, const Real* dx);

  void FORT_SCGRD2(Real* R, intS, intS, const int& idim,
		   Real* kappa_r, intS, Real* Er, 
#if (BL_SPACEDIM >= 2)
		   intS, Real* dtmp1, 
#endif
#if (BL_SPACEDIM == 3)
		   Real* dtmp2,
#endif
                   const Real* dx);

  void FORT_SCGRD3(Real* R, intS, intS, const int& idim,
		   Real* kappa_r, intS, Real* Er, 
#if (BL_SPACEDIM >= 2)
		   intS, Real* dtmp1, 
#endif
#if (BL_SPACEDIM == 3)
		   Real* dtmp2,
#endif
                   const Real* dx);

  void FORT_LRHS(Real* rhs, intS, intS,
		 Real* temp, Real* fkp, Real* eta, Real* etainv,
		 Real* rhoem, Real* rhoes,
		 Real* dflux_old, Real* Er_old, intS, Real* Edot,
		 Real* r, Real* s, const Real& delta_t, const Real& sigma,
		 const Real& c, const Real& theta);

  void FORT_ANATW2(Real* test, intS,
		   Real* temp,
		   const Real& p, const Real& xf, const Real& Tc,
		   const Real* dx, const Real* xlo, const int* lo);

  void FORT_CFRHOE(intS, Real* frhoe, intS, Real* state, intS);

  void FORT_GTEMP(intS,
		  Real* temp, intS,
		  Real* con, Real* em, Real* en,
                  Real* state, intS);

  void FORT_GCV(intS,
                Real* c_v, intS, const Real* temp, intS,
		const Real* con, const Real* em, const Real* en, const Real* tf,
		const Real* state, intS);

  void FORT_CEXCH(intS,
		  Real* exch, intS,
		  Real* er, intS, 
		  Real* fkp, intS,
		  const Real& sigma, const Real& c);

  void FORT_CETA2(intS,
		  Real* eta, Real* etainv, intS,
		  Real* frho, intS, 
		  Real* temp, intS,
		  Real* c_v, intS,
		  Real* fkp, intS,
		  Real* er, intS,
		  const Real& dT, const Real& delta_t, const Real& sigma,
		  const Real& c, const Real& underr, const int& lagpla);

  void FORT_CEUP(intS, Real& relres, Real& absres,
		 Real* frhoes, intS,
		 Real* frhoem, Real* eta, Real* etainv,
		 Real* dflux_old, Real* dflux_new,
		 Real* exch,
		 const Real& delta_t, const Real& theta);

  void FORT_CEUPDTERM(intS, Real& relres, Real& absres,
		      Real* frhoes, intS,
		      Real* frhoem, Real* eta, Real* etainv,
		      Real* dflux_old, Real* dflux_new,
		      Real* exch, Real* dterm, 
		      const Real& delta_t, const Real& theta);

  void FORT_NCEUP(intS, Real& relres, Real& absres,
		  Real* frhoes, intS,
		  Real* frhoem, Real* eta, Real* etainv,
		  Real* er, intS, Real* dflux_old, Real* dflux_new,
		  Real* temp, Real* fkp, Real* cv, Real* state, intS,
		  const Real& sigma, const Real& c,
		  const Real& delta_t, const Real& theta);

  void FORT_CETOT(intS, Real* state, intS, Real* rhoe, intS);

  void FORT_FKPN(intS,
		 Real* fkp, intS,
		 Real* con, Real* em, Real* en,
		 Real* ep, const Real& xnu, Real* tf,
		 Real* temp, intS, 
		 Real* state, intS);

  void FORT_ROSSE1(intS,
		   Real* kappar, intS,
		   Real* con, Real* em, Real* en,
		   Real* ep, const Real& xnu,
		   Real* tf, const Real& kfloor,
		   Real* temp, intS, Real* state, intS);

  void FORT_ROSSE1S(intS,
		    Real* kappar, intS,
		    Real* con, Real* em, Real* en,
		    Real* ep, 
		    Real* scon, Real* sem, Real* sen,
		    Real* sep, 
		    const Real& xnu,
		    Real* tf, const Real& kfloor,
		    Real* temp, intS, Real* state, intS);

  void FORT_NFLOOR(Real* dest, intS, intS,
		   int& nflr, const Real& flr, const int& nvar);

  void FORT_LACOEFMGFLD(Real* a, intS, intS,
			Real* kappa, intS,
			Real* r, Real* s,
			const Real& delta_t, const Real& c);

  void FORT_RFFACE(Real* fine, intS, Real* cres, intS,
                   const int& idim, const int* irat);

  void FORT_RADBNDRY(Real* bf, intS, intS,
		     const Real* dx, const Real* xlo, const Real& time,
                     const int& dir, const int& face);

  void FORT_RADBNDRY2(Real* bf, intS, int* tfab, intS,
		      const Real* dx, const Real* xlo, const Real& time);

  void FORT_BEXTRP(Real* f, intS, intS);

  // neutrino routines below

  void FORT_INIT_OPACITY_TABLE(const int& iverb);

  void FORT_LBCOEFNA(Real* bcoef, Real* bcgrp, intS, intS,
                     Real* spec, intS,
                     const int& idim);

  void FORT_LJUPNA(Real* jnew, intS, intS,
                   Real* spec, intS,
                   Real* accel, intS,
		   const int& nTotal);
};

#endif
