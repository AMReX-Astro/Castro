#ifndef _RAD_F_H_
#define _RAD_F_H_

#include <AMReX_ArrayLim.H>

#ifdef BL_FORT_USE_UNDERSCORE
#  define   FORT_KAVG    kavg_
#  define   FORT_RADBNDRY  rbndry_
#  define   FORT_RADBNDRY2 rbndr2_
#  define   FORT_INIT_OPACITY_TABLE init_opacity_table_

#elif defined(BL_FORT_USE_UPPERCASE)

#  define   FORT_KAVG    KAVG
#  define   FORT_RADBNDRY  RBNDRY
#  define   FORT_RADBNDRY2 RBNDR2
#  define   FORT_INIT_OPACITY_TABLE INIT_OPACITY_TABLE

#endif

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
void ca_init_radhydro_pars
    (const int & fspace_advection_type, const int & do_inelastic_scattering,
     const int & comoving, const amrex::Real & flatten_pp_threshold);

void ca_initradconstants
  (const amrex::Real& pi, amrex::Real& clight, amrex::Real& hplanck,
   amrex::Real& kboltz, amrex::Real& stefbol, amrex::Real& avogadro,
   amrex::Real& mev2erg, int& J_is_used);

void ca_initsinglegroup
  (const int& ngroups);

#ifdef __cplusplus
}
#endif


BL_FORT_PROC_DECL(CA_INITGROUPS,ca_initgroups)
  (const amrex::Real* nugroup, const amrex::Real* dnugroup,
   const int& ngroups, const int& ng0, const int& ng1);

BL_FORT_PROC_DECL(CA_INITGROUPS2,ca_initgroups2)
  (const amrex::Real* nugroup, const amrex::Real* dnugroup,
   const amrex::Real* xnu, const int& ngroups);

BL_FORT_PROC_DECL(CA_INITGROUPS3,ca_initgroups3)
  (const amrex::Real* nugroup, const amrex::Real* dnugroup, const amrex::Real* dlognugroup,
   const amrex::Real* xnu, const int& ngroups, const int& ng0, const int& ng1);

#ifdef NEUTRINO
#ifdef __cplusplus
extern "C"
{
#endif
  void ca_derneut
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
     const int* level, const int* grid_no);

  void ca_derrhoyl
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
     const int* level, const int* grid_no);

  void ca_deryl
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
     const int* level, const int* grid_no);

  void ca_derynue
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
     const int* level, const int* grid_no);

  void ca_derynuae
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
     const int* level, const int* grid_no);
#ifdef __cplusplus
}
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif
  void ca_compute_powerlaw_kappa
     (const int lo[], const int hi[],
      BL_FORT_FAB_ARG(kappa), const BL_FORT_FAB_ARG(u), 
      const amrex::Real* kappa0, const amrex::Real* m, const amrex::Real* n,  const amrex::Real* p, 
      const amrex::Real* Tfloor, const amrex::Real* kfloor);

  void ca_compute_powerlaw_kappa_s
     (const int lo[], const int hi[],
      BL_FORT_FAB_ARG(kappa), const BL_FORT_FAB_ARG(u), 
      const amrex::Real* kappa0, const amrex::Real* m, const amrex::Real* n, const amrex::Real* p, 
      const amrex::Real* s0, const amrex::Real* sm, const amrex::Real* sn, const amrex::Real* sp, 
      const amrex::Real* Tfloor, const amrex::Real* kfloor);
#ifdef __cplusplus
}
#endif

BL_FORT_PROC_DECL(CA_COMPUTE_KAPKAP, ca_compute_kapkap)
     (BL_FORT_FAB_ARG(kapkap), const BL_FORT_FAB_ARG(kap_r)); 

BL_FORT_PROC_DECL(CA_COMPUTE_LAMBORDER, ca_compute_lamborder)
   (const BL_FORT_FAB_ARG(Er), const BL_FORT_FAB_ARG(kap), 
    BL_FORT_FAB_ARG(lam),  const amrex::Real* dx, const int* ngrow, const int* limiter,
    const int* filter_lambda_T, const int* filter_lambda_S);

#ifdef __cplusplus
extern "C" {
#endif
  void ca_get_v_dcf
   (const int lo[], const int hi[],
    const BL_FORT_FAB_ARG(er), 
    const BL_FORT_FAB_ARG(state), 
    const BL_FORT_FAB_ARG(temp),
    const BL_FORT_FAB_ARG(c_v),
    const BL_FORT_FAB_ARG(kr),
    const BL_FORT_FAB_ARG(kp),
    const BL_FORT_FAB_ARG(kp2),
    const amrex::Real* dtemp, const amrex::Real* dtime, const amrex::Real* sigma, const amrex::Real* c,
    BL_FORT_FAB_ARG(v),
    BL_FORT_FAB_ARG(dcf));

  void ca_update_dcf
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(dcf),
    const BL_FORT_FAB_ARG(etainv),
    const BL_FORT_FAB_ARG(kp),
    const BL_FORT_FAB_ARG(kr));

  void ca_compute_dcoefs
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(d), 
    const BL_FORT_FAB_ARG(lam), 
    const BL_FORT_FAB_ARG(v), const BL_FORT_FAB_ARG(dcf),
    const amrex::Real* r, const int* idir);

#ifdef __cplusplus
}
#endif
 
BL_FORT_PROC_DECL(CA_EST_GPR0, ca_est_gpr0)
   (const BL_FORT_FAB_ARG(Er),
    BL_FORT_FAB_ARG(gPr));

BL_FORT_PROC_DECL(CA_EST_GPR2, ca_est_gpr2)
   (const BL_FORT_FAB_ARG(kap),
    const BL_FORT_FAB_ARG(Er),
    BL_FORT_FAB_ARG(gPr), 
    const int lo[], const int hi[],
    const amrex::Real* dx, const int* limiter, const int* comoving); 

#ifdef __cplusplus
extern "C" {
#endif
  void ca_set_dterm_face
   (const int lo[], const int hi[],
    const BL_FORT_FAB_ARG(Er),
    const BL_FORT_FAB_ARG(dcoef), 
    BL_FORT_FAB_ARG(dterm_face), 
    const amrex::Real* dx, const int* idir);

  void ca_face2center
    (const int lo[], const int hi[], 
     const int& scomp, const int& dcomp, const int& ncomp, const int& nf, const int& nc,
     D_DECL(const BL_FORT_FAB_ARG(foox),
	    const BL_FORT_FAB_ARG(fooy),
	    const BL_FORT_FAB_ARG(fooz)),
     BL_FORT_FAB_ARG(fooc));

  void ca_correct_dterm
    (D_DECL(BL_FORT_FAB_ARG(dfx),
	    BL_FORT_FAB_ARG(dfy),
	    BL_FORT_FAB_ARG(dfz)),
     const amrex::Real* re, const amrex::Real* rc);
#ifdef __cplusplus
}
#endif

// <MGFLD>
#ifdef __cplusplus
extern "C" {
#endif
void ca_initfluxlimiter
   (const int* limiter, const int* closure);
#ifdef __cplusplus
}
#endif

#ifdef NEUTRINO
BL_FORT_PROC_DECL(CA_ACCEL_ACOE_NEUT, ca_accel_acoe_neut)
   (const int lo[], const int hi[],
    const BL_FORT_FAB_ARG(eta1),
    const BL_FORT_FAB_ARG(thetaT),
    const BL_FORT_FAB_ARG(thetaY),
    const BL_FORT_FAB_ARG(spc),
    const BL_FORT_FAB_ARG(kap),
    BL_FORT_FAB_ARG(aco),
    const amrex::Real* dt, const amrex::Real* ptc_tau);
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
    const amrex::Real* dt);
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
    const amrex::Real* dt, const amrex::Real* ptc_tau);
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
    amrex::Real* rel, amrex::Real* abs, amrex::Real *errr, const amrex::Real* dt);
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
    amrex::Real* rel_rhoe, amrex::Real* abs_rhoe, 
    amrex::Real* rel_FT  , amrex::Real* abs_FT, 
    amrex::Real* rel_T   , amrex::Real* abs_T, 
    amrex::Real* rel_FY  , amrex::Real* abs_FY, 
    amrex::Real* rel_Ye  , amrex::Real* abs_Ye,
    const amrex::Real* delta_t);
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
     const amrex::Real* delta_t, const amrex::Real* tau);

#ifdef __cplusplus
extern "C" {
#endif
  void ca_compute_rhs_neut
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
    const amrex::Real* r,
    const amrex::Real* dt, const int* igroup, const amrex::Real* ptc_tau);
#ifdef __cplusplus
}
#endif

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
    const amrex::Real* dt, const amrex::Real* ptc_tau);
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
    amrex::Real* derat, amrex::Real* dTrat, amrex::Real* dye);
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
    const amrex::Real* dt, const amrex::Real* ptc_tau);
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
    const amrex::Real* dt);

#else   // photons

BL_FORT_PROC_DECL(CA_ACCEL_ACOE, ca_accel_acoe)
   (const int lo[], const int hi[],
    const BL_FORT_FAB_ARG(eta1),
    const BL_FORT_FAB_ARG(spc),
    const BL_FORT_FAB_ARG(kap),
    BL_FORT_FAB_ARG(aco),
    const amrex::Real* dt, const amrex::Real* ptc_tau);
BL_FORT_PROC_DECL(CA_ACCEL_RHS, ca_accel_rhs)
   (const int lo[], const int hi[],
    const BL_FORT_FAB_ARG(Ern),
    const BL_FORT_FAB_ARG(Erl),
    const BL_FORT_FAB_ARG(kap),
    const BL_FORT_FAB_ARG(etaT),
    BL_FORT_FAB_ARG(rhs),
    const amrex::Real* dt);
BL_FORT_PROC_DECL(CA_ACCEL_SPEC, ca_accel_spec)
   (const int lo[], const int hi[],
    const BL_FORT_FAB_ARG(kap),
    const BL_FORT_FAB_ARG(mugT),
    BL_FORT_FAB_ARG(spec), 
    const amrex::Real* dt, const amrex::Real* ptc_tau);
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
    amrex::Real* rel_rhoe, amrex::Real* abs_rhoe, 
    amrex::Real* rel_FT  , amrex::Real* abs_FT, 
    amrex::Real* rel_T   , amrex::Real* abs_T, 
    const amrex::Real* delta_t);
BL_FORT_PROC_DECL(CA_CHECK_CONV_ER, ca_check_conv_er)
   (const int lo[], const int hi[],
    const BL_FORT_FAB_ARG(Ern),
    const BL_FORT_FAB_ARG(Erl),
    const BL_FORT_FAB_ARG(kap),
    const BL_FORT_FAB_ARG(etTz),
    const BL_FORT_FAB_ARG(temp),
    amrex::Real* rel, amrex::Real* abs, amrex::Real *errr, const amrex::Real* dt);
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
     const amrex::Real* delta_t, const amrex::Real* tau);

#ifdef __cplusplus
extern "C" {
#endif
  void ca_compute_rhs
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
    const amrex::Real* r,
    const amrex::Real* dt, const int* igroup, const amrex::Real* ptc_tau);
#ifdef __cplusplus
}
#endif

BL_FORT_PROC_DECL(CA_LOCAL_ACCEL, ca_local_accel)
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(Ern),
    const BL_FORT_FAB_ARG(Erl),
    const BL_FORT_FAB_ARG(kap),
    const BL_FORT_FAB_ARG(etaT),
    const BL_FORT_FAB_ARG(mugT),
    const amrex::Real* dt, const amrex::Real* ptc_tau);
BL_FORT_PROC_DECL(CA_STATE_UPDATE, ca_state_update)
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(state),
    const BL_FORT_FAB_ARG(rhoe),
    const BL_FORT_FAB_ARG(temp),
    const BL_FORT_FAB_ARG(msk),
    amrex::Real* derat, amrex::Real* dTrat);
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
    const amrex::Real* dt, const amrex::Real* ptc_tau);
BL_FORT_PROC_DECL(CA_NCUPDATE_MATTER, ca_ncupdate_matter)
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(temp_new),
    const BL_FORT_FAB_ARG(Er_new),
    const BL_FORT_FAB_ARG(rhoe_star),
    const BL_FORT_FAB_ARG(rhoe_step),
    const BL_FORT_FAB_ARG(etaTz),
    const BL_FORT_FAB_ARG(kappa_p),
    const BL_FORT_FAB_ARG(jg),
    const amrex::Real* dt);

#endif

#ifdef __cplusplus
extern "C" {
#endif
  void ca_compute_rosseland
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(kpr),
     const BL_FORT_FAB_ARG_3D(state));
#ifdef __cplusplus
}
#endif

BL_FORT_PROC_DECL(CA_COMPUTE_SCATTERING, ca_compute_scattering)
   (const int * lo, const int* hi,
    BL_FORT_FAB_ARG_3D(kps),
    const BL_FORT_FAB_ARG_3D(state));
BL_FORT_PROC_DECL(CA_COMPUTE_SCATTERING_2, ca_compute_scattering_2)
   (const int * lo, const int* hi,
    BL_FORT_FAB_ARG_3D(kps),
    const BL_FORT_FAB_ARG_3D(state),
    const amrex::Real* k0_p, const amrex::Real* m_p, const amrex::Real* n_p,
    const amrex::Real* k0_r, const amrex::Real* m_r, const amrex::Real* n_r,
    const amrex::Real* Tfloor, const amrex::Real* kfloor);

#ifdef __cplusplus
extern "C" {
#endif
  void ca_compute_planck
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(kpp),
     const BL_FORT_FAB_ARG_3D(state),
     amrex::Real temp_offset);

  void ca_filt_prim
    (const int lo[], const int hi[],
     BL_FORT_FAB_ARG(Stmp),
     BL_FORT_FAB_ARG(Snew),
     BL_FORT_FAB_ARG(mask),
     const int* filter_prim_T, const int* filter_prim_S,
     const int* domlo, const int* domhi,
     const amrex::Real* delta, const amrex::Real* xlo, const amrex::Real* problo,
     const amrex::Real* time, const int* level);
#ifdef __cplusplus
}
#endif

// photon problems solved by the neutrino solver need these
BL_FORT_PROC_DECL(CA_COMPUTE_EMISSIVITY, ca_compute_emissivity)
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(jg),
    BL_FORT_FAB_ARG(djdT),
    const BL_FORT_FAB_ARG(temp),
    const BL_FORT_FAB_ARG(kappa_p),
    const BL_FORT_FAB_ARG(dkdT),
    const amrex::Real* PFcoef, 
    const int& use_WiensLaw, const int& integrate_Planck,
    const amrex::Real& Tf);
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
    const amrex::Real* const_kappa_p,
    const amrex::Real* kappa_p_exp_m,
    const amrex::Real* kappa_p_exp_n,
    const amrex::Real* kappa_p_exp_p,
    const amrex::Real* const_kappa_r,
    const amrex::Real* kappa_r_exp_m,
    const amrex::Real* kappa_r_exp_n,
    const amrex::Real* kappa_r_exp_p,
    const amrex::Real* const_scattering,
    const amrex::Real* scattering_exp_m,
    const amrex::Real* scattering_exp_n,
    const amrex::Real* scattering_exp_p,
    const amrex::Real* prop_temp_floor);
#ifdef MG_SU_OLSON
#ifdef __cplusplus
extern "C" {
#endif
  void ca_compute_rhs_so
   (const int* lo, const int *hi, 
    BL_FORT_FAB_ARG(rhs),
    const BL_FORT_FAB_ARG(jg),
    const BL_FORT_FAB_ARG(mugT),
    const BL_FORT_FAB_ARG(coupT),
    const BL_FORT_FAB_ARG(etaT),
    const BL_FORT_FAB_ARG(Er_step),
    const BL_FORT_FAB_ARG(rhoe_step),
    const BL_FORT_FAB_ARG(rhoe_star),    
    const amrex::Real* x,
    const amrex::Real* t, const amrex::Real* dt, const int* igroup);
#ifdef __cplusplus
}
#endif
#endif

// other functions that work for both neutrino and photons
BL_FORT_PROC_DECL(CA_ACCEL_CCOE, ca_accel_ccoe)
   (const int* lo, const int *hi, 
    const BL_FORT_FAB_ARG(bcgr),
    const BL_FORT_FAB_ARG(spec),
    BL_FORT_FAB_ARG(ccoe),
    const amrex::Real* dx, const int* idim, const int* igroup);


#ifdef __cplusplus
extern "C"
{
#endif

  void ca_flux_face2center
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(test), 
    const BL_FORT_FAB_ARG(flx),
    const amrex::Real* x, const int* xlo, const int* xhi, 
    const int* ntest, const int* idim, const int* iflx);

  void ca_spalpha
   (const int* lo, const int* hi,
    BL_FORT_FAB_ARG(SPa),
    D_DECL(BL_FORT_FAB_ARG(lamx),
	   BL_FORT_FAB_ARG(lamy),
	   BL_FORT_FAB_ARG(lamz)),
    const int* igroup);
#ifdef __cplusplus
}
#endif

BL_FORT_PROC_DECL(CA_RHSTOER, ca_rhstoer)
   (const int* lo, const int* hi,
    BL_FORT_FAB_ARG(rhs),
    const amrex::Real* r, const amrex::Real* dt);

// in RadDerive_?d.f90
#ifdef __cplusplus
extern "C"
{
#endif
  void ca_derertot
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);
  void ca_derenue
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);
  void ca_derenuae
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const amrex::Real* delta, const amrex::Real* xlo,
     const amrex::Real* time, const amrex::Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);
#ifdef __cplusplus
}
#endif
// </ MGFLD>

// RadEOS
#ifdef __cplusplus
extern "C" 
{
#endif
  void ca_compute_c_v
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(cv),
     const BL_FORT_FAB_ARG_3D(temp),
     const BL_FORT_FAB_ARG_3D(state));

  void ca_get_rhoe
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(rhoe),
     const BL_FORT_FAB_ARG_3D(temp),
     const BL_FORT_FAB_ARG_3D(state));

  void ca_compute_temp_given_rhoe
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(temp),
     const BL_FORT_FAB_ARG_3D(state),
     int update_state);

  void ca_compute_temp_given_cv
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(temp), 
     const BL_FORT_FAB_ARG_3D(state),
     const amrex::Real const_c_v,
     const amrex::Real c_v_exp_m,
     const amrex::Real c_v_exp_n);

#ifdef NEUTRINO
  ca_compute_temp_given_reye
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(temp),
     const BL_FORT_FAB_ARG_3D(rhoe),
     const BL_FORT_FAB_ARG_3D(ye),
     const BL_FORT_FAB_ARG_3D(state));

  ca_compute_reye_given_ty
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(rhoe),
     BL_FORT_FAB_ARG_3D(rhoYe),
     const BL_FORT_FAB_ARG_3D(temp),
     const BL_FORT_FAB_ARG_3D(ye),
     const BL_FORT_FAB_ARG_3D(state));
#endif

#ifdef __cplusplus
}
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
    (const int lo[], const int hi[], const amrex::Real& flag,
     const BL_FORT_FAB_ARG(Snew),
     const BL_FORT_FAB_ARG(f),
     const BL_FORT_FAB_ARG(Er), 
     const BL_FORT_FAB_ARG(Fi), int& ifi, int& nfi, 
     BL_FORT_FAB_ARG(Fo), int& ifo, int& nfo);

#ifdef __cplusplus
extern "C" {
#endif
  void multrs(BL_FORT_FAB_ARG(cc), 
	      ARLIM_P(reglo), ARLIM_P(reghi), amrex::Real* r, amrex::Real* s);

  void sphc(amrex::Real* r, amrex::Real* s, ARLIM_P(reglo), ARLIM_P(reghi), const amrex::Real* dx);

  void sphe(amrex::Real* r, amrex::Real* s, const int& idim, 
	    ARLIM_P(edgeboxlo), ARLIM_P(edgeboxhi), const amrex::Real* dx);

  void lacoef(const int* lo, const int* hi,
              BL_FORT_FAB_ARG_3D(acoefs), 
	      BL_FORT_FAB_ARG_3D(fkp),
              BL_FORT_FAB_ARG_3D(eta),
              BL_FORT_FAB_ARG_3D(etainv),
              const amrex::Real* dx,
	      amrex::Real c, amrex::Real dt, amrex::Real theta);

  amrex::Real FORT_KAVG(const amrex::Real& a, const amrex::Real& b, const amrex::Real& d, const int& iopt);

  void bclim(amrex::Real* bcoefs, 
	     BL_FORT_FAB_ARG(lambda), 
	     ARLIM_P(reglo), ARLIM_P(reghi), const int& idim,
	     BL_FORT_FAB_ARG(kappa_r), 
	     amrex::Real* r, amrex::Real* s,
	     const amrex::Real& c, const amrex::Real* dx);

  void flxlim(BL_FORT_FAB_ARG(lambda), 
	      ARLIM_P(reglo), ARLIM_P(reghi), const int& limiter);

  void eddfac(amrex::Real* efact, ARLIM_P(rlo), ARLIM_P(rhi),
	      ARLIM_P(reglo), ARLIM_P(reghi), const int& limiter, const int& idim);

  void scgrd(const int* lo, const int* hi,
             BL_FORT_FAB_ARG_3D(R), 
             int idim, const amrex::Real* dx,
             BL_FORT_FAB_ARG_3D(kappa_r), 
             BL_FORT_FAB_ARG_3D(Er),
             int include_cross_terms);

  void lrhs(BL_FORT_FAB_ARG(rhs), 
	    ARLIM_P(reglo), ARLIM_P(reghi),
	    amrex::Real* temp, amrex::Real* fkp, amrex::Real* eta, amrex::Real* etainv,
	    amrex::Real* rhoem, amrex::Real* rhoes,
	    amrex::Real* dflux_old, 
	    BL_FORT_FAB_ARG(Er_old), 
	    amrex::Real* Edot,
	    amrex::Real* r, amrex::Real* s, const amrex::Real& delta_t, const amrex::Real& sigma,
	    const amrex::Real& c, const amrex::Real& theta);

  void anatw2(amrex::Real* test, ARLIM_P(reglo), ARLIM_P(reghi),
	      amrex::Real* temp,
	      const amrex::Real& p, const amrex::Real& xf, const amrex::Real& Tc,
	      const amrex::Real* dx, const amrex::Real* xlo, const int* lo);

  void cfrhoe(const int* lo, const int* hi, 
	      BL_FORT_FAB_ARG_3D(frhoe), 
	      BL_FORT_FAB_ARG_3D(state));

  void gtemp(const int* lo, const int* hi,
	     BL_FORT_FAB_ARG_3D(temp), 
	     amrex::Real con, amrex::Real em, amrex::Real en,
	     BL_FORT_FAB_ARG_3D(state),
             int update_state);

  void gcv(const int* lo, const int* hi,
	   BL_FORT_FAB_ARG_3D(c_v), 
	   const BL_FORT_FAB_ARG_3D(temp), 
	   amrex::Real con, amrex::Real em, amrex::Real en,
           amrex::Real tf,
	   const BL_FORT_FAB_ARG_3D(state));

  void cexch(ARLIM_P(reglo), ARLIM_P(reghi),
	     BL_FORT_FAB_ARG(exch), 
	     BL_FORT_FAB_ARG(Er), 
	     BL_FORT_FAB_ARG(fkp), 
	     const amrex::Real& sigma, const amrex::Real& c);

  void ceta2(const int* lo, const int* hi,
	     BL_FORT_FAB_ARG_3D(eta),
             BL_FORT_FAB_ARG_3D(etainv),
	     BL_FORT_FAB_ARG_3D(frho), 
	     BL_FORT_FAB_ARG_3D(temp), 
	     BL_FORT_FAB_ARG_3D(c_v), 
	     BL_FORT_FAB_ARG_3D(fkp), 
	     BL_FORT_FAB_ARG_3D(er), 
	     amrex::Real dT, amrex::Real delta_t, amrex::Real sigma,
	     amrex::Real c, amrex::Real underr, int lagpla);

  void ceup(ARLIM_P(reglo), ARLIM_P(reghi), amrex::Real& relres, amrex::Real& absres,
	    BL_FORT_FAB_ARG(frhoes), 
	    amrex::Real* frhoem, amrex::Real* eta, amrex::Real* etainv,
	    amrex::Real* dflux_old, amrex::Real* dflux_new,
	    amrex::Real* exch,
	    const amrex::Real& delta_t, const amrex::Real& theta);

  void ceupdterm(ARLIM_P(reglo), ARLIM_P(reghi), 
		 amrex::Real& relres, amrex::Real& absres,
		 BL_FORT_FAB_ARG(frhoes), 
		 amrex::Real* frhoem, amrex::Real* eta, amrex::Real* etainv,
		 amrex::Real* dflux_old, amrex::Real* dflux_new,
		 amrex::Real* exch, amrex::Real* dterm, 
		 const amrex::Real& delta_t, const amrex::Real& theta);

  void nceup(ARLIM_P(reglo), ARLIM_P(reghi), 
	     amrex::Real& relres, amrex::Real& absres,
	     BL_FORT_FAB_ARG(frhoes), 
	     amrex::Real* frhoem, amrex::Real* eta, amrex::Real* etainv,
	     BL_FORT_FAB_ARG(er), 
	     amrex::Real* dflux_old, amrex::Real* dflux_new,
	     amrex::Real* temp, amrex::Real* fkp, amrex::Real* cv, 
	     BL_FORT_FAB_ARG(state), 
	     const amrex::Real& sigma, const amrex::Real& c,
	     const amrex::Real& delta_t, const amrex::Real& theta);

  void cetot(ARLIM_P(reglo), ARLIM_P(reghi), 
	     BL_FORT_FAB_ARG(state), 
	     BL_FORT_FAB_ARG(rhoe));

  void fkpn(const int* lo, const int* hi,
	    BL_FORT_FAB_ARG_3D(fkp),	   
	    amrex::Real con, amrex::Real em, amrex::Real en,
	    amrex::Real ep, amrex::Real xnu, amrex::Real tf,
	    BL_FORT_FAB_ARG_3D(temp), 
	    BL_FORT_FAB_ARG_3D(state),
            amrex::Real temp_offset);

  void rosse1(const int* lo, const int* hi,
	      amrex::Real con, amrex::Real em, amrex::Real en,
	      amrex::Real ep, amrex::Real nu,
              amrex::Real tf, amrex::Real kfloor,
	      BL_FORT_FAB_ARG_3D(state),
      	      BL_FORT_FAB_ARG_3D(kappar));

  void rosse1s(const int* lo, const int* hi,
	       amrex::Real con, amrex::Real em, amrex::Real en,
	       amrex::Real ep,  amrex::Real scon,
               amrex::Real sem, amrex::Real sen,
	       amrex::Real sep, amrex::Real nu,
	       amrex::Real tf, amrex::Real kfloor,
	       BL_FORT_FAB_ARG_3D(state),
      	       BL_FORT_FAB_ARG_3D(kappar));
  
  void nfloor(const int* lo, const int* hi,
              BL_FORT_FAB_ARG_3D(dest), 
	      amrex::Real flr, int nvar);

  void lacoefmgfld(BL_FORT_FAB_ARG(a), 
		   ARLIM_P(reglo), ARLIM_P(reghi),
		   amrex::Real* kappa, ARLIM_P(klo), ARLIM_P(khi),
		   amrex::Real* r, amrex::Real* s,
		   const amrex::Real& delta_t, const amrex::Real& c);

  void rfface(BL_FORT_FAB_ARG(fine), 
	      BL_FORT_FAB_ARG(crse), 
	      const int& idim, const int* irat);

  void FORT_RADBNDRY(amrex::Real* bf, ARLIM_P(blo), ARLIM_P(bhi),
		     ARLIM_P(dlo), ARLIM_P(dhi), 
		     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real& time,
                     const int& dir, const int& face);

  void FORT_RADBNDRY2(amrex::Real* bf, ARLIM_P(blo), ARLIM_P(bhi), 
		      int* tfab, ARLIM_P(dlo), ARLIM_P(dhi),
		      const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real& time);

  void bextrp(BL_FORT_FAB_ARG(f), 
	      ARLIM_P(reglo), ARLIM_P(reghi));

  // neutrino routines below

  void FORT_INIT_OPACITY_TABLE(const int& iverb);

  void lbcoefna(amrex::Real* bcoef, amrex::Real* bcgrp, 
		ARLIM_P(blo), ARLIM_P(bhi), 
		ARLIM_P(bxlo), ARLIM_P(bxhi),
		BL_FORT_FAB_ARG(spec),
		const int& idim);

  void ljupna(BL_FORT_FAB_ARG(jnew), 
	      ARLIM_P(reglo), ARLIM_P(reghi), 
	      BL_FORT_FAB_ARG(spec),
	      BL_FORT_FAB_ARG(accel), 
	      const int& nTotal);
#ifdef __cplusplus
};
#endif

#endif
