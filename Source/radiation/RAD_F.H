#ifndef CASTRO_RAD_F_H
#define CASTRO_RAD_F_H

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
     const int & comoving);

void ca_initradconstants
  (const amrex::Real& pi, amrex::Real& clight, amrex::Real& hplanck,
   amrex::Real& kboltz, amrex::Real& stefbol, amrex::Real& avogadro,
   amrex::Real& mev2erg, int& J_is_used);

void ca_initsinglegroup
  (const int& ngroups);

void ca_get_dlognu
  (amrex::Real* dlognu_out);

void ca_get_nugroup
  (amrex::Real* nugroup);

void ca_get_dnugroup
  (amrex::Real* dnugroup);

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

  void ca_compute_dcoefs
   (const int* lo, const int* hi,
    BL_FORT_FAB_ARG_3D(d), 
    const BL_FORT_FAB_ARG_3D(lam), 
    const BL_FORT_FAB_ARG_3D(v),
    const BL_FORT_FAB_ARG_3D(dcf),
    const amrex::Real* dx, int idir);

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

extern "C"
{
  void ca_accel_acoe
   (const int* lo, const int* hi,
    const BL_FORT_FAB_ARG_3D(eta1),
    const BL_FORT_FAB_ARG_3D(spc),
    const BL_FORT_FAB_ARG_3D(kap),
    BL_FORT_FAB_ARG_3D(aco),
    amrex::Real dt, amrex::Real ptc_tau);

  void ca_accel_rhs
   (const int* lo, const int* hi,
    const BL_FORT_FAB_ARG_3D(Ern),
    const BL_FORT_FAB_ARG_3D(Erl),
    const BL_FORT_FAB_ARG_3D(kap),
    const BL_FORT_FAB_ARG_3D(etaT),
    BL_FORT_FAB_ARG_3D(rhs),
    amrex::Real dt);

  void ca_accel_spec
   (const int* lo, const int* hi,
    const BL_FORT_FAB_ARG_3D(kap),
    const BL_FORT_FAB_ARG_3D(mugT),
    BL_FORT_FAB_ARG_3D(spec), 
    amrex::Real dt, amrex::Real ptc_tau);

  void ca_check_conv_er
   (const int* lo, const int* hi,
    const BL_FORT_FAB_ARG_3D(Ern),
    const BL_FORT_FAB_ARG_3D(Erl),
    const BL_FORT_FAB_ARG_3D(kap),
    const BL_FORT_FAB_ARG_3D(etTz),
    const BL_FORT_FAB_ARG_3D(temp),
    amrex::Real* rel, amrex::Real* abs, amrex::Real *errr,
    amrex::Real dt);

  void ca_compute_etat
    (const int* lo, const int* hi, 
     BL_FORT_FAB_ARG_3D(etat), 
     BL_FORT_FAB_ARG_3D(etaTz), 
     BL_FORT_FAB_ARG_3D(eta1), 
     BL_FORT_FAB_ARG_3D(djdT),
     const BL_FORT_FAB_ARG_3D(dkdT),
     const BL_FORT_FAB_ARG_3D(dedT),
     const BL_FORT_FAB_ARG_3D(Er_star),
     const BL_FORT_FAB_ARG_3D(rho),
     amrex::Real delta_t, amrex::Real tau);

  void ca_compute_rhs
    (const int* lo, const int* hi, 
     BL_FORT_FAB_ARG_3D(rhs),
     const BL_FORT_FAB_ARG_3D(jg),
     const BL_FORT_FAB_ARG_3D(mugT),
     const BL_FORT_FAB_ARG_3D(coupT),
     const BL_FORT_FAB_ARG_3D(etaT),
     const BL_FORT_FAB_ARG_3D(Er_step),
     const BL_FORT_FAB_ARG_3D(rhoe_step),
     const BL_FORT_FAB_ARG_3D(Er_star),
     const BL_FORT_FAB_ARG_3D(rhoe_star),    
     const amrex::Real* dx, amrex::Real dt,
     int igroup, amrex::Real ptc_tau);

  void ca_rad_source
    (const int* lo, const int* hi, 
     BL_FORT_FAB_ARG_3D(rhs),
     const amrex::Real* dx, amrex::Real dt,
     amrex::Real time, int igroup);

  void ca_local_accel
   (const int* lo, const int* hi,
    BL_FORT_FAB_ARG_3D(Ern),
    const BL_FORT_FAB_ARG_3D(Erl),
    const BL_FORT_FAB_ARG_3D(kap),
    const BL_FORT_FAB_ARG_3D(etaT),
    const BL_FORT_FAB_ARG_3D(mugT),
    amrex::Real dt, amrex::Real ptc_tau);

}

BL_FORT_PROC_DECL(CA_STATE_UPDATE, ca_state_update)
   (const int lo[], const int hi[],
    BL_FORT_FAB_ARG(state),
    const BL_FORT_FAB_ARG(rhoe),
    const BL_FORT_FAB_ARG(temp),
    const BL_FORT_FAB_ARG(msk),
    amrex::Real* derat, amrex::Real* dTrat);

extern "C"
{
  void ca_update_matter
   (const int* lo, const int* hi,
    BL_FORT_FAB_ARG_3D(rhoe_new),
    const BL_FORT_FAB_ARG_3D(Er_new),
    const BL_FORT_FAB_ARG_3D(Er_pi),
    const BL_FORT_FAB_ARG_3D(rhoe_star),
    const BL_FORT_FAB_ARG_3D(rhoe_step),
    const BL_FORT_FAB_ARG_3D(eta1),
    const BL_FORT_FAB_ARG_3D(coupT),
    const BL_FORT_FAB_ARG_3D(kappa_p),
    amrex::Real dt, amrex::Real ptc_tau);
}

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

#ifdef __cplusplus
extern "C" {
#endif
  void ca_compute_rosseland
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(kpr),
     const BL_FORT_FAB_ARG_3D(state),
     int first_group, int last_group, int num_groups);

  void ca_compute_scattering
   (const int * lo, const int* hi,
    BL_FORT_FAB_ARG_3D(kps),
    const BL_FORT_FAB_ARG_3D(state));

  void ca_compute_planck
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(kpp),
     const BL_FORT_FAB_ARG_3D(state),
     int first_group, int last_group, int num_groups,
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

extern "C"
{
  void ca_opacs
   (const int* lo, const int* hi,
    const BL_FORT_FAB_ARG_3D(Snew),
    const BL_FORT_FAB_ARG_3D(T),
    const BL_FORT_FAB_ARG_3D(Ts),
    BL_FORT_FAB_ARG_3D(kpp),
    BL_FORT_FAB_ARG_3D(kpr),
    BL_FORT_FAB_ARG_3D(dkdT),
    int use_dkdT, int validStar, int lag_opac);

  void ca_compute_emissivity
   (const int* lo, const int* hi,
    BL_FORT_FAB_ARG_3D(jg),
    BL_FORT_FAB_ARG_3D(djdT),
    const BL_FORT_FAB_ARG_3D(temp),
    const BL_FORT_FAB_ARG_3D(kappa_p),
    const BL_FORT_FAB_ARG_3D(dkdT));
}

extern "C"
{
  void ca_accel_ccoe
   (const int* lo, const int* hi, 
    const BL_FORT_FAB_ARG_3D(bcgr),
    const BL_FORT_FAB_ARG_3D(spec),
    BL_FORT_FAB_ARG_3D(ccoe),
    const amrex::Real* dx, int idim, int igroup);
}

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

  void ca_rhstoer
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(rhs),
     const amrex::Real* dx, const amrex::Real dt);

  void ca_derertot
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
  void multrs(const int* lo, const int* hi,
              BL_FORT_FAB_ARG_3D(cc),
              const amrex::Real* dx);

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

  void bclim(const int* lo, const int* hi,
             BL_FORT_FAB_ARG_3D(bcoefs), 
             BL_FORT_FAB_ARG_3D(lambda), 
             int idim,
             BL_FORT_FAB_ARG_3D(kappa_r), 
             amrex::Real c, const amrex::Real* dx);

  void flxlim(const int* lo, const int* hi,
              BL_FORT_FAB_ARG_3D(lambda), 
              int limiter);

  void scgrd(const int* lo, const int* hi,
             BL_FORT_FAB_ARG_3D(R), 
             int idim, const amrex::Real* dx,
             BL_FORT_FAB_ARG_3D(kappa_r), 
             BL_FORT_FAB_ARG_3D(Er),
             int include_cross_terms);

  void lrhs(const int* lo, const int* hi,
            BL_FORT_FAB_ARG_3D(rhs),
            BL_FORT_FAB_ARG_3D(temp),
            BL_FORT_FAB_ARG_3D(fkp),
            BL_FORT_FAB_ARG_3D(eta),
            BL_FORT_FAB_ARG_3D(etainv),
            BL_FORT_FAB_ARG_3D(rhoem),
            BL_FORT_FAB_ARG_3D(rhoes),
            BL_FORT_FAB_ARG_3D(dflux_old),
            BL_FORT_FAB_ARG_3D(Er_old),
            BL_FORT_FAB_ARG_3D(Edot), 
            amrex::Real delta_t, const amrex::Real* dx,
            amrex::Real sigma, amrex::Real c, amrex::Real theta);

  void cfrhoe(const int* lo, const int* hi, 
              BL_FORT_FAB_ARG_3D(frhoe), 
              BL_FORT_FAB_ARG_3D(state));

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

  void ceup(const int* lo, const int* hi,
            amrex::Real* relres, amrex::Real* absres,
            BL_FORT_FAB_ARG_3D(frhoes), 
            const BL_FORT_FAB_ARG_3D(frhoem),
            const BL_FORT_FAB_ARG_3D(eta),
            const BL_FORT_FAB_ARG_3D(etainv),
            const BL_FORT_FAB_ARG_3D(dflux_old),
            const BL_FORT_FAB_ARG_3D(dflux_new),
            const BL_FORT_FAB_ARG_3D(exch),
            amrex::Real delta_t, amrex::Real theta);

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

  void cetot(const int* lo, const int* hi,
             BL_FORT_FAB_ARG_3D(state), 
             BL_FORT_FAB_ARG_3D(rhoe));

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
              int nvar);

  void lacoefmgfld(const int* lo, const int* hi,
                   BL_FORT_FAB_ARG_3D(a), 
                   BL_FORT_FAB_ARG_3D(kappa),
                   const amrex::Real* dx, amrex::Real delta_t, amrex::Real c);

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

  void lbcoefna(const int* lo, const int* hi,
                BL_FORT_FAB_ARG_3D(bcoef),
                BL_FORT_FAB_ARG_3D(bcgrp),
                BL_FORT_FAB_ARG_3D(spec),
                int idim);

  void ljupna(const int* lo, const int* hi,
              BL_FORT_FAB_ARG_3D(jnew), 
              BL_FORT_FAB_ARG_3D(spec),
              BL_FORT_FAB_ARG_3D(accel), 
              int nTotal);
#ifdef __cplusplus
}
#endif

#endif
