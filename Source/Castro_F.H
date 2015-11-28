#ifndef _Castro_F_H_
#define _Castro_F_H_
#include <BLFort.H>

BL_FORT_PROC_DECL(FILCC,filcc)
                (Real * q, ARLIM_P(q_lo), ARLIM_P(q_hi),
                 const int * domlo, const int * domhi,
                 const Real * dx_crse, const Real * xlo,
                 const int * bc);

BL_FORT_PROC_DECL(CA_NETWORK_INIT,ca_network_init)();

BL_FORT_PROC_DECL(CA_EXTERN_INIT,ca_extern_init)(const int* name, const int* namlen);

BL_FORT_PROC_DECL(GET_NUM_SPEC,get_num_spec)(int* nspec);

BL_FORT_PROC_DECL(GET_NUM_AUX,get_num_aux)(int* naux);

BL_FORT_PROC_DECL(GET_SPEC_NAMES,get_spec_names)(int* spec_names, int* ispec, int* len);

BL_FORT_PROC_DECL(GET_SPEC_AZ,get_spec_az)(int* ispec, Real* A, Real* Z);

BL_FORT_PROC_DECL(GET_AUX_NAMES,get_aux_names)(int* aux_names, int* iaux, int* len);

BL_FORT_PROC_DECL(GET_METHOD_PARAMS,get_method_params)(int* HYP_GROW);

#if (BL_SPACEDIM > 1)
BL_FORT_PROC_DECL(GET_NUMPTS_1D,get_numpts_1d)(int* numpts_1d);
BL_FORT_PROC_DECL(SET_NUMPTS_1D,set_numpts_1d)(int* numpts_1d);
BL_FORT_PROC_DECL(ALLOCATE_OUTFLOW_DATA,allocate_outflow_data)(const int* n1d, const int* nc);
BL_FORT_PROC_DECL(SET_OLD_OUTFLOW_DATA,set_old_outflow_data)(Real* radial, Real* time, 
                                                             int* numpts, int* nc);
BL_FORT_PROC_DECL(SET_NEW_OUTFLOW_DATA,set_new_outflow_data)(Real* radial, Real* time, 
                                                             int* numpts, int* nc);
BL_FORT_PROC_DECL(SWAP_OUTFLOW_DATA,swap_outflow_data)();
#endif

BL_FORT_PROC_DECL(SET_METHOD_PARAMS,set_method_params)
    (const int& dm,
     const int& Density, const int& Xmom, 
     const int& Eden,    const int& Eint, 
     const int& Temp     , const int& FirstAdv, 
     const int& FirstSpec, const int& FirstAux, 
     const int& NumAdv, 
     const int* gravity_type, const int* gravity_type_len,
     const int& get_g_from_phi,
     const int& use_sgs,
     const Real& diffuse_cutoff_density, 
     const Real& const_grav);

#include <castro_set_meth.H>

BL_FORT_PROC_DECL(SET_PROBLEM_PARAMS,set_problem_params)
    (const int& dm,
     const int* physbc_lo, const int* physbc_hi,
     const int& Outflow_value, const int& Symmetry_value, 
     const int& SlipWall_value, const int& NoSlipWall_value,
     const int& coord_type,
     const Real* problo, const Real* probhi,
     const Real* center);

BL_FORT_PROC_DECL(SET_REFINEMENT_PARAMS,set_refinement_params)
    (const int& max_level, const Real* dx_level);

BL_FORT_PROC_DECL(GET_TAGGING_PARAMS,get_tagging_params)(const int* name, const int* namlen);
BL_FORT_PROC_DECL(GET_SPONGE_PARAMS,get_sponge_params)(const int* name, const int* namlen);
BL_FORT_PROC_DECL(UPDATE_SPONGE_PARAMS,update_sponge_params)(const Real* time);

#ifdef DIMENSION_AGNOSTIC

BL_FORT_PROC_DECL(CA_INITDATA,ca_initdata)
    (const int& level, const Real& time, 
     const int* lo, const int* hi,
     const int& num_state,
     BL_FORT_FAB_ARG_3D(state),
     const Real* dx, const Real* xlo, const Real* xhi);

#else

BL_FORT_PROC_DECL(CA_INITDATA,ca_initdata)
    (const int& level, const Real& time, 
     const int* lo, const int* hi,
     const int& num_state,
     BL_FORT_FAB_ARG(state),
     const Real* dx, const Real* xlo, const Real* xhi);

#endif

BL_FORT_PROC_DECL(CA_CHECK_INITIAL_SPECIES, ca_check_initial_species)
    (const int* lo, const int* hi, BL_FORT_FAB_ARG(state));

BL_FORT_PROC_DECL(CA_GET_OMB,ca_get_omb)(Real* frac);

BL_FORT_PROC_DECL(CA_ENFORCE_CONSISTENT_E,ca_enforce_consistent_e)
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state));

#ifdef MAESTRO_INIT
BL_FORT_PROC_DECL(CA_INITDATA_MAESTRO,ca_initdata_maestro)
  (const int* lo, const int* hi, const int& MAESTRO_init_type,
     BL_FORT_FAB_ARG(state),
     const Real dx[], const Real& dr,
     const Real xlo[], const Real xhi[],
     const Real p0[], const int& MAESTRO_npts_model, const int& level);

BL_FORT_PROC_DECL(CA_INITDATA_MAKEMODEL,ca_initdata_makemodel)
  (Real model[], const int& model_size, const int& MAESTRO_npts_model,
   const Real rho0[], const Real tempbar[], 
   const Real dx[], const Real& dr, const int& r_model_start);

BL_FORT_PROC_DECL(CA_INITDATA_OVERWRITE,ca_initdata_overwrite)
  (const int* lo, const int* hi, 
   BL_FORT_FAB_ARG(state),
   const Real model[], const int& model_size, const int& MAESTRO_npts_model, 
   const Real dx[], const Real& dr, const Real xlo[], const Real xhi[],
   const int& r_model_start);
#endif

#ifdef LEVELSET
BL_FORT_PROC_DECL(CA_INITPHI,ca_initphi)
    (const int& level, const Real& time, 
     const int* lo, const int* hi,
     const int& num_state,
     BL_FORT_FAB_ARG(state),
     int* type, ARLIM_P(type_lo), ARLIM_P(type_hi),
     const Real dx[], const Real xlo[], const Real xhi[]);
#endif

#ifdef RADIATION
BL_FORT_PROC_DECL(CA_INITRAD,ca_initrad)
  (const int& level, const Real& time,
   const int* lo, const int* hi,
   const int& num_state,
   BL_FORT_FAB_ARG(state),
   const Real dx[], const Real xlo[], const Real xhi[]);
#endif

#ifdef DIMENSION_AGNOSTIC

BL_FORT_PROC_DECL(CA_HYPFILL,ca_hypfill)
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo, 
     const Real* time, const int* bc);

BL_FORT_PROC_DECL(CA_DENFILL,ca_denfill)
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo, 
     const Real* time, const int* bc);

BL_FORT_PROC_DECL(CA_XMOMFILL,ca_xmomfill)
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo, 
     const Real* time, const int* bc);

BL_FORT_PROC_DECL(CA_YMOMFILL,ca_ymomfill)
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo, 
     const Real* time, const int* bc);

BL_FORT_PROC_DECL(CA_ZMOMFILL,ca_zmomfill)
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo, 
     const Real* time, const int* bc);

BL_FORT_PROC_DECL(CA_ZEROFILL,ca_zerofill)
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo, 
     const Real* time, const int* bc);

BL_FORT_PROC_DECL(CA_DEFAULT_FILL,ca_default_fill)
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo,
     const Real* time, const int* bc);

#ifdef RADIATION
BL_FORT_PROC_DECL(CA_RADFILL,ca_radfill)
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo, 
     const Real* time, const int* bc);
#endif

#ifdef LEVELSET
BL_FORT_PROC_DECL(CA_PHIFILL,ca_phifill)
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo, 
     const Real* time, const int* bc);
#endif

#ifdef REACTIONS
BL_FORT_PROC_DECL(CA_REACTFILL,ca_reactfill)
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo, 
     const Real* time, const int* bc);
#endif

#ifdef SGS
BL_FORT_PROC_DECL(CA_SGSFILL,ca_sgsfill)
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo, 
     const Real* time, const int* bc);
#endif

#if defined(GRAVITY) || defined(ROTATION)
BL_FORT_PROC_DECL(CA_PHIGRAVFILL,ca_phigravfill)
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo, 
     const Real* time, const int* bc);

BL_FORT_PROC_DECL(CA_GRAVXFILL,ca_gravxfill)
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo, 
     const Real* time, const int* bc);

BL_FORT_PROC_DECL(CA_GRAVYFILL,ca_gravyfill)
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo, 
     const Real* time, const int* bc);

BL_FORT_PROC_DECL(CA_GRAVZFILL,ca_gravzfill)
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo, 
     const Real* time, const int* bc);
#endif

#else

BL_FORT_PROC_DECL(CA_HYPFILL,ca_hypfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(CA_DENFILL,ca_denfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(CA_XMOMFILL,ca_xmomfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(CA_YMOMFILL,ca_ymomfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(CA_ZMOMFILL,ca_zmomfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(CA_ZEROFILL,ca_zerofill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(CA_DEFAULT_FILL,ca_default_fill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);

#ifdef RADIATION
BL_FORT_PROC_DECL(CA_RADFILL,ca_radfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);
#endif

#ifdef LEVELSET
BL_FORT_PROC_DECL(CA_PHIFILL,ca_phifill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);
#endif

#ifdef REACTIONS
BL_FORT_PROC_DECL(CA_REACTFILL,ca_reactfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);
#endif

#ifdef SGS
BL_FORT_PROC_DECL(CA_SGSFILL,ca_sgsfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);
#endif

#if defined(GRAVITY) || defined(ROTATION)
BL_FORT_PROC_DECL(CA_PHIGRAVFILL,ca_phigravfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(CA_GRAVXFILL,ca_gravxfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(CA_GRAVYFILL,ca_gravyfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);

BL_FORT_PROC_DECL(CA_GRAVZFILL,ca_gravzfill)
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[], 
     const Real* time, const int bc[]);
#endif

#endif

BL_FORT_PROC_DECL(CA_ESTDT,ca_estdt)
    (const int* lo, const int* hi, 
     const BL_FORT_FAB_ARG_3D(state),
     const Real* dx, Real* dt);

#ifdef DIFFUSION
BL_FORT_PROC_DECL(CA_ESTDT_DIFFUSION,ca_estdt_diffusion)
    (const int* lo, const int* hi, 
     const BL_FORT_FAB_ARG_3D(state),
     const Real* dx, Real* dt);
#endif

#ifdef RADIATION
BL_FORT_PROC_DECL(CA_ESTDT_RAD,ca_estdt_rad)
    (const BL_FORT_FAB_ARG(state),
     const BL_FORT_FAB_ARG(gpr),
     const int lo[], const int hi[],
     const Real dx[], Real* dt);
#endif

#ifdef REACTIONS
BL_FORT_PROC_DECL(CA_ESTDT_BURNING,ca_estdt_burning)
    (const BL_FORT_FAB_ARG_3D(state),
     const BL_FORT_FAB_ARG_3D(reactions),
     const int* lo, const int* hi,
     const Real* dx, Real* dt);
#endif

BL_FORT_PROC_DECL(CA_UMDRV,ca_umdrv)
    (const int* is_finest_level,
     const Real* time, 
     const int    lo[], const int    hi[],
     const int domlo[], const int domhi[],
     const BL_FORT_FAB_ARG(state),
     BL_FORT_FAB_ARG(stateout),
     BL_FORT_FAB_ARG(ugdnvx),
#if (BL_SPACEDIM >= 2)
     BL_FORT_FAB_ARG(ugdnvy),
#endif
#if (BL_SPACEDIM == 3)
     BL_FORT_FAB_ARG(ugdnvz),
#endif
     const BL_FORT_FAB_ARG(src),
     const Real dx[], const Real* dt,
     D_DECL(BL_FORT_FAB_ARG(xflux),
            BL_FORT_FAB_ARG(yflux),
            BL_FORT_FAB_ARG(zflux)),
     D_DECL(BL_FORT_FAB_ARG(xarea),
            BL_FORT_FAB_ARG(yarea),
            BL_FORT_FAB_ARG(zarea)),
#if (BL_SPACEDIM < 3)
     const BL_FORT_FAB_ARG(dloga),
#endif
     const BL_FORT_FAB_ARG(volume),
     Real* cflLoc,
     const int&  verbose,
     Real& mass_added,
     Real& eint_added,
     Real& eden_added,
     Real& xmom_added_flux,
     Real& ymom_added_flux,
     Real& zmom_added_flux,
     Real& E_added_flux);

#ifdef RADIATION
BL_FORT_PROC_DECL(CA_UMDRV_RAD,ca_umdrv_rad)
    (const int* is_finest_level,
     const Real* time, 
     const int    lo[], const int    hi[],
     const int domlo[], const int domhi[],
     const BL_FORT_FAB_ARG(state),
     BL_FORT_FAB_ARG(stateout),
     const BL_FORT_FAB_ARG(Er),
     const BL_FORT_FAB_ARG(lam),
     BL_FORT_FAB_ARG(Erout),
     BL_FORT_FAB_ARG(ugdnvx),
#if (BL_SPACEDIM >= 2)
     BL_FORT_FAB_ARG(ugdnvy),
#endif
#if (BL_SPACEDIM == 3)
     BL_FORT_FAB_ARG(ugdnvz),
#endif
     const BL_FORT_FAB_ARG(src),
     const Real dx[], const Real* dt,
     D_DECL(BL_FORT_FAB_ARG(xflux),
            BL_FORT_FAB_ARG(yflux),
            BL_FORT_FAB_ARG(zflux)),
     D_DECL(BL_FORT_FAB_ARG(xradflux),
            BL_FORT_FAB_ARG(yradflux),
            BL_FORT_FAB_ARG(zradflux)),
     D_DECL(const BL_FORT_FAB_ARG(xarea),
            const BL_FORT_FAB_ARG(yarea),
            const BL_FORT_FAB_ARG(zarea)),
#if (BL_SPACEDIM < 3)
     const BL_FORT_FAB_ARG(dloga),
#endif
     const BL_FORT_FAB_ARG(volume),
     Real* cflLoc,
     const int& verbose, 
     int* nstep_fsp); 

BL_FORT_PROC_DECL(CA_INELASTIC_SCT, ca_inelastic_sct)
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(stateout),
     BL_FORT_FAB_ARG_3D(Erout),
     const BL_FORT_FAB_ARG_3D(kappa_s),
     const Real& dt);
#endif

BL_FORT_PROC_DECL(CA_SPONGE,ca_sponge)
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const Real* dx, const Real& dt, const Real* time, 
     const Real& E_added, const Real* mom_added);

#ifdef GRAVITY

BL_FORT_PROC_DECL(CA_GSRC,ca_gsrc)
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const BL_FORT_FAB_ARG_3D(phi),
     const BL_FORT_FAB_ARG_3D(grav),
     const BL_FORT_FAB_ARG_3D(s_old),
     BL_FORT_FAB_ARG_3D(s_new),
     const Real* dx, const Real& dt, const Real* time,
     const Real& E_added, const Real* mom_added);

BL_FORT_PROC_DECL(CA_CORRGSRC,ca_corrgsrc)
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const BL_FORT_FAB_ARG_3D(phi_old),
     const BL_FORT_FAB_ARG_3D(phi_new),
     const BL_FORT_FAB_ARG_3D(grav_old),
     const BL_FORT_FAB_ARG_3D(grav_new),
     const BL_FORT_FAB_ARG_3D(S_old),
     BL_FORT_FAB_ARG_3D(S_new),
     BL_FORT_FAB_ARG_3D(xflux),
     BL_FORT_FAB_ARG_3D(yflux),
     BL_FORT_FAB_ARG_3D(zflux),
     const Real* dx, const Real& dt, const Real* time,
     const BL_FORT_FAB_ARG_3D(volume), 
     const Real& E_added, const Real* mom_added);

#endif

#ifdef ROTATION

BL_FORT_PROC_DECL(CA_FILL_ROTATIONAL_POTENTIAL,ca_fill_rotational_potential)
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(phi),
     const Real* dx, const Real& time);

BL_FORT_PROC_DECL(CA_FILL_ROTATIONAL_ACCELERATION,ca_fill_rotational_acceleration)
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(rot),
     const BL_FORT_FAB_ARG_3D(state),
     const Real* dx, const Real& time);

BL_FORT_PROC_DECL(CA_RSRC,ca_rsrc)
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const BL_FORT_FAB_ARG_3D(phi),
     const BL_FORT_FAB_ARG_3D(rot),
     const BL_FORT_FAB_ARG_3D(s_old),
     BL_FORT_FAB_ARG_3D(s_new),
     const Real* dx, const Real& dt, const Real* time, 
     const Real& E_added, const Real* mom_added);

BL_FORT_PROC_DECL(CA_CORRRSRC,ca_corrrsrc)
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const BL_FORT_FAB_ARG_3D(phi_old),
     const BL_FORT_FAB_ARG_3D(phi_new),
     const BL_FORT_FAB_ARG_3D(rot_old),
     const BL_FORT_FAB_ARG_3D(rot_new),
     const BL_FORT_FAB_ARG_3D(S_old),
     BL_FORT_FAB_ARG_3D(S_new),
     BL_FORT_FAB_ARG_3D(xflux),
     BL_FORT_FAB_ARG_3D(yflux),
     BL_FORT_FAB_ARG_3D(zflux),
     const Real* dx, const Real& dt, const Real* time, 
     const BL_FORT_FAB_ARG_3D(volume), 
     const Real& E_added, const Real* mom_added);

#endif

BL_FORT_PROC_DECL(RESET_INTERNAL_E,reset_internal_e)
    (const int* lo, const int* hi, 
     BL_FORT_FAB_ARG_3D(S_new), 
     const int& verbose);

#ifdef SGS
BL_FORT_PROC_DECL(CA_RESET_SGS,ca_reset_sgs)
    (BL_FORT_FAB_ARG(S_new), 
     BL_FORT_FAB_ARG(sgs_old), 
     BL_FORT_FAB_ARG(sgs_new), 
     const int lo[], const int hi[], 
     const int& verbose, const int& is_old, const Real& dt);
#endif

BL_FORT_PROC_DECL(CA_SYNCGSRC,ca_syncgsrc)
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(gPhi),
     const BL_FORT_FAB_ARG_3D(gdPhi),
     const BL_FORT_FAB_ARG_3D(S),
     const BL_FORT_FAB_ARG_3D(dS),
     BL_FORT_FAB_ARG_3D(src),
     const Real& dt);

BL_FORT_PROC_DECL(CA_SUMMASS,ca_summass)
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(rho),
     const Real* dx, BL_FORT_FAB_ARG_3D(vol), Real* mass);
BL_FORT_PROC_DECL(CA_SUMSQUARED,ca_sumsquared)
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(rho),
     const Real* dx, BL_FORT_FAB_ARG_3D(vol), Real* mass);
BL_FORT_PROC_DECL(CA_SUMLOCMASS,ca_sumlocmass)
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(rho),
     const Real* dx, BL_FORT_FAB_ARG_3D(vol), Real* mass, const int& idir);
BL_FORT_PROC_DECL(CA_SUMLOCMASS2D,ca_sumlocmass2d)
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(rho),
     const Real* dx, BL_FORT_FAB_ARG_3D(vol), Real* mass, const int& idir1, const int& idir2);
BL_FORT_PROC_DECL(CA_SUMLOCSQUAREDMASS,ca_sumlocsquaredmass)
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(rho),
     const Real* dx, BL_FORT_FAB_ARG_3D(vol), Real* mass, const int& idir); 
BL_FORT_PROC_DECL(CA_SUMPRODUCT,ca_sumproduct)
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(f1), BL_FORT_FAB_ARG_3D(f2),
     const Real* dx, BL_FORT_FAB_ARG_3D(vol), Real* s);

BL_FORT_PROC_DECL(CA_ENFORCE_NONNEGATIVE_SPECIES,ca_enforce_nonnegative_species)
    (BL_FORT_FAB_ARG(S_new), const int lo[], const int hi[]);

#ifdef SGS
BL_FORT_PROC_DECL(CA_EXT_SRC,ca_ext_src)
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(old_state),
     BL_FORT_FAB_ARG(fluxx),
     BL_FORT_FAB_ARG(fluxy),
     BL_FORT_FAB_ARG(fluxz),
     BL_FORT_FAB_ARG(ext_src),
     BL_FORT_FAB_ARG(prod_sgs),
     BL_FORT_FAB_ARG(diff_sgs),
     BL_FORT_FAB_ARG(turb_src),
     const Real* dx, const Real* time, const Real* dt);
#else
#ifdef DIMENSION_AGNOSTIC
BL_FORT_PROC_DECL(CA_EXT_SRC,ca_ext_src)
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(old_state),
     const BL_FORT_FAB_ARG_3D(new_state),
     BL_FORT_FAB_ARG_3D(ext_src),
     const Real* prob_lo, const Real* dx, 
     const Real* time, const Real* dt);
#else
BL_FORT_PROC_DECL(CA_EXT_SRC,ca_ext_src)
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG(old_state),
     const BL_FORT_FAB_ARG(new_state),
     BL_FORT_FAB_ARG(ext_src),
     const Real* prob_lo, const Real* dx, 
     const Real* time, const Real* dt);
#endif
#endif

BL_FORT_PROC_DECL(CA_REACT_STATE,ca_react_state)
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(reactions),
#ifdef TAU
     const BL_FORT_FAB_ARG_3D(Tau),
#endif
     const Real& time, const Real& dt_react);

BL_FORT_PROC_DECL(CA_TEMPDIFFEXTRAP,ca_tempdiffextrap)
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(TempDiffTerm));

BL_FORT_PROC_DECL(COMPUTE_TEMP,compute_temp)
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state));

BL_FORT_PROC_DECL(CA_FILL_TEMP_COND,ca_fill_temp_cond)
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
#ifdef TAU
     const BL_FORT_FAB_ARG_3D(tau_diff),
#endif
     const BL_FORT_FAB_ARG_3D(xcoeffs),
     const BL_FORT_FAB_ARG_3D(ycoeffs),
     const BL_FORT_FAB_ARG_3D(zcoeffs),
     const Real* dx);

BL_FORT_PROC_DECL(CA_FILL_SPEC_COEFF,ca_fill_spec_coeff)
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const BL_FORT_FAB_ARG_3D(xcoeffs),
     const BL_FORT_FAB_ARG_3D(ycoeffs),
     const BL_FORT_FAB_ARG_3D(zcoeffs),
     const Real* dx);

BL_FORT_PROC_DECL(CA_FILL_FIRST_VISC_COEFF,ca_fill_first_visc_coeff)
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const BL_FORT_FAB_ARG_3D(xcoeffs),
     const BL_FORT_FAB_ARG_3D(ycoeffs),
     const BL_FORT_FAB_ARG_3D(zcoeffs),
     const Real* dx);

BL_FORT_PROC_DECL(CA_FILL_SECND_VISC_COEFF,ca_fill_secnd_visc_coeff)
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const BL_FORT_FAB_ARG_3D(xcoeffs),
     const BL_FORT_FAB_ARG_3D(ycoeffs),
     const BL_FORT_FAB_ARG_3D(zcoeffs),
     const Real* dx);

BL_FORT_PROC_DECL(CA_COMPUTE_DIV_TAU_U,ca_compute_div_tau_u)
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(div_tau_u),
     const BL_FORT_FAB_ARG_3D(state),
     const Real* dx, const int* coord_type);

#ifdef TAU
BL_FORT_PROC_DECL(CA_DEFINE_TAU,ca_define_tau)
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(tau_diff),
     const BL_FORT_FAB_ARG(state),
     const BL_FORT_FAB_ARG(grav),
     const Real dx[]);
#endif

BL_FORT_PROC_DECL(CA_SET_SPECIAL_TAGGING_FLAG,ca_set_special_tagging_flag)
    (const Real& max_den, const int* flag);

#ifdef POINTMASS
BL_FORT_PROC_DECL(PM_COMPUTE_DELTA_MASS,pm_compute_delta_mass)
    (const Real* pm_delta_mass, const int* lo, const int* hi,
     const BL_FORT_FAB_ARG(state),
     BL_FORT_FAB_ARG(stateout),
     const BL_FORT_FAB_ARG(volume),
     const Real* problo, const Real dx[], 
     const Real* time, const Real* dt);

BL_FORT_PROC_DECL(PM_FIX_SOLUTION,pm_fix_solution)
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(stateold),
     BL_FORT_FAB_ARG(statenew),
     const Real* problo, const Real dx[], 
     const Real* time, const Real* dt);
#endif

#ifdef AUX_UPDATE
BL_FORT_PROC_DECL(CA_AUXUPDATE,ca_auxupdate)
    (BL_FORT_FAB_ARG(state_old),
     BL_FORT_FAB_ARG(state_new),
     const int* lo, const int* hi,
     const Real * dt);
#endif

#if (BL_SPACEDIM > 1)
BL_FORT_PROC_DECL(CA_COMPUTE_AVGSTATE,ca_compute_avgstate)
    (const int lo[], const int hi[], 
     const Real* dx, const Real* dr, const int* nc,
     const BL_FORT_FAB_ARG(S  ), const Real* avgden,
     const BL_FORT_FAB_ARG(Vol), const Real* avgvol,
     const Real* problo, const int* numpts_1d);
#endif

BL_FORT_PROC_DECL( GET_CENTER, get_center)(Real* center);
BL_FORT_PROC_DECL( SET_CENTER, set_center)(Real* center);
BL_FORT_PROC_DECL(FIND_CENTER,find_center)(Real* data, Real* center, int* icen,
                                           const Real* dx, const Real* problo);

#ifdef NEUTRINO
BL_FORT_PROC_DECL(CA_SETGROUP,ca_setgroup)
    (const int& igroup);
#endif

#ifdef ROTATION
BL_FORT_PROC_DECL(CA_ROTATE,ca_rotate)
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(state),
     BL_FORT_FAB_ARG(rotation_terms),
     const Real* prob_lo, const Real* dx);
#endif

/* problem-specific includes */
#include <Problem_F.H>

#endif
