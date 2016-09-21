#ifndef _Castro_F_H_
#define _Castro_F_H_
#include <BLFort.H>

#include <BoxLib.H>

extern "C"
{

  void ca_network_init();

  void burner_init();

  void ca_extern_init(const int* name, const int* namlen);

  void get_num_spec(int* nspec);

  void get_num_aux(int* naux);

  void get_spec_names(int* spec_names, int* ispec, int* len);

  void get_spec_az(int* ispec, Real* A, Real* Z);

  void get_aux_names(int* aux_names, int* iaux, int* len);

  void get_qvar(int* qvar);

  void get_nqaux(int* nqaux);

  void ca_set_special_tagging_flag(const Real& max_den, const int* flag);

  void set_amr_info(const int& level, const int& iteration, const int& ncycle,
		    const Real& time, const Real& dt);

  void get_method_params(int* HYP_GROW);

#if (BL_SPACEDIM > 1)
  void get_numpts_1d(int* numpts_1d);
  void set_numpts_1d(int* numpts_1d);
  void allocate_outflow_data(const int* n1d, const int* nc);
  void set_old_outflow_data(Real* radial, Real* time,
			    int* numpts, int* nc);
  void set_new_outflow_data(Real* radial, Real* time,
			    int* numpts, int* nc);
  void swap_outflow_data();
#endif

  void set_method_params
    (const int& dm,
     const int& Density, const int& Xmom, 
     const int& Eden,    const int& Eint, 
     const int& Temp     , const int& FirstAdv, 
     const int& FirstSpec, const int& FirstAux, 
     const int& NumAdv,
#ifdef SHOCK_VAR
     const int& Shock,
#endif
     const int* gravity_type, const int& gravity_type_len,
     const Real& diffuse_cutoff_density);

  void set_castro_method_params();

#ifdef RADIATION
  void init_godunov_indices_rad();
#else
  void init_godunov_indices();
#endif

  void set_problem_params
    (const int& dm,
     const int* physbc_lo, const int* physbc_hi,
     const int& Interior_value, const int& Inflow_value,
     const int& Outflow_value, const int& Symmetry_value,
     const int& SlipWall_value, const int& NoSlipWall_value,
     const int& coord_type,
     const Real* problo, const Real* probhi,
     const Real* center);

  void set_grid_info
    (const int& max_level, const Real* dx_level,
     const int* domlo_level, const int* domhi_level,
     const int* ref_ratio, const int* n_error_buf,
     const int* blocking_factor);

  void get_tagging_params(const int* name, const int* namlen);
  void get_sponge_params(const int* name, const int* namlen);

#if (BL_SPACEDIM > 1)
  void ca_compute_avgstate
    (const int* lo, const int* hi, 
     const Real* dx, const Real* dr, const int* nc,
     const BL_FORT_FAB_ARG_3D(S  ), const Real* avgden,
     const BL_FORT_FAB_ARG_3D(Vol), const Real* avgvol,
     const Real* problo, const int* numpts_1d);
#endif

#ifdef HYBRID_MOMENTUM
  void init_hybrid_momentum
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state));

  void ca_hybrid_hydro_source
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(sources));

  void hybrid_update
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state));
#endif

  void ctoprim(const int* lo, const int* hi,
	       const Real* u, const int* u_lo, const int* u_hi,
	       const Real* q, const int* q_lo, const int* q_hi,
               const Real* qaux, const int* qa_lo, const int* qa_hi);

  void srctoprim(const int* lo,   const int* hi,
		 const Real* q,   const int* q_lo,   const int* q_hi,
		 const Real* qaux, const int* qa_lo, const int* qa_hi,
		 const Real* src, const int* src_lo, const int* src_hi,
		 const Real* srQ, const int* srQ_lo, const int* srQ_hi);

  void ca_check_initial_species
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state));

  void enforce_minimum_density
    (const Real* S_old, const int* s_old_lo, const int* s_old_hi,
     const Real* S_new, const int* s_new_lo, const int* s_new_hi,
     const Real* vol, const int* vol_lo, const int* vol_hi,
     const int* lo, const int* hi,
     const Real* mass_added, const Real* e_added,
     const Real* E_added, const Real* frac_change,
     const int* verbose);

  void ca_normalize_species
    (BL_FORT_FAB_ARG_3D(S_new), const int* lo, const int* hi);

  void get_center(Real* center);
  void set_center(Real* center);
  void find_center(Real* data, Real* center, int* icen,
		   const Real* dx, const Real* problo);

  void compute_temp
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state));

  void ca_enforce_consistent_e
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state));

  void reset_internal_e
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(S_new),
     const int& verbose);

#ifdef DIMENSION_AGNOSTIC

  void ca_hypfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo, 
     const Real* time, const int* bc);

  void ca_denfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo, 
     const Real* time, const int* bc);

  void ca_xmomfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo,
     const Real* time, const int* bc);

  void ca_ymomfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo,
     const Real* time, const int* bc);

  void ca_zmomfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo,
     const Real* time, const int* bc);

  void ca_zerofill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo,
     const Real* time, const int* bc);

  void ca_default_fill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo,
     const Real* time, const int* bc);

#ifdef RADIATION
  void ca_radfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo,
     const Real* time, const int* bc);
#endif

#ifdef LEVELSET
  void ca_phifill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo,
     const Real* time, const int* bc);
#endif

#ifdef REACTIONS
  void ca_reactfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo, 
     const Real* time, const int* bc);
#endif

#ifdef GRAVITY
  void ca_phigravfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo, 
     const Real* time, const int* bc);

  void ca_gravxfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo,
     const Real* time, const int* bc);

  void ca_gravyfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo,
     const Real* time, const int* bc);

  void ca_gravzfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo,
     const Real* time, const int* bc);
#endif

#ifdef ROTATION
  void ca_phirotfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo,
     const Real* time, const int* bc);

  void ca_rotxfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo,
     const Real* time, const int* bc);

  void ca_rotyfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo,
     const Real* time, const int* bc);

  void ca_rotzfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo,
     const Real* time, const int* bc);
#endif

  inline void ca_nullfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo, 
     const Real* time, const int* bc)
  { BoxLib::Abort("Not supposed to call ca_nullfill"); }

#else

  void ca_hypfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);

  void ca_denfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);

  void ca_xmomfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);

  void ca_ymomfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);

  void ca_zmomfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);

  void ca_zerofill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);

  void ca_default_fill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);

#ifdef RADIATION
  void ca_radfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);
#endif

#ifdef LEVELSET
  void ca_phifill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);
#endif

#ifdef REACTIONS
  void ca_reactfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);
#endif

#ifdef GRAVITY
  void ca_phigravfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);

  void ca_gravxfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);

  void ca_gravyfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);

  void ca_gravzfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);
#endif

#ifdef ROTATION
  void ca_phirotfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);

  void ca_rotxfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);

  void ca_rotyfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);

  void ca_rotzfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[]);
#endif

  inline void ca_nullfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const Real dx[], const Real glo[],
     const Real* time, const int bc[])
  { BoxLib::Abort("Not supposed to call ca_nullfill"); }

#endif

  void ca_estdt
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const Real* dx, Real* dt);

#ifdef DIFFUSION
  void ca_estdt_temp_diffusion
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const Real* dx, Real* dt);
  void ca_estdt_enth_diffusion
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const Real* dx, Real* dt);
#endif

#ifdef RADIATION
  void ca_estdt_rad
    (const BL_FORT_FAB_ARG(state),
     const BL_FORT_FAB_ARG(gpr),
     const int lo[], const int hi[],
     const Real dx[], Real* dt);
#endif

#ifdef REACTIONS
    void ca_estdt_burning
    (const BL_FORT_FAB_ARG_3D(state_old),
     const BL_FORT_FAB_ARG_3D(state_new),
     const BL_FORT_FAB_ARG_3D(reactions_old),
     const BL_FORT_FAB_ARG_3D(reactions_new),
     const int* lo, const int* hi,
     const Real* dx, const Real* dt_old, Real* dt);
#endif

  void ca_check_timestep
    (const BL_FORT_FAB_ARG_3D(state_old),
     const BL_FORT_FAB_ARG_3D(state_new),
#ifdef REACTIONS
     const BL_FORT_FAB_ARG_3D(reactions_old),
     const BL_FORT_FAB_ARG_3D(reactions_new),
#endif
     const int* lo, const int* hi,
     const Real* dx, const Real* dt_old, Real* dt);

  void ca_umdrv
    (const int* is_finest_level,
     const Real* time,
     const int    lo[], const int    hi[],
     const int domlo[], const int domhi[],
     const BL_FORT_FAB_ARG(state),
     BL_FORT_FAB_ARG(stateout),
     BL_FORT_FAB_ARG(q),
     BL_FORT_FAB_ARG(qaux),
     BL_FORT_FAB_ARG(src_q),
     BL_FORT_FAB_ARG(update),
     BL_FORT_FAB_ARG(ugdnvx),
#if (BL_SPACEDIM >= 2)
     BL_FORT_FAB_ARG(ugdnvy),
#endif
#if (BL_SPACEDIM == 3)
     BL_FORT_FAB_ARG(ugdnvz),
#endif
     const Real dx[], const Real* dt,
     D_DECL(BL_FORT_FAB_ARG(xflux),
            BL_FORT_FAB_ARG(yflux),
            BL_FORT_FAB_ARG(zflux)),
#if (BL_SPACEDIM < 3)
     BL_FORT_FAB_ARG(pradial),
#endif
     D_DECL(BL_FORT_FAB_ARG(xarea),
            BL_FORT_FAB_ARG(yarea),
            BL_FORT_FAB_ARG(zarea)),
#if (BL_SPACEDIM < 3)
     const BL_FORT_FAB_ARG(dloga),
#endif
     const BL_FORT_FAB_ARG(volume),
     Real* cflLoc,
     const int&  verbose,
     Real& mass_added_flux,
     Real& xmom_added_flux,
     Real& ymom_added_flux,
     Real& zmom_added_flux,
     Real& E_added_flux,
     Real& mass_lost,
     Real& xmom_lost,
     Real& ymom_lost,
     Real& zmom_lost,
     Real& eden_lost,
     Real& xang_lost,
     Real& yang_lost,
     Real& zang_lost);

#ifdef RADIATION
  void ca_umdrv_rad
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
#if (BL_SPACEDIM < 3)
     BL_FORT_FAB_ARG(pradial),
#endif
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

  void ca_inelastic_sct
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(stateout),
     BL_FORT_FAB_ARG_3D(Erout),
     const BL_FORT_FAB_ARG_3D(kappa_s),
     const Real& dt);
#endif

  void ca_sponge
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(source),
     BL_FORT_FAB_ARG_3D(vol),
     const Real* dx, const Real& dt, const Real* time,
     const Real& E_added, const Real* mom_added);

  void update_sponge_params(const Real* time);

  void flush_output();

#ifdef POINTMASS
  void pm_compute_delta_mass
    (const Real* pm_delta_mass, const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(stateout),
     const BL_FORT_FAB_ARG_3D(volume),
     const Real* problo, const Real* dx,
     const Real* time, const Real* dt);

  void pm_fix_solution
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(stateold),
     BL_FORT_FAB_ARG_3D(statenew),
     const Real* problo, const Real* dx,
     const Real* time, const Real* dt);

  void set_pointmass(const Real* pointmass);
#endif

#ifdef GRAVITY

  void ca_gsrc
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const BL_FORT_FAB_ARG_3D(phi),
     const BL_FORT_FAB_ARG_3D(grav),
     const BL_FORT_FAB_ARG_3D(s_old),
     BL_FORT_FAB_ARG_3D(source),
     const BL_FORT_FAB_ARG_3D(vol),
     const Real* dx, const Real& dt, const Real* time,
     const Real& E_added, const Real* mom_added);

  void ca_corrgsrc
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const BL_FORT_FAB_ARG_3D(phi_old),
     const BL_FORT_FAB_ARG_3D(phi_new),
     const BL_FORT_FAB_ARG_3D(grav_old),
     const BL_FORT_FAB_ARG_3D(grav_new),
     const BL_FORT_FAB_ARG_3D(S_old),
     const BL_FORT_FAB_ARG_3D(S_new),
     BL_FORT_FAB_ARG_3D(source),
     BL_FORT_FAB_ARG_3D(xflux),
     BL_FORT_FAB_ARG_3D(yflux),
     BL_FORT_FAB_ARG_3D(zflux),
     const Real* dx, const Real& dt, const Real* time,
     const BL_FORT_FAB_ARG_3D(volume), 
     const Real& E_added, const Real* mom_added);

  void ca_syncgsrc
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(gPhi),
     const BL_FORT_FAB_ARG_3D(gdPhi),
     const BL_FORT_FAB_ARG_3D(S),
     const BL_FORT_FAB_ARG_3D(dS),
     BL_FORT_FAB_ARG_3D(src), const int* ncomp,
     const Real& dt);

#endif

#ifdef ROTATION

  void ca_fill_rotational_potential
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(phi),
     const Real* dx, const Real& time);

  void ca_fill_rotational_acceleration
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(rot),
     const BL_FORT_FAB_ARG_3D(state),
     const Real* dx, const Real& time);

  void ca_rsrc
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const BL_FORT_FAB_ARG_3D(phi),
     const BL_FORT_FAB_ARG_3D(rot),
     const BL_FORT_FAB_ARG_3D(s_old),
     BL_FORT_FAB_ARG_3D(source),
     const BL_FORT_FAB_ARG_3D(vol),
     const Real* dx, const Real& dt, const Real* time,
     const Real& E_added, const Real* mom_added);

  void ca_corrrsrc
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const BL_FORT_FAB_ARG_3D(phi_old),
     const BL_FORT_FAB_ARG_3D(phi_new),
     const BL_FORT_FAB_ARG_3D(rot_old),
     const BL_FORT_FAB_ARG_3D(rot_new),
     const BL_FORT_FAB_ARG_3D(S_old),
     const BL_FORT_FAB_ARG_3D(S_new),
     BL_FORT_FAB_ARG_3D(source),
     const BL_FORT_FAB_ARG_3D(xflux),
     const BL_FORT_FAB_ARG_3D(yflux),
     const BL_FORT_FAB_ARG_3D(zflux),
     const Real* dx, const Real& dt, const Real* time,
     const BL_FORT_FAB_ARG_3D(volume),
     const Real& E_added, const Real* mom_added);

#endif

  void ca_summass
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(rho),
     const Real* dx, BL_FORT_FAB_ARG_3D(vol), Real* mass);

  void ca_sumsquared
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(rho),
     const Real* dx, BL_FORT_FAB_ARG_3D(vol), Real* mass);

  void ca_sumlocmass
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(rho),
     const Real* dx, BL_FORT_FAB_ARG_3D(vol), Real* mass, const int& idir);

  void ca_sumlocmass2d
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(rho),
     const Real* dx, BL_FORT_FAB_ARG_3D(vol), Real* mass, const int& idir1, const int& idir2);
  
  void ca_sumlocsquaredmass
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(rho),
     const Real* dx, BL_FORT_FAB_ARG_3D(vol), Real* mass, const int& idir);
  
  void ca_sumproduct
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(f1), BL_FORT_FAB_ARG_3D(f2),
     const Real* dx, BL_FORT_FAB_ARG_3D(vol), Real* s);

#ifdef REACTIONS
#ifdef SDC
  void ca_react_state
    (const int* lo, const int* hi,
     const Real* uold, const int* uo_lo, const int* uo_hi,
     const Real* unew, const int* un_lo, const int* un_hi,
     const Real* asrc, const int* as_lo, const int* as_hi,
     const Real* reactions, const int* r_lo, const int* r_hi,
     const int* mask, const int* m_lo, const int* m_hi,
     const Real& time, const Real& dt_react);
#else
  void ca_react_state
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(reactions),
     BL_FORT_FAB_ARG_3D(weights),
     const BL_FORT_IFAB_ARG_3D(mask),
     const Real& time, const Real& dt_react);
#endif
#endif

  void ca_tempdiffextrap
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(TempDiffTerm));

  void make_enthalpy
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
           BL_FORT_FAB_ARG_3D(enth));

  void ca_fill_temp_cond
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const BL_FORT_FAB_ARG_3D(xcoeffs),
     const BL_FORT_FAB_ARG_3D(ycoeffs),
     const BL_FORT_FAB_ARG_3D(zcoeffs));

  void ca_fill_enth_cond
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const BL_FORT_FAB_ARG_3D(xcoeffs),
     const BL_FORT_FAB_ARG_3D(ycoeffs),
     const BL_FORT_FAB_ARG_3D(zcoeffs));
  
  void ca_fill_spec_coeff
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const BL_FORT_FAB_ARG_3D(xcoeffs),
     const BL_FORT_FAB_ARG_3D(ycoeffs),
     const BL_FORT_FAB_ARG_3D(zcoeffs));

  void ca_fill_first_visc_coeff
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const BL_FORT_FAB_ARG_3D(xcoeffs),
     const BL_FORT_FAB_ARG_3D(ycoeffs),
     const BL_FORT_FAB_ARG_3D(zcoeffs));

  void ca_fill_secnd_visc_coeff
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const BL_FORT_FAB_ARG_3D(xcoeffs),
     const BL_FORT_FAB_ARG_3D(ycoeffs),
     const BL_FORT_FAB_ARG_3D(zcoeffs));

  void ca_compute_div_tau_u
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(div_tau_u),
     const BL_FORT_FAB_ARG_3D(state),
     const Real* dx, const int* coord_type);


#ifdef AUX_UPDATE
  void ca_auxupdate
    (BL_FORT_FAB_ARG(state_old),
     BL_FORT_FAB_ARG(state_new),
     const int* lo, const int* hi,
     const Real * dt);
#endif

#ifdef NEUTRINO
  void ca_setgroup(const int& igroup);
#endif

}

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

/* problem-specific includes */
#include <Problem_F.H>

#endif
