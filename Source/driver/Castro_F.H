#ifndef _Castro_F_H_
#define _Castro_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void ca_network_init();

  void ca_network_finalize();

  void burner_init();

  void ca_extern_init(const int* name, const int* namlen);

  void ca_get_num_spec(int* nspec);

  void ca_get_num_aux(int* naux);

  void ca_get_spec_names(int* spec_names, int* ispec, int* len);

  void ca_get_spec_az(int* ispec, amrex::Real* A, amrex::Real* Z);

  void ca_get_aux_names(int* aux_names, int* iaux, int* len);

  void ca_get_qvar(int* qvar);

#ifdef RADIATION
  void ca_get_qradvar(int* qradvar);
#endif

  void ca_get_nqaux(int* nqaux);

  void ca_set_special_tagging_flag(const amrex::Real& max_den, const int* flag);

  void ca_set_amr_info(const int& level, const int& iteration, const int& ncycle,
		       const amrex::Real& time, const amrex::Real& dt);

  void ca_get_method_params(int* HYP_GROW);

#ifdef SELF_GRAVITY
#if (BL_SPACEDIM > 1)
  void get_numpts_1d(int* numpts_1d);
  void set_numpts_1d(int* numpts_1d);
  void allocate_outflow_data(const int* n1d, const int* nc);
  void set_old_outflow_data(amrex::Real* radial, amrex::Real* time,
			    int* numpts, int* nc);
  void set_new_outflow_data(amrex::Real* radial, amrex::Real* time,
			    int* numpts, int* nc);
  void swap_outflow_data();
#endif
#endif

  void ca_set_method_params
    (const int& dm,
     const int& Density, const int& Xmom, 
     const int& Eden,    const int& Eint, 
     const int& Temp     , const int& FirstAdv, 
     const int& FirstSpec, const int& FirstAux, 
     const int& NumAdv,
#ifdef SHOCK_VAR
     const int& Shock,
#endif
     const int* gravity_type, const int& gravity_type_len);

  void ca_finalize_meth_params();

  void ca_set_castro_method_params();

#ifdef RADIATION
  void ca_init_godunov_indices_rad();
#else
  void ca_init_godunov_indices();
#endif

  void ca_set_problem_params
    (const int& dm,
     const int* physbc_lo, const int* physbc_hi,
     const int& Interior_value, const int& Inflow_value,
     const int& Outflow_value, const int& Symmetry_value,
     const int& SlipWall_value, const int& NoSlipWall_value,
     const int& coord_type,
     const amrex::Real* problo, const amrex::Real* probhi,
     const amrex::Real* center);

  void ca_set_grid_info
    (const int& max_level, const amrex::Real* dx_level,
     const int* domlo_level, const int* domhi_level,
     const int* ref_ratio, const int* n_error_buf,
     const int* blocking_factor);

  void ca_get_tagging_params(const int* name, const int* namlen);
#ifdef SPONGE
  void ca_get_sponge_params(const int* name, const int* namlen);
#endif

#ifdef SELF_GRAVITY
#if (BL_SPACEDIM > 1)
  void ca_compute_avgstate
    (const int* lo, const int* hi, 
     const amrex::Real* dx, const amrex::Real* dr, const int* nc,
     const BL_FORT_FAB_ARG_3D(S  ), const amrex::Real* avgden,
     const BL_FORT_FAB_ARG_3D(Vol), const amrex::Real* avgvol,
     const amrex::Real* problo, const int* numpts_1d);
#endif
#endif

  void ca_apply_metric
    (const int* lo, const int* hi,
     D_DECL(const int* xlo, const int *ylo, const int *zlo),
     D_DECL(const int* xhi, const int *yhi, const int *zhi),
     BL_FORT_FAB_ARG(cc),
     D_DECL(const BL_FORT_FAB_ARG(xedge),
            const BL_FORT_FAB_ARG(yedge),
            const BL_FORT_FAB_ARG(zedge)), 
            const amrex::Real* dx, const int* coord_type);

  void ca_weight_cc
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(cc),
     const amrex::Real* dx, const int* coord_type);

  void ca_unweight_cc
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(cc),
     const amrex::Real* dx, const int* coord_type);

  void ca_unweight_edges
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG(edge),
     const amrex::Real* dx, const int* coord_type, const int* idir);


#ifdef HYBRID_MOMENTUM
  void ca_init_hybrid_momentum
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state));

  void ca_hybrid_hydro_source
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(sources));

  void ca_hybrid_update
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state));
#endif

  void ca_ctoprim(const int* lo, const int* hi,
	          const amrex::Real* u, const int* u_lo, const int* u_hi,
#ifdef RADIATION
		  const amrex::Real* Er, const int* Er_lo, const int* Er_hi,
		  const amrex::Real* lam, const int* lam_lo, const int* lam_hi,
#endif
		  const amrex::Real* q, const int* q_lo, const int* q_hi,
		  const amrex::Real* qaux, const int* qa_lo, const int* qa_hi, const int* idx);

  void ca_srctoprim(const int* lo,   const int* hi,
		    const amrex::Real* q,   const int* q_lo,   const int* q_hi,
		    const amrex::Real* qaux, const int* qa_lo, const int* qa_hi,
		    const amrex::Real* src, const int* src_lo, const int* src_hi,
		    const amrex::Real* srQ, const int* srQ_lo, const int* srQ_hi, const int* idx);

  void ca_check_initial_species
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state), const int* idx);

  void ca_enforce_minimum_density
    (const amrex::Real* S_old, const int* s_old_lo, const int* s_old_hi,
     const amrex::Real* S_new, const int* s_new_lo, const int* s_new_hi,
     const amrex::Real* vol, const int* vol_lo, const int* vol_hi,
     const int* lo, const int* hi,
     const amrex::Real* frac_change, const int* verbose, const int* idx);

  void ca_normalize_species
    (BL_FORT_FAB_ARG_3D(S_new), const int* lo, const int* hi, const int* idx);

  void ca_get_center(amrex::Real* center);
  void ca_set_center(amrex::Real* center);
  void ca_find_center(amrex::Real* data, amrex::Real* center, int* icen,
		      const amrex::Real* dx, const amrex::Real* problo);

  void ca_compute_temp
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(state), const int* idx);

  void ca_enforce_consistent_e
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state), const int* idx);

  void ca_reset_internal_e
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(S_new),
     const int& verbose, const int* idx);

#ifdef DIMENSION_AGNOSTIC

  void ca_hypfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo, 
     const amrex::Real* time, const int* bc);

  void ca_denfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo, 
     const amrex::Real* time, const int* bc);

  void ca_xmomfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

  void ca_ymomfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

  void ca_zmomfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

  void ca_zerofill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

  void ca_default_fill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

#ifdef RADIATION
  void ca_radfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);
#endif

#ifdef REACTIONS
  void ca_reactfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo, 
     const amrex::Real* time, const int* bc);
#endif

#ifdef SELF_GRAVITY
  void ca_phigravfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo, 
     const amrex::Real* time, const int* bc);

  void ca_gravxfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

  void ca_gravyfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

  void ca_gravzfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);
#endif

#ifdef ROTATION
  void ca_phirotfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

  void ca_rotxfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

  void ca_rotyfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

  void ca_rotzfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);
#endif

  inline void ca_nullfill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo, 
     const amrex::Real* time, const int* bc)
#ifdef __cplusplus
  { amrex::Abort("Not supposed to call ca_nullfill"); }
#else
    ;
#endif

#else

  void ca_hypfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void ca_denfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void ca_xmomfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void ca_ymomfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void ca_zmomfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void ca_zerofill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void ca_default_fill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

#ifdef RADIATION
  void ca_radfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);
#endif

#ifdef REACTIONS
  void ca_reactfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);
#endif

#ifdef SELF_GRAVITY
  void ca_phigravfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void ca_gravxfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void ca_gravyfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void ca_gravzfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);
#endif

#ifdef ROTATION
  void ca_phirotfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void ca_rotxfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void ca_rotyfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);

  void ca_rotzfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[]);
#endif

  inline void ca_nullfill
    (BL_FORT_FAB_ARG(state),
     const int dlo[], const int dhi[],
     const amrex::Real dx[], const amrex::Real glo[],
     const amrex::Real* time, const int bc[])
#ifdef __cplusplus
  { amrex::Abort("Not supposed to call ca_nullfill"); }
#else
    ;
#endif

#endif

  void ca_estdt
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, amrex::Real* dt);

#ifdef DIFFUSION
  void ca_estdt_temp_diffusion
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, amrex::Real* dt);
  void ca_estdt_enth_diffusion
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, amrex::Real* dt);
#endif

#ifdef RADIATION
  void ca_estdt_rad
    (const BL_FORT_FAB_ARG(state),
     const BL_FORT_FAB_ARG(gpr),
     const int lo[], const int hi[],
     const amrex::Real dx[], amrex::Real* dt);
#endif

#ifdef REACTIONS
    void ca_estdt_burning
    (const BL_FORT_FAB_ARG_3D(state_old),
     const BL_FORT_FAB_ARG_3D(state_new),
     const BL_FORT_FAB_ARG_3D(reactions_old),
     const BL_FORT_FAB_ARG_3D(reactions_new),
     const int* lo, const int* hi,
     const amrex::Real* dx, const amrex::Real* dt_old, amrex::Real* dt);
#endif

  void ca_check_timestep
    (const BL_FORT_FAB_ARG_3D(state_old),
     const BL_FORT_FAB_ARG_3D(state_new),
#ifdef REACTIONS
     const BL_FORT_FAB_ARG_3D(reactions_old),
     const BL_FORT_FAB_ARG_3D(reactions_new),
#endif
     const int* lo, const int* hi,
     const amrex::Real* dx, const amrex::Real* dt_old, amrex::Real* dt);

  void ca_ctu_update
    (const int* is_finest_level,
     const amrex::Real* time,
     const int    lo[], const int    hi[],
     const int domlo[], const int domhi[],
     const BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(stateout),
#ifdef RADIATION
     BL_FORT_FAB_ARG_3D(Er), 
     BL_FORT_FAB_ARG_3D(Erout),
#endif
     BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(qaux),
     BL_FORT_FAB_ARG_3D(src_q),
     BL_FORT_FAB_ARG_3D(update),
     const amrex::Real dx[], const amrex::Real* dt,
     D_DECL(BL_FORT_FAB_ARG_3D(xflux),
            BL_FORT_FAB_ARG_3D(yflux),
            BL_FORT_FAB_ARG_3D(zflux)),
#ifdef RADIATION
     D_DECL(BL_FORT_FAB_ARG_3D(rxflux),
            BL_FORT_FAB_ARG_3D(ryflux),
            BL_FORT_FAB_ARG_3D(rzflux)),
#endif
#if (BL_SPACEDIM < 3)
     BL_FORT_FAB_ARG_3D(pradial),
#endif
     D_DECL(const BL_FORT_FAB_ARG_3D(xarea),
            const BL_FORT_FAB_ARG_3D(yarea),
            const BL_FORT_FAB_ARG_3D(zarea)),
#if (BL_SPACEDIM < 3)
     const BL_FORT_FAB_ARG_3D(dloga),
#endif
     const BL_FORT_FAB_ARG_3D(volume),
     amrex::Real* cflLoc,
     const int&  verbose,
#ifdef RADIATION
     const int* priv_nstep_fsp,
#endif
     amrex::Real& mass_lost,
     amrex::Real& xmom_lost,
     amrex::Real& ymom_lost,
     amrex::Real& zmom_lost,
     amrex::Real& eden_lost,
     amrex::Real& xang_lost,
     amrex::Real& yang_lost,
     amrex::Real& zang_lost);


  void ca_mol_single_stage
    (const amrex::Real* time,
     const int    lo[], const int    hi[],
     const int domlo[], const int domhi[],
     const amrex::Real *stage_weight,
     const BL_FORT_FAB_ARG_3D(statein),
     BL_FORT_FAB_ARG_3D(stateout),
     BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(qaux), 
     BL_FORT_FAB_ARG_3D(srcU),
     BL_FORT_FAB_ARG_3D(update),
     BL_FORT_FAB_ARG_3D(update_flux),
     const amrex::Real dx[], const amrex::Real* dt,
     D_DECL(BL_FORT_FAB_ARG_3D(flux1),
	    BL_FORT_FAB_ARG_3D(flux2),
	    BL_FORT_FAB_ARG_3D(flux3)),
#if (BL_SPACEDIM < 3)
     BL_FORT_FAB_ARG_3D(pradial), 
#endif
     D_DECL(const BL_FORT_FAB_ARG_3D(area1),
	    const BL_FORT_FAB_ARG_3D(area2),
	    const BL_FORT_FAB_ARG_3D(area3)),
#if (BL_SPACEDIM < 3)
     BL_FORT_FAB_ARG_3D(dloga),
#endif
     const BL_FORT_FAB_ARG_3D(volume),
     amrex::Real* cflLoc,
     const int&  verbose);


#ifdef RADIATION
  void ca_inelastic_sct
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(stateout),
     BL_FORT_FAB_ARG_3D(Erout),
     const BL_FORT_FAB_ARG_3D(kappa_s),
     const amrex::Real& dt);
#endif

#ifdef SPONGE
  void ca_sponge
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(source),
     BL_FORT_FAB_ARG_3D(vol),
     const amrex::Real* dx, const amrex::Real& dt, const amrex::Real* time);

  void update_sponge_params(const amrex::Real* time);
#endif

  void flush_output();

#ifdef POINTMASS
  void pm_compute_delta_mass
    (const amrex::Real* pm_delta_mass, const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(stateout),
     const BL_FORT_FAB_ARG_3D(volume),
     const amrex::Real* problo, const amrex::Real* dx,
     const amrex::Real* time, const amrex::Real* dt);

  void pm_fix_solution
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(stateold),
     BL_FORT_FAB_ARG_3D(statenew),
     const amrex::Real* problo, const amrex::Real* dx,
     const amrex::Real* time, const amrex::Real* dt);

  void set_pointmass(const amrex::Real* pointmass);
#endif

#ifdef GRAVITY
  void ca_gsrc
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const BL_FORT_FAB_ARG_3D(s_old),
#ifdef SELF_GRAVITY
     const BL_FORT_FAB_ARG_3D(phi),
     const BL_FORT_FAB_ARG_3D(grav),
#endif
     BL_FORT_FAB_ARG_3D(source),
     const amrex::Real* dx, const amrex::Real& dt, const amrex::Real* time);

  void ca_corrgsrc
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const BL_FORT_FAB_ARG_3D(S_old),
     const BL_FORT_FAB_ARG_3D(S_new),
#ifdef SELF_GRAVITY
     const BL_FORT_FAB_ARG_3D(phi_old),
     const BL_FORT_FAB_ARG_3D(phi_new),
     const BL_FORT_FAB_ARG_3D(grav_old),
     const BL_FORT_FAB_ARG_3D(grav_new),
#endif
     const BL_FORT_FAB_ARG_3D(volume),
     BL_FORT_FAB_ARG_3D(xflux),
     BL_FORT_FAB_ARG_3D(yflux),
     BL_FORT_FAB_ARG_3D(zflux),
     BL_FORT_FAB_ARG_3D(source),
     const amrex::Real* dx, const amrex::Real& dt, const amrex::Real* time);
#endif

#ifdef ROTATION

  void ca_fill_rotational_potential
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(phi),
     const amrex::Real* dx, const amrex::Real& time);

  void ca_fill_rotational_acceleration
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(rot),
     const BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, const amrex::Real& time);

  void ca_rsrc
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const BL_FORT_FAB_ARG_3D(phi),
     const BL_FORT_FAB_ARG_3D(rot),
     const BL_FORT_FAB_ARG_3D(s_old),
     BL_FORT_FAB_ARG_3D(source),
     const BL_FORT_FAB_ARG_3D(vol),
     const amrex::Real* dx, const amrex::Real& dt, const amrex::Real* time);

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
     const amrex::Real* dx, const amrex::Real& dt, const amrex::Real* time,
     const BL_FORT_FAB_ARG_3D(volume));

#endif

  void ca_summass
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(rho),
     const amrex::Real* dx, BL_FORT_FAB_ARG_3D(vol), amrex::Real* mass);

  void ca_sumsquared
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(rho),
     const amrex::Real* dx, BL_FORT_FAB_ARG_3D(vol), amrex::Real* mass);

  void ca_sumlocmass
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(rho),
     const amrex::Real* dx, BL_FORT_FAB_ARG_3D(vol), amrex::Real* mass, const int& idir);

  void ca_sumlocmass2d
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(rho),
     const amrex::Real* dx, BL_FORT_FAB_ARG_3D(vol), amrex::Real* mass, const int& idir1, const int& idir2);
  
  void ca_sumlocsquaredmass
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(rho),
     const amrex::Real* dx, BL_FORT_FAB_ARG_3D(vol), amrex::Real* mass, const int& idir);
  
  void ca_sumproduct
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(f1), BL_FORT_FAB_ARG_3D(f2),
     const amrex::Real* dx, BL_FORT_FAB_ARG_3D(vol), amrex::Real* s);

#ifdef REACTIONS
#ifdef SDC
  void ca_react_state
    (const int* lo, const int* hi,
     const amrex::Real* uold, const int* uo_lo, const int* uo_hi,
     const amrex::Real* unew, const int* un_lo, const int* un_hi,
     const amrex::Real* asrc, const int* as_lo, const int* as_hi,
     const amrex::Real* reactions, const int* r_lo, const int* r_hi,
     const int* mask, const int* m_lo, const int* m_hi,
     const amrex::Real& time, const amrex::Real& dt_react, const int& sdc_iter);
#else
  void ca_react_state
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(reactions),
     BL_FORT_FAB_ARG_3D(weights),
     const BL_FORT_IFAB_ARG_3D(mask),
     const amrex::Real& time, const amrex::Real& dt_react, const int& strang_half);
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
     const amrex::Real* dx, const int* coord_type);


#ifdef AUX_UPDATE
  void ca_auxupdate
    (BL_FORT_FAB_ARG(state_old),
     BL_FORT_FAB_ARG(state_new),
     const int* lo, const int* hi,
     const amrex::Real * dt);
#endif

#ifdef NEUTRINO
  void ca_setgroup(const int& igroup);
#endif

#ifdef __cplusplus
}
#endif

#ifdef DIMENSION_AGNOSTIC

BL_FORT_PROC_DECL(CA_INITDATA,ca_initdata)
    (const int& level, const amrex::Real& time,
     const int* lo, const int* hi,
     const int& num_state,
     BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* xhi);

#else

BL_FORT_PROC_DECL(CA_INITDATA,ca_initdata)
    (const int& level, const amrex::Real& time,
     const int* lo, const int* hi,
     const int& num_state,
     BL_FORT_FAB_ARG(state),
     const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* xhi);

#endif

#ifdef MAESTRO_INIT
BL_FORT_PROC_DECL(CA_INITDATA_MAESTRO,ca_initdata_maestro)
  (const int* lo, const int* hi, const int& MAESTRO_init_type,
     BL_FORT_FAB_ARG(state),
     const amrex::Real dx[], const amrex::Real& dr,
     const amrex::Real xlo[], const amrex::Real xhi[],
     const amrex::Real p0[], const int& MAESTRO_npts_model, const int& level);

BL_FORT_PROC_DECL(CA_INITDATA_MAKEMODEL,ca_initdata_makemodel)
  (Real model[], const int& model_size, const int& MAESTRO_npts_model,
   const amrex::Real rho0[], const amrex::Real tempbar[], 
   const amrex::Real dx[], const amrex::Real& dr, const int& r_model_start);

BL_FORT_PROC_DECL(CA_INITDATA_OVERWRITE,ca_initdata_overwrite)
  (const int* lo, const int* hi,
   BL_FORT_FAB_ARG(state),
   const amrex::Real model[], const int& model_size, const int& MAESTRO_npts_model,
   const amrex::Real dx[], const amrex::Real& dr, const amrex::Real xlo[], const amrex::Real xhi[],
   const int& r_model_start);
#endif

#ifdef RADIATION
BL_FORT_PROC_DECL(CA_INITRAD,ca_initrad)
  (const int& level, const amrex::Real& time,
   const int* lo, const int* hi,
   const int& num_state,
   BL_FORT_FAB_ARG(state),
   const amrex::Real dx[], const amrex::Real xlo[], const amrex::Real xhi[]);
#endif

#ifdef DIMENSION_AGNOSTIC
BL_FORT_PROC_DECL(CA_EXT_SRC,ca_ext_src)
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(old_state),
     const BL_FORT_FAB_ARG_3D(new_state),
     BL_FORT_FAB_ARG_3D(ext_src),
     const amrex::Real* prob_lo, const amrex::Real* dx,
     const amrex::Real* time, const amrex::Real* dt);
#else
BL_FORT_PROC_DECL(CA_EXT_SRC,ca_ext_src)
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG(old_state),
     const BL_FORT_FAB_ARG(new_state),
     BL_FORT_FAB_ARG(ext_src),
     const amrex::Real* prob_lo, const amrex::Real* dx,
     const amrex::Real* time, const amrex::Real* dt);
#endif

/* problem-specific includes */
#include <Problem_F.H>

#endif
