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

  void ca_eos_finalize();

  void burner_init();

  void ca_extern_init(const int* name, const int* namlen);

  void ca_get_num_spec(int* nspec);

  void ca_get_num_aux(int* naux);

  void ca_get_num_adv(int* nadv);

  void ca_get_spec_names(int* spec_names, const int* ispec, int* len);

  void ca_get_spec_az(const int* ispec, amrex::Real* A, amrex::Real* Z);

  void ca_get_aux_names(int* aux_names,  const int* iaux, int* len);

  void ca_get_qvar(int* qvar);

  void ca_get_nq(int* nq);

  void ca_get_nqaux(int* nqaux);

  void ca_get_ngdnv(int* ngdnv);

  void ca_amrinfo_init();
  void ca_amrinfo_finalize();

  void ca_set_amr_info(const int& level, const int& iteration, const int& ncycle,
		       const amrex::Real& time, const amrex::Real& dt);

  void ca_get_method_params(int* HYP_GROW);

  void runtime_pretty_print(int* jobinfo_file_name, const int* jobinfo_file_length);

#ifdef SELF_GRAVITY
#if (BL_SPACEDIM > 1)
  void get_numpts_1d(int* numpts_1d);
  void set_numpts_1d(int* numpts_1d);
  void allocate_outflow_data(const int* n1d, const int* nc);
  void set_old_outflow_data(const amrex::Real* radial, const amrex::Real* time,
			    const int* numpts, const int* nc);
  void set_new_outflow_data(const amrex::Real* radial, const amrex::Real* time,
			    const int* numpts, const int* nc);
  void swap_outflow_data();
#endif
#endif

  void ca_set_method_params
    (const int& dm,
     const int& Density, const int& Xmom,
#ifdef HYBRID_MOMENTUM
     const int& Rmom,
#endif
     const int& Eden,    const int& Eint,
     const int& Temp     , const int& FirstAdv,
     const int& FirstSpec, const int& FirstAux,
#ifdef SHOCK_VAR
     const int& Shock,
#endif
#ifdef MHD
     const int& QMAGX_in, const int& QMAGY_in, const int& QMAGZ_in,
#endif
#ifdef RADIATION
     const int& QPTOT_in, const int& QREITOT_in, const int& QRAD_in,
#endif
     const int& QRHO_in,
     const int& QU_in, const int& QV_in, const int& QW_in,
     const int& QGAME_in, const int& QPRES_in, const int& QREINT_in,
     const int& QTEMP_in,
     const int& QFA_in, const int& QFS_in, const int& QFX_in);

  void ca_finalize_meth_params();

  void ca_set_castro_method_params();

  void ca_init_godunov_indices();

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
  void ca_allocate_sponge_params();
  void ca_deallocate_sponge_params();
#endif

#ifdef SELF_GRAVITY
#if (BL_SPACEDIM > 1)
  void ca_compute_avgstate
    (const int* lo, const int* hi,
     const amrex::Real* dx, const amrex::Real* dr, const int* nc,
     const BL_FORT_FAB_ARG_3D(S  ),  amrex::Real* avgden,
     const BL_FORT_FAB_ARG_3D(Vol),  amrex::Real* avgvol,
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
     BL_FORT_FAB_ARG_3D(cc),
     const amrex::Real* dx, const int coord_type);

  void ca_unweight_cc
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(cc),
     const amrex::Real* dx, const int coord_type);

  void ca_unweight_edges
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(edge),
     const amrex::Real* dx, const int coord_type, const int idir);


#ifdef HYBRID_MOMENTUM
  void ca_init_hybrid_momentum
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state));

  void ca_hybrid_hydro_source
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(sources),
     const amrex::Real mult_factor);

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
		  amrex::Real* q, const int* q_lo, const int* q_hi,
		  amrex::Real* qaux, const int* qa_lo, const int* qa_hi);

  void ca_srctoprim(const int* lo,   const int* hi,
		    const amrex::Real* q,   const int* q_lo,   const int* q_hi,
		    const amrex::Real* qaux, const int* qa_lo, const int* qa_hi,
		    const amrex::Real* src, const int* src_lo, const int* src_hi,
		          amrex::Real* srQ, const int* srQ_lo, const int* srQ_hi);

  void ca_compute_cfl(const int* lo, const int* hi,
                      const BL_FORT_FAB_ARG_ANYD(q),
                      const BL_FORT_FAB_ARG_ANYD(qaux),
                      const amrex::Real dt, const amrex::Real* dx,
                      amrex::Real* courno, const int verbose);

  void ca_check_initial_species
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(state));

  void ca_enforce_minimum_density
    (const int* lo, const int* hi,
     const amrex::Real* S_old, const int* s_old_lo, const int* s_old_hi,
           amrex::Real* S_new, const int* s_new_lo, const int* s_new_hi,
     const amrex::Real* vol, const int* vol_lo, const int* vol_hi,
           amrex::Real* frac_change, const int verbose);

  void ca_normalize_species
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(S_new));

  void ca_get_center(amrex::Real* center);
  void ca_set_center(const amrex::Real* center);
  void ca_find_center(amrex::Real* data, amrex::Real* center, const int* icen,
		      const amrex::Real* dx, const amrex::Real* problo);

  void ca_compute_temp
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state));

  void ca_enforce_consistent_e
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state));

  void ca_reset_internal_e
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(S_new),
     const int verbose);

  void ca_generic_single_fill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

  void ca_generic_multi_fill
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real* time, const int* bc);

#ifdef AMREX_DIMENSION_AGNOSTIC

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
    (const int lo[], const int hi[],
     const BL_FORT_FAB_ARG(state),
     const BL_FORT_FAB_ARG(gpr),
     const amrex::Real dx[], amrex::Real* dt);
#endif

#ifdef REACTIONS
    void ca_estdt_burning
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state_old),
     const BL_FORT_FAB_ARG_3D(state_new),
     const BL_FORT_FAB_ARG_3D(reactions_old),
     const BL_FORT_FAB_ARG_3D(reactions_new),
     const amrex::Real* dx, const amrex::Real* dt_old, amrex::Real* dt);
#endif

  void ca_check_timestep
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state_old),
     const BL_FORT_FAB_ARG_3D(state_new),
#ifdef REACTIONS
     const BL_FORT_FAB_ARG_3D(reactions_old),
     const BL_FORT_FAB_ARG_3D(reactions_new),
#endif
     const amrex::Real* dx, const amrex::Real* dt_old, amrex::Real* dt);

  void divu
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(q),
     const amrex::Real* dx,
     BL_FORT_FAB_ARG_3D(div));

  void ca_uflatten
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(flatn), int pres_comp);

  void ca_rad_flatten
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(flatn),
     BL_FORT_FAB_ARG_3D(flatg));

  void ctu_clean_fluxes
    (const int* lo, const int* hi,
     const int idir_f,
     const BL_FORT_FAB_ARG_3D(sborder),
     const BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(flux),
#ifdef RADIATION
     const BL_FORT_FAB_ARG_3D(Er),
     BL_FORT_FAB_ARG_3D(radflux),
#endif
     const BL_FORT_FAB_ARG_3D(area),
     const BL_FORT_FAB_ARG_3D(volume),
     const BL_FORT_FAB_ARG_3D(div),
     const amrex::Real* dx, const amrex::Real dt);

  void ctu_consup
    (const int*lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(Sborder),
     const BL_FORT_FAB_ARG_3D(q),
     const BL_FORT_FAB_ARG_3D(shk),
     BL_FORT_FAB_ARG_3D(S_new),
     BL_FORT_FAB_ARG_3D(hydro_source),
     D_DECL(const BL_FORT_FAB_ARG_3D(flux0),
            const BL_FORT_FAB_ARG_3D(flux1),
            const BL_FORT_FAB_ARG_3D(flux2)),
#ifdef RADIATION
     const BL_FORT_FAB_ARG_3D(Erborder),
     BL_FORT_FAB_ARG_3D(Er_out),
     D_DECL(const BL_FORT_FAB_ARG_3D(rad_flux0),
            const BL_FORT_FAB_ARG_3D(rad_flux1),
            const BL_FORT_FAB_ARG_3D(rad_flux2)),
     int* nstep_fsp,
#endif
     D_DECL(const BL_FORT_FAB_ARG_3D(qe0),
            const BL_FORT_FAB_ARG_3D(qe1),
            const BL_FORT_FAB_ARG_3D(qe2)),
     D_DECL(const BL_FORT_FAB_ARG_3D(area0),
            const BL_FORT_FAB_ARG_3D(area1),
            const BL_FORT_FAB_ARG_3D(area2)),
     const BL_FORT_FAB_ARG_3D(volume),
     const BL_FORT_FAB_ARG_3D(pdivu),
     const amrex::Real* dx, const amrex::Real dt);

  void scale_flux
    (const int* lo, const int* hi,
#if AMREX_SPACEDIM == 1
     const BL_FORT_FAB_ARG_3D(qe0),
#endif
     BL_FORT_FAB_ARG_3D(flux),
     BL_FORT_FAB_ARG_3D(area),
     const amrex::Real dt);

  void scale_rad_flux
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(rad_flux),
     BL_FORT_FAB_ARG_3D(area),
     const amrex::Real dt);

  void store_pradial
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(qe),
     BL_FORT_FAB_ARG_3D(pradial),
     const amrex::Real dt);

  void ca_track_grid_losses
    (const int* lo, const int* hi,
     D_DECL(const BL_FORT_FAB_ARG_3D(flux0),
            const BL_FORT_FAB_ARG_3D(flux1),
            const BL_FORT_FAB_ARG_3D(flux2)),
     amrex::Real& mass_lost, amrex::Real& xmom_lost, amrex::Real& ymom_lost, amrex::Real& zmom_lost,
     amrex::Real& eden_lost, amrex::Real& xang_lost, amrex::Real& yang_lost, amrex::Real& zang_lost);

  void cmpflx_plus_godunov
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(ql),
     const BL_FORT_FAB_ARG_3D(qr), const int nc, const int comp,
     BL_FORT_FAB_ARG_3D(flux),
     BL_FORT_FAB_ARG_3D(q_int),
#ifdef RADIATION
     BL_FORT_FAB_ARG_3D(radflux),
     BL_FORT_FAB_ARG_3D(lambda_int),
#endif
     BL_FORT_FAB_ARG_3D(qe),
     const BL_FORT_FAB_ARG_3D(qaux),
     const BL_FORT_FAB_ARG_3D(shk),
     const int idir, const int* domlo, const int* domhi);


#ifdef AMREX_USE_CUDA
  void ca_ppm_reconstruct
    (const int* lo, const int* hi, int put_on_edges,
     BL_FORT_FAB_ARG_3D(q), int nc, int nstart, int nend,
     BL_FORT_FAB_ARG_3D(flatn),
     BL_FORT_FAB_ARG_3D(qm),
     BL_FORT_FAB_ARG_3D(qp), int nqc, int nqstart, int nqend);

  void ca_shock
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(shk),
     const amrex::Real* dx);

  void ca_store_godunov_state
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(qint),
#ifdef RADIATION
     BL_FORT_FAB_ARG_3D(lambda_int),
#endif
     BL_FORT_FAB_ARG_3D(qgdnv));

  void ca_construct_flux_cuda
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const amrex::Real* dx, const amrex::Real dt,
     const int idir,
     const BL_FORT_FAB_ARG_3D(statein),
     BL_FORT_FAB_ARG_3D(div),
     BL_FORT_FAB_ARG_3D(qaux),
     BL_FORT_FAB_ARG_3D(shk),
     BL_FORT_FAB_ARG_3D(qm),
     BL_FORT_FAB_ARG_3D(qp),
     BL_FORT_FAB_ARG_3D(qe),
     BL_FORT_FAB_ARG_3D(flux),
     BL_FORT_FAB_ARG_3D(area));

  void ca_construct_hydro_update_cuda
    (const int* lo, const int* hi,
     const amrex::Real* dx, const amrex::Real dt,
     BL_FORT_FAB_ARG_3D(qx),
     BL_FORT_FAB_ARG_3D(qy),
     BL_FORT_FAB_ARG_3D(qz),
     BL_FORT_FAB_ARG_3D(xflux),
     BL_FORT_FAB_ARG_3D(yflux),
     BL_FORT_FAB_ARG_3D(zflux),
     BL_FORT_FAB_ARG_3D(xarea),
     BL_FORT_FAB_ARG_3D(yarea),
     BL_FORT_FAB_ARG_3D(zarea),
     const BL_FORT_FAB_ARG_3D(volume),
     const BL_FORT_FAB_ARG_3D(srcU),
     BL_FORT_FAB_ARG_3D(update));
#endif

  void ca_ctu_update
    (const int    lo[], const int    hi[],
     const int* is_finest_level,
     const amrex::Real time,
     const int domlo[], const int domhi[],
     const BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(stateout),
#ifdef RADIATION
     BL_FORT_FAB_ARG_3D(Er),
     BL_FORT_FAB_ARG_3D(Erout),
#endif
     BL_FORT_FAB_ARG_3D(q),
     const BL_FORT_FAB_ARG_3D(qaux),
     const BL_FORT_FAB_ARG_3D(src_q),
     BL_FORT_FAB_ARG_3D(update),
     const amrex::Real dx[], const amrex::Real dt,
     D_DECL(BL_FORT_FAB_ARG_3D(xflux),
            BL_FORT_FAB_ARG_3D(yflux),
            BL_FORT_FAB_ARG_3D(zflux)),
#ifdef RADIATION
     D_DECL(BL_FORT_FAB_ARG_3D(rxflux),
            BL_FORT_FAB_ARG_3D(ryflux),
            BL_FORT_FAB_ARG_3D(rzflux)),
#endif
     D_DECL(const BL_FORT_FAB_ARG_3D(xarea),
            const BL_FORT_FAB_ARG_3D(yarea),
            const BL_FORT_FAB_ARG_3D(zarea)),
#if (BL_SPACEDIM < 3)
     BL_FORT_FAB_ARG_3D(pradial),
     const BL_FORT_FAB_ARG_3D(dloga),
#endif
     const BL_FORT_FAB_ARG_3D(volume),
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

  void ctu_normal_states
     (const int* lo, const int* hi,
      const int* vlo, const int* vhi,
      const BL_FORT_FAB_ARG_3D(q),
      const BL_FORT_FAB_ARG_3D(flatn),
      const BL_FORT_FAB_ARG_3D(qaux),
      const BL_FORT_FAB_ARG_3D(srcQ),
      const BL_FORT_FAB_ARG_3D(shk),
      BL_FORT_FAB_ARG_3D(Ip),
      BL_FORT_FAB_ARG_3D(Im),
      BL_FORT_FAB_ARG_3D(Ip_src),
      BL_FORT_FAB_ARG_3D(Im_src),
      BL_FORT_FAB_ARG_3D(Ip_gc),
      BL_FORT_FAB_ARG_3D(Im_gc),
      BL_FORT_FAB_ARG_3D(dq),
      BL_FORT_FAB_ARG_3D(sm),
      BL_FORT_FAB_ARG_3D(sp),
      BL_FORT_FAB_ARG_3D(qxm),
      BL_FORT_FAB_ARG_3D(qxp),
#if AMREX_SPACEDIM >= 2
      BL_FORT_FAB_ARG_3D(qym),
      BL_FORT_FAB_ARG_3D(qyp),
#endif
#if AMREX_SPACEDIM == 3
      BL_FORT_FAB_ARG_3D(qzm),
      BL_FORT_FAB_ARG_3D(qzp),
#endif
      const amrex::Real* dx, const amrex::Real dt,
#if AMREX_SPACEDIM < 3
      const BL_FORT_FAB_ARG_3D(dloga),
#endif
      const int* domlo, const int* domhi);


  void transx_on_ystates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qym),
                         BL_FORT_FAB_ARG_3D(qymo),
                         const BL_FORT_FAB_ARG_3D(qyp),
                         BL_FORT_FAB_ARG_3D(qypo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fx),
#ifdef RADIATION
                         const BL_FORT_FAB_ARG_3D(rfx),
#endif
                         const BL_FORT_FAB_ARG_3D(qx),
#if AMREX_SPACEDIM == 2
                         const BL_FORT_FAB_ARG_3D(area1),
                         const BL_FORT_FAB_ARG_3D(vol),
#endif
                         const amrex::Real hdt, const amrex::Real cdtdx);

  void transx_on_zstates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qzm),
                         BL_FORT_FAB_ARG_3D(qzmo),
                         const BL_FORT_FAB_ARG_3D(qzp),
                         BL_FORT_FAB_ARG_3D(qzpo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fx),
#ifdef RADIATION
                         const BL_FORT_FAB_ARG_3D(rfx),
#endif
                         const BL_FORT_FAB_ARG_3D(qx),
                         const amrex::Real hdt, const amrex::Real cdtdx);

  void transy_on_xstates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qxm),
                         BL_FORT_FAB_ARG_3D(qxmo),
                         const BL_FORT_FAB_ARG_3D(qxp),
                         BL_FORT_FAB_ARG_3D(qxpo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fy),
#ifdef RADIATION
                         const BL_FORT_FAB_ARG_3D(rfy),
#endif
                         const BL_FORT_FAB_ARG_3D(qy),
                         const amrex::Real cdtdy);

  void transy_on_zstates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qzm),
                         BL_FORT_FAB_ARG_3D(qzmo),
                         const BL_FORT_FAB_ARG_3D(qzp),
                         BL_FORT_FAB_ARG_3D(qzpo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fy),
#ifdef RADIATION
                         const BL_FORT_FAB_ARG_3D(rfy),
#endif
                         const BL_FORT_FAB_ARG_3D(qy),
                         const amrex::Real cdtdy);

  void transz_on_xstates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qxm),
                         BL_FORT_FAB_ARG_3D(qxmo),
                         const BL_FORT_FAB_ARG_3D(qxp),
                         BL_FORT_FAB_ARG_3D(qxpo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fz),
#ifdef RADIATION
                         const BL_FORT_FAB_ARG_3D(rfz),
#endif
                         const BL_FORT_FAB_ARG_3D(qz),
                         const amrex::Real cdtdz);

  void transz_on_ystates(const int* lo, const int* hi,
                         const BL_FORT_FAB_ARG_3D(qym),
                         BL_FORT_FAB_ARG_3D(qymo),
                         const BL_FORT_FAB_ARG_3D(qyp),
                         BL_FORT_FAB_ARG_3D(qypo),
                         const BL_FORT_FAB_ARG_3D(qaux),
                         const BL_FORT_FAB_ARG_3D(fz),
#ifdef RADIATION
                         const BL_FORT_FAB_ARG_3D(rfz),
#endif
                         const BL_FORT_FAB_ARG_3D(qz),
                         const amrex::Real cdtdz);

  void transyz(const int* lo, const int* hi,
               const BL_FORT_FAB_ARG_3D(qm),
               BL_FORT_FAB_ARG_3D(qmo),
               const BL_FORT_FAB_ARG_3D(qp),
               BL_FORT_FAB_ARG_3D(qpo),
               const BL_FORT_FAB_ARG_3D(qaux),
               const BL_FORT_FAB_ARG_3D(fyz),
#ifdef RADIATION
               const BL_FORT_FAB_ARG_3D(rfyz),
#endif
               const BL_FORT_FAB_ARG_3D(fzy),
#ifdef RADIATION
               const BL_FORT_FAB_ARG_3D(rfzy),
#endif
               const BL_FORT_FAB_ARG_3D(qy),
               const BL_FORT_FAB_ARG_3D(qz),
               const amrex::Real hdt, const amrex::Real cdtdy, const amrex::Real cdtdz);

  void transxz(const int* lo, const int* hi,
               const BL_FORT_FAB_ARG_3D(qm),
               BL_FORT_FAB_ARG_3D(qmo),
               const BL_FORT_FAB_ARG_3D(qp),
               BL_FORT_FAB_ARG_3D(qpo),
               const BL_FORT_FAB_ARG_3D(qaux),
               const BL_FORT_FAB_ARG_3D(fxz),
#ifdef RADIATION
               const BL_FORT_FAB_ARG_3D(rfxz),
#endif
               const BL_FORT_FAB_ARG_3D(fzx),
#ifdef RADIATION
               const BL_FORT_FAB_ARG_3D(rfzx),
#endif
               const BL_FORT_FAB_ARG_3D(qx),
               const BL_FORT_FAB_ARG_3D(qz),
               const amrex::Real hdt, const amrex::Real cdtdx, const amrex::Real cdtdz);

  void transxy(const int* lo, const int* hi,
               const BL_FORT_FAB_ARG_3D(qm),
               BL_FORT_FAB_ARG_3D(qmo),
               const BL_FORT_FAB_ARG_3D(qp),
               BL_FORT_FAB_ARG_3D(qpo),
               const BL_FORT_FAB_ARG_3D(qaux),
               const BL_FORT_FAB_ARG_3D(fxy),
#ifdef RADIATION
               const BL_FORT_FAB_ARG_3D(rfxy),
#endif
               const BL_FORT_FAB_ARG_3D(fyx),
#ifdef RADIATION
               const BL_FORT_FAB_ARG_3D(rfyx),
#endif
               const BL_FORT_FAB_ARG_3D(qx),
               const BL_FORT_FAB_ARG_3D(qy),
               const amrex::Real hdt, const amrex::Real cdtdx, const amrex::Real cdtdy);


  void ca_make_fourth_average(const int lo[], const int hi[],
                              BL_FORT_FAB_ARG_3D(q), const int* nc,
                              BL_FORT_FAB_ARG_3D(q_bar), const int* nc_bar);

  void ca_make_fourth_in_place(const int lo[], const int hi[],
                               BL_FORT_FAB_ARG_3D(q), const int* nc);

  void ca_make_cell_center(const int lo[], const int hi[],
                           const BL_FORT_FAB_ARG_3D(U), const int* nc,
                           BL_FORT_FAB_ARG_3D(U_cc), const int* nc_cc);

  void ca_mol_single_stage
    (const int    lo[], const int    hi[],
     const amrex::Real* time,
     const int domlo[], const int domhi[],
     const amrex::Real *stage_weight,
     const BL_FORT_FAB_ARG_3D(statein),
     BL_FORT_FAB_ARG_3D(stateout),
     BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(qaux),
     const BL_FORT_FAB_ARG_3D(srcU),
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
     const int&  verbose);


  void ca_fourth_single_stage
    (const int    lo[], const int    hi[],
     const amrex::Real* time,
     const int domlo[], const int domhi[],
     const amrex::Real *stage_weight,
     const BL_FORT_FAB_ARG_3D(statein),
     BL_FORT_FAB_ARG_3D(stateout),
     BL_FORT_FAB_ARG_3D(q),
     BL_FORT_FAB_ARG_3D(q_bar),
     BL_FORT_FAB_ARG_3D(qaux),
     BL_FORT_FAB_ARG_3D(qaux_bar),
     const BL_FORT_FAB_ARG_3D(srcU),
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
     const amrex::Real* dx, const amrex::Real dt,
     const amrex::Real time, const amrex::Real mult_factor);

  void update_sponge_params(const amrex::Real* time);
#endif


  void flush_output();

#ifdef POINTMASS
  void pm_compute_delta_mass
    (const int* lo, const int* hi, amrex::Real* pm_delta_mass,
     const BL_FORT_FAB_ARG_3D(state),
     const BL_FORT_FAB_ARG_3D(stateout),
     const BL_FORT_FAB_ARG_3D(volume),
     const amrex::Real* problo, const amrex::Real* dx,
     const amrex::Real* time, const amrex::Real* dt);

  void pm_fix_solution
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(stateold),
     BL_FORT_FAB_ARG_3D(statenew),
     const amrex::Real* problo, const amrex::Real* dx,
     const amrex::Real time, const amrex::Real dt);

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
     const amrex::Real* dx, const amrex::Real dt, const amrex::Real time);

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
     const BL_FORT_FAB_ARG_3D(xflux),
     const BL_FORT_FAB_ARG_3D(yflux),
     const BL_FORT_FAB_ARG_3D(zflux),
     const BL_FORT_FAB_ARG_3D(source),
     const amrex::Real* dx, const amrex::Real& dt, const amrex::Real* time);
#endif

#ifdef ROTATION

  void ca_fill_rotational_potential
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(phi),
     const amrex::Real* dx, const amrex::Real time);

  void ca_fill_rotational_acceleration
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(rot),
     const BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, const amrex::Real time);

  void ca_rsrc
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const BL_FORT_FAB_ARG_3D(phi),
     const BL_FORT_FAB_ARG_3D(rot),
     const BL_FORT_FAB_ARG_3D(s_old),
     BL_FORT_FAB_ARG_3D(source),
     const BL_FORT_FAB_ARG_3D(vol),
     const amrex::Real* dx, const amrex::Real dt, const amrex::Real time);

  void ca_corrrsrc
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const BL_FORT_FAB_ARG_3D(phi),
     const BL_FORT_FAB_ARG_3D(rot_old),
     const BL_FORT_FAB_ARG_3D(rot_new),
     const BL_FORT_FAB_ARG_3D(S_old),
     const BL_FORT_FAB_ARG_3D(S_new),
     BL_FORT_FAB_ARG_3D(source),
     const BL_FORT_FAB_ARG_3D(xflux),
     const BL_FORT_FAB_ARG_3D(yflux),
     const BL_FORT_FAB_ARG_3D(zflux),
     const amrex::Real* dx, const amrex::Real dt, const amrex::Real time,
     const BL_FORT_FAB_ARG_3D(volume));

#endif

  void ca_summass
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(rho),
     const amrex::Real* dx, const BL_FORT_FAB_ARG_3D(vol), amrex::Real* mass);

  void ca_sumsquared
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(rho),
     const amrex::Real* dx, const BL_FORT_FAB_ARG_3D(vol), amrex::Real* mass);

  void ca_sumlocmass
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(rho),
     const amrex::Real* dx, const BL_FORT_FAB_ARG_3D(vol), amrex::Real* mass, const int& idir);

  void ca_sumlocmass2d
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(rho),
     const amrex::Real* dx, const BL_FORT_FAB_ARG_3D(vol), amrex::Real* mass, const int& idir1, const int& idir2);

  void ca_sumlocsquaredmass
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(rho),
     const amrex::Real* dx, const BL_FORT_FAB_ARG_3D(vol), amrex::Real* mass, const int& idir);

  void ca_sumproduct
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(f1), const BL_FORT_FAB_ARG_3D(f2),
     const amrex::Real* dx, const BL_FORT_FAB_ARG_3D(vol), amrex::Real* s);

#ifdef REACTIONS
#ifdef SDC
  void ca_react_state
    (const int* lo, const int* hi,
     const amrex::Real* uold, const int* uo_lo, const int* uo_hi,
     amrex::Real* unew, const int* un_lo, const int* un_hi,
     const amrex::Real* asrc, const int* as_lo, const int* as_hi,
     amrex::Real* reactions, const int* r_lo, const int* r_hi,
     const int* mask, const int* m_lo, const int* m_hi,
     const amrex::Real& time, const amrex::Real& dt_react, const int& sdc_iter);
#else
  void ca_react_state
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(reactions),
     BL_FORT_FAB_ARG_3D(weights),
     const BL_FORT_IFAB_ARG_3D(mask),
     const amrex::Real time, const amrex::Real dt_react, const int strang_half,
     amrex::Real* burn_failure);
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
     BL_FORT_FAB_ARG_3D(xcoeffs),
     BL_FORT_FAB_ARG_3D(ycoeffs),
     BL_FORT_FAB_ARG_3D(zcoeffs));

  void ca_fill_enth_cond
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(xcoeffs),
     BL_FORT_FAB_ARG_3D(ycoeffs),
     BL_FORT_FAB_ARG_3D(zcoeffs));

  void ca_fill_spec_coeff
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
      BL_FORT_FAB_ARG_3D(xcoeffs),
      BL_FORT_FAB_ARG_3D(ycoeffs),
      BL_FORT_FAB_ARG_3D(zcoeffs));

  void ca_fill_first_visc_coeff
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(xcoeffs),
     BL_FORT_FAB_ARG_3D(ycoeffs),
     BL_FORT_FAB_ARG_3D(zcoeffs));

  void ca_fill_secnd_visc_coeff
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(xcoeffs),
     BL_FORT_FAB_ARG_3D(ycoeffs),
     BL_FORT_FAB_ARG_3D(zcoeffs));

  void ca_compute_div_tau_u
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(div_tau_u),
     const BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, const int* coord_type);

  void ca_thermo_src
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state_old),
     const BL_FORT_FAB_ARG_3D(state_new),
     BL_FORT_FAB_ARG_3D(thermo_src),
     const amrex::Real* prob_lo, const amrex::Real* dx,
     const amrex::Real time, const amrex::Real dt);

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

#ifdef AMREX_DIMENSION_AGNOSTIC

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
#ifdef AMREX_DIMENSION_AGNOSTIC
BL_FORT_PROC_DECL(CA_INITRAD,ca_initrad)
  (const int& level, const amrex::Real& time,
   const int* lo, const int* hi,
   const int& num_state,
   BL_FORT_FAB_ARG_3D(state),
   const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* xhi);
#else
BL_FORT_PROC_DECL(CA_INITRAD,ca_initrad)
  (const int& level, const amrex::Real& time,
   const int* lo, const int* hi,
   const int& num_state,
   BL_FORT_FAB_ARG(state),
   const amrex::Real dx[], const amrex::Real xlo[], const amrex::Real xhi[]);
#endif
#endif

#ifdef AMREX_DIMENSION_AGNOSTIC
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
