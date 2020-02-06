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

  void ca_get_spec_names(int* spec_names, const int* ispec, int* len);

  void ca_get_spec_az(const int* ispec, amrex::Real* A, amrex::Real* Z);

  void ca_get_aux_names(int* aux_names,  const int* iaux, int* len);

#ifdef REACTIONS
  void ca_get_abort_on_failure(int * abort_on_failure);

  void ca_set_abort_on_failure(int * abort_on_failure);
#endif

  void ca_amrinfo_init();
  void ca_amrinfo_finalize();

  void ca_set_amr_info(const int& level, const int& iteration, const int& ncycle,
		       const amrex::Real& time, const amrex::Real& dt);

  void ca_get_method_params(int* HYP_GROW);

  void runtime_pretty_print(int* jobinfo_file_name, const int* jobinfo_file_length);
#ifdef PROB_PARAMS
  void prob_params_pretty_print(int* jobinfo_file_name, const int* jobinfo_file_length);
#endif

#ifdef GRAVITY
#if (BL_SPACEDIM > 1)
  void allocate_outflow_data(const int* n1d, const int* nc);
  void set_old_outflow_data(const amrex::Real* radial, const amrex::Real* time,
			    const int* numpts, const int* nc);
  void set_new_outflow_data(const amrex::Real* radial, const amrex::Real* time,
			    const int* numpts, const int* nc);
  void swap_outflow_data();
#endif
#endif

  void ca_set_method_params
    (const int& dm);

  void ca_finalize_meth_params();

  void ca_set_castro_method_params();

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
  void ca_get_ambient_params(const int* name, const int* namlen);

#ifdef GRAVITY
#if (BL_SPACEDIM > 1)
  void ca_compute_avgstate
    (const int* lo, const int* hi,
     const amrex::Real* dx, const amrex::Real* dr, const int* nc,
     const BL_FORT_FAB_ARG_3D(S  ),  amrex::Real* avgden,
     const BL_FORT_FAB_ARG_3D(Vol),  amrex::Real* avgvol,
     const amrex::Real* problo, const int* numpts_1d);
#endif
#endif

#ifdef GPU_COMPATIBLE_PROBLEM
  void ca_initdata(const int* lo, const int* hi,
                   BL_FORT_FAB_ARG_3D(state),
                   const amrex::Real* dx, const amrex::Real* problo);

  void ca_initrad(const int* lo, const int* hi,
                  BL_FORT_FAB_ARG_3D(rad_state),
                  const amrex::Real* dx, const amrex::Real* problo);
#endif

  void ca_apply_metric
    (const int* lo, const int* hi,
     const int* xlo, const int* xhi,
#if AMREX_SPACEDIM >= 2
     const int* ylo, const int* yhi,
#endif
     BL_FORT_FAB_ARG_3D(cc),
     const BL_FORT_FAB_ARG_3D(xedge),
#if AMREX_SPACEDIM >= 2
     const BL_FORT_FAB_ARG_3D(yedge),
#endif
     const amrex::Real* dx, const int coord_type);

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
  void ca_linear_to_hybrid_momentum
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state));

  void ca_hybrid_to_linear_momentum
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state));

  void ca_hybrid_hydro_source
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     BL_FORT_FAB_ARG_3D(sources),
     const amrex::Real mult_factor);
#endif

  void ca_check_initial_species
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(state));

  void ca_enforce_minimum_density
    (const int* lo, const int* hi,
     const amrex::Real* state, const int* s_lo, const int* s_hi,
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

  void ca_recompute_energetics
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(state));

  void ca_reset_internal_e
    (const int* lo, const int* hi, BL_FORT_FAB_ARG_3D(S_new),
     const int verbose);

  void ca_estdt
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state),
     const amrex::Real* dx, amrex::Real* dt);

#ifdef DIFFUSION
  void ca_estdt_temp_diffusion
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
     const BL_FORT_FAB_ARG_3D(state_new),
     const BL_FORT_FAB_ARG_3D(reactions_new),
     const amrex::Real* dx, amrex::Real* dt);
#endif

  void ca_check_timestep
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(state_old),
     const BL_FORT_FAB_ARG_3D(state_new),
#ifdef REACTIONS
     const BL_FORT_FAB_ARG_3D(reactions_old),
     const BL_FORT_FAB_ARG_3D(reactions_new),
#endif
     const amrex::Real* dx, const amrex::Real dt_old, amrex::Real* dt);

  void ca_make_fourth_average(const int lo[], const int hi[],
                              BL_FORT_FAB_ARG_3D(q), const int* nc,
                              const BL_FORT_FAB_ARG_3D(q_bar), const int* nc_bar,
                              const int* domlo, const int* domhi);

  void ca_make_fourth_in_place(const int lo[], const int hi[],
                               BL_FORT_FAB_ARG_3D(q), const int* nc,
                               const int* domlo, const int* domhi);

  void ca_make_fourth_in_place_n(const int lo[], const int hi[],
                                 BL_FORT_FAB_ARG_3D(q), const int* nc, const int* ncomp,
                                 const int* domlo, const int* domhi);

  void ca_compute_lap_term(const int lo[], const int hi[],
                           const BL_FORT_FAB_ARG_3D(U), const int* nc,
                           BL_FORT_FAB_ARG_ANYD(lap), const int* ncomp,
                           const int* domlo, const int* domhi);

  void ca_make_cell_center(const int lo[], const int hi[],
                           const BL_FORT_FAB_ARG_3D(U), const int* nc,
                           BL_FORT_FAB_ARG_3D(U_cc), const int* nc_cc,
                           const int* domlo, const int* domhi);

  void ca_make_cell_center_in_place(const int lo[], const int hi[],
                                    BL_FORT_FAB_ARG_3D(q), const int* nc,
                                    const int* domlo, const int* domhi);

#ifdef RADIATION
  void ca_inelastic_sct
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(stateout),
     BL_FORT_FAB_ARG_3D(Erout),
     const BL_FORT_FAB_ARG_3D(kappa_s),
     const amrex::Real& dt);
#endif



  void flush_output();


  void ca_summass
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(rho),
     const amrex::Real* dx, const BL_FORT_FAB_ARG_3D(vol), amrex::Real* mass);

  void ca_sumsquared
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(rho),
     const amrex::Real* dx, const BL_FORT_FAB_ARG_3D(vol), amrex::Real* mass);

  void ca_sumlocmass
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(rho),
     const amrex::Real* dx, const BL_FORT_FAB_ARG_3D(vol), amrex::Real* mass, const int idir);

  void ca_sumlocmass2d
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(rho),
     const amrex::Real* dx, const BL_FORT_FAB_ARG_3D(vol), amrex::Real* mass, const int idir1, const int idir2);

  void ca_sumlocsquaredmass
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(rho),
     const amrex::Real* dx, const BL_FORT_FAB_ARG_3D(vol), amrex::Real* mass, const int idir);

  void ca_sumproduct
    (const int* lo, const int* hi, const BL_FORT_FAB_ARG_3D(f1), const BL_FORT_FAB_ARG_3D(f2),
     const amrex::Real* dx, const BL_FORT_FAB_ARG_3D(vol), amrex::Real* s);


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


/* hydrodynamics */
#include <Castro_hydro_F.H>


/* reactions */
#ifdef REACTIONS
#include <Castro_react_F.H>
#endif


/* sources */
#include <Castro_sources_F.H>


/* gravity */
#ifdef GRAVITY
#include <Castro_gravity_F.H>
#endif

/* rotation */
#ifdef ROTATION
#include <Castro_rotation_F.H>
#endif


/* SDC */
#ifdef TRUE_SDC
#include <Castro_sdc_F.H>
#endif


/* diffusion */
#ifdef DIFFUSION
#include <Castro_diffusion_F.H>
#endif


/* SCF */
#include <scf_relax_F.H>

#ifdef __cplusplus
}
#endif


#ifndef GPU_COMPATIBLE_PROBLEM
BL_FORT_PROC_DECL(CA_INITDATA,ca_initdata)
    (const int& level, const amrex::Real& time,
     const int* lo, const int* hi,
     const int& num_state,
     BL_FORT_FAB_ARG_3D(state),
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
#ifndef GPU_COMPATIBLE_PROBLEM
BL_FORT_PROC_DECL(CA_INITRAD,ca_initrad)
  (const int& level, const amrex::Real& time,
   const int* lo, const int* hi,
   const int& num_state,
   BL_FORT_FAB_ARG_3D(state),
   const amrex::Real* dx, const amrex::Real* xlo, const amrex::Real* xhi);
#endif
#endif


/* problem-specific includes */
#include <Problem_F.H>

#endif
