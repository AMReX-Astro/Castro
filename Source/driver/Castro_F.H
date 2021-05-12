#ifndef CASTRO_F_H
#define CASTRO_F_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void ca_microphysics_init();

  void ca_extern_init(const int* name, const int* namlen);

#ifdef REACTIONS
  void ca_set_abort_on_failure(int * abort_on_failure);
#endif

  void prob_params_pretty_print(int* jobinfo_file_name, const int* jobinfo_file_length);

  void ca_set_method_params
    (const int& dm);

  void ca_finalize_meth_params();

  void ca_set_castro_method_params();

  void ca_set_problem_params
    (const int& dm,
     const int& coord_type,
     const amrex::Real* problo, const amrex::Real* probhi);

  void ca_get_tagging_params(const int* name, const int* namlen);
#ifdef SPONGE
  void ca_read_sponge_params(const int* name, const int* namlen);

  void ca_get_sponge_params(amrex::Real& sponge_lower_factor,
                            amrex::Real& sponge_upper_factor,
                            amrex::Real& sponge_lower_radius,
                            amrex::Real& sponge_upper_radius,
                            amrex::Real& sponge_lower_density,
                            amrex::Real& sponge_upper_density,
                            amrex::Real& sponge_lower_pressure,
                            amrex::Real& sponge_upper_pressure,
                            amrex::Real* sponge_target_velocity,
                            amrex::Real& sponge_timescale);

  void ca_allocate_sponge_params();

  void ca_deallocate_sponge_params();
#endif
  void ca_get_ambient_params(const int* name, const int* namlen);

  void get_ambient_data(Real* ambient_state);

#ifdef GPU_COMPATIBLE_PROBLEM
  void ca_initdata(const int* lo, const int* hi,
                   BL_FORT_FAB_ARG_3D(state),
                   const amrex::Real* dx, const amrex::Real* problo);

  void ca_initrad(const int* lo, const int* hi,
                  BL_FORT_FAB_ARG_3D(rad_state),
                  const amrex::Real* dx, const amrex::Real* problo);
#endif


#ifdef RADIATION
  void ca_inelastic_sct
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(stateout),
     BL_FORT_FAB_ARG_3D(Erout),
     const BL_FORT_FAB_ARG_3D(kappa_s),
     const amrex::Real& dt);
#endif

#ifdef AUX_UPDATE
  void ca_auxupdate
    (BL_FORT_FAB_ARG(state_old),
     BL_FORT_FAB_ARG(state_new),
     const int* lo, const int* hi,
     const amrex::Real * dt);
#endif




/* sources */
#include <Castro_sources_F.H>



/* SCF */
#ifdef GRAVITY
#include <scf_relax_F.H>
#endif

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

#ifdef MHD
BL_FORT_PROC_DECL(CA_INITMAG,ca_initmag)
  (const int& level, const amrex::Real& time,
   const int* lo, const int* hi,
   const int& nx, BL_FORT_FAB_ARG_3D(magx),
   const int& ny, BL_FORT_FAB_ARG_3D(magy),
   const int& nz, BL_FORT_FAB_ARG_3D(magz),
   const amrex::Real dx[], const amrex::Real xlo[], const amrex::Real xhi[]);
#endif

#endif
