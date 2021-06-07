#ifndef CASTRO_F_H
#define CASTRO_F_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void ca_microphysics_init(const Real small_dens_in, const Real small_temp_in);

  void ca_extern_init(const int* name, const int* namlen);

#ifdef REACTIONS
  void ca_set_abort_on_failure(int * abort_on_failure);
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

#ifdef __cplusplus
}
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
