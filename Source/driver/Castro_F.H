#ifndef CASTRO_F_H
#define CASTRO_F_H

#include <AMReX_BLFort.H>

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

#endif
