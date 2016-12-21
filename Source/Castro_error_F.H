#ifndef _Castro_error_F_H_
#define _Castro_error_F_H_
#include <AMReX_BLFort.H>

extern "C"
{

  void ca_laplac_error
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(var),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);

  void ca_denerror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(den),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);

  void ca_velerror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(vel),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);

  void ca_vorterror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(vel),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);

  void ca_temperror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(temp),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);

  void ca_presserror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(press),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);

  void ca_nucerror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(state),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);

  void ca_mfracerror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(press),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);

#ifdef NEUTRINO
  void ca_yeerror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(Ye),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);
#endif

#ifdef SELF_GRAVITY
  void ca_graverror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(grav),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);
#endif

#ifdef RADIATION
  void ca_raderror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(rad),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);
#endif
  
  void ca_enterror
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(den),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);

  void ca_state_error
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(state),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);

#ifdef DIMENSION_AGNOSTIC
  void set_problem_tags
    (int* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* tagval, const int* clearval,
     const int* lo, const int* hi,
     const Real* dx, const Real* problo,
     const Real* time, const int* level);
#else
  void set_problem_tags
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     BL_FORT_FAB_ARG(state),
     const int* tagval, const int* clearval,
     const int lo[], const int hi[],
     const Real dx[], const Real problo[],
     const Real* time, const int* level);
#endif

}

#endif
