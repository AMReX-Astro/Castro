#ifndef _Castro_error_F_H_
#define _Castro_error_F_H_
#include <BLFort.H>

BL_FORT_PROC_DECL(CA_LAPLAC_ERROR,ca_laplac_error)
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(var),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CA_DENERROR,ca_denerror)
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(den),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CA_VELERROR,ca_velerror)
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(vel),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CA_VORTERROR,ca_vorterror)
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(vel),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CA_TEMPERROR,ca_temperror)
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(temp),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CA_PRESSERROR,ca_presserror)
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(press),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);

#ifdef PARTICLES
BL_FORT_PROC_DECL(CA_PART_CNT_ERR,ca_part_cnt_err)
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(ls),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);
#endif

#ifdef LEVELSET
BL_FORT_PROC_DECL(CA_LSERROR,ca_lserror)
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(ls),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);
#endif

BL_FORT_PROC_DECL(CA_MFRACERROR,ca_mfracerror)
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(press),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);

#ifdef NEUTRINO
BL_FORT_PROC_DECL(CA_YEERROR,ca_yeerror)
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(Ye),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);
#endif

#ifdef GRAVITY
BL_FORT_PROC_DECL(CA_GRAVERROR,ca_graverror)
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(grav),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);
#endif

BL_FORT_PROC_DECL(CA_RADERROR,ca_raderror)
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(rad),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CA_ENTERROR,ca_enterror)
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(den),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);

BL_FORT_PROC_DECL(CA_STATE_ERROR,ca_state_error)
    (int* tag, const int* tag_lo, const int* tag_hi,
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG_3D(state),
     const int* lo, const int* hi,
     const int* ncomp,
     const int* domlo, const int* domhi,
     const Real* dx, const Real* xlo, const Real* problo,
     const Real* time, const int* level);

#ifdef DIMENSION_AGNOSTIC
BL_FORT_PROC_DECL(SET_PROBLEM_TAGS,set_problem_tags)
    (int* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* tagval, const int* clearval,
     const int* lo, const int* hi,
     const Real* dx, const Real* problo,
     const Real* time, const int* level);
#else
BL_FORT_PROC_DECL(SET_PROBLEM_TAGS,set_problem_tags)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     BL_FORT_FAB_ARG(state),
     const int* tagval, const int* clearval,
     const int lo[], const int hi[],
     const Real dx[], const Real problo[],
     const Real* time, const int* level);
#endif
#endif
