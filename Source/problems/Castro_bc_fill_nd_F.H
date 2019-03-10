#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void hypfill
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real time, const int* bc);

  void denfill
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real time, const int* bc);

#ifdef GRAVITY
  void ca_phigravfill
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real time, const int* bc);

  void ca_gravxfill
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real time, const int* bc);

  void ca_gravyfill
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real time, const int* bc);

  void ca_gravzfill
    (const int* lo, const int* hi
     BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real time, const int* bc);
#endif

#ifdef ROTATION
  void ca_phirotfill
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real time, const int* bc);

  void ca_rotxfill
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real time, const int* bc);

  void ca_rotyfill
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real time, const int* bc);

  void ca_rotzfill
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real time, const int* bc);
#endif

#ifdef REACTIONS
  void ca_reactfill
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real time, const int* bc);
#endif

#ifdef RADIATION
  void ca_radfill
    (const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const amrex::Real* dx, const amrex::Real* glo,
     const amrex::Real time, const int* bc);
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

#ifdef __cplusplus
}
#endif
