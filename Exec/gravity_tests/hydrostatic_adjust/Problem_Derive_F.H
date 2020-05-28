#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

void ca_derpi(BL_FORT_FAB_ARG_3D(der), const int* nvar,
              const BL_FORT_FAB_ARG_3D(data), const int* ncomp, const int* lo,
              const int* hi, const int* domain_lo, const int* domain_hi,
              const amrex::Real* delta, const amrex::Real* xlo,
              const amrex::Real* time, const amrex::Real* dt, const int* bcrec,
              const int* level, const int* grid_no);

void ca_derpioverp0(BL_FORT_FAB_ARG_3D(der), const int* nvar,
                    const BL_FORT_FAB_ARG_3D(data), const int* ncomp,
                    const int* lo, const int* hi, const int* domain_lo,
                    const int* domain_hi, const amrex::Real* delta,
                    const amrex::Real* xlo, const amrex::Real* time,
                    const amrex::Real* dt, const int* bcrec, const int* level,
                    const int* grid_no);

#ifdef __cplusplus
}
#endif
