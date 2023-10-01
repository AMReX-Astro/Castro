#ifndef CASTRO_HABEC_F_H
#define CASTRO_HABEC_F_H

#include <AMReX_ArrayLim.H>
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#else
#define RadBoundCond int
#endif

  void hdterm(BL_FORT_FAB_ARG(dterm),
              BL_FORT_FAB_ARG(soln),
              ARLIM_P(reglo), ARLIM_P(reghi),
              const int& cdir, const RadBoundCond& bct,
              const amrex::Real& bcl,
              const BL_FORT_FAB_ARG(bcval),
              const BL_FORT_IFAB_ARG(mask),
              BL_FORT_FAB_ARG(dcoefs),
              const amrex::Real* dx);

  void hdterm3(BL_FORT_FAB_ARG(dterm),
               BL_FORT_FAB_ARG(soln),
               ARLIM_P(reglo), ARLIM_P(reghi),
               const int& cdir, const int& bctype, const int* tf,
               const amrex::Real& bcl,
               const BL_FORT_FAB_ARG(bcval),
               const BL_FORT_IFAB_ARG(mask),
               BL_FORT_FAB_ARG(dcoefs),
               const amrex::Real* dx);

#ifdef __cplusplus
};
#endif

#endif
