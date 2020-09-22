#ifndef CASTRO_HYDRO_F_H
#define CASTRO_HYDRO_F_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

#ifdef RADIATION
  void update_one_species(const int n, Real* u, Real* a, const Real* dx, const Real tend, int& nstepmax);
#endif

#ifdef __cplusplus
}
#endif

#endif
