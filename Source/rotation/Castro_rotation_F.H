#ifndef CASTRO_ROTATION_F_H
#define CASTRO_ROTATION_F_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void set_rot_period(const amrex::Real* period);

  void get_omega(amrex::Real* omega);

#ifdef __cplusplus
}
#endif

#endif
