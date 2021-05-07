#ifndef CASTRO_RADDERIVE_H
#define CASTRO_RADDERIVE_H

#include <AMReX_BLFort.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void ca_derertot(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                   const FArrayBox& datfab, const Geometry& /*geomdata*/,
                   Real /*time*/, const int* /*bcrec*/, int /*level*/);


#ifdef __cplusplus
}
#endif

#endif
