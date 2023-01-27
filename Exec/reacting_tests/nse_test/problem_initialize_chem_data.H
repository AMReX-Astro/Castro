#ifndef problem_initialize_chem_data_H
#define problem_initialize_chem_data_H

#include <prob_parameters.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_chem_data (int i, int j, int k,
				   Array4<Real> const& state,
				   const GeometryData& geomdata)
{
  state(i,j,k,CMUP) = -3.0_rt;
  state(i,j,k,CMUN) = -12.0_rt;
}

#endif
