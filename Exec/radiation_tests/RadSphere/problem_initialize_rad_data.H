#ifndef problem_initialize_rad_data_H
#define problem_initialize_rad_data_H

#include <blackbody.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_rad_data (int i, int j, int k,
                                  Array4<Real> const& rad,
                                  GpuArray<Real, NGROUPS+1> const& xnu,
                                  GpuArray<Real, NGROUPS> const& nugroup,
                                  GpuArray<Real, NGROUPS> const& dnugroup,
                                  const GeometryData& geomdata)
{

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real xcen = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

    for (int igroup = 0; igroup < NGROUPS; igroup++) {

        Real xx = C::hplanck * nugroup[igroup] / (C::k_B * problem::T_0);
        if (xx > 708.e0_rt) {
            // exp(708+eps) will give 1.e308 -- the limits of IEEE floating point,
            // and generate an overflow exception
            rad(i,j,k,igroup) = 0.0_rt;
        } else {
            rad(i,j,k,igroup) =
                (8.0_rt * M_PI * C::hplanck * std::pow(nugroup[igroup], 3) / std::pow(C::c_light, 3) ) /
                (std::exp(xx) - 1.0_rt) * dnugroup[igroup];
        }
    }
}
#endif

