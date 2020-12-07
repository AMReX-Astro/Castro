#ifndef problem_initialize_rad_data_H
#define problem_initialize_rad_data_H

#include <fundamental_constants.H>
#include <blackbody.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_rad_data (int i, int j, int k,
                                  Array4<Real> const& rad,
                                  GpuArray<Real, NGROUPS+1> const& xnu,
                                  const GeometryData& geomdata)
{

    constexpr Real Tindex = 0.5;
    constexpr int nsub = 16;

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real dx_sub = dx[0] / static_cast<Real>(nsub);

    Real rlast = model::profile(0).r(model::npts-1);
    Real rbase = rlast * problem::rbasefac;

    Real Twind0 = problem::Twind1 * std::pow(problem::rwind1 / problem::rwind0, Tindex);
    Real Tlast = model::profile(0).state(model::npts-1, model::itemp);

    Real xcl = problo[0] + dx[0] * static_cast<Real>(i);

    Real vtot = 0.0_rt;
    Real T = 0.0_rt;

    for (int ii = 0; ii < nsub; ii++) {
        Real xx = xcl + (static_cast<Real>(ii) + 0.5_rt) * dx_sub;
        Real vsub = xx * xx;
        vtot += vsub;

        if (xx >= model::profile(0).r(model::npts-1)) {
            Real T_tmp;
            if (xx >= problem::rwind0) {
                T_tmp = problem::Twind1 * std::pow(problem::rwind1 / xx, Tindex);
            } else {
                T_tmp = Tlast * std::pow(Twind0 / Tlast,
                                         ((std::log(xx - rbase) - std::log(rlast - rbase)) /
                                          (std::log(problem::rwind0 - rbase) - std::log(rlast - rbase))));
            }
            T += vsub * T_tmp;

        } else if (xx <= model::profile(0).r(0)) {
            T += model::profile(0).state(0, model::itemp) * vsub;

        } else {
            T += interpolate(xx, model::itemp) * vsub;
        }
    }

    T /= vtot;

    // set radiation energy density to a T**4

    if (NGROUPS == 1) {
        rad(i,j,k,0) = C::a_rad * std::pow(T, 4);

    } else {
        for (int g = 0; g < NGROUPS; g++) {
            rad(i,j,k,g) = BGroup(T, xnu[g], xnu[g+1]);
        }
    }
}
#endif
