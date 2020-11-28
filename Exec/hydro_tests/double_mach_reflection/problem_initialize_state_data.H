#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>
#include <model_parser.H>
#include <initial_model.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

    constexpr Real pi_over_3 = M_PI / 3.0_rt;
    constexpr Real ff = 0.25_rt;

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real xcen = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

    Real ycen = 0.0;
#if AMREX_SPACEDIM >= 2
    ycen = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
#endif

    Real zcen = 0.0;
#if AMREX_SPACEDIM == 3
    zcen = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
#endif

    state(i,j,k,URHO) = 0.0_rt;
    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;
    state(i,j,k,UEDEN) = 0.0_rt;
    state(i,j,k,UEINT) = 0.0_rt;

    for (int jj = -1; jj <= 1; ++jj) {
        if (jj == 0) {
            continue;
        }
        Real y = ycen + 0.5_rt * dx[1] * (static_cast<Real>(jj) / std::sqrt(3.0_rt));

        for (int ii = -1; ii <= 1; ++ii) {
            if (ii == 0) {
                continue;
            }
            Real x = xcen + 0.5_rt * dx[0] * (static_cast<Real>(ii) / std::sqrt(3.0_rt));

            Real shockfront = std::tan(pi_over_3) * (x - 1.0_rt/6.0_rt); // initial shock front

            if (y >= shockfront ) {
                state(i,j,k,URHO) += ff * problem::rho_l;
                state(i,j,k,UMX) += ff * problem::rho_l * problem::u_l;
                state(i,j,k,UMY) += ff * problem::rho_l * problem::v_l;

                state(i,j,k,UEDEN) += ff *
                    (problem::rhoe_l +
                     0.5_rt * problem::rho_l * (problem::u_l * problem::u_l +
                                                problem::v_l * problem::v_l));
                state(i,j,k,UEINT) += ff * problem::rhoe_l;
                state(i,j,k,UTEMP) += ff * problem::T_l;
            } else {
                state(i,j,k,URHO) += ff * problem::rho_r;
                state(i,j,k,UMX) += ff * problem::rho_r * problem::u_r;
                state(i,j,k,UMY) += ff * problem::rho_r * problem::v_r;

                state(i,j,k,UEDEN) += ff *
                    (problem::rhoe_r +
                     0.5_rt * problem::rho_r * (problem::u_r * problem::u_r +
                                                problem::v_r * problem::v_r));
                state(i,j,k,UEINT) += ff * problem::rhoe_r;
                state(i,j,k,UTEMP) += ff * problem::T_r;
            }
        }
    }

    state(i,j,k,UFS) = state(i,j,k,URHO);
    for (int n = 1; n < NumSpec; ++n) {
        state(i,j,k,UFS+n) = 0.0_rt;
    }
}
#endif
