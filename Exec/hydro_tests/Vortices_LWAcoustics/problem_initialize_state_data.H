#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

    const Real* problo = geomdata.ProbLo();
    const Real* probhi = geomdata.ProbHi();

    const Real* dx = geomdata.CellSize();

    // density

    Real rho = problem::rho_0;

    Real y = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1];
    Real x = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0];

    Real r_1 = std::sqrt((x - problem::x_c1) * (x - problem::x_c1) +
                         (y - problem::y_c1) * (y - problem::y_c1));
    Real r_2 = std::sqrt((x - problem::x_c2) * (x - problem::x_c2) +
                         (y - problem::y_c2) * (y - problem::y_c2));

    Real vel_theta_1 = problem::circ * r_1 /
        (2.0_rt * M_PI * (problem::r_c * problem::r_c + r_1 * r_1));
    Real vel_theta_2 = problem::circ * r_2 /
        (2.0_rt * M_PI * (problem::r_c * problem::r_c + r_2 * r_2));

    Real sin_theta_1 = (y - problem::y_c1) / r_1;
    Real cos_theta_1 = (x - problem::x_c1) / r_1;

    Real sin_theta_2 = (y - problem::y_c2) / r_2;
    Real cos_theta_2 = (x - problem::x_c2) / r_2;

    Real u = vel_theta_1 * sin_theta_1 + vel_theta_2 * sin_theta_2;
    Real v = -vel_theta_1 * cos_theta_1 - vel_theta_2 * cos_theta_2;

    // single species for all zones

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = 0.0_rt;
    }
    state(i,j,k,UFS) = 1.0_rt;

    // momentum field

    state(i,j,k,UMX) = rho * u;
    state(i,j,k,UMY) = rho * v;
    state(i,j,k,UMZ) = 0.0_rt;

    // density

    state(i,j,k,URHO) = rho;

    // internal energy

    state(i,j,k,UEINT) = problem::p_ref / (eos_rp::eos_gamma - 1.0_rt);

    // Total energy

    state(i,j,k,UEDEN) = state(i,j,k,UEINT) + 0.5_rt *
        (state(i,j,k,UMX) * state(i,j,k,UMX) +
         state(i,j,k,UMY) * state(i,j,k,UMY) +
         state(i,j,k,UMZ) * state(i,j,k,UMZ)) / state(i,j,k,URHO);

    // Convert mass fractions to conserved quantity

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) *= rho;
    }
}
#endif
