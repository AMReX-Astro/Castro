#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{
    int coord_type = geomdata.Coord();

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

    Real y = 0.0;
#if AMREX_SPACEDIM >= 2
    y = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
#endif

    Real z = 0.0;
#if AMREX_SPACEDIM == 3
    z = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
#endif


    Real pert = 1.0e-5_rt * std::sin(2.0_rt * M_PI * (problem::k_x * x + problem::k_y * y));

    // initialize using MM eq. 55
    state(i,j,k,URHO) = problem::rho_0;
    state(i,j,k,UMX) = (problem::u_x + (std::sqrt(2.0_rt) * problem::cw * problem::cw *
                                        problem::u_by - problem::c_s * problem::u_ky) *
                        pert / problem::cw) * problem::rho_0;
    state(i,j,k,UMY) = (problem::u_y + (problem::c_s * problem::u_kx - std::sqrt(2.0_rt) *
                                        problem::cw * problem::cw * problem::u_bx) *
                        pert / problem::cw) * problem::rho_0;
    state(i,j,k,UMZ) = 0.0_rt;

    state(i,j,k,UFS) = state(i,j,k,URHO);
    for (int n = 1; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = 0.0_rt;
    }

    Real pressure = problem::p_0 + problem::rho_0 * problem::c_s * pert;

    // compute the internal energy (erg/cc)

    eos_t eos_state;
    eos_state.rho = problem::rho_0;
    eos_state.p = pressure;
    eos_state.T = 100.0_rt;  // initial guess
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state(i,j,k,UFS+n) / state(i,j,k,URHO);
    }

    eos(eos_input_rp, eos_state);

    Real rhoe_0 = problem::rho_0 * eos_state.e;
    Real T_0 = eos_state.T;

    Real u = state(i,j,k,UMX) / problem::rho_0;
    Real v = state(i,j,k,UMY) / problem::rho_0;
    Real w = state(i,j,k,UMZ) / problem::rho_0;

    state(i,j,k,UEDEN) = rhoe_0 + 0.5_rt * problem::rho_0 * (u * u + v * v +w * w);
    state(i,j,k,UEINT) = rhoe_0;
    state(i,j,k,UTEMP) = T_0;

}
#endif
