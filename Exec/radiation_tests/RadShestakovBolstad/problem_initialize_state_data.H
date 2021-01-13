#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real Tcgs = problem::T_0 * 1.e3_rt * C::ev2erg / C::k_B;

    Real B0 = 8.0_rt * M_PI * C::hplanck / std::pow(C::c_light, 3);
    Real nu0 = C::k_B * Tcgs / C::hplanck;
    Real l0 = std::pow(nu0, 3) / problem::kappa_0;
    Real x0 = l0 / std::sqrt(3.0_rt);
    Real u0 = B0 * std::pow(nu0, 3);

    // rhoe_0 = R * k_B * Tcgs * u0 / hplanck
    Real rhoe_0 = 99968636.6828e0_rt * Tcgs * problem::rho_0;
    //                cv               * T        * rho

    Real xcen = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

    state(i,j,k,URHO) = problem::rho_0;
    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

    // set the composition to be all in the first species

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = 0.0_rt;
    }
    state(i,j,k,UFS) = state(i,j,k,URHO);

#if NAUX_NET > 0
    for (int n = 0; n < NumAux; n++) {
        state(i,j,k,UFX+n) = 0.0_rt;
    }
    state(i,j,k,UFX) = state(i,j,k,URHO);
#endif

    if (std::abs(xcen)/x0 < problem::x_jump) {
        state(i,j,k,UEDEN) = rhoe_0;
        state(i,j,k,UEINT) = rhoe_0;
        state(i,j,k,UTEMP) = Tcgs;
    } else {
        state(i,j,k,UEDEN) = rhoe_0 * 1.e-12_rt;
        state(i,j,k,UEINT) = rhoe_0 * 1.e-12_rt;
        state(i,j,k,UTEMP) = Tcgs * 1.e-12_rt;
    }
}
#endif
