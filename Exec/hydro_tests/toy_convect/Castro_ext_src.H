#ifndef CASTRO_EXT_SRC_H
#define CASTRO_EXT_SRC_H

#include <prob_parameters.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void do_ext_src(int i, int j, int k,
                GeometryData const& geomdata, 
                Array4<const Real> const& state,
                Array4<Real> const& src,
                GpuArray<Real, 3> const& center,
                const Real dt, const Real time)
{

    for (int n = 0; n < NSRC; ++n) {
        src(i,j,k,n) = 0.0_rt;
    }

    auto rho = state(i,j,k,URHO);

    auto T6 = state(i,j,k,UTEMP) / 1.0e6_rt;
    auto T613 = std::pow(T6, 1.0_rt/3.0_rt);

    // CNO abundance
    auto X_CNO = (state(i,j,k,UFS-1+Species::C12) +
                    state(i,j,k,UFS-1+Species::N14) +
                    state(i,j,k,UFS-1+Species::O16)) / rho;

    // H abundance
    auto X_1 = state(i,j,k,UFS-1+Species::H1) / rho;

    // CNO heating from Kippenhahn & Weigert, Eq. 18.65
    auto g14 = 1.0_rt + 2.7e-3_rt*T613 - 7.78e-3_rt*T613*T613 - 1.49e-4_rt*T6;
    auto eps_CNO = 8.67e27_rt * g14 * X_CNO * X_1 * rho * std::exp(-152.28_rt/T613) / (T613*T613);

    // source terms
    src(i,j,k,UEDEN) = rho * eps_CNO;
    src(i,j,k,UEINT) = rho * eps_CNO;
}

#endif
