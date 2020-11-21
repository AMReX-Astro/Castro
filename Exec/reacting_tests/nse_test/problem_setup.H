#ifndef problem_setup_H
#define problem_setup_H

#include <prob_parameters.H>
#include <eos.H>
#include <network.H>

AMREX_INLINE
void problem_initialize ()
{
    const Geometry& dgeom = DefaultGeometry();

    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
        problem::center[idir] = 0.5_rt * (problo[idir] + probhi[idir]);
    }

    for (int idir = AMREX_SPACEDIM; idir < 3; idir++) {
        problem::center[idir] = 0.0_rt;
    }

    // we only work in NSE mode, so the mass fractions are not used by
    // the EOS

#ifndef NSE_THERMO
    amrex::Error("Error: this problem requires USE_NSE=TRUE");
#endif

}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{
    int coord_type = geomdata.Coord();

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();
    const Real* probhi = geomdata.ProbHi();

    Real xx = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt) - problem::center[0];

    Real yy = 0.0;
#if AMREX_SPACEDIM >= 2
    yy = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt) - problem::center[1];
#endif

    Real zz = 0.0;
#if AMREX_SPACEDIM == 3
    zz = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt) - problem::center[2];
#endif

    constexpr Real ye0 = 0.5_rt;
    constexpr Real dye = -0.05_rt;

    Real dist = std::sqrt(xx * xx + yy * yy + zz * zz);

    Real T;
    Real ye;

    if (dist <= problem::center[0]) {
        T = problem::T0 *
            (1.0_rt + problem::dT_fact *
             std::exp(-(dist/problem::L_pert) * (dist/problem::L_pert)) *
             std::pow(std::cos(M_PI*(dist/(probhi[0]-problo[0]))), 6));
             ye = ye0 *
             (1.0_rt + dye * std::exp(-(dist/problem::L_pert) * (dist/problem::L_pert) ) *
              std::pow(std::cos(M_PI*(dist/(probhi[0]-problo[0]))), 6));
    } else {
        T = problem::T0;
        ye = ye0;
    }

    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

    // we are isentropic, so find rho

    eos_t eos_state;
    eos_state.T = T;
    eos_state.rho = problem::rho0;

    Real xn[NumSpec];
    Real abar;
    Real dq;
    Real dyedt;

    nse_interp(T, problem::rho0, ye, abar, dq, dyedt, xn);

    // since the species are interpolated, normalize them

    Real sumX = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
        sumX += xn[n];
    }
    for (int n = 0; n < NumSpec; n++) {
        xn[n] /= sumX;
    }

    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = xn[n];
    }
    eos_state.aux[AuxZero::iye] = ye;
    eos_state.aux[AuxZero::iabar] = abar;
    eos_state.aux[AuxZero::ibea] = dq;

    eos(eos_input_rt, eos_state);

    state(i,j,k,URHO) = eos_state.rho;

    state(i,j,k,UEDEN) = eos_state.rho * eos_state.e;
    state(i,j,k,UEINT) = eos_state.rho * eos_state.e;

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * xn[n];
    }

    for (int n = 0; n < NumAux; n++) {
        state(i,j,k,UFX+n) = state(i,j,k,URHO) * eos_state.aux[n];
    }

    state(i,j,k,UTEMP) = eos_state.T;

}
#endif
