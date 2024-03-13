#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>
#include <Castro_util.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k, Array4<Real> const& state, const GeometryData& geomdata)
{
    Real volinv = 1.0_rt;
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        volinv /= static_cast<Real>(problem::nsub);
    }

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real dx_sub = dx[0] / static_cast<Real>(problem::nsub);

#if AMREX_SPACEDIM >= 2
    Real dy_sub = dx[1] / static_cast<Real>(problem::nsub);
#else
    Real dy_sub = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
    Real dz_sub = dx[2] / static_cast<Real>(problem::nsub);
#else
    Real dz_sub = 0.0_rt;
#endif

    Real xl = problo[0] + static_cast<Real>(i) * dx[0];

#if AMREX_SPACEDIM >= 2
    Real yl = problo[1] + static_cast<Real>(j) * dx[1];
#else
    Real yl = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
    Real zl = problo[2] + static_cast<Real>(k) * dx[2];
#else
    Real zl = 0.0_rt;
#endif

    Real avg_rho = 0.0_rt;

    for (int kk = 0; kk <= dg2 * (problem::nsub - 1); ++kk) {
        Real zz = zl + (static_cast<Real>(kk) + 0.5_rt) * dz_sub;

        for (int jj = 0; jj <= dg1 * (problem::nsub - 1); ++jj) {
            Real yy = yl + (static_cast<Real>(jj) + 0.5_rt) * dy_sub;

            for (int ii = 0; ii <= (problem::nsub - 1); ++ii) {
                Real xx = xl + (static_cast<Real>(ii) + 0.5_rt) * dx_sub;

                Real dist = std::sqrt((xx - problem::center[0]) * (xx - problem::center[0]) +
                                      (yy - problem::center[1]) * (yy - problem::center[1]) +
                                      (zz - problem::center[2]) * (zz - problem::center[2]));

                // use a tanh profile to smooth the transition between rho_0
                // and rho_ambient
                Real rho_n = problem::rho_0 - 0.5_rt * (problem::rho_0 - problem::rho_ambient) *
                                              (1.0_rt + std::tanh((dist - problem::r_0) /
                                                                  problem::smooth_delta));

                // allow for the center to be empty
                if (problem::r_offset > 0.0_rt) {
                    rho_n -= 0.5_rt * (rho_n - problem::rho_ambient) *
                             (1.0_rt + std::tanh((problem::r_offset - dist) /
                                                 problem::offset_smooth_delta));
                }

                avg_rho += rho_n;
            }
        }
    }

    state(i,j,k,URHO) = avg_rho * volinv;

    eos_t eos_state;

    eos_state.rho = state(i,j,k,URHO);
    eos_state.p   = problem::p_0;
    eos_state.T   = castro::small_temp; // Initial guess for the EOS
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = problem::X_0[n];
    }

    eos(eos_input_rp, eos_state);

    Real temp = eos_state.T;
    Real eint = eos_state.e;

    state(i,j,k,UTEMP) = temp;
    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;
    state(i,j,k,UEDEN) = state(i,j,k,URHO) * eint;
    state(i,j,k,UEINT) = state(i,j,k,URHO) * eint;
    for (int n = 0; n < NumSpec; ++n) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * problem::X_0[n];
    }
}
#endif
