// This sets up the Gresho vortex problem as described in
// Miczek, Roeple, and Edelmann 2015
//
// By choosing the reference pressure, p0, we can specify the
// Mach number

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

    Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);
    Real xl = problo[0] + dx[0] * (static_cast<Real>(i));

    Real y = 0.0;
    Real yl = 0.0;

#if AMREX_SPACEDIM >= 2
    y = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
    yl = problo[1] + dx[1] * (static_cast<Real>(j));
#endif

    Real z = 0.0;
    Real zl = 0.0;

#if AMREX_SPACEDIM == 3
    z = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
    zl = problo[2] + dx[2] * (static_cast<Real>(k));
#endif

    Real reint = 0.0;
    Real u_tot = 0.0;

    for (int kk = 0; kk < problem::nsub; kk++) {
        Real zz = zl + dx[2] * (static_cast<Real>(kk) + 0.5_rt) / problem::nsub;

        for (int jj = 0; jj < problem::nsub; jj++) {
            Real yy = yl + dx[1] * (static_cast<Real>(jj) + 0.5_rt) / problem::nsub;

            for (int ii = 0; ii < problem::nsub; ii++) {
                Real xx = xl + dx[0] * (static_cast<Real>(ii) + 0.5_rt) / problem::nsub;

                Real r = std::sqrt((xx - problem::center[0]) * (xx - problem::center[0]) +
                                   (yy - problem::center[1]) * (yy - problem::center[1]));

                Real u_phi;
                Real p;

                if (r < 0.2_rt) {
                    u_phi = 5.0_rt * r;
                    p = problem::p0 + 12.5_rt * r * r;

                } else if (r < 0.4_rt) {
                    u_phi = 2.0_rt - 5.0_rt * r;
                    p = problem::p0 + 12.5_rt * r * r + 4.0_rt *
                        (1.0_rt - 5.0_rt * r - std::log(0.2_rt) + std::log(r));

                } else {
                    u_phi = 0.0_rt;
                    p = problem::p0 - 2.0_rt + 4.0_rt * std::log(2.0_rt);
                }

                u_tot += u_phi;
                reint += p/(gamma_const - 1.0_rt);
            }
        }
    }

    Real u_phi = u_tot / (problem::nsub * problem::nsub * problem::nsub);
    reint = reint / (problem::nsub * problem::nsub * problem::nsub);

    state(i,j,k,URHO) = problem::rho0;

    // phi unit vector: \hat{\phi} = -sin(phi) \hat{x} + cos(phi) \hat{y}
    // with cos(phi) = x/r; sin(phi) = y/r
    Real r = std::sqrt((x - problem::center[0]) * (x - problem::center[0]) +
                       (y - problem::center[1]) * (y - problem::center[1]));

    //  -sin(phi) = y/r
    state(i,j,k,UMX) = -problem::rho0 * problem::q_r * u_phi * ((y - problem::center[1]) / r);
    //  cos(phi) = x/r
    state(i,j,k,UMY) =  problem::rho0 * problem::q_r * u_phi * ((x - problem::center[0]) / r);
    state(i,j,k,UMZ) = 0.0_rt;

    state(i,j,k,UEDEN) = reint +
        0.5_rt * (state(i,j,k,UMX) * state(i,j,k,UMX) / state(i,j,k,URHO) +
                  state(i,j,k,UMY) * state(i,j,k,UMY) / state(i,j,k,URHO) +
                  state(i,j,k,UMZ) * state(i,j,k,UMZ) / state(i,j,k,URHO));

    state(i,j,k,UEINT) = reint;

    state(i,j,k,UFS) = state(i,j,k,URHO);
    for (int n = 1; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = 0.0_rt;
    }
}

#endif
