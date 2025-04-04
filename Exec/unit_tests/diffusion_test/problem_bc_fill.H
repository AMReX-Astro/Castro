#ifndef problem_bc_fill_H
#define problem_bc_fill_H

#include <prob_util.H>
#include <eos.H>

///
/// Problem-specific boundary condition fill
///
/// @param i         x-index
/// @param j         y-index
/// @param k         z-index
/// @param state     simulation state (Fab)
/// @param time      simulation time
/// @param bcs       boundary conditions
/// @param geomdata  geometry data
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_bc_fill(int i, int j, int k,
                      Array4<Real> const& state,
                      Real time, const Array1D<BCRec, 0, NUM_STATE-1>& bcs,
                      const GeometryData& geomdata) {

    // set Dirchilet boundary condition using analytic solution for diffusion problem.

    amrex::ignore_unused(j);
    amrex::ignore_unused(k);

    const auto coord = geomdata.Coord();
    const int* domlo = geomdata.Domain().loVect();
    const int* domhi = geomdata.Domain().hiVect();
    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real r[3] = {0.0_rt};
    r[0] = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);
#if AMREX_SPACEDIM >= 2
    r[1] = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
#endif
#if AMREX_SPACEDIM == 3
    r[2] = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
#endif
    if (coord_type <= 1) {
        for (int i; i < AMREX_SPACEDIM; ++i) {
            r[i] -= problem::center[i];
        }
    }

    if (diffuse_use_amrex_mlmg) {
        if (   ((bcs(UTEMP).lo(0) == amrex::BCType::ext_dir) && (i < domlo[0]))
            || ((bcs(UTEMP).hi(0) == amrex::BCType::ext_dir) && (i > domhi[0]))) {
            r[0] += dx[0] * 0.5_rt;
        }
#if AMREX_SPACEDIM >= 2
        if (   ((bcs(UTEMP).lo(1) == amrex::BCType::ext_dir) && (j < domlo[1]))
            || ((bcs(UTEMP).hi(1) == amrex::BCType::ext_dir) && (j > domhi[1]))) {
            r[1] += dx[1] * 0.5_rt;
        }
#endif
#if AMREX_SPACEDIM == 3
        if (   ((bcs(UTEMP).lo(2) == amrex::BCType::ext_dir) && (k < domlo[2]))
            || ((bcs(UTEMP).hi(2) == amrex::BCType::ext_dir) && (k > domhi[2]))) {
            r[2] += dx[2] * 0.5_rt;
        }
#endif
    }

    if (   ((bcs(UTEMP).lo(0) == amrex::BCType::ext_dir) && (i < domlo[0]))
        || ((bcs(UTEMP).hi(0) == amrex::BCType::ext_dir) && (i > domhi[0]))
#if AMREX_SPACEDIM >= 2
        || ((bcs(UTEMP).lo(1) == amrex::BCType::ext_dir) && (j < domlo[1]))
        || ((bcs(UTEMP).hi(1) == amrex::BCType::ext_dir) && (j > domhi[1]))
#endif
#if AMREX_SPACEDIM == 3
        || ((bcs(UTEMP).lo(2) == amrex::BCType::ext_dir) && (k < domlo[2]))
        || ((bcs(UTEMP).hi(2) == amrex::BCType::ext_dir) && (k > domhi[2]))
#endif
        ) {

        state(i,j,k,URHO) = problem::rho0;
        state(i,j,k,UMX) = 0.0_rt;
        state(i,j,k,UMY) = 0.0_rt;
        state(i,j,k,UMZ) = 0.0_rt;
        state(i,j,k,UTEMP) = analytic(r, time, coord);

        // compute internal energy
        eos_t eos_state;
        eos_state.T = state(i,j,k,UTEMP);
        eos_state.rho = state(i,j,k,URHO);
        Real X[NumSpec] = {0.0_rt};
        X[0] = 1.0_rt;
        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = X[n];
        }
        eos(eos_input_rt, eos_state);

        state(i,j,k,UEDEN) = problem::rho0 * eos_state.e;
        state(i,j,k,UEINT) = problem::rho0 * eos_state.e;
        for (int n = 0; n < NumSpec; n++) {
            state(i,j,k,UFS+n) = problem::rho0 * X[n];
        }
    }
}

#endif
