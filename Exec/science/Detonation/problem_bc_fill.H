#ifndef problem_bc_fill_H
#define problem_bc_fill_H

// Given a zone state, fill it with ambient material.

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_fill_ambient(int i, int j, int k,
                          Array4<Real> const& state,
                          Real x, Real time,
                          const GeometryData& geomdata)
{
    const Real* problo = geomdata.ProbLo();
    const Real* probhi = geomdata.ProbHi();

    Real c_T = problo[0] + problem::center_T * (probhi[0] - problo[0]);

    state(i,j,k,URHO) = problem::ambient_dens;

    // This should be problem::ambient_comp but we currently don't have
    // arrays working in C++. This should be synced when we do.

    Real ambient_xn[NumSpec] = { problem::smallx };
    ambient_xn[problem::ic12-1] = amrex::max(problem::cfrac, problem::smallx);
    ambient_xn[problem::io16-1] = amrex::max(problem::ofrac, problem::smallx);
    ambient_xn[problem::ihe4-1] = 1.e0_rt - problem::cfrac - problem::ofrac - (NumSpec - 2) * problem::smallx;

    for (int n = 0; n < NumSpec; ++n) {
        state(i,j,k,UFS+n) = problem::ambient_dens * ambient_xn[n];
    }

    if (x < c_T) {
        state(i,j,k,UTEMP) = problem::T_l;
        state(i,j,k,UEINT) = state(i,j,k,URHO) * problem::ambient_e_l;
        state(i,j,k,UMX  ) = state(i,j,k,URHO) * (problem::vel + problem::grav_acceleration * time);
    }
    else {
        state(i,j,k,UTEMP) = problem::T_r;
        state(i,j,k,UEINT) = state(i,j,k,URHO) * problem::ambient_e_r;
        state(i,j,k,UMX  ) = -state(i,j,k,URHO) * (problem::vel + problem::grav_acceleration * time);
    }

    state(i,j,k,UMY) = 0.0;
    state(i,j,k,UMZ) = 0.0;
    state(i,j,k,UEDEN) = state(i,j,k,UEINT) + 0.5_rt * (state(i,j,k,UMX) * state(i,j,k,UMX) +
                                                        state(i,j,k,UMY) * state(i,j,k,UMY) +
                                                        state(i,j,k,UMZ) * state(i,j,k,UMZ)) /
                                                       state(i,j,k,URHO);
}



AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_bc_fill(int i, int j, int k,
                     Array4<Real> const& state,
                     Real time,
                     Array1D<BCRec, 0, NUM_STATE-1> bcs,
                     GeometryData geomdata)
{
    // Override the generic routine at the physical boundaries by
    // setting the material to the ambient state. Note that we
    // don't want to do this for interior/periodic boundaries
    // or for reflecting boundaries.

    if (castro::fill_ambient_bc == 1) {

        const int* domlo = geomdata.Domain().loVect();
        const int* domhi = geomdata.Domain().hiVect();

        bool do_ambient_fill = false;

        // Note: in these checks we're only looking at the boundary
        // conditions on the first state component (density).

        if (i < domlo[0] && (bcs(0).lo(0) != REFLECT_ODD &&
                             bcs(0).lo(0) != INT_DIR &&
                             bcs(0).lo(0) != REFLECT_EVEN)) {
            do_ambient_fill = true;
        }
        else if (i > domhi[0] && (bcs(0).hi(0) != REFLECT_ODD &&
                                  bcs(0).hi(0) != INT_DIR &&
                                  bcs(0).hi(0) != REFLECT_ODD)) {
            do_ambient_fill = true;
        }

        if (AMREX_SPACEDIM >= 2) {
            if (j < domlo[1] && (bcs(0).lo(1) != REFLECT_ODD &&
                                 bcs(0).lo(1) != INT_DIR &&
                                 bcs(0).lo(1) != REFLECT_ODD)) {
                do_ambient_fill = true;
            }
            else if (j > domhi[1] && (bcs(0).hi(1) != REFLECT_ODD &&
                                      bcs(0).hi(1) != INT_DIR &&
                                      bcs(0).hi(1) != REFLECT_ODD)) {
                do_ambient_fill = true;
            }
        }
                
        if (AMREX_SPACEDIM == 3) {
            if (k < domlo[2] && (bcs(0).lo(2) != REFLECT_ODD &&
                                 bcs(0).lo(2) != INT_DIR &&
                                 bcs(0).lo(2) != REFLECT_ODD)) {
                do_ambient_fill = true;
            }
            else if (k > domhi[2] && (bcs(0).hi(2) != REFLECT_ODD &&
                                      bcs(0).hi(2) != INT_DIR &&
                                      bcs(0).hi(2) != REFLECT_ODD)) {
                do_ambient_fill = true;
            }
        }

        if (do_ambient_fill) {
            const Real* dx = geomdata.CellSize();
            const Real* problo = geomdata.ProbLo();

            Real x = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0];

            problem_fill_ambient(i, j, k, state, x, time, geomdata);
        }

    }
}

#endif
