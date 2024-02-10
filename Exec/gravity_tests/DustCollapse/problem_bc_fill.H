#ifndef problem_bc_fill_H
#define problem_bc_fill_H

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_bc_fill(int i, int j, int k,
                     Array4<Real> const& state,
                     Real time,
                     const Array1D<BCRec, 0, NUM_STATE-1>& bcs,
                     const GeometryData& geomdata)
{
#if AMREX_SPACEDIM == 3
    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    const int* domlo = geomdata.Domain().loVect();
    const int* domhi = geomdata.Domain().hiVect();

    Real x = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0];

    Real y = 0.0;
#if AMREX_SPACEDIM >= 2
    y = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1];
#endif

    Real z = 0.0;
#if AMREX_SPACEDIM == 3
    z = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - problem::center[2];
#endif

    Real r = std::sqrt(x * x + y * y + z * z);

    // XLO
    if (bcs(0).lo(0) == amrex::BCType::foextrap && i < domlo[0]) {

        int ic = domlo[0];
        Real xc = problo[0] + (static_cast<Real>(ic) + 0.5_rt) * dx[0] - problem::center[0];
        Real rc = std::sqrt(xc * xc + y * y + z * z);

        Real momc = std::sqrt(state(ic,j,k,UMX) * state(ic,j,k,UMX) +
                              state(ic,j,k,UMY) * state(ic,j,k,UMY) +
                              state(ic,j,k,UMZ) * state(ic,j,k,UMZ));

        Real mom  = momc * (rc / r) * (rc / r);

        // Project along the normal
        state(i,j,k,UMX) = std::copysign(1.e0_rt, state(ic,j,k,UMX)) * mom * (x / r);
        state(i,j,k,UMY) = std::copysign(1.e0_rt, state(ic,j,k,UMY)) * mom * (y / r);
        state(i,j,k,UMZ) = std::copysign(1.e0_rt, state(ic,j,k,UMZ)) * mom * (z / r);

    }

    // XHI
    if (bcs(0).hi(0) == amrex::BCType::foextrap && i > domhi[0]) {

        int ic = domhi[0];
        Real xc = problo[0] + (static_cast<Real>(ic) + 0.5_rt) * dx[0] - problem::center[0];
        Real rc = std::sqrt(xc * xc + y * y + z * z);

        Real momc = std::sqrt(state(ic,j,k,UMX) * state(ic,j,k,UMX) +
                              state(ic,j,k,UMY) * state(ic,j,k,UMY) +
                              state(ic,j,k,UMZ) * state(ic,j,k,UMZ));

        Real mom  = momc * (rc / r) * (rc / r);

        // Project along the normal
        state(i,j,k,UMX) = std::copysign(1.e0_rt, state(ic,j,k,UMX)) * mom * (x / r);
        state(i,j,k,UMY) = std::copysign(1.e0_rt, state(ic,j,k,UMY)) * mom * (y / r);
        state(i,j,k,UMZ) = std::copysign(1.e0_rt, state(ic,j,k,UMZ)) * mom * (z / r);

    }

#if AMREX_SPACEDIM >= 2
    // YLO
    if (bcs(0).lo(1) == amrex::BCType::foextrap && j < domlo[1]) {

        int jc = domlo[1];
        Real yc = problo[1] + (static_cast<Real>(jc) + 0.5_rt) * dx[1] - problem::center[1];
        Real rc = std::sqrt(x * x + yc * yc + z * z);

        Real momc = std::sqrt(state(i,jc,k,UMX) * state(i,jc,k,UMX) +
                              state(i,jc,k,UMY) * state(i,jc,k,UMY) +
                              state(i,jc,k,UMZ) * state(i,jc,k,UMZ));

        Real mom  = momc * (rc / r) * (rc / r);

        // Project along the normal
        state(i,j,k,UMX) = std::copysign(1.e0_rt, state(i,jc,k,UMX)) * mom * (x / r);
        state(i,j,k,UMY) = std::copysign(1.e0_rt, state(i,jc,k,UMY)) * mom * (y / r);
        state(i,j,k,UMZ) = std::copysign(1.e0_rt, state(i,jc,k,UMZ)) * mom * (z / r);

    }

    // YHI
    if (bcs(0).hi(1) == amrex::BCType::foextrap && j > domhi[1]) {

        int jc = domhi[1];
        Real yc = problo[1] + (static_cast<Real>(jc) + 0.5_rt) * dx[1] - problem::center[1];
        Real rc = std::sqrt(x * x + yc * yc + z * z);

        Real momc = std::sqrt(state(i,jc,k,UMX) * state(i,jc,k,UMX) +
                              state(i,jc,k,UMY) * state(i,jc,k,UMY) +
                              state(i,jc,k,UMZ) * state(i,jc,k,UMZ));

        Real mom  = momc * (rc / r) * (rc / r);

        // Project along the normal
        state(i,j,k,UMX) = std::copysign(1.e0_rt, state(i,jc,k,UMX)) * mom * (x / r);
        state(i,j,k,UMY) = std::copysign(1.e0_rt, state(i,jc,k,UMY)) * mom * (y / r);
        state(i,j,k,UMZ) = std::copysign(1.e0_rt, state(i,jc,k,UMZ)) * mom * (z / r);
    }
#endif

#if AMREX_SPACEDIM == 3
    // ZLO
    if (bcs(0).lo(2) == amrex::BCType::foextrap && k < domlo[2]) {

        int kc = domlo[2];
        Real zc = problo[2] + (static_cast<Real>(kc) + 0.5_rt) * dx[2] - problem::center[2];
        Real rc = std::sqrt(x * x + y * y + zc * zc);

        Real momc = std::sqrt(state(i,j,kc,UMX) * state(i,j,kc,UMX) +
                              state(i,j,kc,UMY) * state(i,j,kc,UMY) +
                              state(i,j,kc,UMZ) * state(i,j,kc,UMZ));

        Real mom  = momc * (rc / r) * (rc / r);

        // Project along the normal
        state(i,j,k,UMX) = std::copysign(1.e0_rt, state(i,j,kc,UMX)) * mom * (x / r);
        state(i,j,k,UMY) = std::copysign(1.e0_rt, state(i,j,kc,UMY)) * mom * (y / r);
        state(i,j,k,UMZ) = std::copysign(1.e0_rt, state(i,j,kc,UMZ)) * mom * (z / r);

    }

    // ZHI
    if (bcs(0).hi(2) == amrex::BCType::foextrap && k > domhi[2]) {

        int kc = domhi[2];
        Real zc = problo[2] + (static_cast<Real>(kc) + 0.5_rt) * dx[2] - problem::center[2];
        Real rc = std::sqrt(x * x + y * y + zc * zc);

        Real momc = std::sqrt(state(i,j,kc,UMX) * state(i,j,kc,UMX) +
                              state(i,j,kc,UMY) * state(i,j,kc,UMY) +
                              state(i,j,kc,UMZ) * state(i,j,kc,UMZ));

        Real mom  = momc * (rc / r) * (rc / r);

        // Project along the normal
        state(i,j,k,UMX) = std::copysign(1.e0_rt, state(i,j,kc,UMX)) * mom * (x / r);
        state(i,j,k,UMY) = std::copysign(1.e0_rt, state(i,j,kc,UMY)) * mom * (y / r);
        state(i,j,k,UMZ) = std::copysign(1.e0_rt, state(i,j,kc,UMZ)) * mom * (z / r);

    }
#endif
#endif
}

#endif
