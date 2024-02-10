#ifndef problem_bc_fill_H
#define problem_bc_fill_H

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_bc_fill(int i, int j, int k,
                     Array4<Real> const& state,
                     Real time,
                     const Array1D<BCRec, 0, NUM_STATE-1>& bcs,
                     const GeometryData& geomdata)
{
    // The strategy here is to set Dirichlet condition for inflow and
    // outflow boundaries, and let the Riemann solver sort out the
    // proper upwinding.

    const int* domlo = geomdata.Domain().loVect();
    const int* domhi = geomdata.Domain().hiVect();

    if (problem::idir == 1) {

        // XLO
        if ((bcs(URHO).lo(0) == amrex::BCType::ext_dir || bcs(URHO).lo(0) == amrex::BCType::foextrap) && i < domlo[0]) {
            state(i,j,k,URHO) = problem::rho0;
            state(i,j,k,UMX) = problem::rho0 * problem::v0;
            state(i,j,k,UMY) = 0.e0_rt;
            state(i,j,k,UMZ) = 0.e0_rt;
            state(i,j,k,UFS) = state(i,j,k,URHO);
            if (NAUX_NET > 0) {
                state(i,j,k,UFX) = state(i,j,k,URHO);
            }
            state(i,j,k,UEINT) = problem::eint0;
            state(i,j,k,UEDEN) = problem::etot0;
            state(i,j,k,UTEMP) = problem::T0;
        }

        // XHI
        if ((bcs(URHO).hi(0) == amrex::BCType::ext_dir || bcs(URHO).hi(0) == amrex::BCType::foextrap) && i > domhi[0]) {
            state(i,j,k,URHO) = problem::rho1;
            state(i,j,k,UMX) = problem::rho1 * problem::v1;
            state(i,j,k,UMY) = 0.e0_rt;
            state(i,j,k,UMZ) = 0.e0_rt;
            state(i,j,k,UFS) = state(i,j,k,URHO);
            if (NAUX_NET > 0) {
                state(i,j,k,UFX) = state(i,j,k,URHO);
            }
            state(i,j,k,UEINT) = problem::eint1;
            state(i,j,k,UEDEN) = problem::etot1;
            state(i,j,k,UTEMP) = problem::T1;
        }

    }
    else if (problem::idir == 2) {

        // YLO
        if ((bcs(URHO).lo(1) == amrex::BCType::ext_dir || bcs(URHO).lo(1) == amrex::BCType::foextrap) && j < domlo[1]) {
            state(i,j,k,URHO) = problem::rho0;
            state(i,j,k,UMX) = 0.e0_rt;
            state(i,j,k,UMY) = problem::rho0 * problem::v0;
            state(i,j,k,UMZ) = 0.e0_rt;
            state(i,j,k,UFS) = state(i,j,k,URHO);
            if (NAUX_NET > 0) {
                state(i,j,k,UFX) = state(i,j,k,URHO);
            }
            state(i,j,k,UEINT) = problem::eint0;
            state(i,j,k,UEDEN) = problem::etot0;
            state(i,j,k,UTEMP) = problem::T0;
        }

        // YHI
        if ((bcs(URHO).hi(1) == amrex::BCType::ext_dir || bcs(URHO).hi(1) == amrex::BCType::foextrap) && j > domhi[1]) {
            state(i,j,k,URHO) = problem::rho1;
            state(i,j,k,UMX) = 0.e0_rt;
            state(i,j,k,UMY) = problem::rho1 * problem::v1;
            state(i,j,k,UMZ) = 0.e0_rt;
            state(i,j,k,UFS) = state(i,j,k,URHO);
            if (NAUX_NET > 0) {
                state(i,j,k,UFX) = state(i,j,k,URHO);
            }
            state(i,j,k,UEINT) = problem::eint1;
            state(i,j,k,UEDEN) = problem::etot1;
            state(i,j,k,UTEMP) = problem::T1;
        }

    }
    else if (problem::idir == 3) {

        // ZLO
        if ((bcs(URHO).lo(2) == amrex::BCType::ext_dir || bcs(URHO).lo(2) == amrex::BCType::foextrap) && k < domlo[2]) {
            state(i,j,k,URHO) = problem::rho0;
            state(i,j,k,UMX) = 0.e0_rt;
            state(i,j,k,UMY) = 0.e0_rt;
            state(i,j,k,UMZ) = problem::rho0 * problem::v0;
            state(i,j,k,UFS) = state(i,j,k,URHO);
            if (NAUX_NET > 0) {
                state(i,j,k,UFX) = state(i,j,k,URHO);
            }
            state(i,j,k,UEINT) = problem::eint0;
            state(i,j,k,UEDEN) = problem::etot0;
            state(i,j,k,UTEMP) = problem::T0;
        }

        // ZHI
        if ((bcs(URHO).hi(2) == amrex::BCType::ext_dir || bcs(URHO).hi(2) == amrex::BCType::foextrap) && k > domhi[2]) {
            state(i,j,k,URHO) = problem::rho1;
            state(i,j,k,UMX) = 0.e0_rt;
            state(i,j,k,UMY) = 0.e0_rt;
            state(i,j,k,UMZ) = problem::rho1 * problem::v1;
            state(i,j,k,UFS) = state(i,j,k,URHO);
            if (NAUX_NET > 0) {
                state(i,j,k,UFX) = state(i,j,k,URHO);
            }
            state(i,j,k,UEINT) = problem::eint1;
            state(i,j,k,UEDEN) = problem::etot1;
            state(i,j,k,UTEMP) = problem::T1;
        }

    }
}

#endif
