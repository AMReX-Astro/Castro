#include <ambient.H>
#include <Castro.H>
#include <Castro_bc_fill_nd.H>

void
ambient_denfill(const Box& bx, Array4<Real> const& state,
                Geometry const& geom, const Vector<BCRec>& bcr)
{
    // Make a copy of the BC that can be passed by value to the ParallelFor.

    BCRec bc = bcr[0];

    const auto domlo = geom.Domain().loVect3d();
    const auto domhi = geom.Domain().hiVect3d();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        bool ambient_x_lo = (castro::ambient_fill_dir == 0 || castro::ambient_fill_dir == -1) &&
                            (bc.lo(0) == FOEXTRAP || bc.lo(0) == HOEXTRAP);
        bool ambient_x_hi = (castro::ambient_fill_dir == 0 || castro::ambient_fill_dir == -1) &&
                            (bc.hi(0) == FOEXTRAP || bc.hi(0) == HOEXTRAP);

#if AMREX_SPACEDIM >= 2
        bool ambient_y_lo = (castro::ambient_fill_dir == 1 || castro::ambient_fill_dir == -1) &&
                            (bc.lo(1) == FOEXTRAP || bc.lo(1) == HOEXTRAP);
        bool ambient_y_hi = (castro::ambient_fill_dir == 1 || castro::ambient_fill_dir == -1) &&
                            (bc.hi(1) == FOEXTRAP || bc.hi(1) == HOEXTRAP);
#endif

#if AMREX_SPACEDIM == 3
        bool ambient_z_lo = (castro::ambient_fill_dir == 2 || castro::ambient_fill_dir == -1) &&
                            (bc.lo(2) == FOEXTRAP || bc.lo(2) == HOEXTRAP);
        bool ambient_z_hi = (castro::ambient_fill_dir == 2 || castro::ambient_fill_dir == -1) &&
                            (bc.hi(2) == FOEXTRAP || bc.hi(2) == HOEXTRAP);
#endif

        if (castro::fill_ambient_bc == 1) {
            if (ambient_x_lo && i < domlo[0] ||
                ambient_x_hi && i > domhi[0]
#if AMREX_SPACEDIM >= 2
                ||
                ambient_y_lo && j < domlo[1] ||
                ambient_y_hi && j > domhi[1]
#endif
#if AMREX_SPACEDIM == 3
                ||
                ambient_z_lo && k < domlo[2] ||
                ambient_z_hi && k > domhi[2]
#endif
                )
            {
                state(i,j,k) = ambient::ambient_state[URHO];
            }
        }
    });
}



void
ambient_fill(const Box& bx, Array4<Real> const& state,
             Geometry const& geom, const Vector<BCRec>& bcr)
{
    // Copy BCs to an Array1D so they can be passed by value to the ParallelFor.

    Array1D<BCRec, 0, NUM_STATE - 1> bcs;
    for (int n = 0; n < NUM_STATE; ++n) {
        bcs(n) = bcr[n];
    }

    const auto domlo = geom.Domain().loVect3d();
    const auto domhi = geom.Domain().hiVect3d();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        bool ambient_x_lo = (castro::ambient_fill_dir == 0 || castro::ambient_fill_dir == -1) &&
                            (bcs(URHO).lo(0) == FOEXTRAP || bcs(URHO).lo(0) == HOEXTRAP);
        bool ambient_x_hi = (castro::ambient_fill_dir == 0 || castro::ambient_fill_dir == -1) &&
                            (bcs(URHO).hi(0) == FOEXTRAP || bcs(URHO).hi(0) == HOEXTRAP);

#if AMREX_SPACEDIM >= 2
        bool ambient_y_lo = (castro::ambient_fill_dir == 1 || castro::ambient_fill_dir == -1) &&
                            (bcs(URHO).lo(1) == FOEXTRAP || bcs(URHO).lo(1) == HOEXTRAP);
        bool ambient_y_hi = (castro::ambient_fill_dir == 1 || castro::ambient_fill_dir == -1) &&
                            (bcs(URHO).hi(1) == FOEXTRAP || bcs(URHO).hi(1) == HOEXTRAP);
#endif

#if AMREX_SPACEDIM == 3
        bool ambient_z_lo = (castro::ambient_fill_dir == 2 || castro::ambient_fill_dir == -1) &&
                            (bcs(URHO).lo(2) == FOEXTRAP || bcs(URHO).lo(2) == HOEXTRAP);
        bool ambient_z_hi = (castro::ambient_fill_dir == 2 || castro::ambient_fill_dir == -1) &&
                            (bcs(URHO).hi(2) == FOEXTRAP || bcs(URHO).hi(2) == HOEXTRAP);
#endif

        if (castro::fill_ambient_bc == 1) {
            if (ambient_x_lo && i < domlo[0] ||
                ambient_x_hi && i > domhi[0]
#if AMREX_SPACEDIM >= 2
                ||
                ambient_y_lo && j < domlo[1] ||
                ambient_y_hi && j > domhi[1]
#endif
#if AMREX_SPACEDIM == 3
                ||
                ambient_z_lo && k < domlo[2] ||
                ambient_z_hi && k > domhi[2]
#endif
                ) {
                for (int n = 0; n < NUM_STATE; ++n) {
                    state(i,j,k,n) = ambient::ambient_state[n];
                }

                if (castro::ambient_outflow_vel == 1) {

                    // extrapolate the normal velocity only if it is outgoing
                    if (i < domlo[0]) {
                        state(i,j,k,UMX) = amrex::min(0.0_rt, state(domlo[0],j,k,UMX));
                        state(i,j,k,UMY) = 0.0_rt;
                        state(i,j,k,UMZ) = 0.0_rt;
                    } else if (i > domhi[0]) {
                        state(i,j,k,UMX) = amrex::max(0.0_rt, state(domhi[0],j,k,UMX));
                        state(i,j,k,UMY) = 0.0_rt;
                        state(i,j,k,UMZ) = 0.0_rt;
                    } else if (j < domlo[1]) {
                        state(i,j,k,UMX) = 0.0_rt;
                        state(i,j,k,UMY) = amrex::min(0.0_rt, state(i,domlo[1],k,UMY));
                        state(i,j,k,UMZ) = 0.0_rt;
                    } else if (j > domhi[1]) {
                        state(i,j,k,UMX) = 0.0_rt;
                        state(i,j,k,UMY) = amrex::max(0.0_rt, state(i,domhi[1],k,UMY));
                        state(i,j,k,UMZ) = 0.0_rt;
                    } else if (k < domlo[2]) {
                        state(i,j,k,UMX) = 0.0_rt;
                        state(i,j,k,UMY) = 0.0_rt;
                        state(i,j,k,UMZ) = amrex::min(0.0_rt, state(i,j,domlo[2],UMZ));
                    } else if (k > domhi[2]) {
                        state(i,j,k,UMX) = 0.0_rt;
                        state(i,j,k,UMY) = 0.0_rt;
                        state(i,j,k,UMZ) = amrex::max(0.0_rt, state(i,j,domhi[2],UMZ));
                    }

                    // now make the energy consistent
                    state(i,j,k,UEDEN) = state(i,j,k,UEINT) + 0.5_rt * (state(i,j,k,UMX) * state(i,j,k,UMX) +
                                                                        state(i,j,k,UMY) * state(i,j,k,UMY) +
                                                                        state(i,j,k,UMZ) * state(i,j,k,UMZ)) /
                                                                       state(i,j,k,URHO);

                }
            }
        }
    });
}

