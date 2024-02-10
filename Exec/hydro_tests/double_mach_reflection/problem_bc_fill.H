#ifndef problem_bc_fill_H
#define problem_bc_fill_H

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_bc_fill(int i, int j, int k,
                     Array4<Real> const& state,
                     Real time,
                     const Array1D<BCRec, 0, NUM_STATE-1>& bcs,
                     const GeometryData& geomdata)
{
    const int* domlo = geomdata.Domain().loVect();
    const int* domhi = geomdata.Domain().hiVect();

    const Real* problo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();

    //-------------------------------------------------------------------------
    // x boundaries
    //-------------------------------------------------------------------------

    if ((bcs(URHO).lo(0) == amrex::BCType::ext_dir || bcs(URHO).lo(0) == amrex::BCType::foextrap) && i < domlo[0]) {
        state(i,j,k,URHO) = problem::rho_l;
        state(i,j,k,UMX) = problem::rho_l * problem::u_l;
        state(i,j,k,UMY) = problem::rho_l * problem::v_l;
        state(i,j,k,UMZ) = 0.0e0_rt;
        state(i,j,k,UEDEN) = problem::rhoe_l + 0.5e0_rt * problem::rho_l *
                                               (problem::u_l * problem::u_l +
                                                problem::v_l * problem::v_l);
        state(i,j,k,UEINT) = problem::rhoe_l;
        state(i,j,k,UTEMP) = problem::T_l;
        for (int n = 0; n < NumSpec; ++n) {
            state(i,j,k,UFS+n) = 0.0e0_rt;
        }
        state(i,j,k,UFS) = state(i,j,k,URHO);
    }

#if AMREX_SPACEDIM >= 2
    //-------------------------------------------------------------------------
    // y boundaries
    //-------------------------------------------------------------------------

    Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5e0_rt);
    Real y = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5e0_rt);

    // YLO
    if ((bcs(URHO).lo(1) == amrex::BCType::ext_dir || bcs(URHO).lo(1) == amrex::BCType::foextrap) && j < domlo[1]) {
        if (x < 1.e0_rt / 6.e0_rt) {
            // ICs
            state(i,j,k,URHO) = problem::rho_l;
            state(i,j,k,UMX) = problem::rho_l * problem::u_l;
            state(i,j,k,UMY) = problem::rho_l * problem::v_l;
            state(i,j,k,UEDEN) = problem::rhoe_l + 0.5e0_rt * problem::rho_l *
                                                   (problem::u_l * problem::u_l +
                                                    problem::v_l * problem::v_l);
            state(i,j,k,UEINT) = problem::rhoe_l;
            state(i,j,k,UTEMP) = problem::T_l;
        } else {
            // reflect
            state(i,j,k,URHO) = state(i,domlo[1],k,URHO);
            state(i,j,k,UMX)  = state(i,domlo[1],k,UMX);
            state(i,j,k,UMY)  = -1.e0_rt * state(i,domlo[1],k,UMY);
            state(i,j,k,UEDEN) = state(i,domlo[1],k,UEDEN);
            state(i,j,k,UEINT) = state(i,domlo[1],k,UEINT);
            state(i,j,k,UTEMP) = state(i,domlo[1],k,UTEMP);
        }
        state(i,j,k,UMZ) = 0.0e0_rt;
        for (int n = 0; n < NumSpec; ++n) {
            state(i,j,k,UFS+n) = 0.0e0_rt;
        }
        state(i,j,k,UFS) = state(i,j,k,URHO);
    }

    // YHI
    if ((bcs(URHO).hi(2) == amrex::BCType::ext_dir || bcs(URHO).hi(2) == amrex::BCType::foextrap) && j > domhi[1]) {
        state(i,j,k,URHO ) = 0.e0_rt;
        state(i,j,k,UMX  ) = 0.e0_rt;
        state(i,j,k,UMY  ) = 0.e0_rt;
        state(i,j,k,UEDEN) = 0.e0_rt;
        state(i,j,k,UEINT) = 0.e0_rt;
        state(i,j,k,UTEMP) = 0.e0_rt;

        for (int jj = -1; jj <= 1; ++jj) {
            if (jj == 0) {
                continue;
            }

            Real y_loc = y + 0.5_rt * dx[1] * (jj / std::sqrt(3.0_rt));

            Real shockfront = (1.0_rt / 6.0_rt) + y_loc / std::tan(M_PI / 3.0_rt) + (10.e0_rt / std::sin(M_PI / 3.0_rt)) * time;

            for (int ii = -1; ii <= 1; ++ii) {
                if (ii == 0) {
                    continue;
                }

                Real x_loc = x + 0.5_rt * dx[0] * (ii / std::sqrt(3.0_rt));

                if (x_loc < shockfront) {
                    // Post shock ICs
                    state(i,j,k,URHO ) = state(i,j,k,URHO ) + (1.0_rt / 4.0_rt) * problem::rho_l;
                    state(i,j,k,UMX  ) = state(i,j,k,UMX  ) + (1.0_rt / 4.0_rt) * problem::rho_l * problem::u_l;
                    state(i,j,k,UMY  ) = state(i,j,k,UMY  ) + (1.0_rt / 4.0_rt) * problem::rho_l * problem::v_l;
                    state(i,j,k,UEDEN) = state(i,j,k,UEDEN) +
                                         (1.0_rt / 4.0_rt) * (problem::rhoe_l + 0.5e0_rt * problem::rho_l *
                                                                                (problem::u_l * problem::u_l +
                                                                                 problem::v_l * problem::v_l));
                    state(i,j,k,UEINT) = state(i,j,k,UEINT) + (1.0_rt / 4.0_rt) * problem::rhoe_l;
                    state(i,j,k,UTEMP) = state(i,j,k,UTEMP) + (1.0_rt / 4.0_rt) * problem::T_l;
                } else {
                    // Pre Shock ICs
                    state(i,j,k,URHO ) = state(i,j,k,URHO ) + (1.0_rt / 4.0_rt) * problem::rho_r;
                    state(i,j,k,UMX  ) = state(i,j,k,UMX  ) + (1.0_rt / 4.0_rt) * problem::rho_r * problem::u_r;
                    state(i,j,k,UMY  ) = state(i,j,k,UMY  ) + (1.0_rt / 4.0_rt) * problem::rho_r * problem::v_r;
                    state(i,j,k,UEDEN) = state(i,j,k,UEDEN) +
                                         (1.0_rt / 4.0_rt) * (problem::rhoe_r + 0.5e0_rt * problem::rho_r *
                                                                                (problem::u_r * problem::u_r +
                                                                                 problem::v_r * problem::v_r));
                    state(i,j,k,UEINT) = state(i,j,k,UEINT) + (1.0_rt / 4.0_rt) * problem::rhoe_r;
                    state(i,j,k,UTEMP) = state(i,j,k,UTEMP) + (1.0_rt / 4.0_rt) * problem::T_r;
                }

            }
        }

        for (int n = 0; n < NumSpec; ++n) {
            state(i,j,k,UFS+n) = 0.0e0_rt;
        }
        state(i,j,k,UFS) = state(i,j,k,URHO);
        state(i,j,k,UMZ) = 0.0e0_rt;
    }
#endif
}

#endif
