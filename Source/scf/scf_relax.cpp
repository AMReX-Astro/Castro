#include <Castro.H>
#include <fundamental_constants.H>
#include <Gravity.H>
#include <Rotation.H>

using namespace amrex;

#ifdef GRAVITY
#ifdef ROTATION
void Castro::scf_relaxation() {

    AMREX_ASSERT(level == 0);

    // First do some sanity checks.

    // Maximum density must be set.

    if (scf_maximum_density <= 0) {
        amrex::Error("castro.scf_maximum_density must be set for SCF relaxation.");
    }

    // Equatorial radius and polar radius must both
    // be positive, and the polar radius cannot be
    // larger than the equatorial radius.

    if (scf_equatorial_radius <= 0) {
        amrex::Error("Equatorial radius must be positive for SCF relaxation.");
    }

    if (scf_polar_radius <= 0) {
        amrex::Error("Polar radius must be positive for SCF relaxation.");
    }

    if (scf_polar_radius > scf_equatorial_radius) {
        amrex::Error("Polar radius must not be larger than equatorial radius for SCF relaxation.");
    }

    // Now do the SCF solve.

    do_hscf_solve();

    // Update the gravitational field. Since we don't need this
    // during the iterations, we only do it once at the end.

    for (int lev = 0; lev <= parent->finestLevel(); ++lev)
    {
        MultiFab& grav_new = getLevel(lev).get_new_data(Gravity_Type);
        gravity->get_new_grav_vector(lev, grav_new, getLevel(lev).state[Gravity_Type].curTime());
    }

    // Perform a final pass to ensure we're consistent with AMR.

    for (int lev = parent->finestLevel()-1; lev >= 0; --lev) {
        getLevel(lev).avgDown();
    }

}

void
Castro::do_hscf_solve()
{

    const int finest_level = parent->finestLevel();
    const int n_levs = finest_level + 1;

    // Do the initial relaxation setup. We need to fix two points
    // to uniquely determine an equilibrium configuration for a
    // rotating star. We are provided the distance from the center,
    // and the axis we're measuring on.

    GpuArray<Real, 3> scf_r_A = {problem::center[0], problem::center[1], problem::center[2]};
    GpuArray<Real, 3> scf_r_B = {problem::center[0], problem::center[1], problem::center[2]};

    // Note that the sign is somewhat arbitrary here.

    scf_r_A[0] += castro::scf_equatorial_radius;
    scf_r_B[2] += castro::scf_polar_radius;

    // To do the relaxation loop, we need to know the target maximum
    // enthalpy on the domain. So we'll loop through the grid and
    // calculate what it should be, based on the density, temperature,
    // and composition. We do not specify the latter two -- we assume
    // that these are filled as desired in initdata. So, to get the
    // maximum enthalpy, we'll loop through the grid, using the composition
    // and temperature in each zone, combined with the requested maximum
    // density, call the EOS to calculate the enthalpy for that combination,
    // and then take the maximum of that. (Note that the current density on
    // the grid is treated as an initial guess and does not need to reach
    // the requested maximum.) The SCF procedure is only intended for a setup
    // with uniform temperature and composition, so in that case this procedure
    // should yield the same result as just picking a zone at random and using
    // its temperature and composition. The reason for doing it this way is so
    // that the user doesn't have to request temperature and composition in their
    // inputs file (which is messy in general) -- they can tell us what composition
    // they want by filling it in the initialization.

    Real target_h_max = 0.0;

    for (int lev = 0; lev <= finest_level; ++lev) {

        MultiFab& state_new = getLevel(lev).get_new_data(State_Type);

        ReduceOps<ReduceOpMax> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(state_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.tilebox();

            auto state_arr = state_new[mfi].array();

            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                eos_t eos_state;

                eos_state.rho = castro::scf_maximum_density;
                eos_state.T   = state_arr(i,j,k,UTEMP);
                for (int n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] = state_arr(i,j,k,UFS+n) / state_arr(i,j,k,URHO);
                }
#if NAUX_NET > 0
                for (int n = 0; n < NumAux; ++n) {
                    eos_state.aux[n] = state_arr(i,j,k,UFX+n) / state_arr(i,j,k,URHO);
                }
#endif

                eos(eos_input_rt, eos_state);

                return {eos_state.h};
            });

        }

        ReduceTuple hv = reduce_data.value();
        target_h_max = amrex::max(target_h_max, amrex::get<0>(hv));

    }

    ParallelDescriptor::ReduceRealMax(target_h_max);

    Vector< std::unique_ptr<MultiFab> > psi(n_levs);
    Vector< std::unique_ptr<MultiFab> > enthalpy(n_levs);
    Vector< std::unique_ptr<MultiFab> > state_vec(n_levs);
    Vector< std::unique_ptr<MultiFab> > phi(n_levs);

    // Iterate until the system is relaxed by filling the level data
    // and then doing a multilevel gravity solve.

    int ctr = 1;

    while (ctr <= scf_max_iterations) {

        Real time = getLevel(0).state[State_Type].curTime();

        // Construct a local MultiFab for the rotational psi.
        // This does not change over the loop iterations, but
        // this (and the data to follow) must be constructed
        // inside the loop on each iteration because of regrids.

        for (int lev = 0; lev <= finest_level; ++lev) {

            psi[lev] = std::make_unique<MultiFab>(getLevel(lev).grids, getLevel(lev).dmap, 1, 0);

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi((*psi[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {

                const Box& bx = mfi.tilebox();

                fill_rotational_psi(bx, (*psi[lev]).array(mfi), time);

            }

        }

        // Construct a local MultiFab for the enthalpy.

        for (int lev = 0; lev <= finest_level; ++lev) {
            enthalpy[lev] = std::make_unique<MultiFab>(getLevel(lev).grids, getLevel(lev).dmap, 1, 0);
        }

        // Construct a local MultiFab for the state data. We do
        // this because we have to mask out the state several times
        // in the below calculation, and it's easiest to have a scratch
        // copy of the data to work with.

        for (int lev = 0; lev <= finest_level; ++lev) {
            state_vec[lev] = std::make_unique<MultiFab>(getLevel(lev).grids, getLevel(lev).dmap, NUM_STATE, 0);
            phi[lev] = std::make_unique<MultiFab>(getLevel(lev).grids, getLevel(lev).dmap, 1, 0);
        }

        // Copy in the state data. Mask it out on coarse levels.

        for (int lev = 0; lev <= finest_level; ++lev) {

            MultiFab::Copy((*state_vec[lev]), getLevel(lev).get_new_data(State_Type), 0, 0, NUM_STATE, 0);
            MultiFab::Copy((*phi[lev]), getLevel(lev).get_new_data(PhiGrav_Type), 0, 0, 1, 0);

            if (lev < finest_level) {
                const MultiFab& mask = getLevel(lev+1).build_fine_mask();

                for (int n = 0; n < NUM_STATE; ++n) {
                    MultiFab::Multiply((*state_vec[lev]), mask, 0, n, 1, 0);
                }

                MultiFab::Multiply((*phi[lev]), mask, 0, 0, 1, 0);
            }

        }

        // First step is to find the rotational frequency.

        Real phi_A = 0.0;
        Real psi_A = 0.0;
        Real phi_B = 0.0;
        Real psi_B = 0.0;

        for (int lev = 0; lev <= finest_level; ++lev) {

            auto geomdata = parent->Geom(lev).data();

            ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum> reduce_op;
            ReduceData<Real, Real, Real, Real> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi((*phi[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {

                const Box& bx = mfi.tilebox();

                auto phi_arr = (*phi[lev])[mfi].array();
                auto psi_arr = (*psi[lev])[mfi].array();

                reduce_op.eval(bx, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {
                    const auto *dx = geomdata.CellSize();
                    const auto *problo = geomdata.ProbLo();
                    const auto coord = geomdata.Coord();

                    // The below assumes we are rotating on the z-axis.

                    Real r[3] = {0.0};

                    r[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0];
#if AMREX_SPACEDIM >= 2
                    r[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1];
#endif
#if AMREX_SPACEDIM == 3
                    r[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - problem::center[2];
#endif

                    // Do a trilinear interpolation to find the contribution from
                    // this grid point. Limit so that only the nearest zone centers
                    // can participate. This implies that the maximum allowable
                    // distance from the target location is 0.5 * dx.

                    Real rr[3] = {0.0};

                    rr[0] = std::abs(r[0] - scf_r_A[0]) / dx[0];
#if AMREX_SPACEDIM >= 2
                    rr[1] = std::abs(r[1] - scf_r_A[1]) / dx[1];
#endif
#if AMREX_SPACEDIM == 3
                    rr[2] = std::abs(r[2] - scf_r_A[2]) / dx[2];
#endif

                    Real scale;

                    if (rr[0] > 1.0_rt || rr[1] > 1.0_rt || rr[2] > 1.0_rt) {
                        scale = 0.0;
                    }
                    else {
                        scale = (1.0_rt - rr[0]) * (1.0_rt - rr[1]) * (1.0_rt - rr[2]);
                    }

                    Real dphi_A = scale * phi_arr(i,j,k);
                    Real dpsi_A = scale * psi_arr(i,j,k);

                    rr[0] = std::abs(r[0] - scf_r_B[0]) / dx[0];
#if AMREX_SPACEDIM >= 2
                    rr[1] = std::abs(r[1] - scf_r_B[1]) / dx[1];
#endif
#if AMREX_SPACEDIM == 3
                    rr[2] = std::abs(r[2] - scf_r_B[2]) / dx[2];
#endif

                    if (rr[0] > 1.0_rt || rr[1] > 1.0_rt || rr[2] > 1.0_rt) {
                        scale = 0.0;
                    }
                    else {
                        scale = (1.0_rt - rr[0]) * (1.0_rt - rr[1]) * (1.0_rt - rr[2]);
                    }

                    Real dphi_B = scale * phi_arr(i,j,k);
                    Real dpsi_B = scale * psi_arr(i,j,k);

                    return {dphi_A, dpsi_A, dphi_B, dpsi_B};
                });
            }

            ReduceTuple hv = reduce_data.value();
            phi_A += amrex::get<0>(hv);
            psi_A += amrex::get<1>(hv);
            phi_B += amrex::get<2>(hv);
            psi_B += amrex::get<3>(hv);

        }

        ParallelDescriptor::ReduceRealSum(phi_A);
        ParallelDescriptor::ReduceRealSum(psi_A);
        ParallelDescriptor::ReduceRealSum(phi_B);
        ParallelDescriptor::ReduceRealSum(psi_B);

        // Now update the square of the rotation frequency, following Hachisu (Equation 16).
        // Deal carefully with the special case where phi_A and phi_B are equal -- we assume
        // this is the non-rotating case, and we want to handle that gracefully. Since the
        // rotation module knows to treat rotational_period = 0 as effectively disabling
        // rotation, we'll use that approach.

        if (std::abs(scf_polar_radius - scf_equatorial_radius) / std::abs(scf_equatorial_radius) < 1.e-6) {

            rotational_period = 0.0;

        } else {

            Real omegasq = -(phi_A - phi_B) / (psi_A - psi_B);

            if (omegasq < 0.0 && ParallelDescriptor::IOProcessor()) {
                std::cerr << "Omega squared = " << omegasq << " is negative in the relaxation step; aborting." << std::endl;
                amrex::Error();
            }

            Real omega = std::sqrt(omegasq);

            // Rotational period is 2 pi / omega.
            // Let's also be sure not to let the period
            // change by too much in a single iteration.

            rotational_period = amrex::min(1.1_rt * rotational_period,
                                           amrex::max(0.9_rt * rotational_period,
                                                      2.0_rt * M_PI / omega));

        }


        // Second step is to evaluate the Bernoulli constant.

        Real bernoulli = 0.0;

        for (int lev = 0; lev <= finest_level; ++lev) {

            auto geomdata = parent->Geom(lev).data();

            ReduceOps<ReduceOpSum> reduce_op;
            ReduceData<Real> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi((*phi[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {

                const Box& bx = mfi.tilebox();

                auto phi_arr = (*phi[lev])[mfi].array();

                reduce_op.eval(bx, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {
                    const auto *dx = geomdata.CellSize();
                    const auto *problo = geomdata.ProbLo();

                    // The below assumes we are rotating on the z-axis.

                    GpuArray<Real, 3> r = {0.0};

                    r[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0];
#if AMREX_SPACEDIM >= 2
                    r[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1];
#endif
#if AMREX_SPACEDIM == 3
                    r[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - problem::center[2];
#endif

                    // Do a trilinear interpolation to find the contribution from
                    // this grid point. Limit so that only the nearest zone centers
                    // can participate. This implies that the maximum allowable
                    // distance from the target location is 0.5 * dx.

                    Real rr[3] = {0.0};

                    rr[0] = std::abs(r[0] - scf_r_A[0]) / dx[0];
#if AMREX_SPACEDIM >= 2
                    rr[1] = std::abs(r[1] - scf_r_A[1]) / dx[1];
#endif
#if AMREX_SPACEDIM == 3
                    rr[2] = std::abs(r[2] - scf_r_A[2]) / dx[2];
#endif
                    Real scale;

                    if (rr[0] > 1.0_rt || rr[1] > 1.0_rt || rr[2] > 1.0_rt) {
                        scale = 0.0;
                    }
                    else {
                        scale = (1.0_rt - rr[0]) * (1.0_rt - rr[1]) * (1.0_rt - rr[2]);
                    }

                    auto omega = get_omega_vec(geomdata, j);
                    Real bernoulli_zone = scale * (phi_arr(i,j,k) + rotational_potential(r, omega, coord));

                    return {bernoulli_zone};
                });

            }

            ReduceTuple hv = reduce_data.value();
            bernoulli += amrex::get<0>(hv);

        }

        ParallelDescriptor::ReduceRealSum(bernoulli);



        // Third step is to construct the enthalpy field and
        // find the maximum enthalpy for the star.

        for (int lev = 0; lev <= finest_level; ++lev) {

            auto geomdata = parent->Geom(lev).data();

            enthalpy[lev]->setVal(0.0);

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi((*phi[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {

                const Box& bx = mfi.tilebox();

                auto enthalpy_arr = (*enthalpy[lev])[mfi].array();
                auto phi_arr = (*phi[lev])[mfi].array();

                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // The Bernoulli equation says that energy is conserved:
                    // enthalpy + gravitational potential + rotational potential = const
                    // We already have the constant, so our goal is to construct the enthalpy field.

                    const auto *dx = geomdata.CellSize();
                    const auto *problo = geomdata.ProbLo();
                    const auto coord = geomdata.Coord();

                    GpuArray<Real, 3> r = {0.0};

                    r[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0];
#if AMREX_SPACEDIM >= 2
                    r[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1];
#endif
#if AMREX_SPACEDIM == 3
                    r[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - problem::center[2];
#endif

                    auto omega = get_omega_vec(geomdata, j);
                    enthalpy_arr(i,j,k) = bernoulli - phi_arr(i,j,k) - rotational_potential(r, omega, coord);
                });

            }

        }

        Real actual_h_max = 0.0;
        Real actual_rho_max = 0.0;

        for (int lev = 0; lev <= finest_level; ++lev) {

            actual_h_max = std::max(actual_h_max, enthalpy[lev]->max(0));
            actual_rho_max = std::max(actual_rho_max, state_vec[lev]->max(URHO));

        }

        Real Linf_norm = 0.0;

        // Finally, update the density using the enthalpy field.

        for (int lev = 0; lev <= finest_level; ++lev) {

            ReduceOps<ReduceOpMax> reduce_op;
            ReduceData<Real> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi((*state_vec[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {

                const Box& bx = mfi.tilebox();

                auto enthalpy_arr = (*enthalpy[lev])[mfi].array();
                auto state_arr = (*state_vec[lev])[mfi].array();

                reduce_op.eval(bx, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {
                    Real norm = 0.0;

                    // We only want to call the EOS for zones with enthalpy > 0.
                    // For distances far enough from the center, the rotation
                    // term can overcome the other terms and make the enthalpy
                    // spuriously negative. If the enthalpy is negative, we just
                    // leave the zone alone -- this should be ambient material.

                    if (enthalpy_arr(i,j,k) > 0.0 && state_arr(i,j,k,URHO) > 0.0) {

                        Real old_rho = state_arr(i,j,k,URHO);

                        // Rescale the enthalpy by the maximum allowed value.

                        enthalpy_arr(i,j,k) = target_h_max * (enthalpy_arr(i,j,k) / actual_h_max);

                        eos_t eos_state;

                        eos_state.rho = state_arr(i,j,k,URHO); // Initial guess for the EOS
                        eos_state.T   = state_arr(i,j,k,UTEMP);
                        for (int n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = state_arr(i,j,k,UFS+n) / state_arr(i,j,k,URHO);
                        }
#if NAUX_NET > 0
                        for (int n = 0; n < NumAux; ++n) {
                            eos_state.aux[n] = state_arr(i,j,k,UFX+n) / state_arr(i,j,k,URHO);
                        }
#endif
                        eos_state.h   = enthalpy_arr(i,j,k);

                        eos(eos_input_th, eos_state);

                        state_arr(i,j,k,URHO)  = eos_state.rho;
                        state_arr(i,j,k,UTEMP) = eos_state.T;
                        state_arr(i,j,k,UEINT) = state_arr(i,j,k,URHO) * eos_state.e;
                        for (int n = 0; n < NumSpec; ++n) {
                            state_arr(i,j,k,UFS+n) = state_arr(i,j,k,URHO) * eos_state.xn[n];
                        }

                        state_arr(i,j,k,UMX) = 0.0;
                        state_arr(i,j,k,UMY) = 0.0;
                        state_arr(i,j,k,UMZ) = 0.0;

                        state_arr(i,j,k,UEDEN) = state_arr(i,j,k,UEINT);

                        // Convergence test

                        // Zones only participate in this test if they have a density
                        // that is above a certain fraction of the peak, to avoid
                        // oscillations in low density zones stalling convergence.

                        Real drho = std::abs(state_arr(i,j,k,URHO) - old_rho) / old_rho;

                        if (state_arr(i,j,k,URHO) / actual_rho_max > 1.0e-3_rt) {
                            norm = drho;
                        }
                    }

                    return {norm};
                });

            }

            ReduceTuple hv = reduce_data.value();
            Linf_norm = amrex::max(Linf_norm, amrex::get<0>(hv));

        }

        ParallelDescriptor::ReduceRealMax(Linf_norm);

        // Copy state data back to its source, and synchronize it on coarser levels.

        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(getLevel(lev).get_new_data(State_Type), (*state_vec[lev]), 0, 0, NUM_STATE, 0);
            MultiFab::Copy(getLevel(lev).get_new_data(PhiGrav_Type), (*phi[lev]), 0, 0, 1, 0);
        }

        for (int lev = finest_level-1; lev >= 0; --lev) {
            getLevel(lev).avgDown();
        }

        // Since we've changed the density distribution on the grid, regrid.

        bool do_io = false;
        if (finest_level > 0) {
            parent->RegridOnly(time, do_io);
        }

        // Update the gravitational field -- only after we've completed cleaning up the state above.

        gravity->multilevel_solve_for_new_phi(0, finest_level);

        // Update diagnostic quantities.

        Real kin_eng = 0.0;
        Real pot_eng = 0.0;
        Real int_eng = 0.0;
        Real mass = 0.0;

        for (int lev = 0; lev <= finest_level; ++lev) {

            auto geomdata = parent->Geom(lev).data();
            const auto dx = parent->Geom(level).CellSizeArray();

            Real dV = dx[0];
#if AMREX_SPACEDIM >= 2
            dV *= dx[1];
#endif
#if AMREX_SPACEDIM == 3
            dV *= dx[2];
#endif

            ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum> reduce_op;
            ReduceData<Real, Real, Real, Real> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi((*state_vec[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {

                const Box& bx = mfi.tilebox();

                auto state_arr = (*state_vec[lev])[mfi].array();
                auto phi_arr = (*phi[lev])[mfi].array();

                reduce_op.eval(bx, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {
                    Real dM = 0.0, dK = 0.0, dU = 0.0, dE = 0.0;

                    const auto* problo = geomdata.ProbLo();
                    const auto coord = geomdata.Coord();

                    GpuArray<Real, 3> r = {0.0};

                    r[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0];
#if AMREX_SPACEDIM >= 2
                    r[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1];
#endif
#if AMREX_SPACEDIM == 3
                    r[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - problem::center[2];
#endif

                    if (state_arr(i,j,k,URHO) > 0.0)
                    {
                        dM = state_arr(i,j,k,URHO) * dV;

                        auto omega = get_omega_vec(geomdata, j);
                        dK = rotational_potential(r, omega, coord) * dM;

                        dU = 0.5_rt * phi_arr(i,j,k) * dM;

                        eos_t eos_state;

                        eos_state.rho = state_arr(i,j,k,URHO);
                        eos_state.T   = state_arr(i,j,k,UTEMP);
                        for (int n = 0; n < NumSpec; ++n) {
                            eos_state.xn[n] = state_arr(i,j,k,UFS+n) / state_arr(i,j,k,URHO);
                        }
#if NAUX_NET > 0
                        for (int n = 0; n < NumAux; ++n) {
                            eos_state.aux[n] = state_arr(i,j,k,UFX+n) / state_arr(i,j,k,URHO);
                        }
#endif

                        eos(eos_input_rt, eos_state);

                        dE = eos_state.p * dV;
                    }

                    return {dM, dK, dU, dE};
                });

            }

            ReduceTuple hv = reduce_data.value();
            mass += amrex::get<0>(hv);
            kin_eng += amrex::get<1>(hv);
            pot_eng += amrex::get<2>(hv);
            int_eng += amrex::get<3>(hv);

        }

        ParallelDescriptor::ReduceRealSum(kin_eng);
        ParallelDescriptor::ReduceRealSum(pot_eng);
        ParallelDescriptor::ReduceRealSum(int_eng);
        ParallelDescriptor::ReduceRealSum(mass);

        Real virial_error = std::abs(2.0 * kin_eng + pot_eng + 3.0 * int_eng) / std::abs(pot_eng);

        // Now check to see if we're converged.

        int is_relaxed = 0;

        if (Linf_norm < scf_relax_tol) {
            is_relaxed = 1;
        }

        if (ParallelDescriptor::IOProcessor()) {

            // Grab the value for the solar mass.

            std::cout << std::endl << std::endl;
            std::cout << "   Relaxation iterations completed: " << ctr << std::endl;
            std::cout << "   L-infinity norm of residual (relative to old state): " << Linf_norm << std::endl;
            std::cout << "   Rotational period (s): " << rotational_period << std::endl;
            std::cout << "   Kinetic energy: " << kin_eng << std::endl;
            std::cout << "   Potential energy: " << pot_eng << std::endl;
            std::cout << "   Internal energy: " << int_eng << std::endl;
            std::cout << "   Virial error: " << virial_error << std::endl;
            std::cout << "   Mass: " << mass / C::M_solar << " solar masses" << std::endl;

            if (is_relaxed == 1) {
                std::cout << "  Relaxation completed!" << std::endl;
            }

            std::cout << std::endl << std::endl;

        }

        if (is_relaxed == 1) {
            break;
        }

        ctr++;

    }

}
#endif
#endif
