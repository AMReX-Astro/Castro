
#include <Castro.H>
#include <advection_util.H>
#ifdef MODEL_PARSER
#include <model_parser.H>
#endif
#include <sdc_cons_to_burn.H>

using std::string;
using namespace amrex;

#ifndef TRUE_SDC

advance_status
Castro::do_old_reactions (Real time, Real dt) {  // NOLINT(readability-convert-member-functions-to-static)

    amrex::ignore_unused(time);
    amrex::ignore_unused(dt);

    advance_status status {};

#ifndef SIMPLIFIED_SDC
    int burn_success{1};

    MultiFab& R_old = get_old_data(Reactions_Type);
    MultiFab& R_new = get_new_data(Reactions_Type);

    if (time_integration_method != SimplifiedSpectralDeferredCorrections) {
        // The result of the reactions is added directly to Sborder.
        burn_success = react_state(Sborder, R_old, time, 0.5 * dt, 0);

        if (burn_success != 1) {
            status.success = false;
            status.reason = "burn unsuccessful";

            return status;
        }

        clean_state(
#ifdef MHD
                    Bx_old_tmp, By_old_tmp, Bz_old_tmp,
#endif
                    Sborder, time, Sborder.nGrow());

        MultiFab::Copy(R_new, R_old, 0, 0, R_new.nComp(), R_new.nGrow());
    }
#endif

    return status;
}

advance_status
Castro::do_new_reactions (Real time, Real dt)
{
    amrex::ignore_unused(time);
    amrex::ignore_unused(dt);

    advance_status status {};

    int burn_success{1};

    MultiFab& R_new = get_new_data(Reactions_Type);
    MultiFab& S_new = get_new_data(State_Type);

#ifdef SIMPLIFIED_SDC
    MultiFab& R_old = get_old_data(Reactions_Type);

    if (time_integration_method == SimplifiedSpectralDeferredCorrections) {

        if (do_react) {

            // Do the ODE integration to capture the reaction source terms.

            burn_success = react_state(time, dt);

            if (burn_success != 1) {
                status.success = false;
                status.reason = "burn unsuccessful";

                return status;
            }

            clean_state(S_new, time + dt, S_new.nGrow());

            // Check for NaN's.

            check_for_nan(S_new);

        }
        else {

            // If we're not burning, just initialize the reactions data to zero.

            MultiFab& SDC_react_new = get_new_data(Simplified_SDC_React_Type);
            SDC_react_new.setVal(0.0, SDC_react_new.nGrow());

            R_old.setVal(0.0, R_old.nGrow());
            R_new.setVal(0.0, R_new.nGrow());

        }

    }

#else // SIMPLIFIED_SDC

    if (time_integration_method != SimplifiedSpectralDeferredCorrections) {

        burn_success = react_state(S_new, R_new, time - 0.5 * dt, 0.5 * dt, 1);

        if (burn_success != 1) {
            status.success = false;
            status.reason = "burn unsuccessful";

            return status;
        }

        clean_state(
#ifdef MHD
                    Bx_new, By_new, Bz_new,
#endif
                    S_new, time, S_new.nGrow());

    }

#endif // SIMPLIFIED_SDC

    return status;
}

// Strang version

int
Castro::react_state(MultiFab& s, MultiFab& r, Real time, Real dt, const int strang_half)
{

    amrex::ignore_unused(time);

    BL_PROFILE("Castro::react_state()");

    // Sanity check: should only be in here if we're doing CTU.

    if (time_integration_method != CornerTransportUpwind) {
        amrex::Error("Strang reactions are only supported for the CTU advance.");
    }

    const Real strt_time = ParallelDescriptor::second();

    // Start off by assuming a successful burn.

    int burn_success = 1;

    if (do_react != 1) {

        // Ensure we always have valid data, even if we don't do the burn.
        r.setVal(0.0, r.nGrow());

        return burn_success;

    }

    // Check if we have any zones to burn.

    if (!valid_zones_to_burn(s)) {

        // Ensure we always have valid data, even if we don't do the burn.
        r.setVal(0.0, r.nGrow());

        return burn_success;

    }

    const int ng = s.nGrow();

    if (verbose) {
        amrex::Print() << "... Entering burner on level " << level << " and doing half-timestep of burning." << std::endl << std::endl;
    }

    // If we're not subcycling, we only need to do the burn on leaf cells.

    bool mask_covered_zones = false;

    if (level < parent->finestLevel() && parent->subcyclingMode() == "None") {
        mask_covered_zones = true;
    }

    MultiFab tmp_mask_mf;
    const MultiFab& mask_mf = mask_covered_zones ? getLevel(level+1).build_fine_mask() : tmp_mask_mf;

#if defined(AMREX_USE_GPU)
    Gpu::Buffer<int> d_num_failed({0});
    auto* p_num_failed = d_num_failed.data();
#endif
    int num_failed = 0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:num_failed)
#endif
    for (MFIter mfi(s, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        const Box& bx = mfi.growntilebox(ng);

        auto U = s.array(mfi);
        auto reactions = r.array(mfi);
        auto weights = store_burn_weights ? burn_weights.array(mfi) : Array4<Real>{};
        Array4<Real> empty_arr{};
        const auto& mask = mask_covered_zones ? mask_mf.array(mfi) : empty_arr;

        const auto dx = geom.CellSizeArray();
#ifdef MODEL_PARSER
        const auto problo = geom.ProbLoArray();
        const auto geomdata = geom.data();
#endif

#if defined(AMREX_USE_GPU)
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
#else
        LoopOnCpu(bx, [&] (int i, int j, int k)
#endif
        {

            burn_t burn_state;
#ifdef NSE_NET
            burn_state.mu_p = U(i,j,k,UMUP);
            burn_state.mu_n = U(i,j,k,UMUN);

            burn_state.y_e = -1.0_rt;
#endif

#if AMREX_SPACEDIM == 1
            burn_state.dx = dx[0];
#else
            burn_state.dx = amrex::min(AMREX_D_DECL(dx[0], dx[1], dx[2]));
#endif

            // Initialize some data for later.

            bool do_burn = true;
            burn_state.success = true;
            int burn_failed = 0;

            // Don't burn on zones inside shock regions, if the relevant option is set.

#ifdef SHOCK_VAR
            if (U(i,j,k,USHK) > 0.0_rt && disable_shock_burning == 1) {
                do_burn = false;
            }
#endif
            // Don't burn on zones that are masked out.

            if (mask_covered_zones && mask.contains(i,j,k)) {
                if (mask(i,j,k) == 0.0_rt) {
                    do_burn = false;
                }
            }

            Real rhoInv = 1.0_rt / U(i,j,k,URHO);

            burn_state.rho = U(i,j,k,URHO);

            // e is used as an input for some NSE solvers

            burn_state.e = U(i,j,k,UEINT) * rhoInv;

            // this T is consistent with UEINT because we did an EOS call before
            // calling this function

            burn_state.T = U(i,j,k,UTEMP);

            burn_state.T_fixed = -1.e30_rt;

#ifdef MODEL_PARSER
            if (drive_initial_convection) {
                GpuArray<Real, 3> rr = {0.0_rt};

                rr[0] = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt) - problem::center[0];
#if AMREX_SPACEDIM >= 2
                rr[1] = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt) - problem::center[1];
#endif
#if AMREX_SPACEDIM == 3
                rr[2] = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt) - problem::center[2];
#endif

#if DIM_MODEL == 1
                Real dist;

                if (domain_is_plane_parallel) {
                    dist = rr[AMREX_SPACEDIM-1];
                } else {
                    dist = distance(geomdata, rr);
                }

                burn_state.T_fixed = interpolate(dist, model::itemp);
#elif DIM_MODEL == 2
                burn_state.T_fix = interpolate(rr[0], rr[1], model::itemp);
            }
#endif

            for (int n = 0; n < NumSpec; ++n) {
                burn_state.xn[n] = U(i,j,k,UFS+n) * rhoInv;
            }

#if NAUX_NET > 0
            for (int n = 0; n < NumAux; ++n) {
                burn_state.aux[n] = U(i,j,k,UFX+n) * rhoInv;
            }
#endif

            // Ensure we start with no RHS or Jacobian calls registered.

            burn_state.n_rhs = 0;
            burn_state.n_jac = 0;

            // for diagnostics

            burn_state.i = i;
            burn_state.j = j;
            burn_state.k = k;

#ifdef NONAKA_PLOT
            burn_state.level = level;
            burn_state.reference_time = time;
#ifdef STRANG
            burn_state.strang_half = strang_half;
#endif
#endif

            // Don't burn if we're outside of the relevant (rho, T) range.

            if (burn_state.T < castro::react_T_min || burn_state.T > castro::react_T_max ||
                burn_state.rho < castro::react_rho_min || burn_state.rho > castro::react_rho_max) {
                do_burn = false;
            }

            if (do_burn) {
                burner(burn_state, dt);

                // If we were unsuccessful, update the failure count.

                if (!burn_state.success) {
                    burn_failed = 1;
                }

                // Add burning rates to reactions MultiFab, but be
                // careful because the reactions and state MFs may
                // not have the same number of ghost cells.

                if (reactions.contains(i,j,k)) {

                    reactions(i,j,k,0) = (U(i,j,k,URHO) * burn_state.e - U(i,j,k,UEINT)) / dt;

                    if (store_omegadot == 1) {
                        if (reactions.contains(i,j,k)) {
                            for (int n = 0; n < NumSpec; ++n) {
                                reactions(i,j,k,1+n) = U(i,j,k,URHO) * (burn_state.xn[n] - U(i,j,k,UFS+n) * rhoInv) / dt;
                            }
#if NAUX_NET > 0
                            for (int n = 0; n < NumAux; ++n) {
                                reactions(i,j,k,1+n+NumSpec) = U(i,j,k,URHO) * (burn_state.aux[n] - U(i,j,k,UFX+n) * rhoInv) / dt;
                            }
#endif
                        }
                    }

                    if (store_burn_weights) {

                        if (integrator_rp::jacobian == 1) {
                            weights(i,j,k,strang_half) = amrex::max(1.0_rt, static_cast<Real>(burn_state.n_rhs + 2 * burn_state.n_jac));
                        } else {
                            // the RHS evals for the numerical differencing in the Jacobian are already accounted for in n_rhs
                            weights(i,j,k,strang_half) = amrex::max(1.0_rt, static_cast<Real>(burn_state.n_rhs));
                        }
                    }
#ifdef NSE
                    if (store_omegadot == 1) {
                        reactions(i,j,k,NumSpec+NumAux+1) = burn_state.nse;
                    }
                    else {
                        reactions(i,j,k,1) = burn_state.nse;
                    }
#endif
                }

                // update the state
#ifdef NSE_NET
                U(i,j,k,UMUP) = burn_state.mu_p;
                U(i,j,k,UMUN) = burn_state.mu_n;
#endif
                for (int n = 0; n < NumSpec; ++n) {
                    U(i,j,k,UFS+n) = U(i,j,k,URHO) * burn_state.xn[n];
                }
#if NAUX_NET > 0
                for (int n = 0; n < NumAux; ++n) {
                    U(i,j,k,UFX+n) = U(i,j,k,URHO) * burn_state.aux[n];
                }
#endif
                Real reint_old = U(i,j,k,UEINT);
                U(i,j,k,UEINT) = U(i,j,k,URHO) * burn_state.e;
                U(i,j,k,UEDEN) += U(i,j,k,UEINT) - reint_old;

            } else {  // do_burn = false

                if (reactions.contains(i,j,k)) {
                    for (int n = 0; n < reactions.nComp(); n++) {
                        reactions(i,j,k,n) = 0.0_rt;
                    }
                }

            }

#if defined(AMREX_USE_GPU)
            if (burn_failed) {
                Gpu::Atomic::Add(p_num_failed, burn_failed);
            }
#else
            num_failed += burn_failed;
#endif
        });

#if defined(AMREX_USE_HIP)
        Gpu::streamSynchronize(); // otherwise HIP may fail to allocate the necessary resources.
#endif

#ifdef ALLOW_GPU_PRINTF
        std::fflush(nullptr);
#endif

    }

#if defined(AMREX_USE_GPU)
    num_failed = *(d_num_failed.copyToHost());
#endif

    burn_success = !num_failed;

    ParallelDescriptor::ReduceIntMin(burn_success);

    if (print_update_diagnostics) {

        Real e_added = r.sum(0);

        if (e_added != 0.0) {
            amrex::Print() << "... (rho e) added from burning: " << e_added * dt << std::endl << std::endl;
        }

    }

    if (verbose) {
        amrex::Print() << "... Leaving burner on level " << level << " after completing half-timestep of burning." << std::endl << std::endl;
    }

    if (verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        amrex::Real run_time = ParallelDescriptor::second() - strt_time;
        amrex::Real llevel = level;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        amrex::Print() << "Castro::react_state() time = " << run_time
                       << " on level " << llevel << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }

    return burn_success;

}

#ifdef SIMPLIFIED_SDC
// Simplified SDC version

int
Castro::react_state(Real time, Real dt)
{

    amrex::ignore_unused(time);

    // The goal is to update S_old to S_new with the effects of both
    // advection and reactions.  We come into this S_new having seen
    // the effects of advection and sources.  We create an advective
    // update of the form: -div{F} + 0.5 (S^n + S^{n+1} and pass this
    // to the reaction integrator where it is applied together with
    // the reactions to update the full state.

    // Note: S_new actually is already updated with just advection, so
    // in the event that we do not react on a zone (e.g., because it
    // doesn't meet the thermodynamic thresholds) we don't have to do
    // anything.  If we do react, then we overwrite what is stored in
    // S_new with the combined effects of advection and reactions.

    BL_PROFILE("Castro::react_state()");

    // Sanity check: should only be in here if we're doing simplified SDC.

    if (time_integration_method != SimplifiedSpectralDeferredCorrections) {
        amrex::Error("This react_state interface is only supported for simplified SDC.");
    }

    const Real strt_time = ParallelDescriptor::second();

    if (verbose) {
        amrex::Print() << "... Entering burner on level " << level << " and doing full timestep of burning." << std::endl << std::endl;
    }

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

#ifdef MHD
    MultiFab& Bx_new = get_new_data(Mag_Type_x);
    MultiFab& By_new = get_new_data(Mag_Type_y);
    MultiFab& Bz_new = get_new_data(Mag_Type_z);
#endif

    const int ng = S_new.nGrow();

    MultiFab& reactions = get_new_data(Reactions_Type);
    MultiFab& SDC_react = get_new_data(Simplified_SDC_React_Type);

    reactions.setVal(0.0, reactions.nGrow());

    // If we're not subcycling, we only need to do the burn on leaf cells.

    bool mask_covered_zones = false;

    if (level < parent->finestLevel() && parent->subcyclingMode() == "None") {
        mask_covered_zones = true;
    }

    MultiFab tmp_mask_mf;
    const MultiFab& mask_mf = mask_covered_zones ? getLevel(level+1).build_fine_mask() : tmp_mask_mf;

    // Start off assuming a successful burn.

    int burn_success = 1;

#if defined(AMREX_USE_GPU)
    Gpu::Buffer<int> d_num_failed({0});
    auto* p_num_failed = d_num_failed.data();
#endif
    int num_failed = 0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:num_failed)
#endif
    for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(ng);

        auto U_old = S_old.array(mfi);
        auto U_new = S_new.array(mfi);
#ifdef MHD
        auto Bx    = Bx_new.array(mfi);
        auto By    = By_new.array(mfi);
        auto Bz    = Bz_new.array(mfi);
#endif
        auto I     = SDC_react.array(mfi);
        auto react_src = reactions.array(mfi);
        auto weights = store_burn_weights ? burn_weights.array(mfi) : Array4<Real>{};
        Array4<Real> empty_arr{};
        const auto& mask = mask_covered_zones ? mask_mf.array(mfi) : empty_arr;

        int lsdc_iteration = sdc_iteration;

        const auto dx = geom.CellSizeArray();
#ifdef MODEL_PARSER
        const auto problo = geom.ProbLoArray();
        const auto geomdata = geom.data();
#endif

#if defined(AMREX_USE_GPU)
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
#else
        LoopOnCpu(bx, [&] (int i, int j, int k)
#endif
        {
            burn_t burn_state;

#if AMREX_SPACEDIM == 1
            burn_state.dx = dx[0];
#else
            burn_state.dx = amrex::min(AMREX_D_DECL(dx[0], dx[1], dx[2]));
#endif

#ifdef NSE_NET
            burn_state.mu_p = U_old(i,j,k,UMUP);
            burn_state.mu_n = U_old(i,j,k,UMUN);

            burn_state.y_e = -1.0_rt;
#endif
            // Initialize some data for later.

            bool do_burn = true;
            burn_state.success = true;
            int burn_failed = 0;

            // Don't burn on zones inside shock regions, if the
            // relevant option is set.

#ifdef SHOCK_VAR
            if (U_new(i,j,k,USHK) > 0.0_rt && disable_shock_burning == 1) {
                do_burn = false;
            }
#endif

            // Don't burn on zones that are masked out.

            if (mask_covered_zones && mask.contains(i,j,k)) {
                if (mask(i,j,k) == 0.0_rt) {
                    do_burn = false;
                }
            }

            // Feed in the old-time state data.
            // this also sets burn_state.{rho,T}

            copy_cons_to_burn_type(i, j, k, U_old, burn_state);

            burn_state.T_fixed = -1.e30_rt;

#ifdef MODEL_PARSER
            if (drive_initial_convection) {
                GpuArray<Real, 3> rr = {0.0_rt};

                rr[0] = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt) - problem::center[0];
#if AMREX_SPACEDIM >= 2
                rr[1] = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt) - problem::center[1];
#endif
#if AMREX_SPACEDIM == 3
                rr[2] = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt) - problem::center[2];
#endif

#if DIM_MODEL == 1
                Real dist;

                if (domain_is_plane_parallel) {
                    dist = rr[AMREX_SPACEDIM-1];
                } else {
                    dist = distance(geomdata, rr);
                }

                burn_state.T_fixed = interpolate(dist, model::itemp);
#elif DIM_MODEL == 2
                burn_state.T_fix = interpolate(rr[0], rr[1], model::itemp);
            }
#endif

            // Don't burn if we're outside of the relevant (rho, T) range.

            if (U_old(i,j,k,UTEMP) < castro::react_T_min || U_old(i,j,k,UTEMP) > castro::react_T_max ||
                U_old(i,j,k,URHO) < castro::react_rho_min || U_old(i,j,k,URHO) > castro::react_rho_max) {
                do_burn = false;
            }

            // Tell the integrator about the non-reacting source terms,
            // which are advective_source + 1/2 (old source + new source).
            //
            // We have the following:
            //
            // U_new:      U^n - dt div{F} + dt/2 (S^n + S^{n+1})
            // U_old:      U^n
            //
            // So we can compute:
            //
            //   (1 / dt) * (U_new - U_old)
            //
            // To get the non-reacting sources:
            //
            //   -div{F} + (1/2) (S^n + S^{n+1})

            Real dtInv = 1.0_rt / dt;

            burn_state.ydot_a[SRHO] = (U_new(i,j,k,URHO) - U_old(i,j,k,URHO)) * dtInv;
            burn_state.ydot_a[SMX] = (U_new(i,j,k,UMX) - U_old(i,j,k,UMX)) * dtInv;
            burn_state.ydot_a[SMY] = (U_new(i,j,k,UMY) - U_old(i,j,k,UMY)) * dtInv;
            burn_state.ydot_a[SMZ] = (U_new(i,j,k,UMZ) - U_old(i,j,k,UMZ)) * dtInv;
            burn_state.ydot_a[SEDEN] = (U_new(i,j,k,UEDEN) - U_old(i,j,k,UEDEN)) * dtInv;
            burn_state.ydot_a[SEINT] = (U_new(i,j,k,UEINT) - U_old(i,j,k,UEINT)) * dtInv;
            for (int n = 0; n < NumSpec; n++) {
                burn_state.ydot_a[SFS+n] = (U_new(i,j,k,UFS+n) - U_old(i,j,k,UFS+n)) * dtInv;
            }
#if NAUX_NET > 0
            for (int n = 0; n < NumAux; n++) {
                burn_state.ydot_a[SFX+n] = (U_new(i,j,k,UFX+n) - U_old(i,j,k,UFX+n)) * dtInv;
            }
#endif

            // Convert the current state to primitive data.
            // This state is U* = U_old + dt A where A = -div U + S_hydro.

            Array1D<Real, 0, NQ-1> q_noreact;
            Array1D<Real, 0, NQAUX-1> qaux_dummy;

            hydro::conservative_to_primitive(i, j, k, U_new,
#ifdef MHD
                                             Bx, By, Bz,
#endif
                                             q_noreact, qaux_dummy, q_noreact.len() == NQ);

            // dual energy formalism: in doing EOS calls in the burn,
            // switch between e and (E - K) depending on (E - K) / E.

            burn_state.i = i;
            burn_state.j = j;
            burn_state.k = k;

#ifdef NONAKA_PLOT
            burn_state.level = level;
            burn_state.reference_time = time;
#endif

            burn_state.sdc_iter = lsdc_iteration;
            burn_state.num_sdc_iters = sdc_iters;

            if (do_burn) {
                burner(burn_state, dt);

                // If we were unsuccessful, update the failure count.

                if (!burn_state.success) {
                    burn_failed = 1;
                }

                // update the state data.
#ifdef NSE_NET
                U_new(i,j,k,UMUP) = burn_state.mu_p;
                U_new(i,j,k,UMUN) = burn_state.mu_n;
#endif
                U_new(i,j,k,UEDEN) = burn_state.y[SEDEN];
                U_new(i,j,k,UEINT) = burn_state.y[SEINT];
                for (int n = 0; n < NumSpec; n++) {
                    U_new(i,j,k,UFS+n) = burn_state.y[SFS+n];
                }
#if NAUX_NET > 0
                for (int n = 0; n < NumAux; n++) {
                    U_new(i,j,k,UFX+n) = burn_state.y[SFX+n];
                }
#endif

                if (react_src.contains(i,j,k)) {
                    // store the reaction data for the plotfile.
                    // Note, we want to just capture the reaction
                    // portion here, so we subtract off the advective
                    // part.

                    // rho enuc
                    react_src(i,j,k,0) = (U_new(i,j,k,UEINT) - U_old(i,j,k,UEINT)) / dt - burn_state.ydot_a[SEINT];

                    if (store_omegadot) {
                        // rho omegadot_k
                        for (int n = 0; n < NumSpec; ++n) {
                            react_src(i,j,k,1+n) = (U_new(i,j,k,UFS+n) - U_old(i,j,k,UFS+n)) / dt - burn_state.ydot_a[SFS+n];
                        }
#if NAUX_NET > 0
                        // rho auxdot_k
                        for (int n = 0; n < NumAux; ++n) {
                            react_src(i,j,k,1+n+NumSpec) = (U_new(i,j,k,UFX+n) - U_old(i,j,k,UFX+n)) / dt - burn_state.ydot_a[SFX+n];
                        }
#endif
                    }

                    // burn weights

                    if (store_burn_weights) {

                         if (integrator_rp::jacobian == 1) {
                             weights(i,j,k,lsdc_iteration) = amrex::max(1.0_rt, static_cast<Real>(burn_state.n_rhs + 2 * burn_state.n_jac));
                         } else {
                             // the RHS evals for the numerical differencing in the Jacobian are already accounted for in n_rhs
                             weights(i,j,k,lsdc_iteration) = amrex::max(1.0_rt, static_cast<Real>(burn_state.n_rhs));
                         }
                     }
#ifdef NSE
                    if (store_omegadot == 1) {
                        react_src(i,j,k,NumSpec+NumAux+1) = burn_state.nse;
                    }
                    else {
                        react_src(i,j,k,1) = burn_state.nse;
                    }
#endif
                 }

            }

            // Convert the updated state (with the contribution from burning) to primitive data.

            Array1D<Real, 0, NQ-1> q_new;

            hydro::conservative_to_primitive(i, j, k, U_new,
#ifdef MHD
                                             Bx, By, Bz,
#endif
                                             q_new, qaux_dummy, q_new.len() == NQ);

            // Compute the reaction source term.

            // I_q = (q^{n+1} - q^n) / dt - A(q)
            //
            // but A(q) = (q* - q^n) / dt -- that's the effect of advection w/o burning
            //
            // so I_q = (q^{n+1} - q*) / dt

            if (castro::add_sdc_react_source_to_advection) {
                for (int n = 0; n < NQ; ++n) {
                    I(i,j,k,n) = (q_new(n) - q_noreact(n)) * dtInv;
                }
            }
            else {
                for (int n = 0; n < NQ; ++n) {
                    I(i,j,k,n) = 0.0_rt;
                }
            }

#if defined(AMREX_USE_GPU)
            if (burn_failed) {
                Gpu::Atomic::Add(p_num_failed, burn_failed);
            }
#else
            num_failed += burn_failed;
#endif
        });

#if defined(AMREX_USE_HIP)
        Gpu::streamSynchronize(); // otherwise HIP may fail to allocate the necessary resources.
#endif

#ifdef ALLOW_GPU_PRINTF
       std::fflush(nullptr);
#endif

    }

#if defined(AMREX_USE_GPU)
    num_failed = *(d_num_failed.copyToHost());
#endif

    burn_success = !num_failed;

    ParallelDescriptor::ReduceIntMin(burn_success);

    if (ng > 0) {
        S_new.FillBoundary(geom.periodicity());
    }

    Real cur_time = get_state_data(Simplified_SDC_React_Type).curTime();
    AmrLevel::FillPatch(*this, SDC_react, SDC_react.nGrow(), cur_time, Simplified_SDC_React_Type, 0, SDC_react.nComp());

    if (print_update_diagnostics) {

        Real e_added = reactions.sum(0);

        if (e_added != 0.0)
            amrex::Print() << "... (rho e) added from burning: " << e_added * dt << std::endl << std::endl;

    }

    if (verbose) {

        amrex::Print() << "... Leaving burner on level " << level << " after completing full timestep of burning." << std::endl << std::endl;

        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        amrex::Real run_time = ParallelDescriptor::second() - strt_time;
        amrex::Real llevel = level;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time, IOProc);

        amrex::Print() << "Castro::react_state() time = " << run_time
                       << " on level " << llevel << std::endl << std::endl;
#ifdef BL_LAZY
        });
#endif

    }

    return burn_success;

}
#endif


bool
Castro::valid_zones_to_burn(MultiFab& State)
{

    // The default values of the limiters are 0 and 1.e200, respectively.

    Real small = 1.e-10;
    Real large = 1.e199;

    // Check whether we are limiting on either rho or T.

    bool limit_small_rho = react_rho_min >= small;
    bool limit_large_rho = react_rho_max <= large;

    bool limit_rho = limit_small_rho || limit_large_rho;

    bool limit_small_T = react_T_min >= small;
    bool limit_large_T = react_T_max <= large;

    bool limit_T = limit_small_T || limit_large_T;

    bool limit = limit_rho || limit_T;

    if (!limit) {
      return true;
    }

    // Now, if we're limiting on rho, collect the
    // minimum and/or maximum and compare.

    amrex::Vector<Real> small_limiters;
    amrex::Vector<Real> large_limiters;

    bool local = true;

    Real smalldens = small;
    Real largedens = large;

    if (limit_small_rho) {
      smalldens = State.min(URHO, 0, local);
      small_limiters.push_back(smalldens);
    }

    if (limit_large_rho) {
      largedens = State.max(URHO, 0, local);
      large_limiters.push_back(largedens);
    }

    Real small_T = small;
    Real large_T = large;

    if (limit_small_T) {
      small_T = State.min(UTEMP, 0, local);
      small_limiters.push_back(small_T);
    }

    if (limit_large_T) {
      large_T = State.max(UTEMP, 0, local);
      large_limiters.push_back(large_T);
    }

    // Now do the reductions. We're being careful here
    // to limit the amount of work and communication,
    // because regularly doing this check only makes sense
    // if it is negligible compared to the amount of work
    // needed to just do the burn as normal.

    int small_size = static_cast<int>(small_limiters.size());

    if (small_size > 0) {
        amrex::ParallelDescriptor::ReduceRealMin(small_limiters.dataPtr(), small_size);

        if (limit_small_rho) {
            smalldens = small_limiters[0];
            if (limit_small_T) {
                small_T = small_limiters[1];
            }
        } else {
            small_T = small_limiters[0];
        }
    }

    int large_size = static_cast<int>(large_limiters.size());

    if (large_size > 0) {
        amrex::ParallelDescriptor::ReduceRealMax(large_limiters.dataPtr(), large_size);

        if (limit_large_rho) {
            largedens = large_limiters[0];
            if (limit_large_T) {
                large_T = large_limiters[1];
            }
        } else {
            large_T = large_limiters[1];
        }
    }

    // Finally check on whether min <= rho <= max
    // and min <= T <= max. The defaults are small
    // and large respectively, so if the limiters
    // are not on, these checks will not be triggered.

    if (largedens >= react_rho_min && smalldens <= react_rho_max &&
        large_T >= react_T_min && small_T <= react_T_max) {
        return true;
    }

    // If we got to this point, we did not survive the limiters,
    // so there are no zones to burn.

    if (verbose > 1) {
        amrex::Print() << "  No valid zones to burn, skipping react_state()." << std::endl;
    }

    return false;

}

#endif
