
#include "Castro.H"
#include "Castro_F.H"

#include "AMReX_DistributionMapping.H"

using std::string;
using namespace amrex;

bool
Castro::strang_react_first_half(Real time, Real dt)
{
    BL_PROFILE("Castro::strang_react_first_half()");

    // Sanity check: should only be in here if we're doing CTU.

    if (time_integration_method != CornerTransportUpwind) {
        amrex::Error("Strang reactions are only supported for the CTU advance.");
    }

    bool burn_success = true;

    // Get the reactions MultiFab to fill in.

    MultiFab& reactions = get_old_data(Reactions_Type);

    if (do_react != 1) {

        // Ensure we always have valid data, even if we don't do the burn.
        reactions.setVal(0.0, reactions.nGrow());

        return burn_success;

    }

    // Get the current state data.

    MultiFab& state_burn = Sborder;

    // Check if we have any zones to burn.

    if (!valid_zones_to_burn(state_burn)) {

        // Ensure we always have valid data, even if we don't do the burn.
        reactions.setVal(0.0, reactions.nGrow());

        return burn_success;

    }

    const int ng = state_burn.nGrow();

    // Reactions are expensive and we would usually rather do a
    // communication step than burn on the ghost zones. So what we
    // will do here is create a mask that indicates that we want to
    // turn on the valid interior zones but NOT on the ghost zones
    // that are interior to the level. However, we DO want to burn on
    // the ghost zones that are on the coarse-fine interfaces, since
    // that is going to be more accurate than interpolating from
    // coarse zones. So we will not mask out those zones, and the
    // subsequent FillBoundary call will not interfere with it.

    iMultiFab& interior_mask = build_interior_boundary_mask(ng);

    if (verbose)
        amrex::Print() << "... Entering burner and doing half-timestep of burning." << std::endl << std::endl;

    burn_success = react_state(state_burn, reactions, interior_mask, time, dt, 1, ng);

    if (verbose)
        amrex::Print() << "... Leaving burner after completing half-timestep of burning." << std::endl << std::endl;

    // Now do a boundary fill so that the masked out zones which were skipped get their valid data.

    state_burn.FillBoundary(geom.periodicity());

    // Ensure consistency in internal energy and recompute temperature.

    clean_state(state_burn, time, state_burn.nGrow());

    return burn_success;

}



bool
Castro::strang_react_second_half(Real time, Real dt)
{
    BL_PROFILE("Castro::strang_react_second_half()");

    // Sanity check: should only be in here if we're doing CTU.

    if (time_integration_method != CornerTransportUpwind) {
        amrex::Error("Strang reactions are only supported for the CTU advance.");
    }

    bool burn_success = true;

    MultiFab& reactions = get_new_data(Reactions_Type);

    if (do_react != 1) {

        // Ensure we always have valid data, even if we don't do the burn.
        reactions.setVal(0.0, reactions.nGrow());

        return burn_success;

    }

    MultiFab& state_burn = get_new_data(State_Type);

    // Check if we have any zones to burn.

    if (!valid_zones_to_burn(state_burn)) {

        // Ensure we always have valid data, even if we don't do the burn.
        reactions.setVal(0.0, reactions.nGrow());

        return burn_success;

    }

    // To be consistent with other source term types,
    // we are only applying this on the interior zones.

    const int ng = 0;

    iMultiFab& interior_mask = build_interior_boundary_mask(ng);

    if (verbose)
        amrex::Print() << "... Entering burner and doing half-timestep of burning." << std::endl << std::endl;

    burn_success = react_state(state_burn, reactions, interior_mask, time, dt, 2, ng);

    if (verbose)
        amrex::Print() << "... Leaving burner after completing half-timestep of burning." << std::endl << std::endl;

    state_burn.FillBoundary(geom.periodicity());

    clean_state(state_burn, time + 0.5 * dt, state_burn.nGrow());

    return burn_success;

}



bool
Castro::react_state(MultiFab& s, MultiFab& r, const iMultiFab& m, Real time, Real dt_react, int strang_half, int ngrow)
{

    BL_PROFILE("Castro::react_state()");

    // Sanity check: should only be in here if we're doing CTU or MOL.

    if (time_integration_method != CornerTransportUpwind) {
        amrex::Error("Strang reactions are only supported for the CTU and MOL advance.");
    }

    const Real strt_time = ParallelDescriptor::second();

    // Start off assuming a successful burn.

    int burn_success = 1;
#ifndef CXX_REACTIONS
    Real burn_failed = 0.0;
#endif

    // If we're not actually doing the burn, do the interpolation now.

    if (level > castro::reactions_max_solve_level) {
        FillCoarsePatch(r, 0, time, Reactions_Type, 0, r.nComp(), r.nGrow());
    }

    FArrayBox U_cc;
    FArrayBox R_cc;

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(s, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        const Box& bx = mfi.tilebox();
        const Box& obx = amrex::grow(bx, 1);

        auto U = s.array(mfi);
        auto reactions = r.array(mfi);
        auto mask = m.array(mfi);

        U_cc.resize(obx, NUM_STATE);
        auto U_cc_arr = U_cc.array();
        auto elix_u_cc = U_cc.elixir();

        R_cc.resize(obx, r.nComp());
        auto R_cc_arr = R_cc.array();
        auto elix_r_cc = R_cc.elixir();

        if (level <= castro::reactions_max_solve_level) {

#ifdef CXX_REACTIONS

          Real lreact_T_min = castro::react_T_min;
          Real lreact_T_max = castro::react_T_max;
          Real lreact_rho_min = castro::react_rho_min;
          Real lreact_rho_max = castro::react_rho_max;

          // convert the cell-averages to centers
          make_cell_center(obx, U, U_cc, domlo, domhi);

          // burn the cell-centers on obx
          reduce_op.eval(obx, reduce_data,
                         [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept -> ReduceTuple
                         {

                           burn_t burn_state;

                           // Initialize some data for later.

                           bool do_burn = true;
                           burn_state.success = true;
                           Real burn_failed = 0.0_rt;

                           // Don't burn on zones that we are intentionally masking out.

                           if (mask(i,j,k) != 1) {
                             do_burn = false;
                           }

                           // Don't burn on zones inside shock regions, if the relevant option is set.

#ifdef SHOCK_VAR
                           if (U(i,j,k,USHK) > 0.0_rt && disable_shock_burning == 1) {
                             do_burn = false;
                           }
#endif

                           Real rhoInv = 1.0_rt / U_cc_arr(i,j,k,URHO);

                           burn_state.rho = U_cc_arr(i,j,k,URHO);
                           burn_state.T   = U_cc_arr(i,j,k,UTEMP);
                           burn_state.e   = 0.0_rt; // Energy generated by the burn

                           for (int n = 0; n < NumSpec; ++n) {
                             burn_state.xn[n] = U_cc_arr(i,j,k,UFS+n) * rhoInv;
                           }

#if naux > 0
                           for (int n = 0; n < NumAux; ++n) {
                             burn_state.aux[n] = U_cc_arr(i,j,k,UFX+n) * rhoInv;
                           }
#endif

                           // Ensure we start with no RHS or Jacobian calls registered.

                           burn_state.n_rhs = 0;
                           burn_state.n_jac = 0;

                           // Don't burn if we're outside of the relevant (rho, T) range.

                           if (burn_state.T < lreact_T_min || burn_state.T > lreact_T_max ||
                               burn_state.rho < lreact_rho_min || burn_state.rho > lreact_rho_max) {
                             do_burn = false;
                           }

                           if (do_burn) {
                             burner(burn_state, dt_react);
                           }

                           // If we were unsuccessful, update the failure count.

                           if (!burn_state.success) {
                             burn_failed = 1.0_rt;
                           }

                           if (do_burn) {

                             // Note that we want to update the total energy by taking
                             // the difference of the old rho*e and the new rho*e. If
                             // the user wants to ensure that rho * E = rho * e + rho *
                             // K, this reset should be enforced through an appropriate
                             // choice for the dual energy formalism parameter
                             // dual_energy_eta2 in reset_internal_energy.

                             Real delta_e     = burn_state.e;
                             Real delta_rho_e = burn_state.rho * delta_e;

                             // Add burning rates to reactions MultiFab, but be
                             // careful because the reactions and state MFs may
                             // not have the same number of ghost cells. Note that
                             // we must do this before we actually update the state
                             // since we have not saved the old state.

                             if (reactions.contains(i,j,k)) {
                               // temporarily store just (rho X)^new here -- we'll correct it later
                               for (int n = 0; n < NumSpec; ++n) {
                                 R_cc_arr(i,j,k,n) = U_cc_arr(i,j,k,URHO) * burn_state.xn[n];
                               }
                               R_cc_arr(i,j,k,NumSpec  ) = delta_e / dt_react;
                               R_cc_arr(i,j,k,NumSpec+1) = delta_rho_e / dt_react;

                               reactions(i,j,k,NumSpec+2) = amrex::max(1.0_rt, static_cast<Real>(burn_state.n_rhs + 2 * burn_state.n_jac));
                             }

#if naux > 0
                             for (int n = 0; n < NumSpec; ++n) {
                               // note: this isn't a 4th order average
                               U(i,j,k,UFX+n)  = U(i,j,k,URHO) * burn_state.aux[n];
                             }
#endif

                           }
                           else {

                             if (reactions.contains(i,j,k)) {
                               for (int n = 0; n < NumSpec+2; ++n) {
                                 reactions(i,j,k,n) = 0.0_rt;
                               }

                               reactions(i,j,k,NumSpec+2) = 1.0_rt;
                             }

                           }

                           return {burn_failed};

                         });

            // now convert R to cell-average
            if (do_burn) {
              make_fourth_average(bx, reactions, R_cc_arr, domlo, domhi);
            }


            // finally update the state with the cell-averages and store in the reaction MF
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {

              U(i,j,k,UEINT) += dt_react * reactions(i,j,k,NumSpec+1);
              U(i,j,k,UEDEN) += dt_react * reactions(i,j,k,NumSpec+1);

              for (int n = 0; n < NumSpec; ++n) {
                Real rhoX_old = U(i,j,k,UFS+n);
                U(i,j,k,UFS+n) = reactions(i,j,k,n);
                // create dX/dt -- this is only 2nd order accurate
                reactions(i,j,k,n) = (reactions(i,j,k,n) - rhoX_old)/(U(i,j,k,URHO) * dt_react);
              }
            });

#else

#pragma gpu box(bx)
            ca_react_state(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                           BL_TO_FORTRAN_ANYD(s[mfi]),
                           BL_TO_FORTRAN_ANYD(r[mfi]),
                           BL_TO_FORTRAN_ANYD(m[mfi]),
                           time, dt_react, strang_half,
                           AMREX_MFITER_REDUCE_SUM(&burn_failed));

#endif

        }
        else {

            // Use the interpolated reactions data to update the state.

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
            {
                if (U.contains(i,j,k) && reactions.contains(i,j,k)) {
                    for (int n = 0; n < NumSpec; ++n) {
                        U(i,j,k,UFS+n) += reactions(i,j,k,n) * U(i,j,k,URHO) * dt_react;
                    }
                    U(i,j,k,UEINT) += reactions(i,j,k,NumSpec) * U(i,j,k,URHO) * dt_react;
                    U(i,j,k,UEDEN) += reactions(i,j,k,NumSpec) * U(i,j,k,URHO) * dt_react;
                }
            });

        }

    }

#ifdef CXX_REACTIONS
    ReduceTuple hv = reduce_data.value();
    Real burn_failed = amrex::get<0>(hv);
#endif

    if (burn_failed != 0.0) burn_success = 0;

    ParallelDescriptor::ReduceIntMin(burn_success);

    if (print_update_diagnostics) {

        Real e_added = r.sum(NumSpec + 1);

        if (e_added != 0.0)
            amrex::Print() << "... (rho e) added from burning: " << e_added << std::endl << std::endl;

    }

    if (verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        amrex::Print() << "Castro::react_state() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }

    // fill ghost cells since we only computed the average in bx
    AmrLevel::FillPatch(*this, s, s.nGrow(), time, State_Type, 0, NUM_STATE);


    if (burn_success)
        return true;
    else
        return false;

}

// Simplified SDC version

bool
Castro::react_state(Real time, Real dt)
{
    BL_PROFILE("Castro::react_state()");

    // Sanity check: should only be in here if we're doing simplified SDC.

    if (time_integration_method != SimplifiedSpectralDeferredCorrections) {
        amrex::Error("This react_state interface is only supported for simplified SDC.");
    }

    const Real strt_time = ParallelDescriptor::second();

    if (verbose)
        amrex::Print() << "... Entering burner and doing full timestep of burning." << std::endl << std::endl;

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    // Build the burning mask, in case the state has ghost zones.

    const int ng = S_new.nGrow();
    const iMultiFab& interior_mask = build_interior_boundary_mask(ng);

    // Create a MultiFab with all of the non-reacting source terms.

    MultiFab A_src(grids, dmap, NUM_STATE, ng);
    sum_of_sources(A_src);

    MultiFab& reactions = get_new_data(Reactions_Type);

    reactions.setVal(0.0, reactions.nGrow());

    // Start off assuming a successful burn.

    int burn_success = 1;
    Real burn_failed = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:burn_failed)
#endif
    for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        const Box& bx = mfi.growntilebox(ng);

        FArrayBox& uold    = S_old[mfi];
        FArrayBox& unew    = S_new[mfi];
        FArrayBox& a       = A_src[mfi];
        FArrayBox& r       = reactions[mfi];
        const IArrayBox& m = interior_mask[mfi];

#pragma gpu box(bx)
        ca_react_state_simplified_sdc(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                                      BL_TO_FORTRAN_ANYD(uold),
                                      BL_TO_FORTRAN_ANYD(unew),
                                      BL_TO_FORTRAN_ANYD(a),
                                      BL_TO_FORTRAN_ANYD(r),
                                      BL_TO_FORTRAN_ANYD(m),
                                      time, dt, sdc_iteration,
                                      AMREX_MFITER_REDUCE_SUM(&burn_failed));

    }

    if (burn_failed != 0.0) burn_success = 0;

    ParallelDescriptor::ReduceIntMin(burn_success);

    if (ng > 0)
        S_new.FillBoundary(geom.periodicity());

    if (print_update_diagnostics) {

        Real e_added = reactions.sum(NumSpec + 1);

        if (e_added != 0.0)
            amrex::Print() << "... (rho e) added from burning: " << e_added << std::endl << std::endl;

    }

    if (verbose) {

        amrex::Print() << "... Leaving burner after completing full timestep of burning." << std::endl << std::endl;

        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time, IOProc);

        amrex::Print() << "Castro::react_state() time = " << run_time << std::endl << std::endl;
#ifdef BL_LAZY
        });
#endif

    }

    // For the ca_check_timestep routine, we need to have both the old
    // and new burn defined, so we simply do a copy here
    MultiFab& R_old = get_old_data(Reactions_Type);
    MultiFab& R_new = get_new_data(Reactions_Type);
    MultiFab::Copy(R_new, R_old, 0, 0, R_new.nComp(), R_new.nGrow());


    if (burn_success)
        return true;
    else
        return false;

}



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

    if (!limit) return true;

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

    int small_size = small_limiters.size();

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

    int large_size = large_limiters.size();

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

    if (verbose > 1)
        amrex::Print() << "  No valid zones to burn, skipping react_state()." << std::endl;

    return false;

}
