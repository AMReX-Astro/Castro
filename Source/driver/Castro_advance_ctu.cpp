
#include <Castro.H>
#include <Castro_F.H>

#ifdef RADIATION
#include <Radiation.H>
#endif

#ifdef GRAVITY
#include <Gravity.H>
#endif

using namespace amrex;

advance_status
Castro::do_advance_ctu(Real time,
                       Real dt,
                       int  amr_iteration,
                       int  amr_ncycle)
{

    amrex::ignore_unused(amr_iteration);
    amrex::ignore_unused(amr_ncycle);

    // this routine will advance the old state data (called S_old here)
    // to the new time, for a single level.  The new data is called
    // S_new here.  The update includes reactions (if we are not doing
    // SDC), hydro, and the source terms.


    BL_PROFILE("Castro::do_advance_ctu()");

    advance_status status;
    status.success = true;
    status.reason = "";

#ifndef TRUE_SDC

    const Real prev_time = state[State_Type].prevTime();
    const Real  cur_time = state[State_Type].curTime();

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

#ifdef MHD
    MultiFab& Bx_old = get_old_data(Mag_Type_x);
    MultiFab& By_old = get_old_data(Mag_Type_y);
    MultiFab& Bz_old = get_old_data(Mag_Type_z);

    MultiFab& Bx_new = get_new_data(Mag_Type_x);
    MultiFab& By_new = get_new_data(Mag_Type_y);
    MultiFab& Bz_new = get_new_data(Mag_Type_z);
#endif 

    // Perform initialization steps.

    initialize_do_advance(time);

    // Create any correctors to the source term data. This must be done
    // before the source term data is overwritten below. Note: we do
    // not create the corrector source if we're currently retrying the
    // step; we will already have done it, and aside from avoiding
    // duplicate work, we have already lost the data needed to do this
    // calculation since we overwrote the data from the previous step.

    if (!in_retry) {
        create_source_corrector();
    }

    // Check for NaN's.

    check_for_nan(S_old);

    // If we're doing a step later than the first on each level, the fluid
    // state might have evolved to the point where the AMR timestep could be
    // significantly too large, but we don't have freedom to adjust the AMR
    // timestep at that point. Trying to evolve with a dt that is too large
    // could result in catastrophic behavior such that we don't even get to
    // the point where we can bail out later in the advance, so let's just
    // go directly into a retry now if we're too far away from the needed dt.

    bool is_first_step_on_this_level = true;

    for (int lev = level; lev >= 0; --lev) {
        if (getLevel(lev).iteration > 1) {
            is_first_step_on_this_level = false;
            break;
        }
    }

    if (castro::check_dt_before_advance && !is_first_step_on_this_level) {

        int is_new = 0;
        Real old_dt = estTimeStep(is_new);

        if (castro::change_max * old_dt < dt) {
            status.success = false;
            status.reason = "pre-advance timestep validity check failed";
            return status;
        }

    }

    // Since we are Strang splitting the reactions, do them now

#ifdef REACTIONS
    bool burn_success = true;

    MultiFab& R_old = get_old_data(Reactions_Type);
    MultiFab& R_new = get_new_data(Reactions_Type);

    if (time_integration_method != SimplifiedSpectralDeferredCorrections) {

        // The result of the reactions is added directly to Sborder.
        burn_success = react_state(Sborder, R_old, prev_time, 0.5 * dt, 0);
        clean_state(
#ifdef MHD
                    Bx_old_tmp, By_old_tmp, Bz_old_tmp,
#endif
                    Sborder, prev_time, Sborder.nGrow());

    }
#endif

    // Initialize the new-time data. This copy needs to come after the
    // reactions.

    MultiFab::Copy(S_new, Sborder, 0, 0, NUM_STATE, S_new.nGrow());

#ifdef REACTIONS
    if (time_integration_method != SimplifiedSpectralDeferredCorrections) {

        // Do this for the reactions as well, in case we cut the timestep
        // short due to it being rejected.

        MultiFab::Copy(R_new, R_old, 0, 0, R_new.nComp(), R_new.nGrow());

        // Skip the rest of the advance if the burn was unsuccessful.

        if (!burn_success) {
            status.success = false;
            status.reason = "first Strang burn unsuccessful";
            return status;
        }

    }
#endif

    // Construct the old-time sources from Sborder.  This will already
    // be applied to S_new (with full dt weighting), to be correctly
    // later.  Note -- this does not affect the prediction of the
    // interface state, an explict source will be traced there as
    // needed.

#ifdef GRAVITY
    construct_old_gravity(amr_iteration, amr_ncycle, prev_time);
#endif

    bool apply_sources_to_state = true;

    MultiFab& old_source = get_old_data(Source_Type);

    if (apply_sources()) {

      do_old_sources(
#ifdef MHD
                      Bx_old, By_old, Bz_old,
#endif                
                      old_source, Sborder, S_new, prev_time, dt, apply_sources_to_state);

      if (do_hydro) {
          // Fill the ghost cells of old_source / Source_Type

          AmrLevel::FillPatch(*this, old_source, old_source.nGrow(), prev_time, Source_Type, 0, NSRC);
      }


    } else {
      old_source.setVal(0.0, NUM_GROW_SRC);

    }


#ifdef SIMPLIFIED_SDC
#ifdef REACTIONS
    // the SDC reactive source ghost cells on coarse levels might not
    // be in sync due to any average down done, so fill them here

    MultiFab& react_src = get_new_data(Simplified_SDC_React_Type);

    AmrLevel::FillPatch(*this, react_src, react_src.nGrow(), cur_time, Simplified_SDC_React_Type, 0, react_src.nComp());
#endif
#endif

    // Do the hydro update.  We build directly off of Sborder, which
    // is the state that has already seen the burn

    if (do_hydro)
    {
#ifndef MHD
      construct_ctu_hydro_source(time, dt);

//      if (print_update_diagnostics) {
//          evaluate_and_print_source_change(hydro_source, dt, "hydro source");
//      }
#else
      construct_ctu_mhd_source(time, dt);
#endif

      // Check for small/negative densities and X > 1 or X < 0.
      // If we detect this, return immediately.

      ReduceOps<ReduceOpMax, ReduceOpMax> reduce_op;
      ReduceData<int, int> reduce_data(reduce_op);
      using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.tilebox();

          auto S_old_arr = S_old.array(mfi);
          auto S_new_arr = S_new.array(mfi);

          reduce_op.eval(bx, reduce_data,
          [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
          {
              int rho_check_failed = 0;
              int X_check_failed = 0;

              Real rho = S_new_arr(i,j,k,URHO);
              Real rhoInv = 1.0_rt / rho;

              // Optionally, the user can ignore this if the starting
              // density is lower than a certain threshold. This is useful
              // if the minimum density occurs in material that is not
              // dynamically important; in that case, a density reset suffices.

              if (S_old_arr(i,j,k,URHO) >= retry_small_density_cutoff && rho < small_dens) {
#ifndef AMREX_USE_GPU
                  std::cout << "Invalid density = " << rho << " at index " << i << ", " << j << ", " << k << "\n";
#endif
                  rho_check_failed = 1;
              }

              if (S_new_arr(i,j,k,URHO) >= castro::abundance_failure_rho_cutoff) {

                  for (int n = 0; n < NumSpec; ++n) {
                      Real X = S_new_arr(i,j,k,UFS+n) * rhoInv;

                      if (X < -castro::abundance_failure_tolerance ||
                          X > 1.0_rt + castro::abundance_failure_tolerance) {
#ifndef AMREX_USE_GPU
                          std::cout << "Invalid X[" << n << "] = " << X << " in zone "
                                    << i << ", " << j << ", " << k
                                    << " with density = " << rho << "\n";
#endif
                          X_check_failed = 1;
                      }
                  }

              }

              return {rho_check_failed, X_check_failed};
          });

      }

      ReduceTuple hv = reduce_data.value();
      int rho_check_failed = amrex::get<0>(hv);
      int X_check_failed = amrex::get<1>(hv);

      ParallelDescriptor::ReduceIntMax(rho_check_failed);
      ParallelDescriptor::ReduceIntMax(X_check_failed);

      if (rho_check_failed == 1) {
          status.success = false;
          status.reason = "invalid density";
          return status;
      }

      if (X_check_failed == 1) {
          status.success = false;
          status.reason = "invalid X";
          return status;
      }
    }


    // Sync up state after old sources and hydro source.
    clean_state(
#ifdef MHD
                Bx_new, By_new, Bz_new,
#endif
                S_new, cur_time, 0);

    // Check for NaN's.

    check_for_nan(S_new);

    // if we are done with the update do the source correction and
    // then the second half of the reactions

#ifdef GRAVITY
    // Must define new value of "center" before we call new gravity
    // solve or external source routine
    if (moving_center == 1)
      define_new_center(S_new, time);
#endif

#ifdef GRAVITY
    // We need to make the new radial data now so that we can use it when we
    // FillPatch in creating the new source.

#if (AMREX_SPACEDIM > 1)
    if ( (level == 0) && (spherical_star == 1) ) {
      int is_new = 1;
      make_radial_data(is_new);
    }
#endif
#endif

    // Construct and apply new-time source terms.

#ifdef GRAVITY
    construct_new_gravity(amr_iteration, amr_ncycle, cur_time);
#endif

    MultiFab& new_source = get_new_data(Source_Type);

    if (apply_sources()) {

      do_new_sources(
#ifdef MHD
                              Bx_new, By_new, Bz_new,
#endif  
                      new_source, Sborder, S_new, cur_time, dt, apply_sources_to_state);

    } else {

      new_source.setVal(0.0, NUM_GROW_SRC);

    }

    // If the state has ghost zones, sync them up now
    // since the hydro source only works on the valid zones.

    if (S_new.nGrow() > 0) {
      clean_state(
#ifdef MHD
                  Bx_new, By_new, Bz_new,
#endif                
                  S_new, cur_time, 0);

      expand_state(S_new, cur_time, S_new.nGrow());
    }

    // Do the second half of the reactions for Strang, or the full burn for simplified SDC.

#ifdef REACTIONS

#ifdef SIMPLIFIED_SDC

    if (time_integration_method == SimplifiedSpectralDeferredCorrections) {

        if (do_react) {

            // Do the ODE integration to capture the reaction source terms.

            burn_success = react_state(time, dt);

            // Skip the rest of the advance if the burn was unsuccessful.

            if (!burn_success) {
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

        burn_success = react_state(S_new, R_new, cur_time - 0.5 * dt, 0.5 * dt, 1);
        clean_state(
#ifdef MHD
                    Bx_new, By_new, Bz_new,
#endif
                    S_new, cur_time, S_new.nGrow());

        // Skip the rest of the advance if the burn was unsuccessful.

        if (!burn_success) {
            status.success = false;
            status.reason = "second Strang burn unsuccessful";
            return status;
        }

    }

#endif // SIMPLIFIED_SDC

#endif // REACTIONS

    // Check if this timestep violated our stability criteria. Our idea is,
    // if the timestep created a velocity v and sound speed at the new time
    // such that (v+c) * dt / dx < CFL / change_max, where CFL is the user's
    // chosen timestep constraint and change_max is the factor that determines
    // how much the timestep can change during an advance, consider the advance
    // to have failed. This prevents the timestep from shrinking too much,
    // whereas in computeNewDt change_max prevents the timestep from growing
    // too much. The same reasoning applies for the other timestep limiters.

    if (castro::check_dt_after_advance) {

        int is_new = 1;
        Real new_dt = estTimeStep(is_new);

        if (castro::change_max * new_dt < dt) {
            status.success = false;
            status.reason = "post-advance timestep validity check failed";
            return status;
        }

    }

    finalize_do_advance();

#endif

    return status;

}





bool
Castro::retry_advance_ctu(Real dt, advance_status status)
{
    BL_PROFILE("Castro::retry_advance_ctu()");

    bool do_retry = false;

    if (!status.success)
        do_retry = true;

    if (do_retry) {

        dt_subcycle = std::min(dt, dt_subcycle) * retry_subcycle_factor;

        if (verbose && ParallelDescriptor::IOProcessor()) {
            std::cout << std::endl;
            std::cout << Font::Bold << FGColor::Red;
            std::cout << "  Timestep " << dt << " rejected at level " << level << "." << std::endl;
            std::cout << "  Performing a retry, with subcycled timesteps of maximum length dt = " << dt_subcycle << std::endl;
            std::cout << ResetDisplay;
            std::cout << std::endl;
        }

        // If we are doing a retry and this is the first attempt
        // at the advance, make a copy of the state data. This will
        // be useful to us at the end of the timestep when we need
        // to restore the original old data.

        save_data_for_retry();

        // Clear the contribution to the fluxes from this step.

        for (int dir = 0; dir < 3; ++dir) {
          fluxes[dir]->setVal(0.0);
        }

        for (int dir = 0; dir < 3; ++dir) {
          mass_fluxes[dir]->setVal(0.0);
        }

#if (AMREX_SPACEDIM <= 2)
        if (!Geom().IsCartesian()) {
          P_radial.setVal(0.0);
        }
#endif

#ifdef RADIATION
        if (Radiation::rad_hydro_combined) {
          for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            rad_fluxes[dir]->setVal(0.0);
          }
        }
#endif

#ifdef REACTIONS
        burn_weights.setVal(0.0);
#endif

        // For simplified SDC, we'll have garbage data if we
        // attempt to use the lagged source terms (both reacting
        // and non-reacting) from the last timestep, since that
        // advance failed and we don't know if we can trust it.
        // So we zero out both source term correctors.

        if (time_integration_method == SimplifiedSpectralDeferredCorrections) {
            source_corrector.setVal(0.0, source_corrector.nGrow());

#ifdef SIMPLIFIED_SDC
#ifdef REACTIONS
            MultiFab& SDC_react_new = get_new_data(Simplified_SDC_React_Type);
            SDC_react_new.setVal(0.0, SDC_react_new.nGrow());
#endif
#endif
        }

    }

    return do_retry;

}



Real
Castro::subcycle_advance_ctu(const Real time, const Real dt, int amr_iteration, int amr_ncycle)
{
    BL_PROFILE("Castro::subcycle_advance_ctu()");

    // Start the subcycle time off with the main dt,
    // unless we already came in here with an estimate
    // that is different from the initial value we assigned,
    // for example from the post-step regrid algorithm.

    if (dt_subcycle == 1.e200)
        dt_subcycle = dt;

    Real subcycle_time = time;

    Real dt_new = 1.e200;

    sub_iteration = 0;

    // Subcycle until we've reached the target time.
    // Compare against a slightly smaller number to
    // avoid roundoff concerns.

    Real eps = 1.0e-14;

    bool do_swap = false;

    Real last_dt_subcycle = 1.e200;

    while (subcycle_time < (1.0 - eps) * (time + dt)) {

        // Save the dt_subcycle before modifying it, we will use it later.

        last_dt_subcycle = dt_subcycle;

        // Shorten the last timestep so that we don't overshoot
        // the ending time. Relatedly, we also don't want to
        // undershoot the ending time, which would cause us to
        // have to take a very short final timestep to get to
        // the desired time. We'll use a slightly larger tolerance
        // factor than we do in the outer while loop, to avoid
        // roundoff issues, under the logic that as long as the
        // tolerance is still small, there is no meaningful harm
        // from the timestep being slightly larger than our recommended
        // safe dt in the subcycling.

        AMREX_ASSERT(dt_cutoff > eps);

        if (subcycle_time + dt_subcycle > (1.0 - dt_cutoff) * (time + dt)) {
          dt_subcycle = (time + dt) - subcycle_time;
        }

        // Determine whether we're below the cutoff timestep.

        if (dt_subcycle <= dt_cutoff * time) {
            if (ParallelDescriptor::IOProcessor()) {
                std::cout << std::endl;
                std::cout << "  The subcycle mechanism requested subcycled timesteps of maximum length dt = " << dt_subcycle << "," << std::endl
                          << "  but this timestep is shorter than the user-defined minimum, " << std::endl
                          << "  castro.dt_cutoff, multiplied by the current time (" << dt_cutoff * time << "). Aborting." << std::endl;
            }
            amrex::Abort("Error: subcycled timesteps too short.");
        }

        // Check on whether we are going to take too many subcycles.

        int num_subcycles_remaining = int(round(((time + dt) - subcycle_time) / dt_subcycle));

        if (num_subcycles_remaining > max_subcycles) {
            amrex::Print() << std::endl
                           << "  The subcycle mechanism requested " << num_subcycles_remaining << " subcycled timesteps, which is larger than the maximum of " << max_subcycles << "." << std::endl
                           << "  If you would like to override this, increase the parameter castro.max_subcycles." << std::endl;
            amrex::Abort("Error: too many subcycles.");
        }

        // If we get to this point, we survived the sanity checks. Print out the current subcycle iteration.

        if (verbose && ParallelDescriptor::IOProcessor()) {
            std::cout << std::endl;
            std::cout << Font::Bold << FGColor::Green << "  Beginning subcycle " << sub_iteration + 1 << " starting at time " << subcycle_time
                      << " with dt = " << dt_subcycle << ResetDisplay << std::endl;
            std::cout << "  Estimated number of subcycles remaining: " << num_subcycles_remaining << std::endl << std::endl;
        }

        // Swap the time levels. Only do this after the first iteration;
        // the first iteration used the swap done in initialize_advance.
        // After that, the only exception will be if we do a retry.

        if (do_swap) {

            swap_state_time_levels(0.0);

#ifdef GRAVITY
            if (do_grav) {
                gravity->swapTimeLevels(level);
            }
#endif

        }
        else {

            do_swap = true;

        }

        // Set the relevant time levels.

        for (int k = 0; k < num_state_type; k++) {
          state[k].setTimeLevel(subcycle_time + dt_subcycle, dt_subcycle, 0.0);
        }

        // Do the advance and construct the relevant source terms. For CTU this
        // will include Strang-split reactions; for simplified SDC, we defer the
        // burn until after the advance.

        int num_sub_iters = 1;

        // Save the number of SDC iterations in case we are about to modify it.

        int sdc_iters_old = sdc_iters;

        if (time_integration_method == SimplifiedSpectralDeferredCorrections) {
            // If we're in a retry we want to add one more iteration. This is
            // because we will have zeroed out the lagged corrector from the
            // last iteration/timestep and so having another iteration will
            // approximately compensate for that.

            if (in_retry) {
                sdc_iters += 1;
                amrex::Print() << "Adding an SDC iteration due to the retry." << std::endl << std::endl;
            }

            num_sub_iters = sdc_iters;
        }

        advance_status status;

        for (int n = 0; n < num_sub_iters; ++n) {

            if (time_integration_method == SimplifiedSpectralDeferredCorrections) {
                sdc_iteration = n;
                amrex::Print() << "Beginning SDC iteration " << n + 1 << " of " << num_sub_iters << "." << std::endl << std::endl;
            }

            // We do the hydro advance here, and record whether we completed it.

            status = do_advance_ctu(subcycle_time, dt_subcycle, amr_iteration, amr_ncycle);

            if (in_retry) {
                in_retry = false;
            }

            if (!status.success) {
                if (use_retry) {
                    amrex::Print() << "Advance was unsuccessful with reason: " << status.reason << "; proceeding to a retry." << std::endl << std::endl;
                } else {
                    amrex::Abort("Advance was unsuccessful.");
                }
                break;
            }

            if (time_integration_method == SimplifiedSpectralDeferredCorrections) {
                amrex::Print() << "Ending SDC iteration " << n + 1 << " of " << num_sub_iters << "." << std::endl << std::endl;
            }

        }

        if (verbose && ParallelDescriptor::IOProcessor()) {
            std::cout << Font::Bold << FGColor::Green << "  Subcycle completed" << ResetDisplay << std::endl << std::endl;
        }

        // Set sdc_iters to its original value, in case we modified it above.

        sdc_iters = sdc_iters_old;

        // If we have hit a CFL violation during this subcycle, we must abort.

        if (cfl_violation && !use_retry) {
          amrex::Abort("CFL is too high at this level; go back to a checkpoint and restart with lower CFL number, or set castro.use_retry = 1");
        }

        // If we're allowing for retries, check for that here.

        if (use_retry) {

            // If we hit a retry, signal that we want to try again.
            // The retry function will handle resetting the state,
            // and updating dt_subcycle.

            if (retry_advance_ctu(dt_subcycle, status)) {
                do_swap = false;
                in_retry = true;

                continue;
            }
            else {
                in_retry = false;
            }

        }

        subcycle_time += dt_subcycle;
        sub_iteration += 1;

        // Continually record the last timestep we took on this level
        // in case we need it later. We only record it if the subcycle
        // was completed successfully (i.e. we got to this point).
        // Note: this is different from last_dt_subcycle. This variable
        // records the actual timestep taken in this subcycle, while
        // the other one records the timestep as if it had not been
        // modified by the constraint of matching the final time.

        lastDt = dt_subcycle;

    }

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "  Subcycling complete" << std::endl << std::endl;

    // Record the number of subcycles we took for diagnostic purposes.

    num_subcycles_taken = sub_iteration;

    if (sub_iteration > 1) {

        // Finally, copy the original data back to the old state
        // data so that externally it appears like we took only
        // a single timestep. We'll do this as a swap so that
        // we still have the last iteration's old data if we need
        // it later.

        for (int k = 0; k < num_state_type; k++) {

            if (prev_state[k]->hasOldData())
                state[k].replaceOldData(*prev_state[k]);

            state[k].setTimeLevel(time + dt, dt, 0.0);
            prev_state[k]->setTimeLevel(time + dt, dt_subcycle, 0.0);

        }

        // If we took more than one step and are going to do a reflux,
        // keep the data past the end of the step.

        if (do_reflux && update_sources_after_reflux) {

            // Note that since we only want to do this if there's actually a
            // reflux immediately following this, skip this if we're on the
            // finest level and this is not the last iteration.

            if (!(amr_iteration < amr_ncycle && level == parent->finestLevel())) {
              keep_prev_state = true;
            }

        }

    }

    // We want to return the subcycled timestep as a suggestion.
    // Let's be sure to return the subcycled timestep that was
    // unmodified by anything we had to do in the last subcycle
    // to reach the target time.

    dt_subcycle = last_dt_subcycle;

    dt_new = std::min(dt_new, dt_subcycle);

    return dt_new;

}
