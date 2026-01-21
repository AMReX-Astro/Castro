
#include <Castro.H>

#ifdef RADIATION
#include <Radiation.H>
#endif

#ifdef GRAVITY
#include <Gravity.H>
#endif

using namespace amrex;

advance_status
Castro::do_advance_ctu (Real time, Real dt)  // NOLINT(readability-convert-member-functions-to-static)
{
    // this routine will advance the old state data (called Sborder here)
    // to the new time, for a single level.  The new data is called
    // S_new here.  The update includes reactions (if we are not doing
    // SDC), hydro, and the source terms.

    amrex::ignore_unused(time);
    amrex::ignore_unused(dt);

    BL_PROFILE("Castro::do_advance_ctu()");

    advance_status status {};

#ifndef TRUE_SDC

    // Advance simultaneously on all levels that are not subcycling
    // relative to this level.

    int max_level_to_advance = level;

    if (parent->subcyclingMode() == "None" && level == 0) {
        max_level_to_advance = parent->finestLevel();
    }

    const Real prev_time = state[State_Type].prevTime();
    const Real  cur_time = state[State_Type].curTime();

    for (int lev = level; lev <= max_level_to_advance; ++lev) {
        // Perform initialization steps.

        status = getLevel(lev).initialize_do_advance(time, dt);

        if (status.success == false) {
            return status;
        }

        // Perform all pre-advance operations and then initialize
        // the new-time state with the output of those operators.

        status = getLevel(lev).pre_advance_operators(prev_time, dt);

        if (status.success == false) {
            return status;
        }

        // Construct the old-time sources from Sborder. This will already
        // be applied to S_new (with full dt weighting), to be corrected
        // later. Note that this does not affect the prediction of the
        // interface state; an explicit source will be traced there as
        // needed.

        status = getLevel(lev).do_old_sources(prev_time, dt);

        if (status.success == false) {
            return status;
        }

        // Perform any operations that occur after the sources but before the hydro.

        status = getLevel(lev).pre_hydro_operators(prev_time, dt);

        if (status.success == false) {
            return status;
        }

        // Do the hydro update. We build directly off of Sborder, which
        // is the state that has already seen the burn.

#ifndef MHD
        status = getLevel(lev).construct_ctu_hydro_source(prev_time, dt);
#else
        status = getLevel(lev).construct_ctu_mhd_source(prev_time, dt);
#endif

        if (status.success == false) {
            return status;
        }
    }

    // We can perform the reflux immediately if there's no subcycling
    // above this level.

    if (do_reflux && level < max_level_to_advance) {
        reflux(level, max_level_to_advance, false);
    }

    for (int lev = level; lev <= max_level_to_advance; ++lev) {
        // Perform any operations that occur after the hydro but before
        // the corrector sources.

        status = getLevel(lev).post_hydro_operators(cur_time, dt);

        if (status.success == false) {
            return status;
        }

        // Construct and apply new-time source terms.

        status = getLevel(lev).do_new_sources(cur_time, dt);

        if (status.success == false) {
            return status;
        }

        // Do the second half of the reactions for Strang, or the full burn for simplified SDC.

        status = getLevel(lev).post_advance_operators(cur_time, dt);

        if (status.success == false) {
            return status;
        }

        // Perform finalization steps.

        status = getLevel(lev).finalize_do_advance(cur_time, dt);

        if (status.success == false) {
            return status;
        }
    }

#endif

    return status;

}





bool
Castro::retry_advance_ctu(Real dt, const advance_status& status)
{
    BL_PROFILE("Castro::retry_advance_ctu()");

    bool do_retry = false;

    if (!status.success) {
        do_retry = true;
    }

    if (do_retry) {

        int max_level_to_advance = level;

        if (parent->subcyclingMode() == "None" && level == 0) {
            max_level_to_advance = parent->finestLevel();
        }

        for (int lev = level; lev <= max_level_to_advance; ++lev) {
            if (status.suggested_dt > 0.0_rt && status.suggested_dt < dt) {
                getLevel(lev).dt_subcycle = status.suggested_dt;
            } else {
                getLevel(lev).dt_subcycle = std::min(dt, getLevel(lev).dt_subcycle) * retry_subcycle_factor;
            }
        }

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

        for (int lev = level; lev <= max_level_to_advance; ++lev) {
            getLevel(lev).save_data_for_retry();

            // Clear the contribution to the fluxes from this step.

            for (int dir = 0; dir < 3; ++dir) {
                getLevel(lev).fluxes[dir]->setVal(0.0);
            }

            for (int dir = 0; dir < 3; ++dir) {
                getLevel(lev).mass_fluxes[dir]->setVal(0.0);
            }

#if (AMREX_SPACEDIM <= 2)
            if (!Geom().IsCartesian()) {
                getLevel(lev).P_radial.setVal(0.0);
            }
#endif

#ifdef RADIATION
            if (Radiation::rad_hydro_combined) {
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    getLevel(lev).rad_fluxes[dir]->setVal(0.0);
                }
            }
#endif

#ifdef REACTIONS
            if (castro::store_burn_weights) {
                getLevel(lev).burn_weights.setVal(0.0);
            }
#endif

            // For simplified SDC, we'll have garbage data if we
            // attempt to use the lagged source terms (both reacting
            // and non-reacting) from the last timestep, since that
            // advance failed and we don't know if we can trust it.
            // So we zero out both source term correctors.

            if (time_integration_method == SimplifiedSpectralDeferredCorrections) {
                getLevel(lev).source_corrector.setVal(0.0, getLevel(lev).source_corrector.nGrow());

#ifdef REACTIONS
                MultiFab& SDC_react_new = getLevel(lev).get_new_data(Simplified_SDC_React_Type);
                SDC_react_new.setVal(0.0, SDC_react_new.nGrow());
#endif
            }

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

    if (dt_subcycle == 1.e200) {
        dt_subcycle = dt;
    }

    Real subcycle_time = time;

    Real dt_new = 1.e200;

    sub_iteration = 0;

    int max_level_to_advance = level;

    if (parent->subcyclingMode() == "None" && level == 0) {
        max_level_to_advance = parent->finestLevel();
    }

    // Subcycle until we've reached the target time.
    // Compare against a slightly smaller number to
    // avoid roundoff concerns.

    Real eps = 1.0e-14;

    bool do_swap = false;

    Real last_dt_subcycle = 1.e200;

#ifdef NSE_NET
    bool do_loosen = true;
    bool old_nse_dx_independent = nse_dx_independent;
    bool old_nse_molar_independent = nse_molar_independent;
    bool old_nse_skip_molar = nse_skip_molar;
#endif
    
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
#ifdef NSE_NET
	  if (loosen_nse_bailout && do_loosen) {

	    // find the smallest allowed dt_subcycle and then loosen nse_net bailout condition

	    dt_subcycle = ((time + dt) - subcycle_time) / (max_subcycles);
	    num_subcycles_remaining = max_subcycles;

	    nse_dx_independent = true;
	    nse_molar_independent = true;
	    nse_skip_molar = true;
	    do_loosen = false;

	    amrex::Print() << std::endl
			   << "  The subcycle mechanism requested " << num_subcycles_remaining << " subcycled timesteps, which is larger than the maximum of " << max_subcycles << "." << std::endl
			   << "  Reperforming subcycle with smallest possible dt_subcycle and loosen nse bailout conditions." << std::endl;
	  } else {
#endif
            amrex::Print() << std::endl
                           << "  The subcycle mechanism requested " << num_subcycles_remaining
                           << " subcycled timesteps, which is larger than the maximum of "
                           << max_subcycles << "." << std::endl
                           << "  If you would like to override this, increase the parameter castro.max_subcycles." << std::endl;
            amrex::Abort("Error: too many subcycles.");
#ifdef NSE_NET
	  }
#endif
        }

        // If we get to this point, we survived the sanity checks. Print out the current subcycle iteration.

        if (verbose && ParallelDescriptor::IOProcessor()) {
            std::cout << std::endl;
            std::cout << Font::Bold << FGColor::Green
                      << "  Beginning subcycle " << sub_iteration + 1
                      << " starting at time " << subcycle_time
                      << " with dt = " << dt_subcycle << ResetDisplay << std::endl;
            std::cout << "  Estimated number of subcycles remaining: "
                      << num_subcycles_remaining << std::endl << std::endl;
        }

        // Swap the time levels. Only do this after the first iteration;
        // the first iteration used the swap done in initialize_advance.
        // After that, the only exception will be if we do a retry.

        if (do_swap) {

            for (int lev = level; lev <= max_level_to_advance; ++lev) {
                getLevel(lev).swap_state_time_levels(0.0);

#ifdef GRAVITY
                if (do_grav) {
                    getLevel(lev).gravity->swapTimeLevels(lev);
                }
#endif
            }

        }  else {
            do_swap = true;
        }

        // Set the relevant time levels.

        for (int lev = level; lev <= max_level_to_advance; ++lev) {
            for (int k = 0; k < num_state_type; k++) {
                getLevel(lev).state[k].setTimeLevel(subcycle_time + dt_subcycle, dt_subcycle, 0.0);
            }
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

        advance_status status {};

        for (int n = 0; n < num_sub_iters; ++n) {

            if (time_integration_method == SimplifiedSpectralDeferredCorrections) {
                for (int lev = level; lev <= max_level_to_advance; ++lev) {
                    getLevel(lev).sdc_iteration = n;
                }

                amrex::Print() << "Beginning SDC iteration " << n + 1 << " of " << num_sub_iters << "." << std::endl << std::endl;
            }

            // We do the hydro advance here, and record whether we completed it.

            status = do_advance_ctu(subcycle_time, dt_subcycle);

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

    if (verbose) {
        amrex::Print() << "  Subcycling complete" << std::endl << std::endl;
    }

    // Record the number of subcycles we took for diagnostic purposes.

    for (int lev = level; lev <= max_level_to_advance; ++lev) {
        getLevel(lev).num_subcycles_taken = sub_iteration;
    }

    if (sub_iteration > 1) {

        // Finally, copy the original data back to the old state
        // data so that externally it appears like we took only
        // a single timestep. We'll do this as a swap so that
        // we still have the last iteration's old data if we need
        // it later.

        for (int lev = level; lev <= max_level_to_advance; ++lev) {
            for (int k = 0; k < num_state_type; k++) {
                if (getLevel(lev).prev_state[k]->hasOldData()) {
                    getLevel(lev).state[k].replaceOldData(*getLevel(lev).prev_state[k]);
                }
                getLevel(lev).state[k].setTimeLevel(time + dt, dt, 0.0);
                getLevel(lev).prev_state[k]->setTimeLevel(time + dt, dt_subcycle, 0.0);
            }
        }

        // If we took more than one step and are going to do a reflux,
        // keep the data past the end of the step.

        if (do_reflux && update_sources_after_reflux) {

            // Note that since we only want to do this if there's actually a
            // reflux immediately following this, skip this if we're on the
            // finest level and this is not the last iteration.

            for (int lev = level; lev <= max_level_to_advance; ++lev) {
                if (!(amr_iteration < amr_ncycle && lev == parent->finestLevel())) {
                    getLevel(lev).keep_prev_state = true;
                }
            }

        }

    }

#ifdef NSE_NET
    // Restore the original nse configuration

    nse_dx_independent = old_nse_dx_independent;
    nse_molar_independent = old_nse_molar_independent;
    nse_skip_molar = old_nse_skip_molar;
#endif
    
    // We want to return the subcycled timestep as a suggestion.
    // Let's be sure to return the subcycled timestep that was
    // unmodified by anything we had to do in the last subcycle
    // to reach the target time.

    dt_subcycle = last_dt_subcycle;

    dt_new = std::min(dt_new, dt_subcycle);

    return dt_new;

}
