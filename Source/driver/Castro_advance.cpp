
#include <Castro.H>

#ifdef RADIATION
#include <Radiation.H>
#endif

#ifdef GRAVITY
#include <Gravity.H>
#endif

#include <cmath>
#include <climits>

#include <problem_initialize_state_data.H>

using std::string;
using namespace amrex;

Real
Castro::advance (Real time,
                 Real dt,
                 int  amr_iteration,
                 int  amr_ncycle)

  // the main driver for a single level.  This will do either the SDC
  // algorithm or the Strang-split reactions algorithm.
  //
  // arguments:
  //    time          : the current simulation time
  //    dt            : the timestep to advance (e.g., go from time to
  //                    time + dt)
  //    amr_iteration : where we are in the current AMR subcycle.  Each
  //                    level will take a number of steps to reach the
  //                    final time of the coarser level below it.  This
  //                    counter starts at 1
  //    amr_ncycle    : the number of subcycles at this level

{
    if (parent->subcyclingMode() == "None" && level > 0) {
        amrex::Print() << "\n  Advance at this level has already been completed.\n\n";
        return dt;
    }

    BL_PROFILE("Castro::advance()");

    // Save the wall time when we started the step.

    wall_time_start = ParallelDescriptor::second();

    MultiFab::RegionTag amrlevel_tag("AmrLevel_Level_" + std::to_string(level));

    Real dt_new = dt;

    int max_level_to_advance = level;

    if (parent->subcyclingMode() == "None" && level == 0) {
        max_level_to_advance = parent->finestLevel();
    }

    for (int lev = level; lev <= max_level_to_advance; ++lev) {
        getLevel(lev).initialize_advance(time, dt, amr_iteration);
    }

    // Do the advance.

    if (time_integration_method == CornerTransportUpwind || time_integration_method == SimplifiedSpectralDeferredCorrections) {

        dt_new = std::min(dt_new, subcycle_advance_ctu(time, dt, amr_iteration, amr_ncycle));

#ifndef MHD
#ifndef AMREX_USE_GPU
#ifdef TRUE_SDC
    } else if (time_integration_method == SpectralDeferredCorrections) {

      for (int iter = 0; iter < sdc_order+sdc_extra; ++iter) {
        sdc_iteration = iter;
        dt_new = do_advance_sdc(time, dt, amr_iteration, amr_ncycle);
      }

#endif // TRUE_SDC
#endif // AMREX_USE_GPU
#endif //MHD
    }

    // If the user requests, indicate that we want a regrid at the end of the step.

    if (use_post_step_regrid == 1) {
        post_step_regrid = 1;
    }

    for (int lev = level; lev <= max_level_to_advance; ++lev) {
#ifdef GRAVITY
        // Update the point mass.
        if (use_point_mass == 1) {
            getLevel(lev).pointmass_update(time, dt);
        }
#endif

#ifdef RADIATION
        MultiFab& S_new = getLevel(lev).get_new_data(State_Type);
        getLevel(lev).final_radiation_call(S_new, amr_iteration, amr_ncycle);
#endif

#ifdef AMREX_PARTICLES
        getLevel(lev).advance_particles(amr_iteration, time, dt);
#endif

        getLevel(lev).finalize_advance();
    }

    return dt_new;
}


advance_status
Castro::initialize_do_advance (Real time, Real dt)
{
    amrex::ignore_unused(time);

    BL_PROFILE("Castro::initialize_do_advance()");

    advance_status status {};

#ifdef RADIATION
    // make sure these are filled to avoid check/plot file errors:
    if (do_radiation) {
      get_old_data(Rad_Type).setBndry(0.0);
      get_new_data(Rad_Type).setBndry(0.0);
    }
    else {
      get_old_data(Rad_Type).setVal(0.0);
      get_new_data(Rad_Type).setVal(0.0);
    }
#endif

#ifdef GRAVITY
    if (moving_center == 1) {
        define_new_center(get_old_data(State_Type), time);
    }
#endif

    // For the hydrodynamics update we need to have NUM_GROW ghost
    // zones available, but the state data does not carry ghost
    // zones. So we use a FillPatch using the state data to give us
    // Sborder, which does have ghost zones.

    if (time_integration_method == CornerTransportUpwind || time_integration_method == SimplifiedSpectralDeferredCorrections) {
#ifdef MHD
      MultiFab& Bx_old = get_old_data(Mag_Type_x);
      MultiFab& By_old = get_old_data(Mag_Type_y);
      MultiFab& Bz_old = get_old_data(Mag_Type_z);

      Bx_old_tmp.define(Bx_old.boxArray(), Bx_old.DistributionMap(), 1, NUM_GROW);
      By_old_tmp.define(By_old.boxArray(), By_old.DistributionMap(), 1, NUM_GROW);
      Bz_old_tmp.define(Bz_old.boxArray(), Bz_old.DistributionMap(), 1, NUM_GROW);

      FillPatch(*this, Bx_old_tmp, NUM_GROW, time, Mag_Type_x, 0, 1);
      FillPatch(*this, By_old_tmp, NUM_GROW, time, Mag_Type_y, 0, 1);
      FillPatch(*this, Bz_old_tmp, NUM_GROW, time, Mag_Type_z, 0, 1);
#endif
      // for the CTU unsplit method, we always start with the old
      // state note: although clean_state has already been done on
      // the old state in initialize_advance, we still need to do
      // another here to ensure the ghost zones are thermodynamically
      // consistent
      Sborder.define(grids, dmap, NUM_STATE, NUM_GROW, MFInfo().SetTag("Sborder"));
      const Real prev_time = state[State_Type].prevTime();
      expand_state(Sborder, prev_time, NUM_GROW);
      clean_state(
#ifdef MHD
                  Bx_old_tmp, By_old_tmp, Bz_old_tmp,
#endif
                  Sborder, prev_time, NUM_GROW);

    } else if (time_integration_method == SpectralDeferredCorrections) {

      // we'll handle the filling inside of do_advance_sdc
      Sborder.define(grids, dmap, NUM_STATE, NUM_GROW, MFInfo().SetTag("Sborder"));

    } else {
      amrex::Abort("invalid time_integration_method");
    }

#ifdef SHOCK_VAR
    // Zero out the shock data, and fill it during the advance.
    // For subcycling cases this will always give the shock
    // variable for the latest subcycle, rather than averaging.

    Sborder.setVal(0.0, USHK, 1, Sborder.nGrow());
#endif

    // Create any correctors to the source term data. This must be done
    // before the source term data is overwritten below. Note: we do
    // not create the corrector source if we're currently retrying the
    // step; we will already have done it, and aside from avoiding
    // duplicate work, we have already lost the data needed to do this
    // calculation since we overwrote the data from the previous step.

    if (castro::time_integration_method == CornerTransportUpwind ||
        castro::time_integration_method == SimplifiedSpectralDeferredCorrections) {
        if (!in_retry) {
            create_source_corrector();
        }
    }

    // Check for NaN's.

    const MultiFab& S_old = get_old_data(State_Type);

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

        // We only display the estTimeStep output if the validity check fails
        std::string estTimeStep_output;

        Real old_dt;

        {
            CoutRedirection redirection;
            old_dt = estTimeStep(is_new);
            estTimeStep_output = redirection.getCapturedOutput();
        }

        if (castro::change_max * old_dt < dt) {
            status.success = false;
            std::cout << estTimeStep_output;
            status.reason = "pre-advance timestep validity check failed";
        }
    }

    return status;
}



advance_status
Castro::finalize_do_advance (Real time, Real dt)
{
    amrex::ignore_unused(time);

    BL_PROFILE("Castro::finalize_do_advance()");

    advance_status status {};

    // Check if this timestep violated our stability criteria. Our idea is,
    // if the timestep created a velocity v and sound speed at the new time
    // such that (v+c) * dt / dx < CFL / change_max, where CFL is the user's
    // chosen timestep constraint and change_max is the factor that determines
    // how much the timestep can change during an advance, consider the advance
    // to have failed. This prevents the timestep from shrinking too much,
    // whereas in computeNewDt change_max prevents the timestep from growing
    // too much. The same reasoning applies for the other timestep limiters.

    if (castro::time_integration_method == CornerTransportUpwind ||
        castro::time_integration_method == SimplifiedSpectralDeferredCorrections) {
        if (castro::check_dt_after_advance) {

            // But don't do this check if we're using simplified SDC and we're not yet
            // on the final SDC iteration, since we're not yet at the final advance.

            bool do_validity_check = true;

            if (castro::time_integration_method == SimplifiedSpectralDeferredCorrections &&
                sdc_iteration < sdc_iters - 1) {
                do_validity_check = false;
            }

            if (do_validity_check) {
                int is_new = 1;

                // We only display the estTimeStep output if the validity check fails
                std::string estTimeStep_output;

                Real new_dt;

                {
                    CoutRedirection redirection;
                    new_dt = estTimeStep(is_new);
                    estTimeStep_output = redirection.getCapturedOutput();
                }

                if (castro::change_max * new_dt < dt) {
                    status.success = false;
                    std::cout << estTimeStep_output;
                    status.reason = "post-advance timestep validity check failed";
                    return status;
                }
            }
        }
    }

#ifdef RADIATION
    if (!do_hydro && Radiation::rad_hydro_combined) {
        MultiFab& Er_old = get_old_data(Rad_Type);
        MultiFab& Er_new = get_new_data(Rad_Type);
        MultiFab::Copy(Er_new, Er_old, 0, 0, Er_old.nComp(), 0);
    }
#endif

    Sborder.clear();

    return status;
}



void
Castro::initialize_advance(Real time, Real dt, int amr_iteration)
{
    BL_PROFILE("Castro::initialize_advance()");

    // Save the current iteration.

    iteration = amr_iteration;

    sub_iteration = 0;
    sub_ncycle = 0;
    dt_subcycle = 1.e200;
    dt_advance = dt;

    keep_prev_state = false;

    // Reset the retry information.

    in_retry = 0;
    num_subcycles_taken = 1;

    if (use_post_step_regrid == 1 && level > 0) {

        if (getLevel(level-1).post_step_regrid == 1 && amr_iteration == 1) {

            // If the level below this just triggered a special regrid,
            // the coarse contribution to this level's FluxRegister
            // is no longer valid because the grids have, in general, changed.
            // Zero it out, and add them back using the saved copy of the fluxes.

            getLevel(level-1).FluxRegCrseInit();

        }

    }

    // The option of whether to do a multilevel initialization is
    // controlled within the radiation class.  This step belongs
    // before the swap.

#ifdef RADIATION
    if (do_radiation) {
        radiation->pre_timestep(level);
    }

    Erborder.define(grids, dmap, Radiation::nGroups, NUM_GROW);
    lamborder.define(grids, dmap, Radiation::nGroups, NUM_GROW);
#endif

#ifdef GRAVITY
    // If we're on level 0, update the maximum density used in the gravity solver
    // for setting the tolerances. This will be used in all level solves to follow.
    // This must be done before the swap because it relies on the new data.

    if (level == 0 && do_grav == 1 && gravity->get_gravity_type() == "PoissonGrav") {
        gravity->update_max_rhs();
    }
#endif

    // This array holds the source term corrector.

    source_corrector.define(grids, dmap, NSRC, NUM_GROW_SRC);
    source_corrector.setVal(0.0, NUM_GROW_SRC);

    // Swap the new data from the last timestep into the old state data.

    swap_state_time_levels(dt);

#ifdef GRAVITY
    if (do_grav == 1) {
        gravity->swapTimeLevels(level);
    }
#endif

    MultiFab& S_old = get_old_data(State_Type);

    // if we are doing drive_initial_convection, check to see if we
    // need to reinitialize the thermodynamic data (while keeping the
    // velocity unchanged)

#ifndef MHD

    const Real cur_time = state[State_Type].curTime();
    const Real dt_level = parent->dtLevel(level);

    if (drive_initial_convection == 1 && cur_time <= drive_initial_convection_tmax) {

        // Calculate the new dt by comparing to the dt needed to get
        // to the next multiple of drive_initial_convection_reinit_period

        const Real dtMod = std::fmod(cur_time, drive_initial_convection_reinit_period);

        Real reinit_dt;

        // Note that if we are just about exactly on a multiple of
        // drive_initial_convection_reinit_period, then we need to be
        // careful to avoid floating point issues.


        if (std::abs(dtMod - drive_initial_convection_reinit_period) <=
            std::numeric_limits<Real>::epsilon() * cur_time) {
                reinit_dt = drive_initial_convection_reinit_period +
                    (drive_initial_convection_reinit_period - dtMod);
        } else {
            reinit_dt = drive_initial_convection_reinit_period - dtMod;
        }

        if (reinit_dt < dt_level) {

          amrex::Print() << "<<<<< drive initial convection reset >>>>" << std::endl;

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
            for (MFIter mfi(S_old, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& box = mfi.tilebox();

                auto s = S_old[mfi].array();
                auto geomdata = geom.data();

#ifdef RNG_STATE_INIT
                amrex::Error("drive initial convection not yet supported for random initialization");
#else
                amrex::ParallelFor(box,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // redo the problem initialization.  We want to preserve
                    // the current velocity though, so save that and then
                    // restore it afterwards.

                    Real vx_orig = s(i,j,k,UMX) / s(i,j,k,URHO);
                    Real vy_orig = s(i,j,k,UMY) / s(i,j,k,URHO);
                    Real vz_orig = s(i,j,k,UMZ) / s(i,j,k,URHO);

                    problem_initialize_state_data(i, j, k, s, geomdata);

                    s(i,j,k,UMX) = s(i,j,k,URHO) * vx_orig;
                    s(i,j,k,UMY) = s(i,j,k,URHO) * vy_orig;
                    s(i,j,k,UMZ) = s(i,j,k,URHO) * vz_orig;

                    s(i,j,k,UEDEN) = s(i,j,k,UEINT) + 0.5_rt * s(i,j,k,URHO) *
                        (vx_orig * vx_orig + vy_orig * vy_orig + vz_orig * vz_orig);

                });
#endif

            }

        }
    }
#endif


    // Ensure data is valid before beginning advance. This addresses
    // the fact that we may have new data on this level that was interpolated
    // from a coarser level, and the interpolation in general cannot be
    // trusted to respect the consistency between certain state variables
    // (e.g. UEINT and UEDEN) that we demand in every zone.

    clean_state(
#ifdef MHD
                 get_old_data(Mag_Type_x),
                 get_old_data(Mag_Type_y),
                 get_old_data(Mag_Type_z),
#endif
                  S_old, time, S_old.nGrow());


    // Initialize the previous state data container now, so that we can
    // always ask if it has valid data.

    for (int k = 0; k < num_state_type; ++k) {
        prev_state[k] = std::make_unique<StateData>();
    }


    // Allocate space for the primitive variables.

#ifdef TRUE_SDC
    q.define(grids, dmap, NQ, NUM_GROW);
    q.setVal(0.0);
    qaux.define(grids, dmap, NQAUX, NUM_GROW);


    if (sdc_order == 4) {
      q_bar.define(grids, dmap, NQ, NUM_GROW);
      qaux_bar.define(grids, dmap, NQAUX, NUM_GROW);
#ifdef DIFFUSION
      T_cc.define(grids, dmap, 1, NUM_GROW);
#endif
    }

    if (time_integration_method == SpectralDeferredCorrections) {

      k_new.resize(SDC_NODES);

      k_new[0] = std::make_unique<MultiFab>(S_old, amrex::make_alias, 0, NUM_STATE);
      for (int n = 1; n < SDC_NODES; ++n) {
        k_new[n] = std::make_unique<MultiFab>(grids, dmap, NUM_STATE, 0);
        k_new[n]->setVal(0.0);
      }

      A_old.resize(SDC_NODES);
      for (int n = 0; n < SDC_NODES; ++n) {
        A_old[n] = std::make_unique<MultiFab>(grids, dmap, NUM_STATE, 0);
        A_old[n]->setVal(0.0);
      }

      A_new.resize(SDC_NODES);
      A_new[0] = std::make_unique<MultiFab>(*A_old[0], amrex::make_alias, 0, NUM_STATE);
      for (int n = 1; n < SDC_NODES; ++n) {
        A_new[n] = std::make_unique<MultiFab>(grids, dmap, NUM_STATE, 0);
        A_new[n]->setVal(0.0);
      }

      // We use Sburn a few ways for the SDC integration.  First, we
      // use it to store the initial guess to the nonlinear solve.
      // Second, at the end of the SDC update, we copy the cell-center
      // reaction source into it, including one ghost cell, for later
      // filling of the plotfile.  Finally, we use it as a temporary
      // buffer for when we convert the state to centers while making the
      // source term
      Sburn.define(grids, dmap, NUM_STATE, 2);

#ifdef REACTIONS
      R_old.resize(SDC_NODES);
      for (int n = 0; n < SDC_NODES; ++n) {
        R_old[n] = std::make_unique<MultiFab>(grids, dmap, NUM_STATE, 0);
        R_old[n]->setVal(0.0);
      }
#endif

    }
#endif

    // Zero out the current fluxes.

    for (int dir = 0; dir < 3; ++dir) {
        fluxes[dir]->setVal(0.0);
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
    if (store_burn_weights == 1) {
        burn_weights.setVal(0.0);
    }
#endif

}



void
Castro::finalize_advance()
{
    BL_PROFILE("Castro::finalize_advance()");

    if (do_reflux == 1 && parent->subcyclingMode() != "None") {
        FluxRegCrseInit();
        FluxRegFineAdd();
    }

    // The mass_fluxes array currently holds only the fluxes
    // from the last subcycle. Override this with the sum of
    // the fluxes from the full timestep (this will be used
    // later during the reflux operation).

    if (do_reflux == 1 && update_sources_after_reflux == 1 && parent->subcyclingMode() != "None") {
        for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {
            MultiFab::Copy(*mass_fluxes[idir], *fluxes[idir], URHO, 0, 1, 0);
        }
    }

#ifdef TRUE_SDC
    q.clear();
    qaux.clear();

    if (sdc_order == 4) {
      q_bar.clear();
      qaux_bar.clear();
#ifdef DIFFUSION
      T_cc.clear();
#endif
    }
#endif

#ifdef RADIATION
    Erborder.clear();
    lamborder.clear();
#endif

    source_corrector.clear();

    if (!keep_prev_state) {
        amrex::FillNull(prev_state);
    }

#ifdef TRUE_SDC
    if (time_integration_method == SpectralDeferredCorrections) {
      k_new.clear();
      A_new.clear();
      A_old.clear();
#ifdef REACTIONS
      R_old.clear();
      Sburn.clear();
#endif
    }
#endif

    // Record how many zones we have advanced.

    int max_level_to_advance = level;

    if (parent->subcyclingMode() == "None") {
        max_level_to_advance = parent->finestLevel();
    }

    long num_pts_advanced = 0;

    for (int lev = level; lev <= max_level_to_advance; ++lev) {
        num_pts_advanced += getLevel(lev).grids.numPts();
    }

    num_zones_advanced += static_cast<Real>(num_pts_advanced) / static_cast<Real>(getLevel(0).grids.numPts());

    Real wall_time = ParallelDescriptor::second() - wall_time_start;

    Real fom_advance = static_cast<Real>(num_pts_advanced) / wall_time / 1.e6;

    if (verbose >= 1) {
        if (max_level_to_advance > 0) {
            if (level == 0) {
                amrex::Print() << "  Zones advanced per microsecond from level " << level << " to level "
                               << max_level_to_advance << ": " << fom_advance << std::endl << std::endl;
            }
        }
        else {
            amrex::Print() << "  Zones advanced per microsecond at this level: "
                           << fom_advance << std::endl << std::endl;
        }
    }
}
