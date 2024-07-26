#include <Castro.H>

#ifdef RADIATION
#include <Radiation.H>
#endif

using namespace amrex;

void
Castro::apply_source_to_state(MultiFab& target_state, MultiFab& source, Real dt, int ng)
{
    BL_PROFILE("Castro::apply_source_to_state()");

    AMREX_ASSERT(source.nGrow() >= ng);
    AMREX_ASSERT(target_state.nGrow() >= ng);

    MultiFab::Saxpy(target_state, dt, source, 0, 0, source.nComp(), ng);
}

void
Castro::time_center_source_terms(MultiFab& S_new, MultiFab& src_old, MultiFab &src_new, Real dt)
{
  BL_PROFILE("Castro::time_center_source_terms()");

  // Subtract off half of the old source term, and add half of the new.

  MultiFab::Saxpy(S_new,-0.5*dt,src_old,0,0,S_new.nComp(),0);
  MultiFab::Saxpy(S_new, 0.5*dt,src_new,0,0,S_new.nComp(),0);
}

bool
Castro::source_flag(int src)
{
    switch(src) {

#ifdef SPONGE
    case sponge_src:
        if (do_sponge) {
            return true;
        } else {
            return false;
        }
#endif

    case ext_src:
        if (add_ext_src) {
            return true;
        } else {
            return false;
        }
#ifndef MHD
    case thermo_src:
        if (time_integration_method == SpectralDeferredCorrections) {
            return true;
        } else {
          return false;
        }
#else
    case thermo_src:
        return true;
#endif

    case geom_src:
        if (geom.Coord() == 1) {
          return true;
        } else {
          return false;
        }

#ifdef DIFFUSION
    case diff_src:
        if (diffuse_temp &&
            !(time_integration_method == SpectralDeferredCorrections)) {
          return true;
        }
        else {
          return false;
        }
#endif

#ifdef HYBRID_MOMENTUM
    case hybrid_src:
        if (castro::hybid_hydro) {
            return true;
        } else {
            return false;
        }
#endif

#ifdef GRAVITY
    case grav_src:
        if (do_grav) {
            return true;
        } else {
            return false;
        }
#endif

#ifdef ROTATION
    case rot_src:
        if (do_rotation) {
            return true;
        } else {
            return false;
        }
#endif

    default:
        return false;

    } // end switch
}

void
Castro::do_old_sources(
#ifdef MHD
                MultiFab& Bx, MultiFab& By, MultiFab& Bz,
#endif
                MultiFab& source, MultiFab& state_old, MultiFab& state_new, Real time, Real dt, bool apply_to_state)
{

    BL_PROFILE("Castro::do_old_sources()");

    const Real strt_time = ParallelDescriptor::second();

    // Construct the old-time sources.

    source.setVal(0.0, source.nGrow());

    if (!apply_sources()) {
        return;
    }

    for (int n = 0; n < num_src; ++n) {
        construct_old_source(n, source, state_old, time, dt);
    }

    if (apply_to_state) {
        apply_source_to_state(state_new, source, dt, 0);
        clean_state(
#ifdef MHD
                     Bx, By, Bz,
#endif
                     state_new, time, 0);
    }

    // Fill the ghost cells (but only for CTU; the time would not be correct for true SDC, so
    // the FillPatch is handled separately in the SDC advance).

    if (castro::time_integration_method == CornerTransportUpwind ||
        castro::time_integration_method == SimplifiedSpectralDeferredCorrections) {
        AmrLevel::FillPatch(*this, source, source.nGrow(), time, Source_Type, 0, NSRC);
    }

    // Optionally print out diagnostic information about how much
    // these source terms changed the state.

    if (apply_to_state && print_update_diagnostics) {
      bool is_new = false;
      print_all_source_changes(dt, is_new);
    }

    if (verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        amrex::Print() << "Castro::do_old_sources() time = " << run_time << " on level " << level << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }
}

advance_status
Castro::do_old_sources (Real time, Real dt, bool apply_to_state)
{
    advance_status status {};

    MultiFab& S_new = get_new_data(State_Type);

    MultiFab& old_source = get_old_data(Source_Type);

#ifdef MHD
    MultiFab& Bx_old = get_old_data(Mag_Type_x);
    MultiFab& By_old = get_old_data(Mag_Type_y);
    MultiFab& Bz_old = get_old_data(Mag_Type_z);
#endif

    do_old_sources(
#ifdef MHD
                   Bx_old, By_old, Bz_old,
#endif
                   old_source, Sborder, S_new, time, dt,
                   apply_to_state);

    return status;
}

void
Castro::do_new_sources(
#ifdef MHD
                MultiFab& Bx, MultiFab& By, MultiFab& Bz,
#endif
                MultiFab& source, MultiFab& state_old, MultiFab& state_new, Real time, Real dt, bool apply_to_state)
{

    BL_PROFILE("Castro::do_new_sources()");

    const Real strt_time = ParallelDescriptor::second();

    source.setVal(0.0, NUM_GROW_SRC);

    if (!apply_sources()) {
        return;
    }

    // Construct the new-time source terms.

    for (int n = 0; n < num_src; ++n) {
        construct_new_source(n, source, state_old, state_new, time, dt);
    }

    if (apply_to_state) {
        apply_source_to_state(state_new, source, dt, 0);
        clean_state(
#ifdef MHD
                     Bx, By, Bz,
#endif
                     state_new, time, 0);
    }


    // Optionally print out diagnostic information about how much
    // these source terms changed the state.

    if (print_update_diagnostics) {
      bool is_new = true;
      print_all_source_changes(dt, is_new);
    }

    if (verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        amrex::Print() << "Castro::do_new_sources() time = " << run_time << " on level " << level << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }

}

advance_status
Castro::do_new_sources (Real time, Real dt)
{
    advance_status status {};

    MultiFab& S_new = get_new_data(State_Type);

    MultiFab& new_source = get_new_data(Source_Type);

#ifdef MHD
    MultiFab& Bx_new = get_new_data(Mag_Type_x);
    MultiFab& By_new = get_new_data(Mag_Type_y);
    MultiFab& Bz_new = get_new_data(Mag_Type_z);
#endif

    do_new_sources(
#ifdef MHD
                   Bx_new, By_new, Bz_new,
#endif
                   new_source, Sborder, S_new, time, dt);

    return status;
}

void
Castro::construct_old_source(int src, MultiFab& source, MultiFab& state_in, Real time, Real dt)
{
    BL_PROFILE("Castro::construct_old_source()");

    BL_ASSERT(src >= 0 && src < num_src);

    switch(src) {

#ifdef SPONGE
    case sponge_src:
        construct_old_sponge_source(source, state_in, time, dt);
        break;
#endif

    case ext_src:
        construct_old_ext_source(source, state_in, time, dt);
        break;

    case thermo_src:
        construct_old_thermo_source(source, state_in, time, dt);
        break;

    case geom_src:
        construct_old_geom_source(source, state_in, time, dt);
        break;

#ifdef DIFFUSION
    case diff_src:
        if (!(time_integration_method == SpectralDeferredCorrections)) {
          // for MOL or SDC, we'll compute a diffusive flux in the MOL routine
          construct_old_diff_source(source, state_in, time, dt);
        }
        break;
#endif

#ifdef HYBRID_MOMENTUM
    case hybrid_src:
        construct_old_hybrid_source(source, state_in, time, dt);
        break;
#endif

#ifdef GRAVITY
    case grav_src:
        construct_old_gravity_source(source, state_in, time, dt);
        break;
#endif

#ifdef ROTATION
    case rot_src:
        construct_old_rotation_source(source, state_in, time, dt);
        break;
#endif

    default:
        break;

    } // end switch
}

void
Castro::construct_new_source(int src, MultiFab& source, MultiFab& state_old, MultiFab& state_new, Real time, Real dt)
{
    BL_PROFILE("Castro::construct_new_source()");

    BL_ASSERT(src >= 0 && src < num_src);

    switch(src) {

#ifdef SPONGE
    case sponge_src:
        construct_new_sponge_source(source, state_old, state_new, time, dt);
        break;
#endif

    case ext_src:
        construct_new_ext_source(source, state_old, state_new, time, dt);
        break;

#ifdef MHD
    case thermo_src:
        construct_new_thermo_source(source, state_old, state_new, time, dt);
        break;
#endif

    case geom_src:
        construct_new_geom_source(source, state_old, state_new, time, dt);
        break;

#ifdef DIFFUSION
    case diff_src:
        construct_new_diff_source(source, state_old, state_new, time, dt);
        break;
#endif

#ifdef HYBRID_MOMENTUM
    case hybrid_src:
        construct_new_hybrid_source(source, state_old, state_new, time, dt);
        break;
#endif

#ifdef GRAVITY
    case grav_src:
        construct_new_gravity_source(source, state_old, state_new, time, dt);
        break;
#endif

#ifdef ROTATION
    case rot_src:
        construct_new_rotation_source(source, state_old, state_new, time, dt);
        break;
#endif

    default:
        break;

    } // end switch
}

// Returns whether any sources are actually applied.

bool
Castro::apply_sources()
{

    for (int n = 0; n < num_src; ++n) {
        if (source_flag(n)) {
            return true;
        }
    }

    return false;

}

// Evaluate diagnostics quantities describing the effect of an
// update on the state. The optional parameter local determines
// whether we want to do this calculation globally over all processes
// or locally just on this processor. The latter is useful if you
// are evaluating the contribution from multiple source changes at once
// and want to amortize the cost of the parallel reductions.
// Note that the resultant output is volume-weighted.

Vector<Real>
Castro::evaluate_source_change(const MultiFab& source, Real dt, bool local)
{

  BL_PROFILE("Castro::evaluate_source_change()");

  Vector<Real> update(source.nComp(), 0.0);

  // Create a temporary array which will hold a single component
  // at a time of the volume-weighted source.

  MultiFab weighted_source(source.boxArray(), source.DistributionMap(), 1, 0);

  for (int n = 0; n < source.nComp(); ++n) {

    weighted_source.setVal(0.0);

    // Fill weighted_source with source x volume.

    MultiFab::AddProduct(weighted_source, source, n, volume, 0, 0, 1, 0);

    update[n] = weighted_source.sum(0, local) * dt;

  }

  return update;

}

// Print the change due to a given source term update.
// We assume here that the input array is lined up with
// the NUM_STATE components of State_Type because we are
// interested in printing changes to energy, mass, etc.

void
Castro::print_source_change(const Vector<Real>& update)
{

  if (ParallelDescriptor::IOProcessor()) {

    std::cout << "       mass added: " << update[URHO] << std::endl;
    std::cout << "       xmom added: " << update[UMX] << std::endl;
    std::cout << "       ymom added: " << update[UMY] << std::endl;
    std::cout << "       zmom added: " << update[UMZ] << std::endl;
    std::cout << "       eint added: " << update[UEINT] << std::endl;
    std::cout << "       ener added: " << update[UEDEN] << std::endl;

    std::cout << std::endl;

  }


}

// Calculate the changes to the state due to a source term,
// and also print the results.

void
Castro::evaluate_and_print_source_change (const MultiFab& source, Real dt, const std::string& source_name)
{
    bool local = true;
    Vector<Real> update = evaluate_source_change(source, dt, local);

#ifdef BL_LAZY
    Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealSum(update.dataPtr(), static_cast<int>(update.size()), ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor()) {
            if (std::abs(update[URHO]) != 0.0 || std::abs(update[UEDEN]) != 0.0) {
                std::cout << std::endl << "  Contributions to the state from " << source_name << ":" << std::endl;

                print_source_change(update);
            }
        }

#ifdef BL_LAZY
    });
#endif
}

// For the old-time or new-time sources update, evaluate the change in the state
// for all source terms, then print the results.

void
Castro::print_all_source_changes(Real dt, bool is_new)
{
    MultiFab& source = is_new ? get_new_data(Source_Type) : get_old_data(Source_Type);

    std::string source_name = is_new? "new-time sources" : "old-time sources";

    evaluate_and_print_source_change(source, dt, source_name);
}

// Perform all operations that occur prior to computing the predictor sources
// and the hydro advance.

advance_status
Castro::pre_advance_operators (Real time, Real dt)  // NOLINT(readability-convert-member-functions-to-static)
{
    amrex::ignore_unused(time);
    amrex::ignore_unused(dt);

    advance_status status {};

    // If we are using gravity, solve for the potential and
    // gravitational field.  note: since reactions don't change
    // density, we can do this before or after the burn.

#ifdef GRAVITY
    construct_old_gravity(time);
#endif

#ifdef SHOCK_VAR
    // we want to compute the shock flag that will be used
    // (optionally) in disabling reactions in shocks.  We compute this
    // only once per timestep using the time-level n data.

    // We need to compute the old sources -- note that we will
    // recompute the old sources after the burn, so this is done here
    // only for evaluating the shock flag.

    const bool apply_to_state{false};
    do_old_sources(time, dt, apply_to_state);

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& old_source = get_old_data(Source_Type);

    FArrayBox shk(The_Async_Arena());
    FArrayBox q(The_Async_Arena()), qaux(The_Async_Arena());

    for (MFIter mfi(Sborder, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const Box& obx = mfi.growntilebox(1);

        shk.resize(bx, 1);
#ifdef RADIATION
        q.resize(obx, NQ);
#else
        q.resize(obx, NQTHERM);
#endif
        qaux.resize(obx, NQAUX);

        Array4<Real> const shk_arr = shk.array();
        Array4<Real> const q_arr = q.array();
        Array4<Real> const qaux_arr = qaux.array();

        Array4<Real const> const Sborder_old_arr = Sborder.array(mfi);
        Array4<Real> const S_old_arr = S_old.array(mfi);
        Array4<Real> const old_src_arr = old_source.array(mfi);

        ctoprim(obx, time, Sborder_old_arr, q_arr, qaux_arr);

        shock(bx, q_arr, old_src_arr, shk_arr);

        // now store it in S_old -- we'll fillpatch into Sborder in a bit

        // Note: we still compute the shock flag in the hydro for the
        // hybrid-Riemann (for now) since that version will have seen the
        // effect of the burning in the first dt), but that version is
        // never stored in State_Type

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            S_old_arr(i,j,k,USHK) = shk_arr(i,j,k);
        });

    }

    // we only computed the shock flag on the interior, but the first
    // burn needs ghost cells, so FillPatch just the shock flag

    if (Sborder.nGrow() > 0) {
      AmrLevel::FillPatch(*this, Sborder, Sborder.nGrow(), time, State_Type, USHK, 1, USHK);
    }

#endif

    // If we are Strang splitting the reactions, do the old-time contribution now.

#ifndef TRUE_SDC
#ifdef REACTIONS
    status = do_old_reactions(time, dt);

    if (status.success == false) {
        return status;
    }
#endif
#endif


    // Initialize the new-time data. This copy needs to come after all Strang-split operators.

    MultiFab& S_new = get_new_data(State_Type);

    MultiFab::Copy(S_new, Sborder, 0, 0, NUM_STATE, S_new.nGrow());

    return status;
}

// Perform all operations that occur after computing the predictor sources
// but before the hydro advance.

advance_status
Castro::pre_hydro_operators (Real time, Real dt)  // NOLINT(readability-convert-member-functions-to-static)
{
    amrex::ignore_unused(time);
    amrex::ignore_unused(dt);

    advance_status status {};

#ifdef SIMPLIFIED_SDC
#ifdef REACTIONS
    // The SDC reactive source ghost cells on coarse levels might not
    // be in sync due to any average down done, so fill them here.

    MultiFab& react_src = get_new_data(Simplified_SDC_React_Type);

    AmrLevel::FillPatch(*this, react_src, react_src.nGrow(), time + dt, Simplified_SDC_React_Type, 0, react_src.nComp());
#endif
#endif

    return status;
}

// Perform all operations that occur after the hydro source
// but before the corrector sources.

advance_status
Castro::post_hydro_operators (Real time, Real dt)  // NOLINT(readability-convert-member-functions-to-static)
{
    amrex::ignore_unused(time);
    amrex::ignore_unused(dt);

    advance_status status {};

#ifdef GRAVITY
    construct_new_gravity(time);
#endif

    return status;
}

// Perform all operations that occur after the corrector sources.

advance_status
Castro::post_advance_operators (Real time, Real dt)  // NOLINT(readability-convert-member-functions-to-static)
{
    amrex::ignore_unused(time);
    amrex::ignore_unused(dt);

    advance_status status {};

#ifndef TRUE_SDC
#ifdef REACTIONS
    status = do_new_reactions(time, dt);

    if (status.success == false) {
        return status;
    }
#endif
#endif

    return status;
}
