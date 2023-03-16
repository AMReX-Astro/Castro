#include <Castro.H>
#include <Castro_F.H>

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
        return true;
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

    for (int n = 0; n < num_src; ++n) {
        construct_old_source(n, source, state_old, time, dt);

        // We can either apply the sources to the state one by one, or we can
        // group them all together at the end.

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

        amrex::Print() << "Castro::do_old_sources() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }

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

    // Construct the new-time source terms.

    for (int n = 0; n < num_src; ++n) {
        construct_new_source(n, source, state_old, state_new, time, dt);

        // We can either apply the sources to the state one by one, or we can
        // group them all together at the end.

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

        amrex::Print() << "Castro::do_new_sources() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }

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
Castro::print_source_change(Vector<Real> update)
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
