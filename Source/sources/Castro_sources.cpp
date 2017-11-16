#include "Castro.H"
#include "Castro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

using namespace amrex;

void
Castro::apply_source_to_state(MultiFab& state, MultiFab& source, Real dt)
{

  MultiFab::Saxpy(state, dt, source, 0, 0, NUM_STATE, 0);

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
	if (do_sponge)
	    return true;
	else
	    return false;
#endif

    case ext_src:
	if (add_ext_src)
	    return true;
	else
	    return false;

#ifdef DIFFUSION
    case diff_src:
	return true;
#endif

#ifdef HYBRID_MOMENTUM
    case hybrid_src:
	return true;
#endif

#ifdef GRAVITY
    case grav_src:
	if (do_grav)
	    return true;
	else
	    return false;
#endif

#ifdef ROTATION
    case rot_src:
	if (do_rotation)
	    return true;
	else
	    return false;
#endif

    default:
	return false;

    } // end switch
}

void
Castro::do_old_sources(Real time, Real dt, int amr_iteration, int amr_ncycle)
{

    // Construct the old-time sources.

    MultiFab& old_sources = get_old_data(Source_Type);

    old_sources.setVal(0.0);

    for (int n = 0; n < num_src; ++n)
        construct_old_source(n, time, dt, amr_iteration, amr_ncycle);

    // Apply the old-time sources directly to the new-time state,
    // S_new -- note that this addition is for full dt, since we
    // will do a predictor-corrector on the sources to allow for
    // state-dependent sources.

    MultiFab& S_new = get_new_data(State_Type);

    // Only do this if at least one source term has a non-zero contribution.

    bool apply_source = false;
    for (int n = 0; n < num_src; ++n) {
        if (source_flag(n)) {
            apply_source = true;
            break;
        }
    }

    if (apply_source)    
        apply_source_to_state(S_new, old_sources, dt);

    // Optionally print out diagnostic information about how much
    // these source terms changed the state.

    if (print_update_diagnostics) {
      bool is_new = false;
      print_all_source_changes(dt, is_new);
    }

}

void
Castro::do_new_sources(Real time, Real dt, int amr_iteration, int amr_ncycle)
{

    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& new_sources = get_new_data(Source_Type);

    // For the new-time source terms, we have an option for how to proceed.
    // We can either construct all of the old-time sources using the same
    // state that comes out of the hydro update, or we can evaluate the sources
    // one by one and apply them as we go.

    new_sources.setVal(0.0);

    // Construct the new-time source terms.

    for (int n = 0; n < num_src; ++n)
        construct_new_source(n, time, dt, amr_iteration, amr_ncycle);

    // Apply the new-time sources to the state.
    // Only do this if at least one source term has a non-zero contribution.

    bool apply_source = false;
    for (int n = 0; n < num_src; ++n) {
        if (source_flag(n)) {
            apply_source = true;
            break;
        }
    }

    if (apply_source) {

        apply_source_to_state(S_new, new_sources, dt);

        clean_state(S_new);

    }

    // Optionally print out diagnostic information about how much
    // these source terms changed the state.

    if (print_update_diagnostics) {
      bool is_new = true;
      print_all_source_changes(dt, is_new);
    }


}

void
Castro::construct_old_source(int src, Real time, Real dt, int amr_iteration, int amr_ncycle)
{
    BL_ASSERT(src >= 0 && src < num_src);

    switch(src) {

#ifdef SPONGE
    case sponge_src:
	construct_old_sponge_source(time, dt);
	break;
#endif

    case ext_src:
	construct_old_ext_source(time, dt);
	break;

#ifdef DIFFUSION
    case diff_src:
	construct_old_diff_source(time, dt);
	break;
#endif

#ifdef HYBRID_MOMENTUM
    case hybrid_src:
	construct_old_hybrid_source(time, dt);
	break;
#endif

#ifdef GRAVITY
    case grav_src:
	construct_old_gravity_source(time, dt);
	break;
#endif

#ifdef ROTATION
    case rot_src:
	construct_old_rotation_source(time, dt);
	break;
#endif

    default:
	break;

    } // end switch
}

void
Castro::construct_new_source(int src, Real time, Real dt, int amr_iteration, int amr_ncycle)
{
    BL_ASSERT(src >= 0 && src < num_src);

    switch(src) {

#ifdef SPONGE
    case sponge_src:
	construct_new_sponge_source(time, dt);
	break;
#endif

    case ext_src:
	construct_new_ext_source(time, dt);
	break;

#ifdef DIFFUSION
    case diff_src:
	construct_new_diff_source(time, dt);
	break;
#endif

#ifdef HYBRID_MOMENTUM
    case hybrid_src:
	construct_new_hybrid_source(time, dt);
	break;
#endif

#ifdef GRAVITY
    case grav_src:
	construct_new_gravity_source(time, dt);
	break;
#endif

#ifdef ROTATION
    case rot_src:
	construct_new_rotation_source(time, dt);
	break;
#endif

    default:
	break;

    } // end switch
}

// Evaluate diagnostics quantities describing the effect of an
// update on the state. The optional parameter local determines
// whether we want to do this calculation globally over all processes
// or locally just on this processor. The latter is useful if you
// are evaluating the contribution from multiple source changes at once
// and want to amortize the cost of the parallel reductions.
// Note that the resultant output is volume-weighted.

Vector<Real>
Castro::evaluate_source_change(MultiFab& source, Real dt, bool local)
{

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

  BL_ASSERT(update.size() == NUM_STATE);

  if (ParallelDescriptor::IOProcessor()) {

    std::cout << "       mass added: " << update[Density] << std::endl;
    std::cout << "       xmom added: " << update[Xmom] << std::endl;
#if (BL_SPACEDIM >= 2)
    std::cout << "       ymom added: " << update[Ymom] << std::endl;
#endif
#if (BL_SPACEDIM == 3)
    std::cout << "       zmom added: " << update[Zmom] << std::endl;
#endif
    std::cout << "       eint added: " << update[Eint] << std::endl;
    std::cout << "       ener added: " << update[Eden] << std::endl;

    std::cout << std::endl;

  }


}

// For the old-time or new-time sources update, evaluate the change in the state
// for all source terms, then pring the results.

void
Castro::print_all_source_changes(Real dt, bool is_new)
{

  Vector<Real> summed_updates;

  bool local = true;

  MultiFab& source = is_new ? get_new_data(Source_Type) : get_old_data(Source_Type);

  summed_updates = evaluate_source_change(source, dt, local);

#ifdef BL_LAZY
  Lazy::QueueReduction( [=] () mutable {
#endif

      ParallelDescriptor::ReduceRealSum(summed_updates.dataPtr(), NUM_STATE, ParallelDescriptor::IOProcessorNumber());

      std::string time = is_new ? "new" : "old";

      if (ParallelDescriptor::IOProcessor())
          std::cout << std::endl << "  Contributions to the state from the " << time << "-time sources:" << std::endl;

      print_source_change(summed_updates);

#ifdef BL_LAZY
    });
#endif

} 

// Obtain the sum of all source terms.

void
Castro::sum_of_sources(MultiFab& source)
{

  // this computes advective_source + 1/2 (old source + new source)
  // 
  // Note: the advective source is defined as -div{F}
  //
  // the time-centering is accomplished since new source is defined
  // to be 1/2 (new source - old source) generally.

  int ng = source.nGrow();

  source.setVal(0.0);

  MultiFab& old_sources = get_old_data(Source_Type);
  MultiFab& new_sources = get_new_data(Source_Type);

  MultiFab::Add(source, old_sources, 0, 0, NUM_STATE, ng);

  MultiFab::Add(source, hydro_source, 0, 0, NUM_STATE, ng);

  MultiFab::Add(source, new_sources, 0, 0, NUM_STATE, ng);

}

// Obtain the effective source term due to reactions on the primitive variables.

#ifdef REACTIONS
#ifdef SDC
void
Castro::get_react_source_prim(MultiFab& react_src, Real dt)
{
    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

    int ng = 0;

    // Carries the contribution of all non-reacting source terms.

    MultiFab A(grids, dmap, NUM_STATE, ng);

    sum_of_sources(A);

    // Compute the state that has effectively only been updated with advection.
    // U* = U_old + dt A
    // where A = -div U + S_hydro
    MultiFab S_noreact(grids, dmap, NUM_STATE, ng);

    MultiFab::Copy(S_noreact, S_old, 0, 0, NUM_STATE, ng);
    MultiFab::Saxpy(S_noreact, dt, A, 0, 0, NUM_STATE, ng);

    clean_state(S_noreact);

    // Compute its primitive counterpart, q*

    MultiFab q_noreact(grids, dmap, QVAR, ng);
    MultiFab qaux_noreact(grids, dmap, NQAUX, ng);

    cons_to_prim(S_noreact, q_noreact, qaux_noreact);

    // Compute the primitive version of the old state, q_old

    MultiFab q_old(grids, dmap, QVAR, ng);
    MultiFab qaux_old(grids, dmap, NQAUX, ng);

    cons_to_prim(S_old, q_old, qaux_old);

    // Compute the effective advective update on the primitive state.
    // A(q) = (q* - q_old)/dt

    MultiFab A_prim(grids, dmap, QVAR, ng);

    A_prim.setVal(0.0);

    if (dt > 0.0) {
        MultiFab::Saxpy(A_prim,  1.0 / dt, q_noreact, 0, 0, QVAR, ng);
	MultiFab::Saxpy(A_prim, -1.0 / dt, q_old,     0, 0, QVAR, ng);
    }

    // Compute the primitive version of the new state.

    MultiFab q_new(grids, dmap, QVAR, ng);
    MultiFab qaux_new(grids, dmap, NQAUX, ng);

    cons_to_prim(S_new, q_new, qaux_new);

    // Compute the reaction source term.

    react_src.setVal(0.0, react_src.nGrow());

    if (dt > 0.0) {
        MultiFab::Saxpy(react_src,  1.0 / dt, q_new, 0, 0, QVAR, ng);
        MultiFab::Saxpy(react_src, -1.0 / dt, q_old, 0, 0, QVAR, ng);
    }

    MultiFab::Saxpy(react_src, -1.0, A_prim, 0, 0, QVAR, ng);

    // Now fill all of the ghost zones.
    Real time = get_state_data(SDC_React_Type).curTime();
    AmrLevel::FillPatch(*this, react_src, react_src.nGrow(), time, SDC_React_Type, 0, NUM_STATE);

}
#endif
#endif
