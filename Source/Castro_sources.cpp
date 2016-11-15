#include "Castro.H"
#include "Castro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

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
Castro::do_old_sources(Real time, Real dt, int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle)
{

    // Construct the old-time sources.

    for (int n = 0; n < num_src; ++n)
	construct_old_source(n, time, dt, amr_iteration, amr_ncycle,
			     sub_iteration, sub_ncycle);

    // Apply the old-time sources directly to the new-time state,
    // S_new -- note that this addition is for full dt, since we
    // will do a predictor-corrector on the sources to allow for
    // state-dependent sources.

    MultiFab& S_new = get_new_data(State_Type);

    for (int n = 0; n < num_src; ++n)
	if (source_flag(n))
	    apply_source_to_state(S_new, old_sources[n], dt);
}

void
Castro::do_new_sources(Real time, Real dt, int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle)
{

    MultiFab& S_new = get_new_data(State_Type);

    // For the new-time source terms, we have an option for how to proceed.
    // We can either construct all of the old-time sources using the same
    // state that comes out of the hydro update, or we can evaluate the sources
    // one by one and apply them as we go.

    if (update_state_between_sources) {

	for (int n = 0; n < num_src; ++n) {
	    construct_new_source(n, time, dt, amr_iteration, amr_ncycle, sub_iteration, sub_ncycle);
	    if (source_flag(n)) {
		apply_source_to_state(S_new, new_sources[n], dt);
		clean_state(S_new);
	    }
	}

    } else {

	// Construct the new-time source terms.

	for (int n = 0; n < num_src; ++n)
	    construct_new_source(n, time, dt, amr_iteration, amr_ncycle, sub_iteration, sub_ncycle);

	// Apply the new-time sources to the state.

	for (int n = 0; n < num_src; ++n)
	    if (source_flag(n))
		apply_source_to_state(S_new, new_sources[n], dt);

	clean_state(S_new);

    }

}

void
Castro::construct_old_source(int src, Real time, Real dt, int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle)
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
Castro::construct_new_source(int src, Real time, Real dt, int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle)
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

// Obtain the sum of all source terms.

void
Castro::sum_of_sources(MultiFab& source)
{

  // this computes advective_source + 1/2 (old source + new source)
  //
  // the time-centering is accomplished since new source is defined
  // to be 1/2 (new source - old source) generally.

  int ng = source.nGrow();

  source.setVal(0.0);

  for (int n = 0; n < num_src; ++n)
      MultiFab::Add(source, old_sources[n], 0, 0, NUM_STATE, ng);

  MultiFab::Add(source, hydro_source, 0, 0, NUM_STATE, ng);

  for (int n = 0; n < num_src; ++n)
      MultiFab::Add(source, new_sources[n], 0, 0, NUM_STATE, ng);

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

    MultiFab A(grids, NUM_STATE, ng, Fab_allocate);

    sum_of_sources(A);

    // Compute the state that has effectively only been updated with advection.

    MultiFab S_noreact(grids, NUM_STATE, ng, Fab_allocate);

    MultiFab::Copy(S_noreact, S_old, 0, 0, NUM_STATE, ng);
    MultiFab::Saxpy(S_noreact, dt, A, 0, 0, NUM_STATE, ng);

    clean_state(S_noreact);

    // Compute its primitive counterpart.

    MultiFab q_noreact(grids, QVAR, ng, Fab_allocate);
    MultiFab qaux_noreact(grids, NQAUX, ng, Fab_allocate);

    cons_to_prim(S_noreact, q_noreact, qaux_noreact);

    // Compute the primitive version of the old state.

    MultiFab q_old(grids, QVAR, ng, Fab_allocate);
    MultiFab qaux_old(grids, NQAUX, ng, Fab_allocate);

    cons_to_prim(S_old, q_old, qaux_old);

    // Compute the effective advective update on the primitive state.

    MultiFab A_prim(grids, QVAR, ng, Fab_allocate);

    A_prim.setVal(0.0);

    if (dt > 0.0) {
        MultiFab::Saxpy(A_prim, -1.0 / dt, q_noreact, 0, 0, QVAR, ng);
	MultiFab::Saxpy(A_prim,  1.0 / dt, q_old,     0, 0, QVAR, ng);
    }

    // Compute the primitive version of the new state.

    MultiFab q_new(grids, QVAR, ng, Fab_allocate);
    MultiFab qaux_new(grids, NQAUX, ng, Fab_allocate);

    cons_to_prim(S_new, q_new, qaux_new);

    // Compute the reaction source term.

    react_src.setVal(0.0, react_src.nGrow());

    if (dt > 0.0) {
        MultiFab::Saxpy(react_src,  1.0 / dt, q_new, 0, 0, QVAR, ng);
        MultiFab::Saxpy(react_src, -1.0 / dt, q_old, 0, 0, QVAR, ng);
    }

    MultiFab::Add(react_src, A_prim, 0, 0, QVAR, ng);

    // Now fill all of the ghost zones.
    react_src.FillBoundary(geom.periodicity());
}
#endif
#endif
