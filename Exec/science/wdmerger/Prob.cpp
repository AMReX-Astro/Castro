#include "Castro.H"
#include "Castro_F.H"
#include "Castro_prob_err_F.H"

#include "Gravity.H"
#include <Gravity_F.H>

#include "AMReX_ParmParse.H"
#include "AMReX_buildInfo.H"

#include <fstream>

using namespace amrex;

int Castro::relaxation_is_done = 0;
int Castro::problem = -1;
int Castro::use_stopping_criterion = 1;
int Castro::use_energy_stopping_criterion = 0;
Real Castro::ts_te_stopping_criterion = 1.e200;

Real Castro::mass_p = 0.0;
Real Castro::mass_s = 0.0;

Real Castro::mdot_p = 0.0;
Real Castro::mdot_s = 0.0;

Real Castro::com_p[3] = { 0.0 };
Real Castro::com_s[3] = { 0.0 };

Real Castro::vel_p[3] = { 0.0 };
Real Castro::vel_s[3] = { 0.0 };

Real Castro::rad_p[7] = { 0.0 };
Real Castro::rad_s[7] = { 0.0 };

Real Castro::vol_p[7] = { 0.0 };
Real Castro::vol_s[7] = { 0.0 };

Real Castro::rho_avg_p = 0.0;
Real Castro::rho_avg_s = 0.0;

Real Castro::t_ff_p = 0.0;
Real Castro::t_ff_s = 0.0;

Real Castro::T_global_max = 0.0;
Real Castro::rho_global_max = 0.0;
Real Castro::ts_te_global_max = 0.0;

Real Castro::T_curr_max = 0.0;
Real Castro::rho_curr_max = 0.0;
Real Castro::ts_te_curr_max = 0.0;

Real Castro::total_ener_array[num_previous_ener_timesteps] = { 0.0 };

int Castro::num_zones_ignited = 0;
int Castro::ignition_level = -1;

#ifdef DO_PROBLEM_POST_TIMESTEP
void
Castro::problem_post_timestep()
{

    if (level != 0) return;

    int finest_level = parent->finestLevel();
    Real time = state[State_Type].curTime();
    Real dt = parent->dtLevel(0);

    if (time == 0.0) dt = 0.0; // dtLevel returns the next timestep for t = 0, so overwrite

    // Update white dwarf masses, positions, velocities, and auxiliary data.

    wd_update(time, dt);

    // If we are doing problem 3, which has an initial relaxation step,
    // perform any post-timestep updates to assist with the relaxation, then
    // determine whether the criterion for terminating the relaxation
    // has been satisfied.

    update_relaxation(time, dt);

    // Update extrema on the domain.

    update_extrema(time);

    // Some of the problems might have stopping conditions that depend on
    // the state of the simulation; those are checked here.

    check_to_stop(time);

}
#endif



//
// This function updates the WD data, including the masses, center-of-mass locations,
// and velocities of the center of masses of the primary and secondary white dwarfs.
//

void
Castro::wd_update (Real time, Real dt)
{
    BL_PROFILE("Castro::wd_update()");

    // Ensure we are either on the coarse level, or on the finest level
    // when we are not doing subcycling. The data should be sychronized
    // in both of these cases.

    BL_ASSERT(level == 0 || (!parent->subCycle() && level == parent->finestLevel()));

    // Get the current stellar data
    get_star_data(com_p, com_s, vel_p, vel_s, &mass_p, &mass_s, &t_ff_p, &t_ff_s);

    // Update the problem center using the system bulk velocity
    update_center(&time);

    for ( int i = 0; i < 3; i++ ) {
      com_p[i] += vel_p[i] * dt;
      com_s[i] += vel_s[i] * dt;
    }

    // Now send this first estimate of the COM to Fortran, and then re-calculate
    // a more accurate result using it as a starting point.

    set_star_data(com_p, com_s, vel_p, vel_s, &mass_p, &mass_s, &t_ff_p, &t_ff_s);

    // Save relevant current data.

    Real old_mass_p = mass_p;
    Real old_mass_s = mass_s;

    Real com_p_x = 0.0;
    Real com_p_y = 0.0;
    Real com_p_z = 0.0;
    Real com_s_x = 0.0;
    Real com_s_y = 0.0;
    Real com_s_z = 0.0;
    Real vel_p_x = 0.0;
    Real vel_p_y = 0.0;
    Real vel_p_z = 0.0;
    Real vel_s_x = 0.0;
    Real vel_s_y = 0.0;
    Real vel_s_z = 0.0;
    Real mp      = 0.0;
    Real ms      = 0.0;

    for (int lev = 0; lev <= parent->finestLevel(); lev++) {

      ca_set_amr_info(lev, -1, -1, -1.0, -1.0);

      Castro& c_lev = getLevel(lev);

      const Real* dx = c_lev.geom.CellSize();

      // Density and momenta

      auto mfrho  = c_lev.derive("density",time,0);
      auto mfxmom = c_lev.derive("xmom",time,0);
      auto mfymom = c_lev.derive("ymom",time,0);
      auto mfzmom = c_lev.derive("zmom",time,0);

      // Masks for the primary and secondary

      auto mfpmask = c_lev.derive("primarymask", time, 0);
      auto mfsmask = c_lev.derive("secondarymask", time, 0);

      BL_ASSERT(mfrho   != nullptr);
      BL_ASSERT(mfxmom  != nullptr);
      BL_ASSERT(mfymom  != nullptr);
      BL_ASSERT(mfzmom  != nullptr);
      BL_ASSERT(mfpmask != nullptr);
      BL_ASSERT(mfsmask != nullptr);

      if (lev < parent->finestLevel())
      {
          const MultiFab& mask = getLevel(lev+1).build_fine_mask();

	  MultiFab::Multiply(*mfrho,   mask, 0, 0, 1, 0);
	  MultiFab::Multiply(*mfxmom,  mask, 0, 0, 1, 0);
	  MultiFab::Multiply(*mfymom,  mask, 0, 0, 1, 0);
	  MultiFab::Multiply(*mfzmom,  mask, 0, 0, 1, 0);
	  MultiFab::Multiply(*mfpmask, mask, 0, 0, 1, 0);
	  MultiFab::Multiply(*mfsmask, mask, 0, 0, 1, 0);
      }

#ifdef _OPENMP
#pragma omp parallel reduction(+:com_p_x,com_p_y,com_p_z,com_s_x,com_s_y,com_s_z) \
                     reduction(+:vel_p_x,vel_p_y,vel_p_z,vel_s_x,vel_s_y,vel_s_z) \
                     reduction(+:mp, ms)
#endif
      for (MFIter mfi(*mfrho,true); mfi.isValid(); ++mfi) {
          FArrayBox& fabrho   = (*mfrho )[mfi];
	  FArrayBox& fabxmom  = (*mfxmom)[mfi];
	  FArrayBox& fabymom  = (*mfymom)[mfi];
	  FArrayBox& fabzmom  = (*mfzmom)[mfi];
	  FArrayBox& fabpmask = (*mfpmask)[mfi];
	  FArrayBox& fabsmask = (*mfsmask)[mfi];
	  FArrayBox& vol      = c_lev.volume[mfi];

	  const Box& box  = mfi.tilebox();
	  const int* lo   = box.loVect();
	  const int* hi   = box.hiVect();

	  wdcom(BL_TO_FORTRAN_3D(fabrho),
		BL_TO_FORTRAN_3D(fabxmom),
		BL_TO_FORTRAN_3D(fabymom),
		BL_TO_FORTRAN_3D(fabzmom),
		BL_TO_FORTRAN_3D(fabpmask),
		BL_TO_FORTRAN_3D(fabsmask),
		BL_TO_FORTRAN_3D(vol),
		ARLIM_3D(lo),ARLIM_3D(hi),
		ZFILL(dx),&time,
		&com_p_x, &com_p_y, &com_p_z,
		&com_s_x, &com_s_y, &com_s_z,
		&vel_p_x, &vel_p_y, &vel_p_z,
		&vel_s_x, &vel_s_y, &vel_s_z,
		&mp, &ms);
	}

    }

    ca_set_amr_info(level, -1, -1, -1.0, -1.0);

    com_p[0] = com_p_x;
    com_p[1] = com_p_y;
    com_p[2] = com_p_z;
    com_s[0] = com_s_x;
    com_s[1] = com_s_y;
    com_s[2] = com_s_z;
    vel_p[0] = vel_p_x;
    vel_p[1] = vel_p_y;
    vel_p[2] = vel_p_z;
    vel_s[0] = vel_s_x;
    vel_s[1] = vel_s_y;
    vel_s[2] = vel_s_z;
    mass_p   = mp;
    mass_s   = ms;

    bool local_flag = true;

    // Compute effective radii of stars at various density cutoffs

    for (int i = 0; i <= 6; ++i)
        Castro::volInBoundary(time, vol_p[i], vol_s[i], pow(10.0,i), local_flag);

    // Do all of the reductions.

    const int nfoo_sum = 28;
    Real foo_sum[nfoo_sum] = { 0.0 };

    for (int i = 0; i <= 6; ++i) {
      foo_sum[i  ] = vol_p[i];
      foo_sum[i+7] = vol_s[i];
    }

    foo_sum[14] = mass_p;
    foo_sum[15] = mass_s;

    for (int i = 0; i <= 2; ++i) {
      foo_sum[i+16] = com_p[i];
      foo_sum[i+19] = com_s[i];
      foo_sum[i+22] = vel_p[i];
      foo_sum[i+25] = vel_s[i];
    }

    amrex::ParallelDescriptor::ReduceRealSum(foo_sum, nfoo_sum);

    for (int i = 0; i <= 6; ++i) {
      vol_p[i] = foo_sum[i  ];
      vol_s[i] = foo_sum[i+7];
    }

    mass_p = foo_sum[14];
    mass_s = foo_sum[15];

    for (int i = 0; i <= 2; ++i) {
      com_p[i] = foo_sum[i+16];
      com_s[i] = foo_sum[i+19];
      vel_p[i] = foo_sum[i+22];
      vel_s[i] = foo_sum[i+25];
    }

    // Compute effective WD radii

    for (int i = 0; i <= 6; ++i) {

        rad_p[i] = std::pow(vol_p[i] * 3.0 / 4.0 / M_PI, 1.0/3.0);
        rad_s[i] = std::pow(vol_s[i] * 3.0 / 4.0 / M_PI, 1.0/3.0);

    }

     // Complete calculations for center of mass quantities

    for ( int i = 0; i < 3; i++ ) {

      if ( mass_p > 0.0 ) {
        com_p[i] = com_p[i] / mass_p;
        vel_p[i] = vel_p[i] / mass_p;
      }

      if ( mass_s > 0.0 ) {
        com_s[i] = com_s[i] / mass_s;
        vel_s[i] = vel_s[i] / mass_s;
      }

    }

    // For 1D we force the masses to remain constant

#if (BL_SPACEDIM == 1)
    mass_p = old_mass_p;
    mass_s = old_mass_s;
#endif

    if (mass_p > 0.0 && dt > 0.0)
      mdot_p = (mass_p - old_mass_p) / dt;
    else
      mdot_p = 0.0;

    if (mass_s > 0.0 && dt > 0.0)
      mdot_s = (mass_s - old_mass_s) / dt;
    else
      mdot_s = 0.0;

    // Free-fall timescale ~ 1 / sqrt(G * rho_avg}

    Real Gconst;

    get_grav_const(&Gconst);

    if (mass_p > 0.0 && vol_p[2] > 0.0) {
      rho_avg_p = mass_p / vol_p[2];
      t_ff_p = sqrt(3.0 * M_PI / (32.0 * Gconst * rho_avg_p));
    }

    if (mass_s > 0.0 && vol_s[2] > 0.0) {
      rho_avg_s = mass_s / vol_s[2];
      t_ff_s = sqrt(3.0 * M_PI / (32.0 * Gconst * rho_avg_s));
    }

    // Send this updated information back to the Fortran probdata module

    set_star_data(com_p, com_s, vel_p, vel_s, &mass_p, &mass_s, &t_ff_p, &t_ff_s);

}



// This function uses the known center of mass of the two white dwarfs,
// and given a density cutoff, computes the total volume of all zones
// whose density is greater or equal to that density cutoff.
// We also impose a distance requirement so that we only look
// at zones that are within twice the original radius of the white dwarf.

void Castro::volInBoundary (Real time, Real& vol_p, Real& vol_s, Real rho_cutoff, bool local)
{
    BL_PROFILE("Castro::volInBoundary()");

    BL_ASSERT(level == 0);

    vol_p = 0.0;
    vol_s = 0.0;

    for (int lev = 0; lev <= parent->finestLevel(); lev++) {

      ca_set_amr_info(lev, -1, -1, -1.0, -1.0);

      Castro& c_lev = getLevel(lev);

      const Real* dx = c_lev.geom.CellSize();
      auto mf = c_lev.derive("density",time,0);

      // Effective potentials of the primary and secondary

      auto mfpmask = c_lev.derive("primarymask", time, 0);
      auto mfsmask = c_lev.derive("secondarymask", time, 0);

      BL_ASSERT(mf      != nullptr);
      BL_ASSERT(mfpmask != nullptr);
      BL_ASSERT(mfsmask != nullptr);

      if (lev < parent->finestLevel())
      {
	  const MultiFab& mask = c_lev.getLevel(lev+1).build_fine_mask();
	  MultiFab::Multiply(*mf,      mask, 0, 0, 1, 0);
	  MultiFab::Multiply(*mfpmask, mask, 0, 0, 1, 0);
	  MultiFab::Multiply(*mfsmask, mask, 0, 0, 1, 0);
      }

      Real vp = 0.0;
      Real vs = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:vp,vs)
#endif
      for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi) {
          FArrayBox& fab      = (*mf)[mfi];
          FArrayBox& fabpmask = (*mfpmask)[mfi];
          FArrayBox& fabsmask = (*mfsmask)[mfi];
          FArrayBox& vol      = c_lev.volume[mfi];

	  Real sp = 0.0;
	  Real ss = 0.0;

	  const Box& box  = mfi.tilebox();
	  const int* lo   = box.loVect();
	  const int* hi   = box.hiVect();

	  ca_volumeindensityboundary(BL_TO_FORTRAN_3D(fab),
		                     BL_TO_FORTRAN_3D(fabpmask),
				     BL_TO_FORTRAN_3D(fabsmask),
				     BL_TO_FORTRAN_3D(vol),
				     ARLIM_3D(lo),ARLIM_3D(hi),
				     ZFILL(dx),&sp,&ss,&rho_cutoff);
	  vp += sp;
	  vs += ss;
      }

      vol_p += vp;
      vol_s += vs;

    }

    if (!local) {

      const int nfoo = 2;
      Real foo[nfoo] = {vol_p, vol_s};

      amrex::ParallelDescriptor::ReduceRealSum(foo, nfoo);

      vol_p = foo[0];
      vol_s = foo[1];

    }

    ca_set_amr_info(level, -1, -1, -1.0, -1.0);

}



//
// Calculate the gravitational wave signal.
//

void
Castro::gwstrain (Real time,
		  Real& h_plus_1, Real& h_cross_1,
		  Real& h_plus_2, Real& h_cross_2,
		  Real& h_plus_3, Real& h_cross_3,
		  bool local) {

    BL_PROFILE("Castro::gwstrain()");

    const Real* dx = geom.CellSize();

    auto mfrho   = derive("density",time,0);
    auto mfxmom  = derive("xmom",time,0);
    auto mfymom  = derive("ymom",time,0);
    auto mfzmom  = derive("zmom",time,0);
    auto mfgravx = derive("grav_x",time,0);
    auto mfgravy = derive("grav_y",time,0);
    auto mfgravz = derive("grav_z",time,0);

    BL_ASSERT(mfrho   != nullptr);
    BL_ASSERT(mfxmom  != nullptr);
    BL_ASSERT(mfymom  != nullptr);
    BL_ASSERT(mfzmom  != nullptr);
    BL_ASSERT(mfgravx != nullptr);
    BL_ASSERT(mfgravy != nullptr);
    BL_ASSERT(mfgravz != nullptr);

    if (level < parent->finestLevel())
    {
	const MultiFab& mask = getLevel(level+1).build_fine_mask();

	MultiFab::Multiply(*mfrho,   mask, 0, 0, 1, 0);
	MultiFab::Multiply(*mfxmom,  mask, 0, 0, 1, 0);
	MultiFab::Multiply(*mfymom,  mask, 0, 0, 1, 0);
	MultiFab::Multiply(*mfzmom,  mask, 0, 0, 1, 0);
	MultiFab::Multiply(*mfgravx, mask, 0, 0, 1, 0);
	MultiFab::Multiply(*mfgravy, mask, 0, 0, 1, 0);
	MultiFab::Multiply(*mfgravz, mask, 0, 0, 1, 0);
    }

    // Qtt stores the second time derivative of the quadrupole moment.
    // We calculate it directly rather than computing the quadrupole moment
    // and differentiating it in time, because the latter method is less accurate
    // and requires the state at other timesteps. See, e.g., Equation 5 of
    // Loren-Aguilar et al. 2005.

    // It is a 3x3 rank-2 tensor, but AMReX expects IntVect() to use BL_SPACEDIM
    // dimensions, so we add a redundant third index in 3D.

    Box bx( IntVect(D_DECL(0, 0, 0)), IntVect(D_DECL(2, 2, 0)) );

    FArrayBox Qtt(bx);

    Qtt.setVal(0.0);

#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    Array< std::unique_ptr<FArrayBox> > priv_Qtt(nthreads);
    for (int i=0; i<nthreads; i++) {
	priv_Qtt[i].reset(new FArrayBox(bx));
    }
#pragma omp parallel
#endif
    {
#ifdef _OPENMP
	int tid = omp_get_thread_num();
	priv_Qtt[tid]->setVal(0.0);
#endif
	for (MFIter mfi(*mfrho,true); mfi.isValid(); ++mfi) {

	    const Box& box  = mfi.tilebox();
	    const int* lo   = box.loVect();
	    const int* hi   = box.hiVect();

	    quadrupole_tensor_double_dot(BL_TO_FORTRAN_3D((*mfrho)[mfi]),
					 BL_TO_FORTRAN_3D((*mfxmom)[mfi]),
					 BL_TO_FORTRAN_3D((*mfymom)[mfi]),
					 BL_TO_FORTRAN_3D((*mfzmom)[mfi]),
					 BL_TO_FORTRAN_3D((*mfgravx)[mfi]),
					 BL_TO_FORTRAN_3D((*mfgravy)[mfi]),
					 BL_TO_FORTRAN_3D((*mfgravz)[mfi]),
					 BL_TO_FORTRAN_3D(volume[mfi]),
					 ARLIM_3D(lo),ARLIM_3D(hi),ZFILL(dx),&time,
#ifdef _OPENMP
					 priv_Qtt[tid]->dataPtr());
#else
	                                 Qtt.dataPtr());
#endif
        }
    }

    // Do an OpenMP reduction on the tensor.

#ifdef _OPENMP
        int n = bx.numPts();
	Real* p = Qtt.dataPtr();
#pragma omp barrier
#pragma omp for nowait
	for (int i=0; i<n; ++i)
	{
	    for (int it=0; it<nthreads; it++) {
		const Real* pq = priv_Qtt[it]->dataPtr();
		p[i] += pq[i];
	    }
	}
#endif

    // Now, do a global reduce over all processes.

    if (!local)
	amrex::ParallelDescriptor::ReduceRealSum(Qtt.dataPtr(),bx.numPts());

    // Now that we have the second time derivative of the quadrupole
    // tensor, we can calculate the transverse-trace gauge strain tensor.

    gw_strain_tensor(&h_plus_1, &h_cross_1,
		     &h_plus_2, &h_cross_2,
		     &h_plus_3, &h_cross_3,
		     Qtt.dataPtr(), &time);

}



// Computes standard dot-product of two three-vectors.

Real Castro::dot_product(const Real a[], const Real b[]) {

  Real c = 0.0;

  c = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

  return c;

}



// Computes standard cross-product of two three-vectors.

void Castro::cross_product(const Real a[], const Real b[], Real c[]) {

  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];

}



// Computes norm of a three-vector.

Real Castro::norm(const Real a[]) {

  Real n = 0.0;

  n = sqrt( dot_product(a, a) );

  return n;

}



#ifdef DO_PROBLEM_POST_INIT

void Castro::problem_post_init() {

  // Read in inputs.

  ParmParse pp("castro");

  pp.query("use_stopping_criterion", use_stopping_criterion);
  pp.query("use_energy_stopping_criterion", use_energy_stopping_criterion);
  pp.query("ts_te_stopping_criterion", ts_te_stopping_criterion);

  // Get the problem number fom Fortran.

  get_problem_number(&problem);

  // Get the relaxation status.

  get_relaxation_status(&relaxation_is_done);

  // If we're doing an initial relaxation step, ensure that we are not subcycling.

  if (problem == 3 && !relaxation_is_done && parent->subCycle())
    amrex::Abort("Error: cannot perform relaxation step if we are sub-cycling in the AMR.");

  // If we're doing an initial relaxation step, ensure that we retain source terms.

  if (problem == 3 && !relaxation_is_done && !keep_sources_until_end)
    amrex::Abort("Error: cannot perform relaxation step if we are not retaining source terms.");

  // Update the rotational period; some problems change this from what's in the inputs parameters.

  get_period(&rotational_period);

  // Initialize the energy storage array.

  for (int i = 0; i < num_previous_ener_timesteps; ++i)
    total_ener_array[i] = -1.e200;

  set_total_ener_array(total_ener_array);

  // Execute the post timestep diagnostics here,
  // so that the results at t = 0 and later are smooth.
  // This should generally be the last operation
  // in this function.

  problem_post_timestep();

}

#endif



#ifdef DO_PROBLEM_POST_RESTART

void Castro::problem_post_restart() {

  // Read in inputs.

  ParmParse pp("castro");

  pp.query("use_stopping_criterion", use_stopping_criterion);
  pp.query("use_energy_stopping_criterion", use_energy_stopping_criterion);
  pp.query("ts_te_stopping_criterion", ts_te_stopping_criterion);

  // Get the problem number from Fortran.

  get_problem_number(&problem);

  // Get the relaxation status.

  get_relaxation_status(&relaxation_is_done);

  // Get the rotational period.

  get_period(&rotational_period);

  // Get the energy data from Fortran.

  get_total_ener_array(total_ener_array);

  // Get the extrema.

  get_extrema(&T_global_max, &rho_global_max, &ts_te_global_max);

  T_curr_max = T_global_max;
  rho_curr_max = rho_global_max;
  ts_te_curr_max = ts_te_global_max;

  // Get the ignition status.

  get_num_zones_ignited(&num_zones_ignited, &ignition_level);

  // If we're restarting from a checkpoint at t = 0 but don't yet
  // have diagnostics, we want to generate the headers and the t = 0
  // data at this time so that future timestep diagnostics make sense.

  Real time = state[State_Type].curTime();

  if (time == 0.0) {

      if (parent->NumDataLogs() > 0) {

          bool do_sums = false;

          const std::string datalogname = parent->DataLogName(0);

          std::ifstream log;
          log.open(datalogname, std::ios::ate);

          if (log.tellg() == 0)
              do_sums = true;
          log.close();

          if (do_sums)
              sum_integrated_quantities();

      }

  }

}

#endif



void Castro::writeGitHashes(std::ostream& log) {

  const char* castro_hash       = buildInfoGetGitHash(1);
  const char* amrex_hash        = buildInfoGetGitHash(2);
  const char* microphysics_hash = buildInfoGetGitHash(3);

  log << "# Castro       git hash: " << castro_hash       << std::endl;
  log << "# AMReX        git hash: " << amrex_hash        << std::endl;
  log << "# Microphysics git hash: " << microphysics_hash << std::endl;

}



void Castro::check_to_stop(Real time) {

    int jobDoneStatus;

    // Get the current job done status.

    get_job_status(&jobDoneStatus);

#if (BL_SPACEDIM > 1)
    if (use_stopping_criterion) {

      if (problem == 0) {

          if (use_energy_stopping_criterion) {

              // For the collision problem, we know we are done when the total energy
              // is positive (indicating that we have become unbound due to nuclear
              // energy release) and when it is decreasing in magnitude (indicating
              // all of the excitement is done and fluid is now just streaming off
              // the grid). We don't need to be super accurate for this, so let's check
              // on the coarse grid only. It is possible that a collision could not
              // generate enough energy to become unbound, so possibly this criterion
              // should be expanded in the future to cover that case.

              Real rho_E = 0.0;
              Real rho_phi = 0.0;

              // Note that we'll define the total energy using only
              // gas energy + gravitational. Rotation is never on
              // for the collision problem so we can ignore it.

              Real E_tot = 0.0;

              Real curTime   = state[State_Type].curTime();

              bool local_flag = true;
              bool fine_mask = false;

              rho_E += volWgtSum("rho_E", curTime,  local_flag, fine_mask);

#ifdef GRAVITY
              if (do_grav) {
                  rho_phi += volWgtSum("rho_phiGrav", curTime,  local_flag, fine_mask);
              }
#endif

              E_tot = rho_E + 0.5 * rho_phi;

              amrex::ParallelDescriptor::ReduceRealSum(E_tot);

              // Put this on the end of the energy array.

              for (int i = num_previous_ener_timesteps - 1; i > 0; --i)
                  total_ener_array[i] = total_ener_array[i - 1];

              total_ener_array[0] = E_tot;

              // Send the data to Fortran.

              set_total_ener_array(total_ener_array);

              bool stop_flag = false;

              int i = 0;

              // Check if energy is positive and has been decreasing for at least the last few steps.

              while (i < num_previous_ener_timesteps - 1) {

                  if (total_ener_array[i] < 0.0)
                      break;
                  else if (total_ener_array[i] > total_ener_array[i + 1])
                      break;

                  ++i;

              }

              if (i == num_previous_ener_timesteps - 1)
                  stop_flag = true;

              if (stop_flag) {

                  jobDoneStatus = 1;

                  set_job_status(&jobDoneStatus);

                  amrex::Print() << std::endl 
                                 << "Ending simulation because total energy is positive and decreasing." 
                                 << std::endl;

              }

          }

          if (ts_te_curr_max >= ts_te_stopping_criterion) {

              jobDoneStatus = 1;

              set_job_status(&jobDoneStatus);

              amrex::Print() << std::endl
                             << "Ending simulation because we are above the threshold for unstable burning."
                             << std::endl;

          }

      } else if (problem == 4) {

	// We can work out the stopping time using the formula
	// t_freefall = rotational_period / (4 * sqrt(2)).
	// We'll stop 90% of the way there because that's about
	// when the stars start coming into contact, and the
	// assumption of spherically symmetric stars breaks down.

	Real stopping_time = 0.90 * rotational_period / (4.0 * std::sqrt(2));

	if (time >= stopping_time) {

	  jobDoneStatus = 1;

	  set_job_status(&jobDoneStatus);

	}

      }

    }
#endif

    // Is the job done? If so, signal this to AMReX.

    get_job_status(&jobDoneStatus);

    if (jobDoneStatus == 1) {

      // Write out a checkpoint. Note that this will
      // only happen if you have amr.message_int = 1.

      if (amrex::ParallelDescriptor::IOProcessor()) {
	std::ofstream dump_file;
	dump_file.open("dump_and_stop", std::ofstream::out);
	dump_file.close();
      }

    }

}



void Castro::update_extrema(Real time) {

    // Compute extrema

    bool local_flag = true;

    T_curr_max     = 0.0;
    rho_curr_max   = 0.0;
    ts_te_curr_max = 0.0;

    int finest_level = parent->finestLevel();

    for (int lev = 0; lev <= finest_level; lev++) {

      MultiFab& S_new = parent->getLevel(lev).get_new_data(State_Type);

      T_curr_max = std::max(T_curr_max, S_new.max(Temp, 0, local_flag));
      rho_curr_max = std::max(rho_curr_max, S_new.max(Density, 0, local_flag));

#ifdef REACTIONS
      if (lev == finest_level) {

        auto ts_te_MF = parent->getLevel(lev).derive("t_sound_t_enuc", time, 0);
	ts_te_curr_max = std::max(ts_te_curr_max, ts_te_MF->max(0,0,local_flag));

      }
#endif

    }

    // Max reductions

    const int nfoo_max = 3;

    Real foo_max[3];

    foo_max[0] = T_curr_max;
    foo_max[1] = rho_curr_max;
    foo_max[2] = ts_te_curr_max;

    amrex::ParallelDescriptor::ReduceRealMax(foo_max, nfoo_max);

    T_curr_max     = foo_max[0];
    rho_curr_max   = foo_max[1];
    ts_te_curr_max = foo_max[2];

    T_global_max     = std::max(T_global_max, T_curr_max);
    rho_global_max   = std::max(rho_global_max, rho_curr_max);
    ts_te_global_max = std::max(ts_te_global_max, ts_te_curr_max);

    // Send extrema data to Fortran

    set_extrema(&T_global_max, &rho_global_max, &ts_te_global_max);

}



void
Castro::update_relaxation(Real time, Real dt) {

    // Check to make sure whether we should be doing the relaxation here.

    if (problem != 3 || relaxation_is_done || mass_p <= 0.0 || mass_s <= 0.0 || dt <= 0.0) return;

    Real old_time = state[State_Type].prevTime();
    Real new_time = state[State_Type].curTime();

    int ng = 0;

    // Construct the desired rotation force data.

    int finest_level = parent->finestLevel();
    int n_levs = finest_level + 1;

    Array< std::unique_ptr<MultiFab> > rot_force(n_levs);

    for (int lev = 0; lev <= finest_level; ++lev) {

	rot_force[lev].reset(new MultiFab(getLevel(lev).grids, getLevel(lev).dmap, NUM_STATE, ng));

	rot_force[lev]->setVal(0.0);

	MultiFab::Add(*rot_force[lev], *(getLevel(lev).old_sources[grav_src]), 0, 0, NUM_STATE, ng);
        MultiFab::Add(*rot_force[lev], *(getLevel(lev).new_sources[grav_src]), 0, 0, NUM_STATE, ng);

	MultiFab::Add(*rot_force[lev], getLevel(lev).hydro_source, 0, 0, NUM_STATE, ng);

	// Mask out regions covered by fine grids.

	if (lev < finest_level) {
	    const MultiFab& mask = getLevel(lev+1).build_fine_mask();
	    MultiFab::Multiply(*rot_force[lev], mask, 0, 0, NUM_STATE, 0);
	}

    }

    // Construct omega from this.

    Real force_p[3] = { 0.0 };
    Real force_s[3] = { 0.0 };

    for (int lev = 0; lev <= finest_level; ++lev) {

	MultiFab& S_new = getLevel(lev).get_new_data(State_Type);

	auto pmask = getLevel(lev).derive("primarymask", time, 0);
	auto smask = getLevel(lev).derive("secondarymask", time, 0);

	MultiFab& vol = getLevel(lev).Volume();

	Real fpx = 0.0;
	Real fpy = 0.0;
	Real fpz = 0.0;
	Real fsx = 0.0;
	Real fsy = 0.0;
	Real fsz = 0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:fpx, fpy, fpz, fsx, fsy, fsz)
#endif
	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {

	    const Box& box = mfi.tilebox();

	    const int* lo  = box.loVect();
	    const int* hi  = box.hiVect();

	    sum_force_on_stars(lo, hi,
			       BL_TO_FORTRAN_3D((*rot_force[lev])[mfi]),
			       BL_TO_FORTRAN_3D(S_new[mfi]),
			       BL_TO_FORTRAN_3D(vol[mfi]),
			       BL_TO_FORTRAN_3D((*pmask)[mfi]),
                               BL_TO_FORTRAN_3D((*smask)[mfi]),
			       &fpx, &fpy, &fpz, &fsx, &fsy, &fsz);

	}

	force_p[0] += fpx;
	force_p[1] += fpy;
	force_p[2] += fpz;

	force_s[0] += fsx;
	force_s[1] += fsy;
	force_s[2] += fsz;

    }

    Real foo[6];

    // Do the reduction over processors, and then
    // normalize by the masses of the stars.

    foo[0] = force_p[0];
    foo[1] = force_p[1];
    foo[2] = force_p[2];

    foo[3] = force_s[0];
    foo[4] = force_s[1];
    foo[5] = force_s[2];

    amrex::ParallelDescriptor::ReduceRealSum(foo, 6);

    force_p[0] = foo[0];
    force_p[1] = foo[1];
    force_p[2] = foo[2];

    force_s[0] = foo[3];
    force_s[1] = foo[4];
    force_s[2] = foo[5];

    // Divide by the mass of the stars to obtain the acceleration, and then get the new rotation frequency.

    Real fp = std::sqrt(std::pow(force_p[0], 2) + std::pow(force_p[1], 2) + std::pow(force_p[2], 2));
    Real fs = std::sqrt(std::pow(force_s[0], 2) + std::pow(force_s[1], 2) + std::pow(force_s[2], 2));

    Real ap = std::sqrt(std::pow(com_p[0], 2) + std::pow(com_p[1], 2) + std::pow(com_p[2], 2));
    Real as = std::sqrt(std::pow(com_s[0], 2) + std::pow(com_s[1], 2) + std::pow(com_s[2], 2));

    Real omega = 0.5 * ( std::sqrt((fp / mass_p) / ap) + std::sqrt((fs / mass_s) / as) );

    Real period = 2.0 * M_PI / omega;

    if (amrex::ParallelDescriptor::IOProcessor()) {
          std::cout << "\n";
          std::cout << "  Updating the rotational period from " << rotational_period << " s to " << period << " s." << "\n";
	  std::cout << "\n";
    }

    rotational_period = period;
    set_period(&period);

    // Update the force applied to the state using the new rotation vector.
    // To do this we'll need to copy the old state into Sborder since that
    // is what the source term constructors rely on for the old time.

    for (int lev = 0; lev <= finest_level; ++lev) {

	MultiFab& S_old = getLevel(lev).get_old_data(State_Type);
	MultiFab& S_new = getLevel(lev).get_new_data(State_Type);

	MultiFab& Sb = getLevel(lev).Sborder;

	// Note that we'll do this exactly the same way as it occurs in the advance
	// since the construction of Sborder includes more than just a FillPatch.
	// The only difference is that we don't need NUM_GROW ghost cells.

	Sb.define(getLevel(lev).grids, getLevel(lev).dmap, NUM_STATE, ng);

	getLevel(lev).expand_state(Sb, old_time, ng);

	MultiFab::Saxpy(S_new, -dt, *(getLevel(lev).old_sources[rot_src]), 0, 0, NUM_STATE, ng);
	MultiFab::Saxpy(S_new, -dt, *(getLevel(lev).new_sources[rot_src]), 0, 0, NUM_STATE, ng);

	getLevel(lev).construct_old_source(rot_src, old_time, dt);
	getLevel(lev).construct_new_source(rot_src, new_time, dt);

	MultiFab::Saxpy(S_new, dt, *(getLevel(lev).old_sources[rot_src]), 0, 0, NUM_STATE, ng);
	MultiFab::Saxpy(S_new, dt, *(getLevel(lev).new_sources[rot_src]), 0, 0, NUM_STATE, ng);

	Sb.clear();

    }

    // Check to see whether the relaxation should be turned off.
    // Note that at present the following check is only done on the
    // coarse grid but if we wanted more accuracy we could do a loop
    // over levels as above.

    Real L1[3] = { -1.0e200 };
    Real L2[3] = { -1.0e200 };
    Real L3[3] = { -1.0e200 };

    // First, calculate the location of the L1 Lagrange point.

    get_lagrange_points(mass_p, mass_s, com_p, com_s, L1, L2, L3);

    // Then, figure out the effective potential corresponding to that
    // Lagrange point.

    Real potential = 0.0;

    auto mfphieff = derive("phiEff", time, 0);

#ifdef _OPENMP
#pragma omp parallel reduction(+:potential)
#endif
    for (MFIter mfi(*mfphieff,true); mfi.isValid(); ++mfi) {

	const Box& box = mfi.tilebox();

	const int* lo  = box.loVect();
	const int* hi  = box.hiVect();

	get_critical_roche_potential(BL_TO_FORTRAN_3D((*mfphieff)[mfi]),
				     lo, hi, L1, &potential);

    }

    amrex::ParallelDescriptor::ReduceRealSum(potential);

    // Now cycle through the grids and determine if any zones
    // have crossed the density threshold outside the critical surface.

    MultiFab& S_new = get_new_data(State_Type);

    int is_done = 0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:is_done)
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi) {

	const Box& box  = mfi.tilebox();

	const int* lo   = box.loVect();
	const int* hi   = box.hiVect();

	check_relaxation(BL_TO_FORTRAN_3D(S_new[mfi]),
			 BL_TO_FORTRAN_3D((*mfphieff)[mfi]),
			 ARLIM_3D(lo),ARLIM_3D(hi),
			 &potential,&is_done);

    }

    amrex::ParallelDescriptor::ReduceIntSum(is_done);

    if (is_done > 0) {
	relaxation_is_done = 1;
	set_relaxation_status(&relaxation_is_done);
    }

    if (relaxation_is_done) {

	turn_off_relaxation(&time);

    }

}
