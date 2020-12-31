#include <cmath>

#include <iomanip>

#include <vector>

#include <Castro.H>
#include <Castro_F.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>

#include <Problem.H>

#include <Gravity.H>
#include <Gravity_F.H>

#include <wdmerger_data.H>

using amrex::Real;

void
Castro::sum_integrated_quantities ()
{
    using namespace wdmerger;
    using namespace problem;

    if (level > 0) return;

    bool local_flag = true;

    int finest_level  = parent->finestLevel();
    Real time         = state[State_Type].curTime();
    Real dt           = parent->dtLevel(0);

    if (time == 0.0) dt = 0.0; // dtLevel returns the next timestep for t = 0, so overwrite

    int timestep = parent->levelSteps(0);

    Real mass                 = 0.0;
    Real momentum[3]          = { 0.0 };
    Real angular_momentum[3]  = { 0.0 };
    Real hybrid_momentum[3]   = { 0.0 };
    Real rho_E                = 0.0;
    Real rho_e                = 0.0;
    Real rho_K                = 0.0;
    Real rho_phi              = 0.0;
    Real rho_phirot           = 0.0;

    // Total energy on the grid, including decomposition
    // into the various components.

    Real gravitational_energy = 0.0;
    Real kinetic_energy       = 0.0;
    Real gas_energy           = 0.0;
    Real rotational_energy    = 0.0;
    Real internal_energy      = 0.0;
    Real total_energy         = 0.0;
    Real total_E_grid         = 0.0;

    // Mass transfer rate

    Real mdot = 0.5 * (std::abs(mdot_P) + std::abs(mdot_S));

    // Center of mass of the system.

    Real com[3]       = { 0.0 };
    Real com_vel[3]   = { 0.0 };

    // Distance between the WDs.

    Real wd_dist[3] = { 0.0 };
    Real wd_dist_init[3] = { 0.0 };

    Real separation = 0.0;
    Real angle = 0.0;

    // Stellar centers of mass and velocities.

    Real com_P_mag = 0.0;
    Real com_S_mag = 0.0;

    Real vel_P_mag = 0.0;
    Real vel_S_mag = 0.0;

    Real vel_P_rad = 0.0;
    Real vel_S_rad = 0.0;

    Real vel_P_phi = 0.0;
    Real vel_S_phi = 0.0;

    std::string name1;
    std::string name2;

    int dataprecision = 16; // Number of digits after the decimal point, for float data

    int datwidth      = 25; // Floating point data in scientific notation
    int fixwidth      = 25; // Floating point data not in scientific notation
    int intwidth      = 12; // Integer data

    wd_dist_init[problem::axis_1 - 1] = 1.0;

    for (int lev = 0; lev <= finest_level; lev++)
    {

      // Update the local level we're on.

      ca_set_amr_info(lev, -1, -1, -1.0, -1.0);

      // Get the current level from Castro

      Castro& ca_lev = getLevel(lev);

      for ( int i = 0; i < 3; i++ ) {
        com[i] += ca_lev.locWgtSum("density", time, i, local_flag);
      }

      // Calculate total mass, momentum, angular momentum, and energy of system.

      mass += ca_lev.volWgtSum("density", time, local_flag);

      momentum[0] += ca_lev.volWgtSum("inertial_momentum_x", time, local_flag);
      momentum[1] += ca_lev.volWgtSum("inertial_momentum_y", time, local_flag);
      momentum[2] += ca_lev.volWgtSum("inertial_momentum_z", time, local_flag);

      angular_momentum[0] += ca_lev.volWgtSum("inertial_angular_momentum_x", time, local_flag);
      angular_momentum[1] += ca_lev.volWgtSum("inertial_angular_momentum_y", time, local_flag);
      angular_momentum[2] += ca_lev.volWgtSum("inertial_angular_momentum_z", time, local_flag);

#ifdef HYBRID_MOMENTUM
      hybrid_momentum[0] += ca_lev.volWgtSum("rmom", time, local_flag);
      hybrid_momentum[1] += ca_lev.volWgtSum("lmom", time, local_flag);
      hybrid_momentum[2] += ca_lev.volWgtSum("pmom", time, local_flag);
#endif

      rho_E += ca_lev.volWgtSum("rho_E", time, local_flag);
      rho_K += ca_lev.volWgtSum("kineng",time, local_flag);
      rho_e += ca_lev.volWgtSum("rho_e", time, local_flag);

#ifdef GRAVITY
      if (do_grav)
        rho_phi += ca_lev.volProductSum("density", "phiGrav", time, local_flag);
#endif

#ifdef ROTATION
      if (do_rotation)
	rho_phirot += ca_lev.volProductSum("density", "phiRot", time, local_flag);
#endif

    }

    // Return to the original level.

    ca_set_amr_info(level, -1, -1, -1.0, -1.0);

    // Do the reductions.

    int nfoo_sum = 18;

    amrex::Vector<Real> foo_sum(nfoo_sum);

    foo_sum[0] = mass;

    for (int i = 0; i < 3; i++) {
      foo_sum[i+1]  = com[i];
      foo_sum[i+4]  = momentum[i];
      foo_sum[i+7]  = angular_momentum[i];
      foo_sum[i+10] = hybrid_momentum[i];
    }

    foo_sum[13] = rho_E;
    foo_sum[14] = rho_K;
    foo_sum[15] = rho_e;
    foo_sum[16] = rho_phi;
    foo_sum[17] = rho_phirot;

    amrex::ParallelDescriptor::ReduceRealSum(foo_sum.dataPtr(), nfoo_sum);

    mass = foo_sum[0];

    for (int i = 0; i < 3; i++) {
      com[i]              = foo_sum[i+1];
      momentum[i]         = foo_sum[i+4];
      angular_momentum[i] = foo_sum[i+7];
      hybrid_momentum[i]  = foo_sum[i+10];
    }

    rho_E      = foo_sum[13];
    rho_K      = foo_sum[14];
    rho_e      = foo_sum[15];
    rho_phi    = foo_sum[16];
    rho_phirot = foo_sum[17];

    // Complete calculations for energy and momenta

    gravitational_energy = rho_phi;
    if (gravity->get_gravity_type() == "PoissonGrav")
      gravitational_energy *= 0.5; // avoids double counting
    internal_energy = rho_e;
    kinetic_energy = rho_K;
    gas_energy = rho_E;
    rotational_energy = rho_phirot;
    total_E_grid = gravitational_energy + rho_E;
    total_energy = total_E_grid + rotational_energy;

    // Complete calculations for center of mass quantities

    for ( int i = 0; i < 3; i++ ) {

      com[i]       = com[i] / mass;
      com_vel[i]   = momentum[i] / mass;

    }

    com_P_mag += std::pow( std::pow(com_P[0],2) + std::pow(com_P[1],2) + std::pow(com_P[2],2), 0.5 );
    com_S_mag += std::pow( std::pow(com_S[0],2) + std::pow(com_S[1],2) + std::pow(com_S[2],2), 0.5 );
    vel_P_mag += std::pow( std::pow(vel_P[0],2) + std::pow(vel_P[1],2) + std::pow(vel_P[2],2), 0.5 );
    vel_S_mag += std::pow( std::pow(vel_S[0],2) + std::pow(vel_S[1],2) + std::pow(vel_S[2],2), 0.5 );

#if (BL_SPACEDIM == 3)
    if (mass_P > 0.0) {
      vel_P_rad = (com_P[problem::axis_1 - 1] / com_P_mag) * vel_P[problem::axis_1 - 1] +
                  (com_P[problem::axis_2 - 1] / com_P_mag) * vel_P[problem::axis_2 - 1];
      vel_P_phi = (com_P[problem::axis_1 - 1] / com_P_mag) * vel_P[problem::axis_2 - 1] -
                  (com_P[problem::axis_2 - 1] / com_P_mag) * vel_P[problem::axis_1 - 1];
    }

    if (mass_S > 0.0) {
      vel_S_rad = (com_S[problem::axis_1 - 1] / com_S_mag) * vel_S[problem::axis_1 - 1] +
                  (com_S[problem::axis_2 - 1] / com_S_mag) * vel_S[problem::axis_2 - 1];
      vel_S_phi = (com_S[problem::axis_1 - 1] / com_S_mag) * vel_S[problem::axis_2 - 1] -
                  (com_S[problem::axis_2 - 1] / com_S_mag) * vel_S[problem::axis_1 - 1];
    }
#else
    if (mass_P > 0.0) {
      vel_P_rad = vel_P[problem::axis_1 - 1];
      vel_P_phi = vel_P[problem::axis_3 - 1];
    }

    if (mass_S > 0.0) {
      vel_S_rad = vel_S[problem::axis_1 - 1];
      vel_S_phi = vel_S[problem::axis_3 - 1];
    }
#endif

    if (mass_P > 0.0 && mass_S > 0.0) {

      // Calculate the distance between the primary and secondary.

      for ( int i = 0; i < 3; i++ )
	wd_dist[i] = com_S[i] - com_P[i];

      separation = norm(wd_dist);

      // Calculate the angle between the initial stellar axis and
      // the line currently joining the two stars. Note that this
      // neglects any motion in the plane perpendicular to the initial orbit.

      angle = std::atan2(wd_dist[problem::axis_2 - 1] - wd_dist_init[problem::axis_2 - 1],
                         wd_dist[problem::axis_1 - 1] - wd_dist_init[problem::axis_1 - 1]) * 180.0 / M_PI;

      // Now let's transform from [-180, 180] to [0, 360].

      if (angle < 0.0) angle += 360.0;

    }

    // Write data out to the log.

    if ( amrex::ParallelDescriptor::IOProcessor() )
    {

      // The data logs are only defined on the IO processor
      // for parallel runs, so the stream should only be opened inside.

      if (parent->NumDataLogs() > 0) {

	 std::ostream& log = parent->DataLog(0);

	 if ( log.good() ) {

	   // Write header row

	   if (time == 0.0) {

	     // Output the git commit hashes used to build the executable.

	     writeGitHashes(log);

             int n = 0;

             std::ostringstream header;

	     header << std::setw(intwidth) << "#   TIMESTEP";              ++n;
	     header << std::setw(fixwidth) << "                     TIME"; ++n;
	     header << std::setw(datwidth) << "             TOTAL ENERGY"; ++n;
	     header << std::setw(datwidth) << "             TOTAL E GRID"; ++n;
	     header << std::setw(datwidth) << "               GAS ENERGY"; ++n;
	     header << std::setw(datwidth) << "              KIN. ENERGY"; ++n;
	     header << std::setw(datwidth) << "              ROT. ENERGY"; ++n;
	     header << std::setw(datwidth) << "             GRAV. ENERGY"; ++n;
	     header << std::setw(datwidth) << "              INT. ENERGY"; ++n;
	     header << std::setw(datwidth) << "                     MASS"; ++n;
#if (BL_SPACEDIM == 3)
	     header << std::setw(datwidth) << "                     XMOM"; ++n;
	     header << std::setw(datwidth) << "                     YMOM"; ++n;
	     header << std::setw(datwidth) << "                     ZMOM"; ++n;
#ifdef HYBRID_MOMENTUM
	     header << std::setw(datwidth) << "              HYB. MOM. R"; ++n;
	     header << std::setw(datwidth) << "              HYB. MOM. L"; ++n;
	     header << std::setw(datwidth) << "              HYB. MOM. P"; ++n;
#endif
	     header << std::setw(datwidth) << "              ANG. MOM. X"; ++n;
	     header << std::setw(datwidth) << "              ANG. MOM. Y"; ++n;
	     header << std::setw(datwidth) << "              ANG. MOM. Z"; ++n;
#else
	     header << std::setw(datwidth) << "                     RMOM"; ++n;
	     header << std::setw(datwidth) << "                     ZMOM"; ++n;
	     header << std::setw(datwidth) << "              ANG. MOM. R"; ++n;
	     header << std::setw(datwidth) << "              ANG. MOM. Z"; ++n;
#endif
#if (BL_SPACEDIM == 3)
	     header << std::setw(datwidth) << "                    X COM"; ++n;
	     header << std::setw(datwidth) << "                    Y COM"; ++n;
	     header << std::setw(datwidth) << "                    Z COM"; ++n;
	     header << std::setw(datwidth) << "                X COM VEL"; ++n;
	     header << std::setw(datwidth) << "                Y COM VEL"; ++n;
	     header << std::setw(datwidth) << "                Z COM VEL"; ++n;
#else
	     header << std::setw(datwidth) << "                    R COM"; ++n;
	     header << std::setw(datwidth) << "                    Z COM"; ++n;
	     header << std::setw(datwidth) << "                R COM VEL"; ++n;
	     header << std::setw(datwidth) << "                Z COM VEL"; ++n;
#endif

	     header << std::endl;

             log << std::setw(intwidth) << "#   COLUMN 1";
             log << std::setw(fixwidth) << "                        2";

             for (int i = 3; i <= n; ++i)
                 log << std::setw(datwidth) << i;

             log << std::endl;

             log << header.str();

	   }

	   // Write data for the present time

	   log << std::fixed;

	   log << std::setw(intwidth)                                     << timestep;
	   log << std::setw(fixwidth) << std::setprecision(dataprecision) << time;

	   log << std::scientific;

	   log << std::setw(datwidth) << std::setprecision(dataprecision) << total_energy;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << total_E_grid;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << gas_energy;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << kinetic_energy;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << rotational_energy;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << gravitational_energy;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << internal_energy;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << mass;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << momentum[0];
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << momentum[1];
#if (BL_SPACEDIM == 3)
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << momentum[2];
#endif
#ifdef HYBRID_MOMENTUM
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << hybrid_momentum[0];
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << hybrid_momentum[1];
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << hybrid_momentum[2];
#endif
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << angular_momentum[0];
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << angular_momentum[1];
#if (BL_SPACEDIM == 3)
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << angular_momentum[2];
#endif
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << com[0];
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << com[1];
#if (BL_SPACEDIM == 3)
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << com[2];
#endif
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << com_vel[0];
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << com_vel[1];
#if (BL_SPACEDIM == 3)
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << com_vel[2];
#endif

	   log << std::endl;
	 }
      }

      if (parent->NumDataLogs() > 1) {

	 std::ostream& log = parent->DataLog(1);

	 if ( log.good() ) {

	   if (time == 0.0) {

	     // Output the git commit hashes used to build the executable.

	     writeGitHashes(log);

             int n = 0;

             std::ostringstream header;

	     header << std::setw(intwidth) << "#   TIMESTEP";              ++n;
	     header << std::setw(fixwidth) << "                     TIME"; ++n;

	     header << std::setw(datwidth) << "              WD DISTANCE"; ++n;
	     header << std::setw(fixwidth) << "                 WD ANGLE"; ++n;
	     header << std::setw(datwidth) << "                     MDOT"; ++n;

	     header << std::endl;

             log << std::setw(intwidth) << "#   COLUMN 1";
             log << std::setw(fixwidth) << "                        2";

             for (int i = 3; i <= n; ++i)
                 log << std::setw(datwidth) << i;

             log << std::endl;

             log << header.str();

	   }

	   log << std::fixed;

	   log << std::setw(intwidth)                                     << timestep;
	   log << std::setw(fixwidth) << std::setprecision(dataprecision) << time;

	   log << std::scientific;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << separation;

	   log << std::fixed;
	   log << std::setw(fixwidth) << std::setprecision(dataprecision) << angle;

	   log << std::scientific;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << mdot;

	   log << std::endl;

	 }
      }

      // Primary star

      if (parent->NumDataLogs() > 4) {

	std::ostream& log = parent->DataLog(4);

	 if ( log.good() ) {

	   if (time == 0.0) {

	     // Output the git commit hashes used to build the executable.

	     writeGitHashes(log);

             int n = 0;

             std::ostringstream header;

	     header << std::setw(intwidth) << "#   TIMESTEP";              ++n;
	     header << std::setw(fixwidth) << "                     TIME"; ++n;
	     header << std::setw(datwidth) << "             PRIMARY MASS"; ++n;
	     header << std::setw(datwidth) << "             PRIMARY MDOT"; ++n;
	     header << std::setw(datwidth) << "          PRIMARY MAG COM"; ++n;
#if (BL_SPACEDIM == 3)
	     header << std::setw(datwidth) << "            PRIMARY X COM"; ++n;
	     header << std::setw(datwidth) << "            PRIMARY Y COM"; ++n;
	     header << std::setw(datwidth) << "            PRIMARY Z COM"; ++n;
#else
	     header << std::setw(datwidth) << "            PRIMARY R COM"; ++n;
	     header << std::setw(datwidth) << "            PRIMARY Z COM"; ++n;
#endif
	     header << std::setw(datwidth) << "          PRIMARY MAG VEL"; ++n;
	     header << std::setw(datwidth) << "          PRIMARY RAD VEL"; ++n;
	     header << std::setw(datwidth) << "          PRIMARY ANG VEL"; ++n;
#if (BL_SPACEDIM == 3)
	     header << std::setw(datwidth) << "            PRIMARY X VEL"; ++n;
	     header << std::setw(datwidth) << "            PRIMARY Y VEL"; ++n;
	     header << std::setw(datwidth) << "            PRIMARY Z VEL"; ++n;
#else
	     header << std::setw(datwidth) << "            PRIMARY R VEL"; ++n;
	     header << std::setw(datwidth) << "            PRIMARY Z VEL"; ++n;
#endif
	     header << std::setw(datwidth) << "       PRIMARY T_FREEFALL"; ++n;
	     for (int i = 0; i <= 6; ++i) {
                 header << "       PRIMARY 1E" << i << " RADIUS";          ++n;
             }

	     header << std::endl;

             log << std::setw(intwidth) << "#   COLUMN 1";
             log << std::setw(fixwidth) << "                        2";

             for (int i = 3; i <= n; ++i)
                 log << std::setw(datwidth) << i;

             log << std::endl;

             log << header.str();

	   }

	   log << std::fixed;

	   log << std::setw(intwidth)                                     << timestep;
	   log << std::setw(fixwidth) << std::setprecision(dataprecision) << time;

	   log << std::scientific;

	   log << std::setw(datwidth) << std::setprecision(dataprecision) << mass_P;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << mdot_P;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << com_P_mag;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << com_P[0];
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << com_P[1];
#if (BL_SPACEDIM == 3)
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << com_P[2];
#endif
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_P_mag;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_P_rad;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_P_phi;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_P[0];
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_P[1];
#if (BL_SPACEDIM == 3)
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_P[2];
#endif
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << t_ff_P;
	   for (int i = 0; i <= 6; ++i)
	       log << std::setw(datwidth) << std::setprecision(dataprecision) << rad_P[i];

	   log << std::endl;

	 }

      }

      // Secondary star

      if (parent->NumDataLogs() > 5) {

	std::ostream& log = parent->DataLog(5);

	 if ( log.good() ) {

	   if (time == 0.0) {

	     // Output the git commit hashes used to build the executable.

	     writeGitHashes(log);

             int n = 0;

             std::ostringstream header;

	     header << std::setw(intwidth) << "#   TIMESTEP";              ++n;
	     header << std::setw(fixwidth) << "                     TIME"; ++n;
	     header << std::setw(datwidth) << "           SECONDARY MASS"; ++n;
	     header << std::setw(datwidth) << "           SECONDARY MDOT"; ++n;
	     header << std::setw(datwidth) << "        SECONDARY MAG COM"; ++n;
#if (BL_SPACEDIM == 3)
	     header << std::setw(datwidth) << "          SECONDARY X COM"; ++n;
	     header << std::setw(datwidth) << "          SECONDARY Y COM"; ++n;
	     header << std::setw(datwidth) << "          SECONDARY Z COM"; ++n;
#else
	     header << std::setw(datwidth) << "          SECONDARY R COM"; ++n;
	     header << std::setw(datwidth) << "          SECONDARY Z COM"; ++n;
#endif
	     header << std::setw(datwidth) << "        SECONDARY MAG VEL"; ++n;
	     header << std::setw(datwidth) << "        SECONDARY RAD VEL"; ++n;
	     header << std::setw(datwidth) << "        SECONDARY ANG VEL"; ++n;
#if (BL_SPACEDIM == 3)
	     header << std::setw(datwidth) << "          SECONDARY X VEL"; ++n;
	     header << std::setw(datwidth) << "          SECONDARY Y VEL"; ++n;
	     header << std::setw(datwidth) << "          SECONDARY Z VEL"; ++n;
#else
	     header << std::setw(datwidth) << "          SECONDARY R VEL"; ++n;
	     header << std::setw(datwidth) << "          SECONDARY Z VEL"; ++n;
#endif
	     header << std::setw(datwidth) << "     SECONDARY T_FREEFALL"; ++n;
	     for (int i = 0; i <= 6; ++i) {
                 header << "     SECONDARY 1E" << i << " RADIUS";          ++n;
             }

	     header << std::endl;

             log << std::setw(intwidth) << "#   COLUMN 1";
             log << std::setw(fixwidth) << "                        2";

             for (int i = 3; i <= n; ++i)
                 log << std::setw(datwidth) << i;

             log << std::endl;

             log << header.str();

	   }

	   log << std::fixed;

	   log << std::setw(intwidth)                                     << timestep;
	   log << std::setw(fixwidth) << std::setprecision(dataprecision) << time;

	   log << std::scientific;

	   log << std::setw(datwidth) << std::setprecision(dataprecision) << mass_S;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << mdot_S;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << com_S_mag;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << com_S[0];
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << com_S[1];
#if (BL_SPACEDIM == 3)
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << com_S[2];
#endif
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_S_mag;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_S_rad;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_S_phi;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_S[0];
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_S[1];
#if (BL_SPACEDIM == 3)
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_S[2];
#endif
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << t_ff_S;
	   for (int i = 0; i <= 6; ++i)
	       log << std::setw(datwidth) << std::setprecision(dataprecision) << rad_S[i];

	   log << std::endl;

	 }

      }

      // Extrema over time of various quantities

      if (parent->NumDataLogs() > 6) {

	 std::ostream& log = parent->DataLog(6);

	 if ( log.good() ) {

	   if (time == 0.0) {

	     // Output the git commit hashes used to build the executable.

	     writeGitHashes(log);

             int n = 0;

             std::ostringstream header;

	     header << std::setw(intwidth) << "#   TIMESTEP";              ++n;
	     header << std::setw(fixwidth) << "                     TIME"; ++n;
	     header << std::setw(datwidth) << "               MAX T CURR"; ++n;
	     header << std::setw(datwidth) << "             MAX RHO CURR"; ++n;
	     header << std::setw(datwidth) << "           MAX TS_TE CURR"; ++n;
	     header << std::setw(datwidth) << "           MAX T ALL TIME"; ++n;
	     header << std::setw(datwidth) << "         MAX RHO ALL TIME"; ++n;
	     header << std::setw(datwidth) << "       MAX TS_TE ALL TIME"; ++n;

	     header << std::endl;

             log << std::setw(intwidth) << "#   COLUMN 1";
             log << std::setw(fixwidth) << "                        2";

             for (int i = 3; i <= n; ++i)
                 log << std::setw(datwidth) << i;

             log << std::endl;

             log << header.str();

	   }

	   log << std::fixed;

	   log << std::setw(intwidth)                                     << timestep;
	   log << std::setw(fixwidth) << std::setprecision(dataprecision) << time;

	   log << std::scientific;

	   log << std::setw(datwidth) << std::setprecision(dataprecision) << T_curr_max;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << rho_curr_max;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << ts_te_curr_max;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << T_global_max;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << rho_global_max;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << ts_te_global_max;

	   log << std::endl;

	 }

      }

      // Rotation period over time

      if (parent->NumDataLogs() > 7) {

	 std::ostream& log = parent->DataLog(7);

	 if ( log.good() ) {

	   if (time == 0.0) {

	     // Output the git commit hashes used to build the executable.

	     writeGitHashes(log);

             int n = 0;

             std::ostringstream header;

	     header << std::setw(intwidth) << "#   TIMESTEP";              ++n;
	     header << std::setw(fixwidth) << "                     TIME"; ++n;
	     header << std::setw(datwidth) << "          ROTATION PERIOD"; ++n;
             header << std::setw(datwidth) << "       ROTATION FREQUENCY"; ++n;

	     header << std::endl;

             log << std::setw(intwidth) << "#   COLUMN 1";
             log << std::setw(fixwidth) << "                        2";

             for (int i = 3; i <= n; ++i)
                 log << std::setw(datwidth) << i;

             log << std::endl;

             log << header.str();

	   }

	   log << std::fixed;

	   log << std::setw(intwidth)                                     << timestep;
	   log << std::setw(fixwidth) << std::setprecision(dataprecision) << time;

	   log << std::scientific;

	   log << std::setw(datwidth) << std::setprecision(dataprecision) << rotational_period;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << (2.0 * M_PI / rotational_period);

	   log << std::endl;

	 }

      }

    }
}

