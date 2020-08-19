#include <cmath>

#include <iomanip>

#include <vector>

#include <Castro.H>
#include <Castro_F.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>

#include <Problem.H>
#include <Problem_F.H>

#include <Gravity.H>
#include <Gravity_F.H>

#include <wdmerger_data.H>

using amrex::Real;

void
Castro::sum_integrated_quantities ()
{
    using namespace wdmerger;

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

    // Rotation frequency.

    Real omega[3] = { 0.0 };

    get_omega_vec(omega);

    // Mass transfer rate

    Real mdot = 0.5 * (std::abs(mdot_p) + std::abs(mdot_s));

    // Center of mass of the system.

    Real com[3]       = { 0.0 };
    Real com_vel[3]   = { 0.0 };

    // Distance between the WDs.

    Real wd_dist[3] = { 0.0 };
    Real wd_dist_init[3] = { 0.0 };

    Real separation = 0.0;
    Real angle = 0.0;

    // Stellar centers of mass and velocities.

    Real com_p_mag = 0.0;
    Real com_s_mag = 0.0;

    Real vel_p_mag = 0.0;
    Real vel_s_mag = 0.0;

    Real vel_p_rad = 0.0;
    Real vel_s_rad = 0.0;

    Real vel_p_phi = 0.0;
    Real vel_s_phi = 0.0;

    // Gravitational wave amplitudes.
    
    Real h_plus_1  = 0.0;
    Real h_cross_1 = 0.0;

    Real h_plus_2  = 0.0;
    Real h_cross_2 = 0.0;

    Real h_plus_3  = 0.0;
    Real h_cross_3 = 0.0;

    // Species names and total masses on the domain.

    const Real M_solar = 1.9884e33;

    std::vector<Real> species_mass(NumSpec);
    std::vector<std::string> species_names(NumSpec);

    std::string name1;
    std::string name2;

    int dataprecision = 16; // Number of digits after the decimal point, for float data

    int datwidth      = 25; // Floating point data in scientific notation
    int fixwidth      = 25; // Floating point data not in scientific notation
    int intwidth      = 12; // Integer data

    int axis_1;
    int axis_2;
    int axis_3;

    // Determine various coordinate axes
    get_axes(&axis_1, &axis_2, &axis_3);

    wd_dist_init[axis_1 - 1] = 1.0;

    // Determine the names of the species in the simulation.

    for (int i = 0; i < NumSpec; i++) {
      species_names[i] = desc_lst[State_Type].name(UFS+i);
      species_names[i] = species_names[i].substr(4,std::string::npos);
      species_mass[i]  = 0.0;
    }

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

#ifdef GRAVITY
#if (BL_SPACEDIM > 1)
      // Gravitational wave signal. This is designed to add to these quantities so we can send them directly.
      ca_lev.gwstrain(time, h_plus_1, h_cross_1, h_plus_2, h_cross_2, h_plus_3, h_cross_3, local_flag);
#endif
#endif

      // Integrated mass of all species on the domain.
      for (int i = 0; i < NumSpec; i++)
	species_mass[i] += ca_lev.volWgtSum("rho_" + species_names[i], time, local_flag) / M_solar;

    }

    // Return to the original level.

    ca_set_amr_info(level, -1, -1, -1.0, -1.0);

    // Do the reductions.

    int nfoo_sum = 24 + NumSpec;

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
    foo_sum[18] = h_plus_1;
    foo_sum[19] = h_cross_1;
    foo_sum[20] = h_plus_2;
    foo_sum[21] = h_cross_2;
    foo_sum[22] = h_plus_3;
    foo_sum[23] = h_cross_3;

    for (int i = 0; i < NumSpec; i++) {
      foo_sum[i + 24] = species_mass[i];
    }

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
    h_plus_1   = foo_sum[18];
    h_cross_1  = foo_sum[19];
    h_plus_2   = foo_sum[20];
    h_cross_2  = foo_sum[21];
    h_plus_3   = foo_sum[22];
    h_cross_3  = foo_sum[23];

    for (int i = 0; i < NumSpec; i++) {
      species_mass[i] = foo_sum[i + 24];
    }

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

    com_p_mag += std::pow( std::pow(com_p[0],2) + std::pow(com_p[1],2) + std::pow(com_p[2],2), 0.5 );
    com_s_mag += std::pow( std::pow(com_s[0],2) + std::pow(com_s[1],2) + std::pow(com_s[2],2), 0.5 );
    vel_p_mag += std::pow( std::pow(vel_p[0],2) + std::pow(vel_p[1],2) + std::pow(vel_p[2],2), 0.5 );
    vel_s_mag += std::pow( std::pow(vel_s[0],2) + std::pow(vel_s[1],2) + std::pow(vel_s[2],2), 0.5 );

#if (BL_SPACEDIM == 3)
    if (mass_p > 0.0) {
      vel_p_rad = (com_p[axis_1 - 1] / com_p_mag) * vel_p[axis_1 - 1] + (com_p[axis_2 - 1] / com_p_mag) * vel_p[axis_2 - 1];
      vel_p_phi = (com_p[axis_1 - 1] / com_p_mag) * vel_p[axis_2 - 1] - (com_p[axis_2 - 1] / com_p_mag) * vel_p[axis_1 - 1];
    }

    if (mass_s > 0.0) {
      vel_s_rad = (com_s[axis_1 - 1] / com_s_mag) * vel_s[axis_1 - 1] + (com_s[axis_2 - 1] / com_s_mag) * vel_s[axis_2 - 1];
      vel_s_phi = (com_s[axis_1 - 1] / com_s_mag) * vel_s[axis_2 - 1] - (com_s[axis_2 - 1] / com_s_mag) * vel_s[axis_1 - 1];
    }
#else
    if (mass_p > 0.0) {
      vel_p_rad = vel_p[axis_1 - 1];
      vel_p_phi = vel_p[axis_3 - 1];
    }

    if (mass_s > 0.0) {
      vel_s_rad = vel_s[axis_1 - 1];
      vel_s_phi = vel_s[axis_3 - 1];
    }
#endif

    if (mass_p > 0.0 && mass_s > 0.0) {

      // Calculate the distance between the primary and secondary.

      for ( int i = 0; i < 3; i++ )
	wd_dist[i] = com_s[i] - com_p[i];

      separation = norm(wd_dist);

      // Calculate the angle between the initial stellar axis and
      // the line currently joining the two stars. Note that this
      // neglects any motion in the plane perpendicular to the initial orbit.

      angle = atan2( wd_dist[axis_2 - 1] - wd_dist_init[axis_2 - 1],
                     wd_dist[axis_1 - 1] - wd_dist_init[axis_1 - 1] ) * 180.0 / M_PI;

      // Now let's transform from [-180, 180] to [0, 360].

      if (angle < 0.0) angle += 360.0;

    }

    // Calculate wall time for the step.

    Real wall_time = 0.0;

    if (time > 0.0)
        wall_time = amrex::ParallelDescriptor::second() - wall_time_start;

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
	     header << std::setw(datwidth) << "             h_+ (axis 1)"; ++n;
	     header << std::setw(datwidth) << "             h_x (axis 1)"; ++n;
	     header << std::setw(datwidth) << "             h_+ (axis 2)"; ++n;
	     header << std::setw(datwidth) << "             h_x (axis 2)"; ++n;
	     header << std::setw(datwidth) << "             h_+ (axis 3)"; ++n;
	     header << std::setw(datwidth) << "             h_x (axis 3)"; ++n;

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
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << h_plus_1;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << h_cross_1;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << h_plus_2;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << h_cross_2;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << h_plus_3;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << h_cross_3;

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

      // Species

      if (parent->NumDataLogs() > 2) {

	 std::ostream& log = parent->DataLog(2);

	 if ( log.good() ) {

	   if (time == 0.0) {

	     // Output the git commit hashes used to build the executable.

	     writeGitHashes(log);

             int n = 0;

             std::ostringstream header;

	     header << std::setw(intwidth) << "#   TIMESTEP";           ++n;
	     header << std::setw(fixwidth) << "                  TIME"; ++n;

	     // We need to be careful here since the species names have differing numbers of characters

	     for (int i = 0; i < NumSpec; i++) {
	       std::string outString  = "";
	       std::string massString = "Mass ";
	       std::string specString = species_names[i];
               while (outString.length() + specString.length() + massString.length() < datwidth) outString += " ";
	       outString += massString;
	       outString += specString;
	       header << std::setw(datwidth) << outString; ++n;
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

	   for (int i = 0; i < NumSpec; i++)
	     log << std::setw(datwidth) << std::setprecision(dataprecision) << species_mass[i];

	   log << std::endl;

	 }
      }

      // Information about the AMR driver.

      if (parent->NumDataLogs() > 3) {

	 std::ostream& log = parent->DataLog(3);

	 if ( log.good() ) {

	   if (time == 0.0) {

	     // Output the git commit hashes used to build the executable.

	     writeGitHashes(log);

             int n = 0;

             std::ostringstream header;

	     header << std::setw(intwidth) << "#   TIMESTEP";              ++n;
	     header << std::setw(fixwidth) << "                     TIME"; ++n;
	     header << std::setw(fixwidth) << "                       DT"; ++n;
	     header << std::setw(intwidth) << "  FINEST LEV";              ++n;
             header << std::setw(fixwidth) << " COARSE TIMESTEP WALLTIME"; ++n;

	     header << std::endl;

             log << std::setw(intwidth) << "#   COLUMN 1";
             log << std::setw(fixwidth) << "                        2";

             for (int i = 3; i < 4; ++i)
                 log << std::setw(datwidth) << i;

             log << std::setw(intwidth) << 4; // Handle the finest lev column

             for (int i = 5; i <= n; ++i)
                 log << std::setw(datwidth) << i;

             log << std::endl;

             log << header.str();

	   }

	   log << std::fixed;

	   log << std::setw(intwidth)                                     << timestep;
	   log << std::setw(fixwidth) << std::setprecision(dataprecision) << time;
	   log << std::setw(fixwidth) << std::setprecision(dataprecision) << dt;
	   log << std::setw(intwidth)                                     << parent->finestLevel();
           log << std::setw(datwidth) << std::setprecision(dataprecision) << wall_time;

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

	   log << std::setw(datwidth) << std::setprecision(dataprecision) << mass_p;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << mdot_p;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << com_p_mag;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << com_p[0];
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << com_p[1];
#if (BL_SPACEDIM == 3)
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << com_p[2];
#endif
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_p_mag;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_p_rad;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_p_phi;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_p[0];
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_p[1];
#if (BL_SPACEDIM == 3)
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_p[2];
#endif
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << t_ff_p;
	   for (int i = 0; i <= 6; ++i)
	       log << std::setw(datwidth) << std::setprecision(dataprecision) << rad_p[i];

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

	   log << std::setw(datwidth) << std::setprecision(dataprecision) << mass_s;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << mdot_s;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << com_s_mag;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << com_s[0];
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << com_s[1];
#if (BL_SPACEDIM == 3)
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << com_s[2];
#endif
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_s_mag;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_s_rad;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_s_phi;
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_s[0];
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_s[1];
#if (BL_SPACEDIM == 3)
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << vel_s[2];
#endif
	   log << std::setw(datwidth) << std::setprecision(dataprecision) << t_ff_s;
	   for (int i = 0; i <= 6; ++i)
	       log << std::setw(datwidth) << std::setprecision(dataprecision) << rad_s[i];

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

