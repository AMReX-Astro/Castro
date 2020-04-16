
// Preprocessor directive for allowing us to do a post-initialization update.

#ifndef DO_PROBLEM_POST_INIT
#define DO_PROBLEM_POST_INIT
#endif

// Same for post-restart.

#ifndef DO_PROBLEM_POST_RESTART
#define DO_PROBLEM_POST_RESTART
#endif

// Preprocessor directive for allowing us to do a post-timestep update.

#ifndef DO_PROBLEM_POST_TIMESTEP
#define DO_PROBLEM_POST_TIMESTEP
#endif

// Routine for determining white dwarf masses, centers of mass and velocities.

void wd_update(amrex::Real time, amrex::Real dt);

// Volume within white dwarfs at a given density cutoff.

void volInBoundary (amrex::Real time, amrex::Real& vol_p, amrex::Real& vol_s, amrex::Real rho_cutoff, bool local = false);

// Calculate gravitational wave signal.

void gwstrain (amrex::Real time,
	       amrex::Real& h_plus_1, amrex::Real& h_cross_1,
	       amrex::Real& h_plus_2, amrex::Real& h_cross_2,
	       amrex::Real& h_plus_3, amrex::Real& h_cross_3,
	       bool local = false);

// Computes standard dot product of two three-vectors.

amrex::Real dot_product(const amrex::Real a[], const amrex::Real b[]);

// Computes standard cross-product of two three-vectors.

void cross_product(const amrex::Real a[], const amrex::Real b[], amrex::Real c[]);

// Computes norm of a three vector.

amrex::Real norm(const amrex::Real a[]);

// Problem post-initialization routine.

void problem_post_init();

// Problem post-restart routine.

void problem_post_restart();

// End of timestep analytics on the white dwarf properties.

void problem_post_timestep();

// Write out the git hashes for the various parts of the code.

void writeGitHashes(std::ostream& log);

// Check whether we should stop the run.

void check_to_stop(amrex::Real time, bool dump = true);

// Calculate extrema on grid.

void update_extrema(amrex::Real time);

// Update relaxation process.

void update_relaxation(amrex::Real time, amrex::Real dt);

// Input parameters

static int use_stopping_criterion;
static int use_energy_stopping_criterion;
static amrex::Real ts_te_stopping_criterion;
static amrex::Real T_stopping_criterion;

// Data we want to store, for interfacing with the Fortran

static int relaxation_is_done;
static int problem;

// Stellar masses

static amrex::Real mass_p;
static amrex::Real mass_s;

// Rate of change of stellar masses

static amrex::Real mdot_p;
static amrex::Real mdot_s;

// Stellar centers of mass and velocities

static amrex::Real com_p[3];
static amrex::Real com_s[3];

static amrex::Real vel_p[3];
static amrex::Real vel_s[3];

// Radii of the WDs at various density thresholds.

static amrex::Real rad_p[7];
static amrex::Real rad_s[7];

// Effective volume of the stars at various density cutoffs.

static amrex::Real vol_p[7];
static amrex::Real vol_s[7];

// Average density of the stars.

static amrex::Real rho_avg_p;
static amrex::Real rho_avg_s;

// Gravitational free-fall timescale of the stars.

static amrex::Real t_ff_p;
static amrex::Real t_ff_s;

// Global extrema of various quantities over the whole simulation

static amrex::Real T_global_max;
static amrex::Real rho_global_max;
static amrex::Real ts_te_global_max;

// Current values of the above extrema

static amrex::Real T_curr_max;
static amrex::Real rho_curr_max;
static amrex::Real ts_te_curr_max;

// Value of the total energy on the domain over the last several timesteps

static const int num_previous_ener_timesteps = 5;
static amrex::Real total_ener_array[num_previous_ener_timesteps];
