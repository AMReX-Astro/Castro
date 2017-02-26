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

void wd_update(Real time, Real dt);

// Volume within white dwarfs at a given density cutoff.

void volInBoundary (Real time, Real& vol_p, Real& vol_s, Real rho_cutoff, bool local = false);

// Calculate gravitational wave signal.

void gwstrain (Real time,
	       Real& h_plus_1, Real& h_cross_1,
	       Real& h_plus_2, Real& h_cross_2,
	       Real& h_plus_3, Real& h_cross_3,
	       bool local = false);

// Computes standard dot product of two three-vectors.

Real dot_product(const Real a[], const Real b[]);

// Computes standard cross-product of two three-vectors.

void cross_product(const Real a[], const Real b[], Real c[]);

// Computes norm of a three vector.

Real norm(const Real a[]);

// Problem post-initialization routine.

void problem_post_init();

// Problem post-restart routine.

void problem_post_restart();

// End of timestep analytics on the white dwarf properties.

void problem_post_timestep();

// Write out the git hashes for the various parts of the code.

void writeGitHashes(std::ostream& log);

// Check whether we should stop the run.

void check_to_stop(Real time);

// Calculate extrema on grid.

void update_extrema(Real time);

// Update relaxation process.

void update_relaxation(Real time, Real dt);

// Input parameters

static int use_stopping_criterion;

// Data we want to store, for interfacing with the Fortran

static int relaxation_is_done;
static int problem;

// Stellar masses

static Real mass_p;
static Real mass_s;

// Rate of change of stellar masses

static Real mdot_p;
static Real mdot_s;

// Stellar centers of mass and velocities

static Real com_p[3];
static Real com_s[3];

static Real vel_p[3];
static Real vel_s[3];

// Radii of the WDs at various density thresholds.

static Real rad_p[7];
static Real rad_s[7];

// Effective volume of the stars at various density cutoffs.

static Real vol_p[7];
static Real vol_s[7];

// Average density of the stars.

static Real rho_avg_p;
static Real rho_avg_s;

// Gravitational free-fall timescale of the stars.

static Real t_ff_p;
static Real t_ff_s;

// Global extrema of various quantities over the whole simulation

static Real T_global_max;
static Real rho_global_max;
static Real ts_te_global_max;

// Current values of the above extrema

static Real T_curr_max;
static Real rho_curr_max;
static Real ts_te_curr_max;

// Value of the total energy on the domain over the last several timesteps

static const int num_previous_ener_timesteps = 5;
static Real total_ener_array[num_previous_ener_timesteps];
