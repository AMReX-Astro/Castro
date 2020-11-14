
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
