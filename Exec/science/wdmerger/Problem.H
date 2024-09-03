
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

// Computes standard dot product of two three-vectors.

static
amrex::Real dot_product(const amrex::Real a[], const amrex::Real b[]);


// Computes norm of a three vector.

static
amrex::Real norm(const amrex::Real a[]);

// Problem post-initialization routine.

void problem_post_init();

// Problem post-restart routine.

void problem_post_restart();

// End of timestep analytics on the white dwarf properties.

void problem_post_timestep();

// Write out the git hashes for the various parts of the code.

static
void writeGitHashes(std::ostream& log);

// Update relaxation process.

void update_relaxation(amrex::Real time, amrex::Real dt);

// Add some integrated quantities.

void problem_sums();
