#ifndef CASTRO_H
#define CASTRO_H

#include <state_indices.H>

// runtime parameters
#include <castro_params.H>
#include <prob_parameters.H>
#include <Castro_util.H>

using namespace castro;

#ifdef AMREX_PARTICLES
#include <particles_params.H>
#endif

#include <fundamental_constants.H>

constexpr int HISTORY_SIZE=40;
constexpr int PSTAR_BISECT_FACTOR = 5;


#include <AMReX_BC_TYPES.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_ErrorList.H>
#include <AMReX_FluxRegister.H>
#include <network.H>
#include <eos.H>
#ifndef TRUE_SDC
#ifdef REACTIONS
#include <burner.H>
#endif
#endif

#ifdef BL_LAZY
#include <AMReX_Lazy.H>
#endif

#ifdef AMREX_PARTICLES
#include <AMReX_AmrParticles.H>
#endif

#ifdef RADIATION
#include <RadSolve.H>
#endif

#include <memory>
#include <iostream>

#include <params_type.H>

using std::istream;
using std::ostream;

enum StateType { State_Type = 0,
#ifdef RADIATION
                 Rad_Type,
#endif
#ifdef GRAVITY
                 PhiGrav_Type,
                 Gravity_Type,
#endif
                 Source_Type,
#ifdef REACTIONS
                 Reactions_Type,
#endif
#ifdef MHD
               Mag_Type_x,
               Mag_Type_y,
               Mag_Type_z,
#endif
#ifdef SIMPLIFIED_SDC
#ifdef REACTIONS
                 Simplified_SDC_React_Type
#endif
#endif
};

// Create storage for all source terms.

enum sources { thermo_src = 0,
#ifdef DIFFUSION
               diff_src,
#endif
#ifdef HYBRID_MOMENTUM
               hybrid_src,
#endif
#ifdef GRAVITY
               grav_src,
#endif
#ifdef ROTATION
               rot_src,
#endif
#ifdef SPONGE
               sponge_src,
#endif
               ext_src,
               geom_src,
               num_src };


// time integration method

enum int_method { CornerTransportUpwind = 0,
                  UnusedTimeIntegration,
                  SpectralDeferredCorrections,
                  SimplifiedSpectralDeferredCorrections
                };

// Struct that returns information about
// why an advance failed.

struct advance_status {
    bool success;
    Real suggested_dt;
    std::string reason;
    advance_status() : success(true), suggested_dt(0.0_rt), reason("") {}
};

///
/// @class Castro
///
/// @brief ``AmrLevel``-derived class for hyperbolic conservation equations for stellar media
///
class Castro
    :
    public amrex::AmrLevel
{
public:

///
/// Default constructor.  Builds invalid object.
///
    Castro ();

///
/// The basic constructor.
///
/// @param papa     parent ``amrex::Amr`` object
/// @param lev      which level are we on
/// @param level_geom level geometry
/// @param bl       ``BoxArray`` object
/// @param dm       the mapping
/// @param time     current time
///
    Castro (amrex::Amr&            papa,
            int             lev,
            const amrex::Geometry& level_geom,
            const amrex::BoxArray& bl,
            const amrex::DistributionMapping& dm,
            amrex::Real            time);

///
/// The destructor.
///
    ~Castro () override;

///
/// Remove copy/move constructors/assignment operators.
///
    Castro (const Castro&) = delete;
    Castro (Castro&&) = delete;
    Castro& operator= (const Castro&) = delete;
    Castro& operator= (Castro&&) = delete;

///
/// Restart from a checkpoint file.
///
/// @param papa             parent ``amrex::Amr`` object
/// @param is
/// @param bReadSpecial
///
    void restart (amrex::Amr&     papa,
                  istream& is,
                  bool bReadSpecial = false) override;

///
/// This is called only when we restart from an old checkpoint.
/// Sets the Vector `state_in_checkpoint` depending on the input
/// version and type.
///
/// @param state_in_checkpoint      Vector of integer flags
///                                 describing the state
///
    void set_state_in_checkpoint (amrex::Vector<int>& state_in_checkpoint) override;

///
/// Call ``amrex::AmrLevel::checkPoint`` and then add radiation info
///
/// @param dir          Directory to store checkpoint in
/// @param os           ``std::ostream`` object
/// @param how          `V`isMF::How`` object
/// @param dump_old     boolean
///
    void checkPoint(const std::string& dir,
                    std::ostream&      os,
                    amrex::VisMF::How         how,
                    bool               dump_old) override;

///
/// A string written as the first item in writePlotFile() at
/// level zero. It is so we can distinguish between different
/// types of plot files. For Castro it has the form: Castro-Vnnn
///
    std::string thePlotFileType () const override;

    void setPlotVariables () override;

///
/// Write a plotfile to specified directory.
///
/// @param dir  directory to write to
/// @param os   output stream
/// @param how  ``VisMF::How`` object
///
    void writePlotFile (const std::string& dir,
                        ostream&       os,
                        amrex::VisMF::How     how) override;

///
/// @param dir  directory to write to
/// @param os   output stream
/// @param how  ``VisMF::How`` object
///
    void writeSmallPlotFile (const std::string& dir,
                             ostream&       os,
                             amrex::VisMF::How     how) override;
  
///
/// @param dir      directory to write to
/// @param os       output stream
/// @param how      ``VisMF::How`` object
/// @param is_small integer
///
    void plotFileOutput(const std::string& dir,
                        ostream& os,
                        amrex::VisMF::How how,
                        const int is_small);


///
/// Write job info to file
///
/// @param dir       directory to write to
/// @param io_time   how long it took to generate the plotfile/checkpoint
///
    void writeJobInfo (const std::string& dir, const amrex::Real io_time);


///
/// Dump build info 
///
    static void writeBuildInfo ();


///
/// Define data descriptors.
///
/// ::variableSetUp is called in the constructor of Amr.cpp, so
/// it should get called every time we start or restart a job
///
    static void variableSetUp ();

///
/// Define tagging functions.
///
    static void ErrorSetUp ();

///
/// Cleanup data descriptors at end of run.
///
/// ::variableSetUp is in ``Castro_setup.cpp``
/// ``variableCleanUp`` is called once at the end of a simulation
///
    static void variableCleanUp ();

///
/// Initialize grid data at problem start-up.
///
    void initData () override;

#ifdef AMREX_PARTICLES

///
/// Read particle-related inputs
///
    static void read_particle_params ();

///
/// Initialize particle locations and velocities (and strengths if relevant)
///
    void init_particles ();

///
/// Write particles in checkpoint directories
///
/// @param dir  directory containing checkpoint
///
    void ParticleCheckPoint (const std::string& dir);

///
/// Write particles in plotfile directories
///
/// @param dir  directory containing plotfile
///
    void ParticlePlotFile (const std::string& dir);

///
/// How to initialize at restart
///
/// @param restart_file     path to restart file
///
    void ParticlePostRestart (const std::string& restart_file);

///
/// Derived quantities associated with particles
///
/// @param name         name of particle variable to derive
/// @param time         current time
/// @param ngrow        Number of ghost cells
///
    std::unique_ptr<amrex::MultiFab> ParticleDerive (const std::string& name,
                                                     amrex::Real           time,
                                                     int            ngrow);

///
/// Timestamp particles
///
/// @param ngrow
///
    void TimestampParticles (int ngrow);

///
/// Advance the particles by dt
///
/// @param iteration    where we are in the current AMR subcycle
/// @param time         current time
/// @param dt           timestep
///
    void advance_particles (int iteration, amrex::Real time, amrex::Real dt);

#endif

#ifdef MAESTRO_INIT
    void MAESTRO_init ();
#endif


///
///
///
    amrex::MultiFab* Area ();

///
/// Area of the multifab with surface normal along ``dir``
///
/// @param dir  direction of normal
///
    amrex::MultiFab& Area (int dir);

///
/// The volume of the multifab.
///
    amrex::MultiFab& Volume ();

///
/// Return the n'th fluxes MultiFab.
///
/// @param dir  direction in which to find the fluxes along
///
    amrex::MultiFab& Fluxes (int dir);

///
/// Set time levels of state data.
///
/// @param time         current time
/// @param dt_old       old timestep
/// @param dt_new       new timestep
///
    void setTimeLevel (amrex::Real time,
                       amrex::Real dt_old,
                       amrex::Real dt_new) override;

///
/// Initialize data on this level from another Castro (during regrid).
///
/// @param old  amrex::AmrLevel object to use to initialize data
///
    void init (amrex::AmrLevel& old) override;

///
/// Initialize data on this level after regridding if old level did not
/// previously exist
///
    void init () override;

///
/// Proceed with next timestep?
///
    int okToContinue () override;

///
/// Advance grids at this level in time.
///
/// The main driver for a single level.  This will do either the SDC
/// algorithm or the Strang-split reactions algorithm.
///
/// @param time     the current simulation time
/// @param dt       the timestep to advance (e.g., go from time to
///                    time + dt)
/// @param iteration    where we are in the current AMR subcycle.  Each
///                    level will take a number of steps to reach the
///                    final time of the coarser level below it.  This
///                    counter starts at 1
/// @param ncycle   the number of subcycles at this level
///
    amrex::Real advance (amrex::Real time,
                         amrex::Real dt,
                         int  iteration,
                         int  ncycle) override;

///
/// This routine will advance the old state data (called ``S_old`` here)
/// to the new time, for a single level.  The new data is called
/// ``S_new`` here.  The update includes reactions (if we are not doing
/// SDC), hydro, and the source terms. The advance is doing using the CTU
/// (corner transport upwind) approach. We also use this interface for the
/// simplified SDC update.
///
/// @param time     the current simulation time
/// @param dt       the timestep to advance (e.g., go from time to
///                    time + dt)
///
    advance_status do_advance_ctu (amrex::Real time, amrex::Real dt);


#ifndef AMREX_USE_GPU
    amrex::Real do_advance_sdc (amrex::Real time,
                                amrex::Real dt,
                                int  amr_iteration,
                                int  amr_ncycle);
#endif

///
/// Save a copy of the old state data in case for the purposes of a retry.
///
    void save_data_for_retry();

///
/// Should we retry advancing the simulation? By default, we
/// don't do a retry unless the criteria are violated. This is
/// used only for the CTU advance.
///
/// @param dt       the timestep to advance (e.g., go from time to
///                    time + dt)
/// @param advance_status  was the advance successful?
///
    bool retry_advance_ctu(amrex::Real dt, const advance_status& status);

///
/// Subcyles until we've reached the target time, ``time`` + ``dt``.
/// The last timestep will be shortened if needed so that
/// it does not overshoot the ending time. This is used only for
/// the CTU advance.
///
/// @param time     the current simulation time
/// @param dt       the timestep to advance (e.g., go from time to
///                    time + dt)
/// @param amr_iteration    where we are in the current AMR subcycle.  Each
///                    level will take a number of steps to reach the
///                    final time of the coarser level below it.  This
///                    counter starts at 1
/// @param amr_ncycle   the number of subcycles at this level
///
    amrex::Real subcycle_advance_ctu(amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle);

///
/// This is run at the start of an ::advance on a single level.
/// It sets up all the necessary arrays and data containers,
/// ensures the data is valid and makes copies of the MultiFabs
/// in the old and new state in case a retry is performed.
///
/// @param time     the current simulation time
/// @param dt       the timestep to advance (e.g., go from time to
///                    time + dt)
/// @param amr_iteration    where we are in the current AMR subcycle.  Each
///                    level will take a number of steps to reach the
///                    final time of the coarser level below it.  This
///                    counter starts at 1
///
    void initialize_advance(amrex::Real time, amrex::Real dt, int amr_iteration);

///
/// If required, adds material lost in timestep to cumulative
/// losses, does refluxing, and clears arrays/data containers
/// initialized at start of the level's advance.
///
    void finalize_advance();

///
/// Performed at the start of ::do_advance, resets flags and
/// ensures required ghost zones are filled.
///
/// @param time     the current simulation time
/// @param dt       the timestep to advance
///
    advance_status initialize_do_advance (amrex::Real time, amrex::Real dt);

///
/// Performed at the end of ::do_advance, cleans up temporary data
/// created by ::initialize_do_advance
///
/// @param time     the current simulation time
/// @param dt       the timestep to advance
///
    advance_status finalize_do_advance (amrex::Real time, amrex::Real dt);

///
/// Check for NaNs in the given MultiFab
///
/// @param state        MultiFab to check
/// @param check_ghost  do we check the ghost cells?
///
    static void check_for_nan(amrex::MultiFab& state, int check_ghost=0);


#include <Castro_sources.H>


#ifdef GRAVITY
#include <Castro_gravity.H>
#endif


#ifdef MHD
#include <Castro_mhd.H>
#endif

///
/// Estimate time step.
///
    amrex::Real estTimeStep (int is_new = 1);


///
/// Compute the CFL timestep
///
    ValLocPair<amrex::Real, IntVect> estdt_cfl (int is_new = 1);


#ifdef MHD
///
/// Compute the CFL timestep
///
    ValLocPair<amrex::Real, IntVect>  estdt_mhd (int is_new = 1);
#endif

///
/// Diffusion-limited timestep
///
    ValLocPair<amrex::Real, IntVect> estdt_temp_diffusion (int is_new = 1);

#ifdef REACTIONS
///
/// Reactions-limited timestep
///
    ValLocPair<amrex::Real, IntVect> estdt_burning (int is_new = 1);
#endif

#ifdef RADIATION
///
/// Radiation hydro timestep
///
    amrex::Real estdt_rad (int is_new = 1);
#endif

///
/// Compute initial time step.
///
    amrex::Real initialTimeStep ();

///
/// Compute initial ``dt``.
///
/// @param finest_level     Index of finest level
/// @param subcycle         unused
/// @param n_cycle          the number of subcycles at this level
/// @param ref_ratio        unused
/// @param dt_level         Real Vector, will contain initial timestep at each level.
/// @param stop_time        End time of simulation
///
    void computeInitialDt (int                   finest_level,
                           int                   /*subcycle*/,
                           amrex::Vector<int>&           n_cycle,
                           const amrex::Vector<amrex::IntVect>& /*ref_ratio*/,
                           amrex::Vector<amrex::Real>&          dt_level,
                           amrex::Real                  stop_time) override;

///
/// Compute new ``dt``.
///
/// @param finest_level     Index of finest level
/// @param subcycle         unused
/// @param n_cycle          the number of subcycles at this level
/// @param ref_ratio        unused
/// @param dt_min
/// @param dt_level
/// @param stop_time        End time of simulation
/// @param post_regrid_flag Have we done regridding yet?
///
    void computeNewDt (int                   finest_level,
                       int                   /*subcycle*/,
                       amrex::Vector<int>&           n_cycle,
                       const amrex::Vector<amrex::IntVect>& /*ref_ratio*/,
                       amrex::Vector<amrex::Real>&          dt_min,
                       amrex::Vector<amrex::Real>&          dt_level,
                       amrex::Real                  stop_time,
                       int                   post_regrid_flag) override;

///
/// Allocate data at old time.
///
    void allocOldData () override;

///
/// Remove data at old time.
///
    void removeOldData () override;


///
/// Print information about energy budget.
///
    void do_energy_diagnostics();

///
/// Do work after timestep().
///
/// @param iteration    where we are in the current AMR subcycle
///
    void post_timestep (int iteration) override;

///
/// Contains operations to be done only after a full coarse timestep.
///
/// @param cumtime
///
    void postCoarseTimeStep (amrex::Real cumtime) override;

///
/// Do work after regrid().
///
/// @param lbase
/// @param new_finest
///
    void post_regrid (int lbase,
                      int new_finest) override;

///
/// Do work after a restart().
///
    void post_restart () override;

///
/// Do work after init().
///
/// @param stop_time    End time of simulation
///
    void post_init (amrex::Real stop_time) override;

#ifdef GRAVITY
#ifdef ROTATION
///
/// Initialize a model using the self-consistent field method.
///
    void scf_relaxation();

///
/// Perform the Hachisu SCF method.
///
    void do_hscf_solve();
#endif
#endif

///
/// Do work after restart with ``grown_factor`` > 1
///
    void post_grown_restart ();

#ifdef DO_PROBLEM_POST_SIMULATION

///
/// Do work at the end of the simulation - before the last outputs
///
/// @param amr_level
///
    static void problem_post_simulation (amrex::Vector<std::unique_ptr<AmrLevel> >& amr_level);
#endif


///
/// Error estimation for regridding.
///
/// @param tags         TagBoxArray of tags
/// @param clearval     integer
/// @param tagval       integer
/// @param time         current time
/// @param n_error_buf  integer
/// @param ngrow        Number of ghost cells
///
    void errorEst (amrex::TagBoxArray& tags,
                   int          clearval,
                   int          tagval,
                   amrex::Real  time,
                   int          n_error_buf = 0,
                   int          ngrow = 0) override;


///
/// Apply the ::problem_tagging routine.
///
/// @param tags         TagBoxArray of tags
/// @param time         current time
///
    void apply_problem_tags (amrex::TagBoxArray& tags, amrex::Real time);


///
/// Apply any tagging restrictions that must be satisfied by all problems.
///
/// @param tags         TagBoxArray of tags
/// @param time         current time
///
    void apply_tagging_restrictions (amrex::TagBoxArray& tags, amrex::Real time);


///
/// Returns a MultiFab containing the derived data for this level.
/// If ngrow>0 the MultiFab is built on the appropriately grown BoxArray.
///
/// @param name     Name of derived data
/// @param time     Current time
/// @param ngrow    Number of ghost cells
///
    std::unique_ptr<amrex::MultiFab> derive (const std::string& name,
                                             amrex::Real        time,
                                             int                ngrow) override;

///
/// This version of derive() fills the dcomp'th component of mf with the derived quantity.
///
/// @param name         Name of quantity to derive
/// @param time         current time
/// @param mf           MultiFab to store derived quantity in
/// @param dcomp        index of component of `mf` to fill with derived quantity
///
    void derive (const std::string& name,
                 amrex::Real               time,
                 amrex::MultiFab&          mf,
                 int                dcomp) override;


#ifdef REACTIONS
#include <Castro_react.H>
#endif

#ifdef ROTATION
#include <Castro_rotation.H>
#endif

#ifdef RADIATION

///
/// @param S_new
/// @param iteration    where we are in the current AMR subcycle
/// @param ncycle       the number of subcycles at this level
///
    void final_radiation_call (amrex::MultiFab& S_new,
                               int iteration,
                               int ncycle);
#endif

///
/// Reset the internal energy
///
/// @param State    Current state (MultiFab)
/// @param ng       number of ghost cells
///
    void reset_internal_energy (
#ifdef MHD
                                amrex::MultiFab& Bx,
                                amrex::MultiFab& By,
                                amrex::MultiFab& Bz,
#endif
                                amrex::MultiFab& State, int ng);


///
/// Reset the internal energy
///
/// @param bx       Box to update
/// @param u        Current state (Fab)
///
    static void reset_internal_energy (const amrex::Box& bx,
#ifdef MHD
                                amrex::Array4<amrex::Real> const Bx,
                                amrex::Array4<amrex::Real> const By,
                                amrex::Array4<amrex::Real> const Bz,
#endif
                                amrex::Array4<amrex::Real> const u);


#ifdef MHD
///
/// Add magnetic contribution to energy density
/// @param Bx       magnetic field in x    
/// @param By       magnetic field in y
/// @param Bz       magnetic field in z
/// @param state    the state to operate on
///
    void add_magnetic_e (
                         amrex::MultiFab& Bx,
                         amrex::MultiFab& By,
                         amrex::MultiFab& Bz,
                         amrex::MultiFab& state);

///
/// Check if divergence of B is zero at initialization
/// @param Bx       magnetic field in x    
/// @param By       magnetic field in y
/// @param Bz       magnetic field in z
/// @param state    the state to operate on
///
    void check_div_B (
                      amrex::MultiFab& Bx,
                      amrex::MultiFab& By,
                      amrex::MultiFab& Bz,
                      amrex::MultiFab& state);

#endif



///
/// Compute the current temperature
///
/// @param state    the state to operate on
/// @param time     current time
/// @param ng       number of ghost cells
///
    void computeTemp (
#ifdef MHD
                      amrex::MultiFab& Bx,
                      amrex::MultiFab& By,
                      amrex::MultiFab& Bz,
#endif
                      amrex::MultiFab& state, amrex::Real time, int ng);


///
/// Add any terms needed to correct the source terms.
/// Currently this is a lagged predictor for CTU and the
/// time-centering correction for simplified SDC.
///
    void create_source_corrector();


///
/// This is a hack to make sure that we only
/// ever have new data for certain state types that only
/// ever need new time data; by doing a swap now, we'll
/// guarantee that allocOldData() does nothing. We do
/// this because we never need the old data, so we
/// don't want to allocate memory for it.
///
/// @param dt   timestep
///
    void swap_state_time_levels (const amrex::Real dt);

#ifdef DIFFUSION
#include <Castro_diffusion.H>
#endif


///
/// Volume weighted sum of given quantity
///
/// @param mf           MultiFab containing quantity
/// @param comp         integer index of quantity in MultiFab
/// @param local        boolean, is sum local (over each patch) or over entire MultiFab?
/// @param finemask     boolean, should we build a mask to exclude finer levels?
///
    amrex::Real volWgtSum (const MultiFab& mf, int comp, bool local=false, bool finemask=true);

///
/// Volume weighted sum of given quantity
///
/// @param name         Name of quantity
/// @param time         current time
/// @param local        boolean, is sum local (over each patch) or over entire MultiFab?
/// @param finemask     boolean, should we build a mask to exclude finer levels?
///
    amrex::Real volWgtSum (const std::string& name, amrex::Real time, bool local=false, bool finemask=true);

///
/// Sum weighted by volume multiplied by distance from center in given direction
///
/// @param mf       MultiFab containing quantity
/// @param comp     integer index of quantity in MultiFab
/// @param idir     Axis along which to compute distance from center
/// @param local    boolean, is sum local (over each patch) or over entire MultiFab?
///
    amrex::Real locWgtSum (const MultiFab& mf, int comp, int idir, bool local=false);

///
/// Sum weighted by volume multiplied by distance from center in given direction
///
/// @param name     Name of quantity
/// @param time     current time
/// @param idir     Axis along which to compute distance from center
/// @param local    boolean, is sum local (over each patch) or over entire MultiFab?
///
    amrex::Real locWgtSum (const std::string& name, amrex::Real time, int idir, bool local=false);

///
/// Volume weighted sum of the product of two quantities
///
/// @param mf1      MultiFab containing first quantity
/// @param mf2      MultiFab containing second quantity
/// @param comp1    integer index of quantity in first MultiFab
/// @param comp2    integer index of quantity in second MultiFab
/// @param local    boolean, is sum local (over each patch) or over entire MultiFab?
///
    amrex::Real volProductSum (const MultiFab& mf1, const MultiFab& mf2, int comp1, int comp2, bool local=false);

///
/// Volume weighted sum of the product of two quantities
///
/// @param name1    Name of first quantity
/// @param name2    Name of second quantity
/// @param time     current time
/// @param local    boolean, is sum local (over each patch) or over entire MultiFab?
///
    amrex::Real volProductSum (const std::string& name1, const std::string& name2, amrex::Real time, bool local=false);


///
/// Location weighted sum of (quantity) squared
///
/// @param name     Name of quantity
/// @param time     current time
/// @param idir     Direction along which to weight
/// @param local    boolean, is sum local (over each patch) or over entire MultiFab?
///
    amrex::Real locSquaredSum (const std::string& name, amrex::Real time, int idir, bool local=false);

#ifdef GRAVITY
///
/// Calculate the gravitational wave signal
///
/// @param time       current time
/// @param h_plus_1   gravitational wave signal in x-direction
/// @param h_cross_1  gravitational wave signal in x-direction
/// @param h_plus_2   gravitational wave signal in y-direction
/// @param h_cross_2  gravitational wave signal in y-direction
/// @param h_plus_3   gravitational wave signal in z-direction
/// @param h_cross_3  gravitational wave signal in z-direction
/// @param local      is sum local

    void gwstrain (Real time,
                   Real& h_plus_1, Real& h_cross_1,
                   Real& h_plus_2, Real& h_cross_2,
                   Real& h_plus_3, Real& h_cross_3,
                   bool local);
#endif

#ifdef GRAVITY
///
/// Are we using point mass gravity?
///
    static int using_point_mass ();

///
/// Return the value of the point mass.
///
    static amrex::Real get_point_mass ();
#endif


///
/// Fill a version of the state with ``ng`` ghost zones from the state data.
/// After we do the fill, optionally clean the data using the iclean parameter.
///
/// @param S        MultiFab to be filled
/// @param time     current time
/// @param ng       number of ghost cells
///
    void expand_state(amrex::MultiFab& S, amrex::Real time, int ng);



// Hydrodynamics
#include <Castro_hydro.H>

// SDC 
#ifdef TRUE_SDC
#include <Castro_sdc.H>
#endif

    static void extern_init ();

    static void read_params ();


///
/// @param lev
///
    Castro& getLevel (int lev);

    void FluxRegCrseInit();
    void FluxRegFineAdd();


///
/// Do refluxing
///
/// @param crse_level
/// @param fine_level
/// @param in_post_timestep
///
    void reflux (int crse_level, int fine_level, bool in_post_timestep);


///
/// Normalize species fractions so they sum to 1
///
/// @param S_new    State
/// @param ng       number of ghost cells
///
    static void normalize_species (amrex::MultiFab& S_new, int ng);


///
/// Enforces
/// \f[
/// \rho E = \rho e + \frac{1}{2} \rho (u^2 + v^2 + w^2)
/// \f]
///
/// @param S  State
///
    static void enforce_consistent_e (
#ifdef MHD
                               amrex::MultiFab& Bx, amrex::MultiFab& By, amrex::MultiFab& Bz,
#endif
                               amrex::MultiFab& S);

///
/// This routine sets the density in ``state`` to be larger than the density floor.
/// Note that it will operate everywhere on ``state``, including ghost zones.
///
/// @param state        state
/// @param ng           number of ghost cells
///
    void enforce_min_density (amrex::MultiFab& state, int ng);

///
/// After a hydro advance, check to see whether any density went too small.
///
    advance_status check_for_negative_density ();

///
/// Ensure the magnitude of the velocity is not larger than ``castro.speed_limit``
/// in any zone.
///
/// @param state        state
/// @param ng           number of ghost cells
///
    static void enforce_speed_limit (amrex::MultiFab& state, int ng);

///
/// Given ``State_Type`` state data, perform a number of cleaning steps to make
/// sure the data is sensible.
///
/// @param state    State data
/// @param time     current time
/// @param ng       number of ghost cells
///
    void clean_state (
#ifdef MHD
                      amrex::MultiFab& Bx, amrex::MultiFab& By, amrex::MultiFab& Bz,
#endif
                      amrex::MultiFab& state, amrex::Real time, int ng);

///
/// Average new state from ``level+1`` down to ``level``
///
    void avgDown ();


///
/// Average given component of new state from ``level+1`` down to ``level``
///
/// @param state_indx   Index of state to average down
///
    void avgDown (int state_indx);

    void buildMetrics ();


///
/// Initialize the MultiFabs and flux registers that live as class members.
///
    void initMFs ();

///
/// Calculates volume weight sums of variables and prints to screen
///
    void sum_integrated_quantities ();

///
/// Problem-specific diagnostics (called by sum_integrated_quantities)
///
    void problem_diagnostics ();

    void write_info ();

///
/// Return BCs
///
    static amrex::BCRec physbc () { return phys_bc; }

#ifdef GRAVITY

///
/// Find the position of the "center" by interpolating from data at cell centers
///
/// @param S        Current state
/// @param time     current time
///
    void define_new_center (amrex::MultiFab& S, amrex::Real time);
    void write_center ();

///
/// @param time     current time
/// @param dt       timestep
///
    void pointmass_update(amrex::Real time, amrex::Real dt);
#endif

    static void stopJob ();

///
/// Returns true if we're on the GPU and the total
/// memory on this level oversubscribes GPU memory.
///
    bool oversubscribing() {    // NOLINT(readability-convert-member-functions-to-static)
#ifdef AMREX_USE_GPU
        return (amrex::MultiFab::queryMemUsage("AmrLevel_Level_" + std::to_string(level)) >= 
                static_cast<Long>(amrex::Gpu::Device::totalGlobalMem()));
#else
        return false;
#endif
    }

#ifdef AMREX_PARTICLES
    static amrex::AmrTracerParticleContainer* theTracerPC () { return TracerPC; }
#endif

    static amrex::Vector<std::string> source_names;

///
/// Vector storing list of error tags.
///
    static amrex::Vector<amrex::AMRErrorTag> error_tags;

///
/// This MultiFab is on the coarser level.  This is useful for the coarser level
///     to mask out the finer level.  We only build this when it is needed.
///     This coarse MultiFab has to live on the fine level because it must be updated
///     even when only the fine level changes.
///
    amrex::MultiFab fine_mask;
    amrex::MultiFab& build_fine_mask();


///
/// A record of how many cells we have advanced throughout the simulation.
///     This is saved as a real because we will be storing the number of zones
///     advanced as a ratio with the number of zones on the coarse grid (to
///     help prevent the number from getting too large), and that may be a
///     non-integer number.
///
    static amrex::Real num_zones_advanced;

///
/// diagnostics
///
    static Vector<std::unique_ptr<std::fstream> > data_logs;
    static Vector<std::unique_ptr<std::fstream> > problem_data_logs;

///
/// runtime parameters
//
    static params_t params;


protected:


///
/// Build a mask that ghost cells overlapping with interior cells in the same multifab
///     are set to 0, whereas others are set to 1.
///
    amrex::Vector<std::unique_ptr<amrex::iMultiFab> > ib_mask;

///
/// @param ng number of ghost cells
///
    amrex::iMultiFab& build_interior_boundary_mask (int ng);

#ifdef GRAVITY
    int get_numpts();


///
/// Difference between composite and level solve for phi.
///
    amrex::MultiFab comp_minus_level_phi;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > comp_minus_level_grad_phi;
#endif


#ifdef REACTIONS
///
/// Storage for the burn_weights
///
    amrex::MultiFab burn_weights;
    static std::vector<std::string> burn_weight_names;
#endif


///
/// A state array with ghost zones.
///
    amrex::MultiFab Sborder;

#ifdef MHD
   amrex::MultiFab Bx_old_tmp;
   amrex::MultiFab By_old_tmp;
   amrex::MultiFab Bz_old_tmp;
#endif

#ifdef RADIATION
    amrex::MultiFab Erborder;
    amrex::MultiFab lamborder;

    std::unique_ptr<RadSolve> rad_solver;
#endif



#ifdef TRUE_SDC
///
/// A state array for the some reaction stuff with SDC
///
    amrex::MultiFab Sburn;


///
/// The primitive variable state array.
///
    amrex::MultiFab q;

///
/// we need a second version when doing fourth order
///
    amrex::MultiFab q_bar;

///
/// for diffusion in fourth order we need the cell-center temperature
//
#ifdef DIFFUSION
    amrex::MultiFab T_cc;
#endif

///
/// The auxiliary primitive variable state array.
///
    amrex::MultiFab qaux;
    amrex::MultiFab qaux_bar;
#endif


///
/// Source term corrector.
///
    amrex::MultiFab source_corrector;


///
/// Hydrodynamic (and radiation) fluxes.
///
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > fluxes;
#if (AMREX_SPACEDIM <= 2)
    amrex::MultiFab         P_radial;
#endif
#ifdef RADIATION
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > rad_fluxes;
#endif

    amrex::Vector<std::unique_ptr<amrex::MultiFab> > mass_fluxes;

    amrex::FluxRegister flux_reg;
#if (AMREX_SPACEDIM <= 2)
    amrex::FluxRegister pres_reg;
#endif
#ifdef RADIATION
    amrex::FluxRegister rad_flux_reg;
#endif
#ifdef GRAVITY
    amrex::FluxRegister phi_reg;
#endif

///
/// Scalings for the flux registers.
///
    amrex::Real flux_crse_scale;
    amrex::Real flux_fine_scale;

#if (AMREX_SPACEDIM <= 2)
    amrex::Real pres_crse_scale;
    amrex::Real pres_fine_scale;
#endif



///
/// State data to hold if we want to do a retry.
///
    amrex::Vector<std::unique_ptr<amrex::StateData> > prev_state;



///
/// Flag for indicating that we want to save prev_state until the reflux.
///
    bool keep_prev_state;


#ifdef TRUE_SDC
    //
    // Storage for the SDC time integration

    // this is the new iterations solution at the time nodes
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > k_new;

    // this is the old value of the advective update at the
    // nodes of the time integration
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > A_old;

    // this is the new value of the advective update at the
    // nodes of the time integration
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > A_new;

    // this is the old value of the reaction source at the
    // nodes of the time integration
#ifdef REACTIONS
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > R_old;
#endif

    static int SDC_NODES;
    static amrex::Vector<amrex::Real> dt_sdc;
    static amrex::Vector<amrex::Real> node_weights;
#endif

///
/// Wall time that we started the timestep
///
    amrex::Real wall_time_start;


///
/// The data.
///
    amrex::MultiFab             volume;
    amrex::Array<amrex::MultiFab,3> area;
    amrex::MultiFab             dLogArea[1];
    amrex::Vector< amrex::Vector<amrex::Real> > radius;


///
/// For initialization of the C++ values of the runtime parameters
  static int do_cxx_prob_initialize;


///
/// Keep track of which AMR iteration we're on.
///
    int iteration;

    //
    // Static data members.
    //

    static bool      signalStopJob;
    static int       radius_grow;
    static std::vector<std::string> err_list_names;
    static std::vector<int> err_list_ng;
    static int              num_err_list_default;
    static amrex::BCRec     phys_bc;
    static int       NUM_GROW;
    static int       NUM_GROW_SRC;

    static int         lastDtPlotLimited;
    static amrex::Real lastDtBeforePlotLimiting;

    int in_retry;
    int num_subcycles_taken;

    amrex::Real lastDt;


///
/// for keeping track of the amount of CPU time used -- this will persist
///     after restarts
///
    static amrex::Real      previousCPUTimeUsed;
    static amrex::Real      startCPUTime;

///
/// Get total CPU time
///
    static amrex::Real getCPUTime();

    bool             FillPatchedOldState_ok;


///
/// There can be only one Gravity object, it covers all levels:
///
    static class Gravity *gravity;


///
/// There can be only one Diffusion object, it covers all levels:
///
    static class Diffusion *diffusion;

#ifdef RADIATION


///
/// There can be only one Radiation object, it covers all levels:
///
    static class Radiation *radiation;

    friend class Radiation;
#endif

#ifdef AMREX_PARTICLES
    static amrex::AmrTracerParticleContainer* TracerPC;
#endif

    static amrex::IntVect hydro_tile_size;
    static amrex::IntVect no_tile_size;

    static int hydro_tile_size_has_been_tuned;
    static Long largest_box_from_hydro_tile_size_tuning;

    static int SDC_Source_Type;
    static int num_state_type;


    // counters for various retries in Castro

///
/// subcycle capability
///
    int sub_iteration;
    int sub_ncycle;
    amrex::Real dt_subcycle;
    amrex::Real dt_advance;


///
/// sdc
///
    int sdc_iteration;

#ifdef TRUE_SDC
    int current_sdc_node;
#endif


/* problem-specific includes */
public:
#include <Problem.H>

};

//
// Inlines.
//

inline
amrex::MultiFab*
Castro::Area ()
{
    return area.data();
}


///
/// @param dir
///
inline
amrex::MultiFab&
Castro::Area (int dir)
{
    return area[dir];
}

inline
amrex::MultiFab&
Castro::Volume ()
{
    return volume;
}


///
/// @param lev
///
inline
Castro&
Castro::getLevel (int lev)
{
    return *dynamic_cast<Castro *>(&parent->getLevel(lev));
}

inline
void
Castro::stopJob()
{
  signalStopJob = true;


}

inline
int
Castro::get_output_at_completion()
{
    return output_at_completion;
}

#ifdef GRAVITY
inline
amrex::Real
Castro::get_point_mass ()
{
    return point_mass;
}

inline
int
Castro::using_point_mass ()
{
    return use_point_mass;
}
#endif

#endif /*_Castro_H_*/
