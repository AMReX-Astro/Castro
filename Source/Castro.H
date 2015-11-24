
#ifndef _Castro_H_
#define _Castro_H_

#include <BC_TYPES.H>
#include <AmrLevel.H>
#include <ErrorList.H>
#include <FluxRegister.H>
#include <BL_CXX11.H>

#ifdef BL_LAZY
#include <Lazy.H>
#endif

#ifdef PARTICLES
#include "AmrParticles.H"
typedef AmrParticleContainer<1+BL_SPACEDIM> DarkMatterParticleContainer;  // mass + velocities
#endif

#include <iostream>

using std::istream;
using std::ostream;

#ifdef LEVELSET
typedef BaseFab<int> IntFab;
typedef FabArray<IntFab> IntMultiFab;
#endif

enum StateType { State_Type = 0,
#ifdef RADIATION
                 Rad_Type,
#endif
#ifdef GRAVITY
                 PhiGrav_Type,
                 Gravity_Type,
#endif
#ifdef ROTATION
		 PhiRot_Type,
		 Rotation_Type,
#endif
		 Source_Type,
#ifdef LEVELSET
                 LS_State_Type,
#endif
#ifdef REACTIONS
                 Reactions_Type,
#endif
#ifdef SGS
                 SGS_Type,
#endif
                 NUM_STATE_TYPE };

//
// AmrLevel-derived class for hyperbolic conservation equations for stellar media
//

class Castro
    :
    public AmrLevel
{
public:
    //
    //Default constructor.  Builds invalid object.
    //
    Castro ();
    //
    //The basic constructor.
    //
    Castro (Amr&            papa,
            int             lev,
            const Geometry& level_geom,
            const BoxArray& bl,
            Real            time);
    //
    //The destructor.
    //
    virtual ~Castro () BL_OVERRIDE;
    //
    //Restart from a checkpoint file.
    //
    virtual void restart (Amr&     papa,
                          istream& is,
			  bool bReadSpecial = false) BL_OVERRIDE;
    //
    //This is called only when we restart from an old checkpoint.
    //
    virtual void set_state_in_checkpoint (Array<int>& state_in_checkpoint) BL_OVERRIDE;
    //
    //Call AmrLevel::checkPoint and then add radiation info
    //
    virtual void checkPoint(const std::string& dir,
                            std::ostream&      os,
                            VisMF::How         how,
                            bool               dump_old) BL_OVERRIDE;

    /*A string written as the first item in writePlotFile() at
               level zero. It is so we can distinguish between different
               types of plot files. For Castro it has the form: Castro-Vnnn
    */
    virtual std::string thePlotFileType () const BL_OVERRIDE;

    virtual void setPlotVariables () BL_OVERRIDE;
    //
    //Write a plotfile to specified directory.
    //
    virtual void writePlotFile (const std::string& dir,
                                ostream&       os,
                                VisMF::How     how) BL_OVERRIDE;
    //
    //Define data descriptors.
    //
    static void variableSetUp ();
    //
    //Define tagging functions.
    //
    static void ErrorSetUp ();
    //
    //Cleanup data descriptors at end of run.
    //
    static void variableCleanUp ();
    //
    //Initialize grid data at problem start-up.
    //
    virtual void initData () BL_OVERRIDE;

#ifdef PARTICLES
    //
    // Read particle-related inputs
    //
    static void read_particle_params ();

    //
    // Write particles in checkpoint directories
    //
    void ParticleCheckPoint (const std::string& dir);

    //
    // Write particles in plotfile directories
    //
    void ParticlePlotFile (const std::string& dir);

    //
    // How to initialize at restart
    //
    void ParticlePostRestart (const std::string& restart_file);

    //
    // Redistribute
    //
    void ParticleRedistribute ();

    //
    // Move randomly 
    //
    void ParticleMoveRandom ();

    //
    //Initialize particle locations and velocities (and strengths if relevant)
    //
    virtual void init_particles () BL_OVERRIDE;

    //
    // Time step control based on particles
    //
    void ParticleEstTimeStep (Real& estdt);

    //
    // Derived quantities associated with particles
    //
    MultiFab* ParticleDerive (const std::string& name,
                              Real           time,
                              int            ngrow);

    //
    // Default verbosity of Particle class
    //
    static int particle_verbose;

    //
    // Default cfl of particles in Particle class
    //
    static Real particle_cfl;

    //
    // Should we write particles into plotfiles?
    //
    static int write_particles_in_plotfile;
#endif

#ifdef MAESTRO_INIT
    void MAESTRO_init ();
#endif

    //
    //Get problem metrics.
    //
    MultiFab* Area ();
    MultiFab& Area (int dir);
    //
    //The volume of the multifab.
    //
    MultiFab& Volume ();
    //
    // Get Edge BoxArrays
    //
    const BoxArray& getEdgeBoxArray (int dir) const;
    //
    //Set time levels of state data.
    //
    virtual void setTimeLevel (Real time,
                               Real dt_old,
                               Real dt_new) BL_OVERRIDE;
    //
    //Initialize data on this level from another Castro (during regrid).
    //
    virtual void init (AmrLevel& old) BL_OVERRIDE;
    //
    //Initialize data on this level after regridding if old level did not previously exist
    //
    virtual void init () BL_OVERRIDE;
    //
    //Proceed with next timestep?
    //
    virtual int okToContinue () BL_OVERRIDE;
    //
    //Advance grids at this level in time.
    //
    virtual Real advance (Real time,
                          Real dt,
                          int  iteration,
                          int  ncycle) BL_OVERRIDE;

    Real advance_hydro (Real time,
                        Real dt,
                        int  iteration,
                        int  ncycle);

    Real advance_no_hydro (Real time,
                           Real dt,
                           int  iteration,
                           int  ncycle);

    void time_center_source_terms (MultiFab& S_new,
                                   MultiFab& ext_src_old,
                                   MultiFab& ext_src_new,
                                   Real dt);

    //
    //Estimate time step.
    //
    Real estTimeStep (Real dt_old);
    //
    //Compute initial time step.
    //
    Real initialTimeStep ();
    //
    //Compute initial `dt'.
    //
    virtual void computeInitialDt (int                   finest_level,
                                   int                   sub_cycle,
                                   Array<int>&           n_cycle,
                                   const Array<IntVect>& ref_ratio,
                                   Array<Real>&          dt_level,
                                   Real                  stop_time) BL_OVERRIDE;
    //
    //Compute new `dt'.
    //
    virtual void computeNewDt (int                   finest_level,
                               int                   sub_cycle,
                               Array<int>&           n_cycle,
                               const Array<IntVect>& ref_ratio,
                               Array<Real>&          dt_min,
                               Array<Real>&          dt_level,
                               Real                  stop_time,
                               int                   post_regrid_flag) BL_OVERRIDE;
    //
    //Allocate data at old time.
    //
    virtual void allocOldData () BL_OVERRIDE;
    //
    //Remove data at old time.
    //
    virtual void removeOldData () BL_OVERRIDE;
    //
    //Print information about energy budget.
    //
    void do_energy_diagnostics();
    //
    //Do work after timestep().
    //
    virtual void post_timestep (int iteration) BL_OVERRIDE;
    //
    // Contains operations to be done only after a full coarse timestep.
    //
    virtual void postCoarseTimeStep (Real cumtime) BL_OVERRIDE;
    //
    //Do work after regrid().
    //
    virtual void post_regrid (int lbase,
                              int new_finest) BL_OVERRIDE;   
    //
    //Do work after a restart().
    //
    virtual void post_restart () BL_OVERRIDE;
    //
    //Do work after init().
    //
    virtual void post_init (Real stop_time) BL_OVERRIDE;
    //
    //Do work after restart with grown_factor > 1
    //
    void post_grown_restart ();
    //
    //Error estimation for regridding.
    //
    virtual void errorEst (TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
                           Real         time,
			   int          n_error_buf = 0,
			   int          ngrow = 0) BL_OVERRIDE;

    // Returns a MultiFab containing the derived data for this level.
    // The user is responsible for deleting this pointer when done
    // with it.  If ngrow>0 the MultiFab is built on the appropriately
    // grown BoxArray.
    virtual MultiFab* derive (const std::string& name,
			      Real               time,
			      int                ngrow) BL_OVERRIDE;

    // This version of derive() fills the dcomp'th component of mf with the derived quantity.
    virtual void derive (const std::string& name,
			 Real               time,
			 MultiFab&          mf,
			 int                dcomp) BL_OVERRIDE;

    static int numGrow();

#ifdef REACTIONS
    void react_half_dt(MultiFab& state,
		       MultiFab& reactions,
#ifdef TAU
		       MultiFab& tau,
#endif
		       Real time, 
		       Real dt,
		       int  ngrow = 0);

    void strang_chem (MultiFab&  state,
                      MultiFab&  react_mf,
#ifdef TAU
                      MultiFab&  tau,
#endif
                      Real       time,
                      Real       dt);
#endif

#ifdef ROTATION
    void fill_rotation_field(MultiFab& phi, MultiFab& rot, MultiFab& state, Real time);
#endif

#ifdef RADIATION
    void final_radiation_call (MultiFab& S_new, 
                               int iteration, 
                               int ncycle);
#endif

#ifdef SGS 
    void reset_old_sgs (Real dt);
    void reset_new_sgs (Real dt);
#endif

    void reset_internal_energy (MultiFab& State);

    void computeTemp (MultiFab& State);

    void update_point_mass (MultiFab& StateIn, MultiFab& StateOut, 
                            PArray<FArrayBox>& flux, MultiFab& volume,
                            Real dt);

#ifdef SGS
    void getSource (Real time, Real dt, MultiFab& state, MultiFab& src, MultiFab& sgs_state, MultiFab* sgs_fluxes);
#else
    void getSource (Real time, Real dt, MultiFab& state_old, MultiFab& state_new, MultiFab& src, int ngrow = 0);
#endif

#ifdef DIFFUSION
    void getTempDiffusionTerm (Real time, MultiFab& TempDiffTerm, MultiFab* tau = 0);
    void getSpecDiffusionTerm (Real time, MultiFab& SpecDiffTerm);
#if (BL_SPACEDIM == 1)
    void getViscousTerm (Real time, MultiFab& ViscousTermforMomentum, MultiFab& ViscousTermforEnergy);
    void getFirstViscousTerm (Real time, MultiFab& ViscousTerm);
    void getSecndViscousTerm (Real time, MultiFab& ViscousTerm);
    void getViscousTermForEnergy (Real time, MultiFab& ViscousTermforEnergy);
#endif
#ifdef TAU
    void define_tau (MultiFab& tau_diff, MultiFab& grav_vector, Real time);
    void add_temp_diffusion_to_source (MultiFab& ext_src, MultiFab& TempDiffTerm, Real t, MultiFab& tau_diff);
    void time_center_temp_diffusion   (MultiFab& S_new, MultiFab& OldTempDiffTerm, Real  cur_time, Real dt, MultiFab& tau_diff);
    void full_temp_diffusion_update   (MultiFab& S_new, Real prev_time, Real cur_time, Real dt, MultiFab& tau_diff);
#else
    void add_temp_diffusion_to_source (MultiFab& ext_src, MultiFab& TempDiffTerm, Real t);
    void time_center_temp_diffusion   (MultiFab& S_new, MultiFab& OldTempDiffTerm, Real  cur_time, Real dt);
    void full_temp_diffusion_update   (MultiFab& S_new, Real prev_time, Real cur_time, Real dt);
#endif
    void add_spec_diffusion_to_source (MultiFab& ext_src, MultiFab& TempDiffTerm, Real t);
    void time_center_spec_diffusion   (MultiFab& S_new, MultiFab& OldTempDiffTerm, Real  cur_time, Real dt);
    void full_spec_diffusion_update   (MultiFab& S_new, Real prev_time, Real cur_time, Real dt);
#if (BL_SPACEDIM == 1)
    void add_viscous_term_to_source   (MultiFab& ext_src, MultiFab& ViscousTermforMomentum, 
                                       MultiFab& ViscousTermforEnergy, Real t);
    void time_center_viscous_term     (MultiFab& S_new, MultiFab& OldViscousTermforMomentum, 
                                       MultiFab& ViscousTermforEnergy, Real  cur_time, Real dt);
#endif
#endif

#ifdef ROTATION
  void add_rotation_to_source(MultiFab& ext_src, MultiFab& RotationTerm, Real old_time);
  void time_center_rotation(MultiFab& S_new, MultiFab& OldRotationTerm, Real cur_time, Real dt);
#endif

    Real volWgtSum (const std::string& name, Real time, bool local=false, bool finemask=true);

    Real volWgtSquaredSum (const std::string& name, Real time, bool local=false);

    Real locWgtSum (const std::string& name, Real time, int idir, bool local=false);

    Real locWgtSum2D (const std::string& name, Real time, int idir1, int idir2, bool local=false);

    Real volWgtSumMF (MultiFab* mf, int comp, bool local=false);

    Real locWgtSumOneSide (const std::string& name, Real time, int idir, int side, int bdir, bool local=false);

    Real volWgtSumOneSide (const std::string& name, Real time, int side, int bdir, bool local=false);

    Real volProductSum (const std::string& name1, const std::string& name2, Real time, bool local=false);

    Real locSquaredSum (const std::string& name, Real time, int idir, bool local=false);

    Real get_point_mass ();

    void set_special_tagging_flag (Real time);

    void make_radial_data (int is_new);

#ifdef LEVELSET
    void advance_levelset(Real time, Real dt);
    void reinit_phi(Real time);
#endif

#ifdef AUX_UPDATE
    void advance_aux(Real time, Real dt);
#endif

#ifdef PARTICLES
    static DarkMatterParticleContainer* theDMPC ();
#endif

    static int       NUM_STATE;
    static int       Density, Xmom, Ymom, Zmom, Eden, Eint, Temp;
#ifdef SGS
    static int       Esgs;
    static Real      sum_turb_src;
#endif
    static int       FirstAdv,  NumAdv;
    static int       FirstSpec, NumSpec;
    static int       FirstAux,  NumAux;

    FluxRegister*        flux_reg;

    FluxRegister& getFluxReg ();
    FluxRegister& getFluxReg (int lev);

    //
    // This MultiFab is on the coarser level.  This is useful for the coarser level
    // to mask out the finer level.  We only build this when it is needed.
    // This coarse MultiFab has to live on the fine level because it must be updated
    // even when only the fine level changes.
    //
    MultiFab* fine_mask;
    MultiFab* build_fine_mask();

protected:

    int get_numpts();

    //
    // Need to pass ugdnv out of cmpflx to do level sets
    //
    MultiFab* u_gdnv;

    //
    //  Call extern/networks/.../network.f90::network_init()
    //
    static void extern_init ();

    static void network_init ();

    static void read_params ();

    Castro& getLevel (int lev);

#ifdef SGS
    FluxRegister& getSGSFluxReg ();
    FluxRegister& getSGSFluxReg (int lev);
#endif

#ifdef RADIATION
    FluxRegister& getRADFluxReg ();
    FluxRegister& getRADFluxReg (int lev);
#endif

    void reflux ();

    void enforce_nonnegative_species (MultiFab& S_new);

    void enforce_consistent_e (MultiFab& S);

    void avgDown ();

    void avgDown (int state_indx);

    void buildMetrics ();

    // integrate derived quantities over domain

    Real sumDerive (const std::string& name, Real time, bool local=false);

    void sum_integrated_quantities ();

    void write_info ();

#ifdef GRAVITY
    void define_new_center (MultiFab& S, Real time);
    void write_center ();
#endif

    //
    // The data.
    //
    MultiFab             volume;
    MultiFab             area[BL_SPACEDIM];
    MultiFab             dLogArea[1];
    Array< Array<Real> > radius;
#ifdef SGS
    FluxRegister*        sgs_flux_reg;
#endif
#ifdef RADIATION
    FluxRegister*        rad_flux_reg;
#endif

    //
    // Static data members.
    //
    static int       checkpoint_version;
    static bool      dump_old;
    static int       radius_grow;
    static int       verbose;
    static Real      cfl;
    static Real      burning_timestep_factor;  
    static Real      init_shrink;
    static Real      change_max;
    static int       do_reflux;
    static ErrorList err_list;
    static BCRec     phys_bc;
    static int       NUM_GROW;

    static int       grown_factor;
    static int       star_at_center;
    static int       moving_center;
    static int       normalize_species;
    static int       fix_mass_flux;
    static int       allow_negative_energy;
    static int       do_special_tagging;

    static int       ppm_type;
    static int       ppm_reference;
    static int       ppm_trace_sources;
    static int       ppm_temp_fix;
    static int       ppm_tau_in_tracing;
    static int       ppm_predict_gammae;
    static int       ppm_reference_edge_limit;
    static int       ppm_flatten_before_integrals;
    static int       ppm_reference_eigenvectors;
    static int       hybrid_riemann;
    static int       use_colglaz;
    static int       riemann_solver;
    static int       use_flattening;
    static int       transverse_use_eos;
    static int       transverse_reset_density;
    static int       transverse_reset_rhoe;
    static int       source_term_predictor;

    static int       dual_energy_update_E_from_e;
    static Real      dual_energy_eta1;
    static Real      dual_energy_eta2;
    static Real      dual_energy_eta3;

    static int       cg_maxiter;
    static Real      cg_tol;
    static int       use_pslope;
    static int       grav_source_type;
    static int       spherical_star;
    static int       hard_cfl_limit;

    static Real      difmag;
    static Real      small_dens;
    static Real      small_temp;
    static Real      small_pres;
    static Real      small_ener;

    static std::string      job_name;

    // for keeping track of the amount of CPU time used -- this will persist
    // after restarts
    static Real      previousCPUTimeUsed;
    static Real      startCPUTime;

    static Real getCPUTime();

    // for keeping track of energy change from burning
    static Real      max_delta_e;

#ifdef POINTMASS
    static Real      point_mass;
#endif

#ifdef DIFFUSION
    static int       diffuse_temp;
    static int       diffuse_spec;
    static int       diffuse_vel;
    static Real      diffuse_cutoff_density;
#endif

#ifdef PARTICLES
    static bool do_dm_particles ;

    // How do we want to initialize the particles?
    // Must be "Random", "Cosmological" or "AsciiFile"
    static  std::string particle_init_type;

    // How do we want to move the particles?
    // Must be "Random" or "Gravitational"
    static  std::string particle_move_type;

    // These control random initialization
    static bool particle_initrandom_serialize;
    static long particle_initrandom_count;
    static Real particle_initrandom_mass;
    static int  particle_initrandom_iseed;
#endif

    bool             FillPatchedOldState_ok;

    // permits hydro to be turned on and off for running pure rad problems:
    static int do_hydro;

    // permits reactions to be turned on and off -- mostly for efficiency's sake
    static int do_react;

    // permits gravity calculation to be turned on and off
    static int do_grav;

    // permits sponge to be turned on and off
    static int do_sponge;

    // determines whether we use accelerators for specific loops
    static int do_acc;

    // permits rotation calculation to be turned on and off
    static int do_rotation;

    // How often to compute integral sums
    static int sum_interval;

    // Display center of mass diagnostics
    static int show_center_of_mass;

    // Display breakdown of energy sources
    static int print_energy_diagnostics;

    // Display warnings in Fortran90 routines 
    static int print_fortran_warnings;

    std::map<std::string,MultiFab*> auxDiag;
    static std::map<std::string,Array<std::string> > auxDiag_names;

    // if true, define an additional source term 
    static int add_ext_src;

    // There can be only one Gravity object, it covers all levels:
    static class Gravity *gravity;

    // There can be only one Diffusion object, it covers all levels:
    static class Diffusion *diffusion;

#ifdef RADIATION
    // permits radiation to be turned on and off without recompiling:
    static int do_radiation;

    // There can be only one Radiation object, it covers all levels:
    static class Radiation *radiation;

    friend class Radiation;
#endif

#ifdef LEVELSET
  IntMultiFab LStype;
  IntMultiFab LSnband;
  IntMultiFab LSmine;
#endif

#ifdef ROTATION
  static Real rotational_period;
  static Real rotational_dPdt;
  static int rot_source_type;
  static int rot_axis;
#endif

    // set it to true to make OMP deterministic
    static int deterministic;  

    // Is boundary function for inflow (i.e., ext_dir) thread safe?
    // Note that other boundary functions must be thread safe too.
    static int bndry_func_thread_safe;

    // Name of the probin file and its length.
    static std::string probin_file;

    static IntVect hydro_tile_size;

/* problem-specific includes */
#include <Problem.H>

};

//
// Inlines.
//

inline
int
Castro::numGrow()
{
    return NUM_GROW;
}

inline
MultiFab*
Castro::Area ()
{
    return area;
}

inline
MultiFab&
Castro::Area (int dir)
{
    return area[dir];
}

inline
MultiFab&
Castro::Volume ()
{
    return volume;
}

inline
const BoxArray&
Castro::getEdgeBoxArray (int dir) const
{
    return area[dir].boxArray();
}

inline
Castro&
Castro::getLevel (int lev)
{
    return *(Castro *) &parent->getLevel(lev);
}

inline
FluxRegister&
Castro::getFluxReg ()
{
    BL_ASSERT(flux_reg);
    return *flux_reg;
}

inline
FluxRegister&
Castro::getFluxReg (int lev)
{
    return getLevel(lev).getFluxReg();
}

#ifdef SGS
inline
FluxRegister&
Castro::getSGSFluxReg ()
{
    BL_ASSERT(sgs_flux_reg);
    return *sgs_flux_reg;
}

inline
FluxRegister&
Castro::getSGSFluxReg (int lev)
{
    return getLevel(lev).getSGSFluxReg();
}
#endif

#ifdef RADIATION
inline
FluxRegister&
Castro::getRADFluxReg ()
{
    BL_ASSERT(rad_flux_reg);
    return *rad_flux_reg;
}

inline
FluxRegister&
Castro::getRADFluxReg (int lev)
{
    return getLevel(lev).getRADFluxReg();
}
#endif

#ifdef POINTMASS
inline
Real
Castro::get_point_mass ()
{
    return point_mass;
}
#endif


#endif /*_Castro_H_*/
