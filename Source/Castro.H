
#ifndef _Castro_H_
#define _Castro_H_

#include <BC_TYPES.H>
#include <AmrLevel.H>
#include <iMultiFab.H>
#include <ErrorList.H>
#include <FluxRegister.H>
#include <BL_CXX11.H>

#ifdef BL_LAZY
#include <Lazy.H>
#endif

#ifdef PARTICLES
#include "AmrParticles.H"
#endif

#include <iostream>

using std::istream;
using std::ostream;


enum StateType { State_Type = 0,
#ifdef RADIATION
                 Rad_Type,
#endif
#ifdef SELF_GRAVITY
                 PhiGrav_Type,
                 Gravity_Type,
#endif
#ifdef ROTATION
		 PhiRot_Type,
		 Rotation_Type,
#endif
		 Source_Type,
#ifdef REACTIONS
                 Reactions_Type,
#endif
#ifdef SDC
		 SDC_Source_Type,
#ifdef REACTIONS
		 SDC_React_Type
#endif
#endif
};

// Create storage for all source terms.

#ifdef SPONGE
enum sources { sponge_src = 0,
               ext_src,
#else
enum sources { ext_src = 0,
#endif
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
	       num_src };

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
    virtual void writeSmallPlotFile (const std::string& dir,
				     ostream&       os,
				     VisMF::How     how) BL_OVERRIDE;
    void writeJobInfo (const std::string& dir);

    //
    // Define data descriptors.
    //
    static void variableSetUp ();
    //
    // Define tagging functions.
    //
    static void ErrorSetUp ();
    //
    // Cleanup data descriptors at end of run.
    //
    static void variableCleanUp ();
    //
    // Initialize grid data at problem start-up.
    //
    virtual void initData () BL_OVERRIDE;

#ifdef PARTICLES
    //
    // Read particle-related inputs
    //
    static void read_particle_params ();
    //
    // Initialize particle locations and velocities (and strengths if relevant)
    //
    void init_particles ();
    //
    // Write particles in checkpoint directories
    //
    void ParticleCheckPoint (const std::string& dir);
    //
    // How to initialize at restart
    //
    void ParticlePostRestart (const std::string& restart_file);
    //
    // Derived quantities associated with particles
    //
    MultiFab* ParticleDerive (const std::string& name,
                              Real           time,
                              int            ngrow);
    //
    // Timestamp particles
    //
    void TimestampParticles (int ngrow); 
    //
    // Advance the particles by dt
    //
    void advance_particles (int iteration, Real time, Real dt);
    //
    // Default verbosity of Particle class
    //
    static int particle_verbose;

#endif

#ifdef MAESTRO_INIT
    void MAESTRO_init ();
#endif

    //
    // Get problem metrics.
    //
    MultiFab* Area ();
    MultiFab& Area (int dir);
    //
    // The volume of the multifab.
    //
    MultiFab& Volume ();
    //
    // Set time levels of state data.
    //
    virtual void setTimeLevel (Real time,
                               Real dt_old,
                               Real dt_new) BL_OVERRIDE;
    //
    // Initialize data on this level from another Castro (during regrid).
    //
    virtual void init (AmrLevel& old) BL_OVERRIDE;
    //
    // Initialize data on this level after regridding if old level did not previously exist
    //
    virtual void init () BL_OVERRIDE;
    //
    // Proceed with next timestep?
    //
    virtual int okToContinue () BL_OVERRIDE;
    //
    // Advance grids at this level in time.
    //
    virtual Real advance (Real time,
                          Real dt,
                          int  iteration,
                          int  ncycle) BL_OVERRIDE;

    Real do_advance (Real time,
                     Real dt,
                     int  amr_iteration,
                     int  amr_ncycle,
                     int  sub_iteration,
                     int  sub_ncycle);

    Real retry_advance (Real time, Real dt, int amr_iteration, int amr_ncycle);

    void initialize_advance(Real time, Real dt, int amr_iteration, int amr_ncycle);

    void finalize_advance(Real time, Real dt, int amr_iteration, int amr_ncycle);

    void initialize_do_advance(Real time, Real dt, int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle);

    void finalize_do_advance(Real time, Real dt, int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle);

    void construct_hydro_source(Real time, Real dt);

    void check_for_nan(MultiFab& state);

#ifdef SDC
    void cons_to_prim(MultiFab& u, MultiFab& q, MultiFab& qaux);
#endif

    void construct_old_source(int src, Real time, Real dt, int amr_iteration = -1, int amr_ncycle = -1, int sub_iteration = -1, int sub_ncycle = -1);

    void construct_new_source(int src, Real time, Real dt, int amr_iteration = -1, int amr_ncycle = -1, int sub_iteration = -1, int sub_ncycle = -1);

    void sum_of_sources(MultiFab& source);

    void time_center_source_terms (MultiFab& S_new,
				   MultiFab& src_old,
				   MultiFab& src_new,
				   Real dt = 1.0);

#ifdef SPONGE
    void construct_old_sponge_source(Real time, Real dt);

    void construct_new_sponge_source(Real time, Real dt);
#endif

    void construct_old_ext_source(Real time, Real dt);

    void construct_new_ext_source(Real time, Real dt);

    void fill_ext_source(Real time, Real dt, MultiFab& S_old, MultiFab& S_new, MultiFab& ext_src, int ng);

#ifdef GRAVITY
#ifdef SELF_GRAVITY
    void construct_old_gravity(int amr_iteration, int amr_ncycle,
			       int sub_iteration, int sub_ncycle,
			       Real time);

    void construct_new_gravity(int amr_iteration, int amr_ncycle,
			       int sub_iteration, int sub_ncycle,
			       Real time);
#endif
    void construct_old_gravity_source(Real time, Real dt);

    void construct_new_gravity_source(Real time, Real dt);
#endif

    //
    // Estimate time step.
    //
    Real estTimeStep (Real dt_old);
    //
    // Compute initial time step.
    //
    Real initialTimeStep ();
    //
    // Compute initial `dt'.
    //
    virtual void computeInitialDt (int                   finest_level,
                                   int                   sub_cycle,
                                   Array<int>&           n_cycle,
                                   const Array<IntVect>& ref_ratio,
                                   Array<Real>&          dt_level,
                                   Real                  stop_time) BL_OVERRIDE;
    //
    // Compute new `dt'.
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
    // Allocate data at old time.
    //
    virtual void allocOldData () BL_OVERRIDE;
    //
    // Remove data at old time.
    //
    virtual void removeOldData () BL_OVERRIDE;
    //
    // Passes some data about the grid to a Fortran module.
    //
    void setGridInfo();
    //
    // Print information about energy budget.
    //
    void do_energy_diagnostics();
    //
    // Do work after timestep().
    //
    virtual void post_timestep (int iteration) BL_OVERRIDE;
    //
    // Contains operations to be done only after a full coarse timestep.
    //
    virtual void postCoarseTimeStep (Real cumtime) BL_OVERRIDE;
    //
    // Do work after regrid().
    //
    virtual void post_regrid (int lbase,
                              int new_finest) BL_OVERRIDE;   
    //
    // Do work after a restart().
    //
    virtual void post_restart () BL_OVERRIDE;
    //
    // Do work after init().
    //
    virtual void post_init (Real stop_time) BL_OVERRIDE;
    //
    // Do work after restart with grown_factor > 1
    //
    void post_grown_restart ();
    //
    // Error estimation for regridding.
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
#ifndef SDC
    void react_state(MultiFab& state,
		     MultiFab& reactions,
		     const iMultiFab& mask,
		     MultiFab& weights,
		     Real time,
		     Real dt_react,
		     int  ngrow = 0);

    void strang_react_first_half(Real time, Real dt);

    void strang_react_second_half(Real time, Real dt);
#else
    void react_state(Real time, Real dt);
    void get_react_source_prim(MultiFab& source, Real dt);
#endif

    void strang_chem (MultiFab&  state,
                      MultiFab&  react_mf,
                      Real       time,
                      Real       dt);
#endif

#ifdef ROTATION
    void construct_old_rotation(int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle, Real time);

    void construct_new_rotation(int amr_iteration, int amr_ncycle, int sub_iteration, int sub_ncycle, Real time);

    void construct_old_rotation_source(Real time, Real dt);

    void construct_new_rotation_source(Real time, Real dt);

    void fill_rotation_field(MultiFab& phi, MultiFab& rot, MultiFab& state, Real time);
#endif

#ifdef RADIATION
    void final_radiation_call (MultiFab& S_new, 
                               int iteration, 
                               int ncycle);
#endif

    void reset_internal_energy (MultiFab& State);

    void computeTemp (MultiFab& State);

    void update_point_mass (MultiFab& StateIn, MultiFab& StateOut, 
                            PArray<FArrayBox>& flux, MultiFab& volume,
                            Real dt);

#ifdef DIFFUSION
    void construct_old_diff_source(Real time, Real dt);
    void construct_new_diff_source(Real time, Real dt);

    void getTempDiffusionTerm (Real time, MultiFab& DiffTerm, int is_old);
    void getEnthDiffusionTerm (Real time, MultiFab& DiffTerm, int is_old);
#if (BL_SPACEDIM == 1)
    void getSpecDiffusionTerm (Real time, MultiFab& DiffTerm, int is_old);
#endif
#if (BL_SPACEDIM == 1)
    void getViscousTerm (Real time, MultiFab& ViscousTermforMomentum, MultiFab& ViscousTermforEnergy);
    void getFirstViscousTerm (Real time, MultiFab& ViscousTerm);
    void getSecndViscousTerm (Real time, MultiFab& ViscousTerm);
    void getViscousTermForEnergy (Real time, MultiFab& ViscousTermforEnergy);
#endif
    void add_temp_diffusion_to_source (MultiFab& ext_src, MultiFab& DiffTerm, Real t, int is_old);
#if (BL_SPACEDIM == 1)
    void add_spec_diffusion_to_source (MultiFab& ext_src, MultiFab& DiffTerm, Real t, int is_old);
#endif
#if (BL_SPACEDIM == 1)
    void add_viscous_term_to_source   (MultiFab& ext_src, MultiFab& ViscousTermforMomentum, 
                                       MultiFab& ViscousTermforEnergy, Real t);
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

#ifdef HYBRID_MOMENTUM
    void construct_old_hybrid_source(Real time, Real dt);

    void construct_new_hybrid_source(Real time, Real dt);

    void fill_hybrid_hydro_source(MultiFab& state, MultiFab& source);

    void hybrid_sync(MultiFab& state);
#endif

    void add_force_to_sources (MultiFab& force, MultiFab& sources, MultiFab& state);

    void apply_source_to_state (MultiFab& state, MultiFab& source, Real dt);

    void expand_state(MultiFab& S, Real time, int ng);

    void make_radial_data (int is_new);

#ifdef AUX_UPDATE
    void advance_aux(Real time, Real dt);
#endif

#ifdef PARTICLES
    static AmrTracerParticleContainer* theTracerPC () { return TracerPC; }
#endif

    static int       NUM_STATE;
    static int       Density, Xmom, Ymom, Zmom, Eden, Eint, Temp;
#ifdef HYBRID_MOMENTUM
    static int       Rmom, Lmom, Pmom;
#endif
#ifdef SHOCK_VAR
    static int       Shock;
#endif
    static int       QVAR;
    static int       QRADVAR;
    static int       NQAUX;
    static int       NQ;

    static int       FirstAdv,  NumAdv;
    static int       FirstSpec, NumSpec;
    static int       FirstAux,  NumAux;

    FluxRegister     flux_reg;

    FluxRegister& getFluxReg ();
    FluxRegister& getFluxReg (int lev);

    FluxRegister     pres_reg;

    FluxRegister& getPresReg ();
    FluxRegister& getPresReg (int lev);

    //
    // This MultiFab is on the coarser level.  This is useful for the coarser level
    // to mask out the finer level.  We only build this when it is needed.
    // This coarse MultiFab has to live on the fine level because it must be updated
    // even when only the fine level changes.
    //
    MultiFab fine_mask;
    MultiFab& build_fine_mask();

protected:

    //
    // Build a mask that ghost cells overlapping with interior cells in the same multifab
    // are set to 0, whereas others are set to 1.
    //
    PArray<iMultiFab> ib_mask;
    iMultiFab& build_interior_boundary_mask (int ng);

    int get_numpts();

    //
    // Difference between composite and level solve for phi.
    //
    MultiFab comp_minus_level_phi;
    PArray<MultiFab> comp_minus_level_grad_phi;

    //
    // A state array with ghost zones.
    //
    MultiFab Sborder;

    //
    // Source terms to the hydrodynamics solve.
    //
    MultiFab sources_for_hydro;

    //
    // Source term representing hydrodynamics update.
    //
    MultiFab hydro_source;

    //
    // Hydrodynamic (and radiation) fluxes.
    //
    PArray<MultiFab> fluxes;
    MultiFab         P_radial;
#ifdef RADIATION
    PArray<MultiFab> rad_fluxes;
#endif

    //
    // Non-hydro source terms.
    //
    PArray<MultiFab> old_sources;
    PArray<MultiFab> new_sources;

    //
    // State data to hold if we want to do a retry.
    //
    PArray<StateData> prev_state;

    //
    //  Call extern/networks/.../network.f90::network_init()
    //
    static void extern_init ();

    static void network_init ();

    static void read_params ();

    Castro& getLevel (int lev);

#ifdef RADIATION
    FluxRegister& getRADFluxReg ();
    FluxRegister& getRADFluxReg (int lev);
#endif

    void reflux ();

    void normalize_species (MultiFab& S_new);

    void enforce_consistent_e (MultiFab& S);

    Real enforce_min_density (MultiFab& S_old, MultiFab& S_new);

    Real clean_state (MultiFab& state);

    Real clean_state (MultiFab& state, MultiFab& state_old);

    void avgDown ();

    void avgDown (int state_indx);

    void buildMetrics ();

    // integrate derived quantities over domain

    Real sumDerive (const std::string& name, Real time, bool local=false);

    void sum_integrated_quantities ();

    void write_info ();

#ifdef SELF_GRAVITY
    void define_new_center (MultiFab& S, Real time);
    void write_center ();
#endif

#ifdef POINTMASS
    void pointmass_update(Real time, Real dt);
#endif

    void stopJob ();

    //
    // The data.
    //
    MultiFab             volume;
    MultiFab             area[BL_SPACEDIM];
    MultiFab             dLogArea[1];
    Array< Array<Real> > radius;
#ifdef RADIATION
    FluxRegister         rad_flux_reg;
#endif

    //
    // Static data members.
    //
#include <castro_params.H>

    static bool      signalStopJob;
    static bool      dump_old;
    static int       radius_grow;
    static int       verbose;
    static ErrorList err_list;
    static BCRec     phys_bc;
    static int       NUM_GROW;

    // for keeping track of mass changes from negative density resets
    static Real      frac_change;

    // For keeping track of fluid quantities lost at physical grid boundaries.
    // This should persist through restarts, but right now only on level 0.
    static const int n_lost = 8;
    Real             material_lost_through_boundary_cumulative[n_lost];
    Real             material_lost_through_boundary_temp[n_lost];

    // for keeping track of the amount of CPU time used -- this will persist
    // after restarts
    static Real      previousCPUTimeUsed;
    static Real      startCPUTime;

    static Real getCPUTime();

#ifdef DIFFUSION
    static int       diffuse_temp;
    static int       diffuse_enth;
    static int       diffuse_spec;
    static int       diffuse_vel;
    static Real      diffuse_cutoff_density;

    MultiFab OldTempDiffTerm;
    MultiFab OldSpecDiffTerm;
    MultiFab OldViscousTermforMomentum;
    MultiFab OldViscousTermforEnergy;

    MultiFab* NewTempDiffTerm;
    MultiFab* NewSpecDiffTerm;
    MultiFab* NewViscousTermforMomentum;
    MultiFab* NewViscousTermforEnergy;
#endif


#ifdef PARTICLES
    static int       do_tracer_particles;
#endif

    bool             FillPatchedOldState_ok;

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

#ifdef PARTICLES
    static AmrTracerParticleContainer* TracerPC;
#endif

    // Name of the probin file and its length.
    static std::string probin_file;

    static IntVect hydro_tile_size;

    static int Knapsack_Weight_Type;
    static int num_state_type;

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
Castro&
Castro::getLevel (int lev)
{
    return *(Castro *) &parent->getLevel(lev);
}

inline
void
Castro::stopJob()
{
  signalStopJob = true;
}

inline
FluxRegister&
Castro::getFluxReg ()
{
    BL_ASSERT(flux_reg.size()>0);
    return flux_reg;
}

inline
FluxRegister&
Castro::getFluxReg (int lev)
{
    return getLevel(lev).getFluxReg();
}

inline
FluxRegister&
Castro::getPresReg ()
{
    BL_ASSERT(pres_reg.size()>0);
    return pres_reg;
}

inline
FluxRegister&
Castro::getPresReg (int lev)
{
    return getLevel(lev).getPresReg();
}

#ifdef RADIATION
inline
FluxRegister&
Castro::getRADFluxReg ()
{
    BL_ASSERT(rad_flux_reg.size()>0);
    return rad_flux_reg;
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
