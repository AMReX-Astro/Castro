#include <winstd.H>

#ifndef WIN32
#include <unistd.h>
#endif

#include <iomanip>

#include <algorithm>
#include <cstdio>
#include <vector>
#include <iostream>
#include <string>
#include <ctime>

using std::cout;
using std::cerr;
using std::endl;
using std::istream;
using std::ostream;
using std::pair;
using std::string;

#include <Utility.H>
#include <CONSTANTS.H>
#include <Castro.H>
#include <Castro_F.H>
#include <Derive_F.H>
#include <VisMF.H>
#include <TagBox.H>
#include <ParmParse.H>
#include <Castro_error_F.H>

#ifdef RADIATION
#include "Radiation.H"
#endif

#ifdef PARTICLES
#include <Particles_F.H>
#endif

#ifdef GRAVITY
#include "Gravity.H"
#endif

#ifdef DIFFUSION
#include "Diffusion.H"
#endif

#ifdef LEVELSET
#include "LevelSet_F.H"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif



static Real fixed_dt     = -1.0;
static Real initial_dt   = -1.0;
static Real dt_cutoff    = 0.0;

int          Castro::checkpoint_version = 1;

bool         Castro::dump_old      = false;

int          Castro::verbose       = 0;
Real         Castro::cfl           = 0.8;
Real         Castro::burning_timestep_factor = 1.e200;
Real         Castro::init_shrink   = 1.0;
Real         Castro::change_max    = 1.1;
ErrorList    Castro::err_list;
int          Castro::radius_grow   = 1;
BCRec        Castro::phys_bc;
int          Castro::NUM_STATE     = -1;
int          Castro::do_reflux     = 1;
int          Castro::NUM_GROW      = -1;

int          Castro::Density       = -1;
int          Castro::Eden          = -1;
int          Castro::Eint          = -1;
#ifdef SGS
int          Castro::Esgs          = -1;
#endif
int          Castro::Temp          = -1;
int          Castro::Xmom          = -1;
int          Castro::Ymom          = -1;
int          Castro::Zmom          = -1;

int          Castro::NumSpec       = 0;
int          Castro::FirstSpec     = -1;

int          Castro::NumAux        = 0;
int          Castro::FirstAux      = -1;

int          Castro::NumAdv        = 0;
int          Castro::FirstAdv      = -1;

Real         Castro::difmag        = 0.1;
Real         Castro::small_dens    = -1.e200;
Real         Castro::small_temp    = -1.e200;
Real         Castro::small_pres    = -1.e200;
Real         Castro::small_ener    = -1.e200;

int          Castro::do_hydro    = -1;
int          Castro::do_react    = -1;
int          Castro::do_grav     = -1;
int          Castro::do_rotation = -1;
int          Castro::add_ext_src = 0;

int          Castro::do_acc = -1;

#ifdef POINTMASS
Real         Castro::point_mass    = 0.0;
#endif

#ifdef GRAVITY
Gravity*     Castro::gravity  = 0;
#endif

#ifdef DIFFUSION
Diffusion*    Castro::diffusion  = 0;
int           Castro::diffuse_temp = 0;
int           Castro::diffuse_spec = 0;
int           Castro::diffuse_vel = 0;
Real          Castro::diffuse_cutoff_density = -1.e200;
#endif

#ifdef RADIATION
int          Castro::do_radiation = -1;
Radiation*   Castro::radiation = 0;
#endif

#ifdef ROTATION
Real         Castro::rotational_period = -1.e200;
Real         Castro::rotational_period_dot = 0.0;
int          Castro::rot_source_type = 1;
int          Castro::rot_axis = 3;
#endif

int          Castro::deterministic = 0;

int          Castro::bndry_func_thread_safe = 1;

int          Castro::grown_factor = 1;
int          Castro::star_at_center = -1;
int          Castro::moving_center = 0;
int          Castro::normalize_species = 0;
int          Castro::fix_mass_flux = 0;
int          Castro::allow_negative_energy = 1;
int          Castro::do_special_tagging = 0;

int          Castro::ppm_type = 1;
int          Castro::ppm_reference = 1;
int          Castro::ppm_trace_sources = 0;
int          Castro::ppm_temp_fix = 0;
int          Castro::ppm_tau_in_tracing = 0;
int          Castro::ppm_predict_gammae = 0;
int          Castro::ppm_reference_edge_limit = 1;
int          Castro::ppm_reference_eigenvectors = 0;

int          Castro::source_term_predictor = 0;

int          Castro::hybrid_riemann = 0;
int          Castro::use_colglaz = 0;
int          Castro::cg_maxiter  = 12;
Real         Castro::cg_tol      = 1.0e-5;

int          Castro::hard_cfl_limit = 1;

int          Castro::use_flattening = 1;
int          Castro::ppm_flatten_before_integrals = 0;

int          Castro::transverse_use_eos = 0;
int          Castro::transverse_reset_density = 0;
int          Castro::transverse_reset_rhoe = 0;

int          Castro::dual_energy_update_E_from_e = 1;
Real         Castro::dual_energy_eta1 = 1.0e0;
Real         Castro::dual_energy_eta2 = 1.0e-4;
Real         Castro::dual_energy_eta3 = 0.0e0;

int          Castro::use_pslope  = 1;
int          Castro::grav_source_type = 2;
int          Castro::spherical_star = 0;
int          Castro::do_sponge  = 0;

#ifdef DEBUG
int          Castro::print_fortran_warnings  = 1;
int          Castro::print_energy_diagnostics  = 1;
#else
int          Castro::print_fortran_warnings  = 0;
int          Castro::print_energy_diagnostics  = 0;
#endif

int          Castro::sum_interval = -1;
int          Castro::show_center_of_mass = 0;

#ifdef SGS
Real         Castro::sum_turb_src = 0.0;
#endif

std::string  Castro::job_name = "";

std::string  Castro::probin_file = "probin";

#if BL_SPACEDIM == 1
IntVect      Castro::hydro_tile_size(1024);
#elif BL_SPACEDIM == 2
IntVect      Castro::hydro_tile_size(1024,16);
#else
IntVect      Castro::hydro_tile_size(1024,16,16);
#endif

// this will be reset upon restart
Real         Castro::previousCPUTimeUsed = 0.0;

Real         Castro::startCPUTime = 0.0;

Real         Castro::max_delta_e = 0.0;

// Note: Castro::variableSetUp is in Castro_setup.cpp

void
Castro::variableCleanUp () 
{
#ifdef GRAVITY
  if (gravity != 0) {
    if (verbose > 1 && ParallelDescriptor::IOProcessor()) {
      cout << "Deleting gravity in variableCleanUp..." << '\n';
    }
    delete gravity;
    gravity = 0;
  }
#endif

#ifdef DIFFUSION
  if (diffusion != 0) {
    if (verbose > 1 && ParallelDescriptor::IOProcessor()) {
      cout << "Deleting diffusion in variableCleanUp..." << '\n';
    }
    delete diffusion;
    diffusion = 0;
  }
#endif

#ifdef RADIATION
  if (radiation != 0) { int report = (verbose || radiation->verbose);
    if (report && ParallelDescriptor::IOProcessor()) {
      cout << "Deleting radiation in variableCleanUp..." << '\n';
    }
    delete radiation;
    radiation = 0;
    if (report && ParallelDescriptor::IOProcessor()) {
      cout << "                                        done" << endl;
    }
  }
#endif

    desc_lst.clear();
}

void
Castro::read_params ()
{
    static bool done = false;

    if (done) return;

    done = true;

    ParmParse pp("castro");   

    pp.query("v",verbose);
//  verbose = (verbose ? 1 : 0);
    pp.query("init_shrink",init_shrink);
    pp.query("cfl",cfl);
    pp.query("burning_timestep_factor",burning_timestep_factor);
    pp.query("change_max",change_max);
    pp.query("fixed_dt",fixed_dt);
    pp.query("initial_dt",initial_dt);
    pp.query("sum_interval",sum_interval);
    pp.query("do_reflux",do_reflux);
    do_reflux = (do_reflux ? 1 : 0);
    pp.query("dt_cutoff",dt_cutoff);

    pp.query("dump_old",dump_old);

    pp.query("difmag",difmag);
    pp.query("small_dens",small_dens);
    pp.query("small_temp",small_temp);
    pp.query("small_pres",small_pres);
    pp.query("small_ener",small_ener);

#ifdef POINTMASS
    pp.get("point_mass",point_mass);
#endif


    // Get boundary conditions
    Array<int> lo_bc(BL_SPACEDIM), hi_bc(BL_SPACEDIM);
    pp.getarr("lo_bc",lo_bc,0,BL_SPACEDIM);
    pp.getarr("hi_bc",hi_bc,0,BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        phys_bc.setLo(i,lo_bc[i]);
        phys_bc.setHi(i,hi_bc[i]);
    }

    //
    // Check phys_bc against possible periodic geometry
    // if periodic, must have internal BC marked.
    //
    if (Geometry::isAnyPeriodic())
    {
        //
        // Do idiot check.  Periodic means interior in those directions.
        //
        for (int dir = 0; dir<BL_SPACEDIM; dir++)
        {
            if (Geometry::isPeriodic(dir))
            {
                if (lo_bc[dir] != Interior)
                {
                    std::cerr << "Castro::read_params:periodic in direction "
                              << dir
                              << " but low BC is not Interior\n";
                    BoxLib::Error();
                }
                if (hi_bc[dir] != Interior)
                {
                    std::cerr << "Castro::read_params:periodic in direction "
                              << dir
                              << " but high BC is not Interior\n";
                    BoxLib::Error();
                }
            }
        }
    }
    else
    {
        //
        // Do idiot check.  If not periodic, should be no interior.
        //
        for (int dir=0; dir<BL_SPACEDIM; dir++)
        {
            if (lo_bc[dir] == Interior)
            {
                std::cerr << "Castro::read_params:interior bc in direction "
                          << dir
                          << " but not periodic\n";
                BoxLib::Error();
            }
            if (hi_bc[dir] == Interior)
            {
                std::cerr << "Castro::read_params:interior bc in direction "
                          << dir
                          << " but not periodic\n";
                BoxLib::Error();
            }
        }
    }

    if ( Geometry::IsRZ() && (lo_bc[0] != Symmetry) ) {
        std::cerr << "ERROR:Castro::read_params: must set r=0 boundary condition to Symmetry for r-z\n";
        BoxLib::Error();
    }

#if (BL_SPACEDIM == 1)
    if ( Geometry::IsSPHERICAL() )
    {
      if ( (lo_bc[0] != Symmetry) && (Geometry::ProbLo(0) == 0.0) ) 
      {
        std::cerr << "ERROR:Castro::read_params: must set r=0 boundary condition to Symmetry for spherical\n";
        BoxLib::Error();
      }
    }
#elif (BL_SPACEDIM == 2)
    if ( Geometry::IsSPHERICAL() )
      {
	BoxLib::Abort("We don't support spherical coordinate systems in 2D");
      }
#elif (BL_SPACEDIM == 3)
    if ( Geometry::IsRZ() )
      {
	BoxLib::Abort("We don't support cylindrical coordinate systems in 3D"); 
      }
    else if ( Geometry::IsSPHERICAL() )
      {
	BoxLib::Abort("We don't support spherical coordinate systems in 3D");
      }
#endif


    pp.get("do_hydro",do_hydro);
    pp.query("add_ext_src",add_ext_src);

    pp.query("do_acc",do_acc);

#ifdef DIFFUSION
    pp.query("diffuse_temp",diffuse_temp);
    pp.query("diffuse_spec",diffuse_spec);
    pp.query("diffuse_vel",diffuse_vel);
    pp.query("diffuse_cutoff_density",diffuse_cutoff_density);
#endif

    pp.query("grown_factor",grown_factor);
    if (grown_factor < 1) 
       BoxLib::Error("grown_factor must be integer >= 1");

    pp.query("star_at_center",star_at_center);

    pp.query("moving_center",moving_center);

    pp.query("normalize_species",normalize_species);
    pp.query("fix_mass_flux",fix_mass_flux);
    pp.query("allow_negative_energy",allow_negative_energy);
    pp.query("do_special_tagging",do_special_tagging);

    pp.query("ppm_type", ppm_type);
    pp.query("ppm_reference", ppm_reference);
    pp.query("ppm_trace_sources", ppm_trace_sources);
    pp.query("ppm_temp_fix", ppm_temp_fix);
    pp.query("ppm_tau_in_tracing", ppm_tau_in_tracing);
    pp.query("ppm_predict_gammae", ppm_predict_gammae);
    pp.query("ppm_reference_edge_limit", ppm_reference_edge_limit);
    pp.query("ppm_flatten_before_integrals", ppm_flatten_before_integrals);
    pp.query("ppm_reference_eigenvectors", ppm_reference_eigenvectors);
    pp.query("source_term_predictor", source_term_predictor);
    pp.query("hybrid_riemann",hybrid_riemann);
    pp.query("use_colglaz",use_colglaz);
    pp.query("use_flattening",use_flattening);
    pp.query("transverse_use_eos",transverse_use_eos);
    pp.query("transverse_reset_density",transverse_reset_density);
    pp.query("transverse_reset_rhoe",transverse_reset_rhoe);

    pp.query("cg_maxiter",cg_maxiter);
    pp.query("cg_tol",cg_tol);
    pp.query("use_pslope",use_pslope);
    pp.query("grav_source_type",grav_source_type);
    pp.query("spherical_star",spherical_star);
    pp.query("do_sponge",do_sponge);
    pp.query("hard_cfl_limit",hard_cfl_limit);

    pp.query("dual_energy_update_E_from_e",dual_energy_update_E_from_e);
    pp.query("dual_energy_eta1",dual_energy_eta1);
    pp.query("dual_energy_eta2",dual_energy_eta2);
    pp.query("dual_energy_eta3",dual_energy_eta3);

    pp.query("show_center_of_mass",show_center_of_mass);
    pp.query("print_energy_diagnostics",print_energy_diagnostics);
    pp.query("print_fortran_warnings",print_fortran_warnings);


    if (ppm_reference > 1 || ppm_reference < 0)
      {
        std::cerr << "invalid ppm_reference\n";
        BoxLib::Error();
      }	

    // for the moment, ppm_type = 0 does not support ppm_trace_sources --
    // we need to add the momentum sources to the states (and not
    // add it in trans_3d
    if (ppm_type == 0 && ppm_trace_sources == 1)
      {
        std::cerr << "ppm_trace_sources = 1 not implemented for ppm_type = 0 \n";
        BoxLib::Error();
      }


    // ppm_flatten_before_integrals is only done for ppm_type != 0
    if (ppm_type == 0 && ppm_flatten_before_integrals > 0)
      {
        std::cerr << "ppm_flatten_before_integrals > 0 not implemented for ppm_type != 0 \n";
        BoxLib::Error();
      }
	

    if (ppm_temp_fix > 0 && BL_SPACEDIM == 1)
      {
        std::cerr << "ppm_temp_fix > 0 not implemented in 1-d \n";
        BoxLib::Error();
      }

    if (ppm_tau_in_tracing == 1 && BL_SPACEDIM == 1)
      {
        std::cerr << "ppm_tau_in_tracing == 1 not implemented in 1-d \n";
        BoxLib::Error();
      }

    if (ppm_predict_gammae == 1 && ppm_tau_in_tracing != 1)
      {
	std::cerr << "ppm_predict_gammae == 1 needs ppm_tau_in_tracing == 1\n";
	BoxLib::Error();
      }

    if (hybrid_riemann == 1 && BL_SPACEDIM == 1)
      {
        std::cerr << "hybrid_riemann only implemented in 2- and 3-d\n";
        BoxLib::Error();
      }


    // Make sure not to call refluxing if we're not actually doing any hydro.
    if (do_hydro == 0) do_reflux = 0;

#ifdef GRAVITY
    pp.get("do_grav",do_grav);

#if (BL_SPACEDIM == 1)
    if (do_grav && !Geometry::IsSPHERICAL()) {
        std::cerr << "ERROR:Castro::Gravity in 1D assumes that the coordinate system is spherical\n";
        BoxLib::Error();
    }
#endif
#endif

#ifdef REACTIONS
    pp.get("do_react",do_react);
#else
    pp.query("do_react",do_react);
#endif

#ifdef PARTICLES
    read_particle_params();
#endif

#ifdef RADIATION
    pp.get("do_radiation",do_radiation);

    // Some radiation parameters are initialized here because they
    // may be used in variableSetUp, well before the call to the
    // Radiation constructor, 

    if (do_radiation) {
      Radiation::read_static_params();
    }
#endif

#ifdef ROTATION
    pp.get("do_rotation",do_rotation);
    if (do_rotation) {
      pp.get("rotational_period",rotational_period);
      if (rotational_period <= 0.0) {
	std::cerr << "Error:Castro::Rotation enabled but rotation period less than zero\n";
	BoxLib::Error();
      }
    }
    else pp.query("rotational_period",rotational_period);
    pp.query("rotational_period_dot",rotational_period_dot);
    pp.query("rot_source_type",rot_source_type);
    if (Geometry::IsRZ())
      rot_axis = 2;
    pp.query("rot_axis",rot_axis);
#if (BL_SPACEDIM == 1)
      if (do_rotation) {
	std::cerr << "ERROR:Castro::Rotation not implemented in 1d\n";
	BoxLib::Error();
      }
#endif
#endif

   pp.query("deterministic", deterministic);

   pp.query("bndry_func_thread_safe", bndry_func_thread_safe);
   StateDescriptor::setBndryFuncThreadSafety(bndry_func_thread_safe);

   pp.query("job_name",job_name);  

   ParmParse ppa("amr");
   ppa.query("probin_file",probin_file);

    Array<int> tilesize(BL_SPACEDIM);
    if (pp.queryarr("hydro_tile_size", tilesize, 0, BL_SPACEDIM))
    {
	for (int i=0; i<BL_SPACEDIM; i++) hydro_tile_size[i] = tilesize[i];
    }
}

Castro::Castro ()
{
    flux_reg = 0;
#ifdef SGS
    sgs_flux_reg = 0;
#endif
#ifdef RADIATION
    rad_flux_reg = 0;
#endif
    fine_mask = 0;
}

Castro::Castro (Amr&            papa,
                int             lev,
                const Geometry& level_geom,
                const BoxArray& bl,
                Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,time) 
{
    buildMetrics();

    flux_reg = 0;
    if (level > 0 && do_reflux)
    {
        flux_reg = new FluxRegister(grids,crse_ratio,level,NUM_STATE);
        flux_reg->setVal(0.0);
    }

#ifdef SGS
    sgs_flux_reg = 0;
    if (level > 0 && do_reflux)
    {
        sgs_flux_reg = new FluxRegister(grids,crse_ratio,level,NUM_STATE);
        sgs_flux_reg->setVal(0.0);
    }
#endif

#ifdef RADIATION    
    rad_flux_reg = 0;
    if (Radiation::rad_hydro_combined && level > 0 && do_reflux) 
    {
      rad_flux_reg = new FluxRegister(grids,crse_ratio,level,Radiation::nGroups);
      rad_flux_reg->setVal(0.0);
    }
#endif

    fine_mask = 0;

#ifdef GRAVITY

   // Initialize to zero here in case we run with do_grav = false.
   MultiFab& new_grav_mf = get_new_data(Gravity_Type);
   new_grav_mf.setVal(0.0);
       
   if (do_grav) {
      // gravity is a static object, only alloc if not already there
      if (gravity == 0) 
	gravity = new Gravity(parent,parent->finestLevel(),&phys_bc,Density);

#if (BL_SPACEDIM > 1)
      // Passing numpts_1d at level 0 
      if (!Geometry::isAllPeriodic() && gravity != 0)
      {
         int numpts_1d = get_numpts();
         gravity->set_numpts_in_gravity(numpts_1d);
      }
#endif

      gravity->install_level(level,this,volume,area);

      if (verbose && level == 0 &&  ParallelDescriptor::IOProcessor()) 
         std::cout << "Setting the gravity type to " << gravity->get_gravity_type() << std::endl;

   } else {
       MultiFab& phi_new = get_new_data(PhiGrav_Type);
       phi_new.setVal(0.0);
   }
#endif

#ifdef ROTATION

   // Initialize rotation data to zero.

   MultiFab& phirot_new = get_new_data(PhiRot_Type);
   phirot_new.setVal(0.0);

   MultiFab& rot_new = get_new_data(Rotation_Type);
   rot_new.setVal(0.0);

#endif

   // Initialize source term data to zero.

   MultiFab& dSdt_new = get_new_data(Source_Type);
   dSdt_new.setVal(0.0);
   
#ifdef REACTIONS

   // Initialize reaction data to zero.

   MultiFab& reactions_new = get_new_data(Reactions_Type);
   reactions_new.setVal(0.0);

#endif
   
#ifdef DIFFUSION
      // diffusion is a static object, only alloc if not already there
      if (diffusion == 0) 
	diffusion = new Diffusion(parent,&phys_bc);

      diffusion->install_level(level,this,volume,area);
#endif

#ifdef RADIATION
    if (do_radiation) {
      if (radiation == 0) {
	// radiation is a static object, only alloc if not already there
	radiation = new Radiation(parent, this);
      }
      radiation->regrid(level, grids);
    }
#endif

#ifdef LEVELSET
    // Build level set narrowband helpers
    LStype.define(bl,1,1,Fab_allocate);
    LSnband.define(bl,BL_SPACEDIM,1,Fab_allocate);
    LSmine.define(bl,BL_SPACEDIM,1,Fab_allocate);
#endif

#ifdef SGS
   MultiFab& new_sgs_mf = get_new_data(SGS_Type);
   new_sgs_mf.setVal(0.0);
#endif
}

Castro::~Castro () 
{
    delete flux_reg;
#ifdef SGS
    delete sgs_flux_reg;
#endif

#ifdef RADIATION
    if (Radiation::rad_hydro_combined) {
      delete rad_flux_reg;
    }
    if (radiation != 0) {
      //radiation->cleanup(level);
      radiation->close(level);
    }
#endif
    delete fine_mask;
}


void
Castro::buildMetrics ()
{
    const int ngrd = grids.size();

    radius.resize(ngrd);

    const Real* dx = geom.CellSize();

    for (int i = 0; i < ngrd; i++)
    {
        const Box& b = grids[i];
        int ilo      = b.smallEnd(0)-radius_grow;
        int ihi      = b.bigEnd(0)+radius_grow;
        int len      = ihi - ilo + 1;

        radius[i].resize(len);

        Real* rad = radius[i].dataPtr();

        if (Geometry::IsCartesian())
        {
            for (int j = 0; j < len; j++)
            {
                rad[j] = 1.0;
            }
        }
        else
        {
            RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());

            const Real xlo = gridloc.lo(0) + (0.5 - radius_grow)*dx[0];

            for (int j = 0; j < len; j++)
            {
                rad[j] = xlo + j*dx[0];
            }
        }
    }
    //
    // Build volume, face area and dLogArea arrays.
    // volume is not PArrayManaged, must manually delete.
    //
    volume.clear();
    //
    // area is not PArrayManaged, must manually delete.
    //
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        area[dir].clear();
    }
    dLogArea[0].clear();
    geom.GetVolume(volume,grids,NUM_GROW);

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        geom.GetFaceArea(area[dir],grids,dir,NUM_GROW);
    }
#if (BL_SPACEDIM <= 2)
    geom.GetDLogA(dLogArea[0],grids,0,NUM_GROW);
#endif
}

void
Castro::setTimeLevel (Real time,
                      Real dt_old,
                      Real dt_new)
{
    AmrLevel::setTimeLevel(time,dt_old,dt_new);
}

void
Castro::initData ()
{
    BL_PROFILE("Castro::initData()");

    //
    // Loop over grids, call FORTRAN function to init with data.
    //
    int ns          = NUM_STATE;
    const Real* dx  = geom.CellSize();
    MultiFab& S_new = get_new_data(State_Type);
    Real cur_time   = state[State_Type].curTime();

    S_new.setVal(0.);

    // make sure dx = dy = dz -- that's all we guarantee to support
#if (BL_SPACEDIM == 2)
    const Real SMALL = 1.e-13;
    if (fabs(dx[0] - dx[1]) > SMALL*dx[0])
      {
	BoxLib::Abort("We don't support dx != dy");
      }
#elif (BL_SPACEDIM == 3)
    const Real SMALL = 1.e-13;
    if ( (fabs(dx[0] - dx[1]) > SMALL*dx[0]) || (fabs(dx[0] - dx[2]) > SMALL*dx[0]) )
      {
	BoxLib::Abort("We don't support dx != dy != dz");
      }
#endif


    if (verbose && ParallelDescriptor::IOProcessor())
       std::cout << "Initializing the data at level " << level << std::endl;

#ifdef RADIATION
    // rad quantities are in the state even if (do_radiation == 0)
    MultiFab &Rad_new = get_new_data(Rad_Type);
    Rad_new.setVal(0.);
#endif

#ifdef REACTIONS
    MultiFab &React_new = get_new_data(Reactions_Type);
    React_new.setVal(0.);
#endif

#ifdef MAESTRO_INIT
    MAESTRO_init();
#else
    {
       for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
       {
	  RealBox gridloc = RealBox(grids[mfi.index()],geom.CellSize(),geom.ProbLo());	 
          const Box& box     = mfi.validbox();
          const int* lo      = box.loVect();
          const int* hi      = box.hiVect();
  
#ifdef DIMENSION_AGNOSTIC
          BL_FORT_PROC_CALL(CA_INITDATA,ca_initdata)
          (level, cur_time, ARLIM_3D(lo), ARLIM_3D(hi), ns,
  	   BL_TO_FORTRAN_3D(S_new[mfi]), ZFILL(dx),
  	   ZFILL(gridloc.lo()), ZFILL(gridloc.hi()));
#else
          BL_FORT_PROC_CALL(CA_INITDATA,ca_initdata)
  	  (level, cur_time, lo, hi, ns,
  	   BL_TO_FORTRAN(S_new[mfi]), dx,
  	   gridloc.lo(), gridloc.hi());
#endif

          // Verify that the sum of (rho X)_i = rho at every cell
          BL_FORT_PROC_CALL(CA_CHECK_INITIAL_SPECIES, ca_check_initial_species)
              (lo, hi, BL_TO_FORTRAN(S_new[mfi]));
       }
       enforce_consistent_e(S_new);
    }

#ifdef RADIATION
    if (do_radiation) {
      for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
          int i = mfi.index();

	  if (radiation->verbose > 2) {
	    cout << "Calling RADINIT at level " << level << ", grid "
		 << i << endl;
	  }

          RealBox    gridloc(grids[mfi.index()],
                             geom.CellSize(), geom.ProbLo());
          const Box& box = mfi.validbox();
          const int* lo  = box.loVect();
          const int* hi  = box.hiVect();

	  Rad_new[mfi].setVal(0.0);

	  BL_FORT_PROC_CALL(CA_INITRAD,ca_initrad)
	      (level, cur_time, lo, hi, Radiation::nGroups,
	       BL_TO_FORTRAN(Rad_new[mfi]),dx,
	       gridloc.lo(),gridloc.hi());

	  if (Radiation::nNeutrinoSpecies > 0 && Radiation::nNeutrinoGroups[0] == 0) {
	      // Hack: running photon radiation through neutrino solver
            Rad_new[mfi].mult(Radiation::Etorad, 
                            0, Radiation::nGroups);
	  }

          if (Rad_new.nComp() > Radiation::nGroups) {
            // Initialize flux components to 0
            Rad_new[mfi].setVal(0.0, box, Radiation::nGroups,
                              Rad_new.nComp() - Radiation::nGroups);
          }
      }
    }
#endif // RADIATION

#endif // MAESTRO_INIT

    set_special_tagging_flag(cur_time);

#if (BL_SPACEDIM > 1)
    if ( (level == 0) && (spherical_star == 1) ) {
       int nc = S_new.nComp();
       int n1d = get_numpts();
       BL_FORT_PROC_CALL(ALLOCATE_OUTFLOW_DATA,allocate_outflow_data)(&n1d,&nc);
       int is_new = 1;
       make_radial_data(is_new);
    }
#endif

#ifdef GRAVITY
    MultiFab& G_new = get_new_data(Gravity_Type);
    G_new.setVal(0.);

    MultiFab& phi_new = get_new_data(PhiGrav_Type);
    phi_new.setVal(0.);
#endif

    MultiFab& dSdt_new = get_new_data(Source_Type);
    dSdt_new.setVal(0.);

#ifdef ROTATION
    MultiFab& rot_new = get_new_data(Rotation_Type);
    rot_new.setVal(0.);
    
    MultiFab& phirot_new = get_new_data(PhiRot_Type);
    phirot_new.setVal(0.);
#endif

#ifdef LEVELSET
    MultiFab& LS_new = get_new_data(LS_State_Type);
    LS_new.setVal(0.);

    for (MFIter mfi(LS_new); mfi.isValid(); ++mfi) 
      {        
        RealBox    gridloc = RealBox(grids[mfi.index()],geom.CellSize(),geom.ProbLo());
	IntFab& type       = LStype[mfi];
	const Box& box     = mfi.validbox();
        const int* lo      = box.loVect();
        const int* hi      = box.hiVect();

        BL_FORT_PROC_CALL(CA_INITPHI,ca_initphi)
	  (level, cur_time, lo, hi, 1,
	   BL_TO_FORTRAN(LS_new[mfi]),
	   BL_TO_FORTRAN(type),
	   dx,gridloc.lo(),gridloc.hi());
      }
#endif

#ifdef PARTICLES
    if (level == 0)
       init_particles();

    // Must redistribute particles before calling init_santa_barbara so that the particles already
    //  live on the higher level when we go to put some of the mass onto the grid.

    if (level > 0) 
       ParticleRedistribute();
#endif

    if (verbose && ParallelDescriptor::IOProcessor())
       std::cout << "Done initializing the level " << level << " data " << std::endl;
}

void
Castro::init (AmrLevel &old)
{
    BL_PROFILE("Castro::init(old)");

    Castro* oldlev = (Castro*) &old;
    //
    // Create new grid data by fillpatching from old.
    //
    Real dt_new    = parent->dtLevel(level);
    Real cur_time  = oldlev->state[State_Type].curTime();
    Real prev_time = oldlev->state[State_Type].prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);

    MultiFab& S_new = get_new_data(State_Type);
    FillPatch(old,S_new,0,cur_time,State_Type,0,NUM_STATE);

    // Set E in terms of e + kinetic energy
    // enforce_consistent_e(S_new);

#ifdef RADIATION
    if (do_radiation) {
      MultiFab& Er_new = get_new_data(Rad_Type);
      int ncomp = Er_new.nComp();

      FillPatch(old,Er_new,0,cur_time,Rad_Type,0,ncomp);
    }
#endif

#ifdef GRAVITY
    if (do_grav) {
	MultiFab& phi_new = get_new_data(PhiGrav_Type);
	FillPatch(old,phi_new,0,cur_time,PhiGrav_Type,0,1);
    }
#endif

#ifdef ROTATION
    if (do_rotation) {
	MultiFab& phirot_new = get_new_data(PhiRot_Type);
	FillPatch(old,phirot_new,0,cur_time,PhiRot_Type,0,1);
    }
#endif

    MultiFab& dSdt_new = get_new_data(Source_Type);
    FillPatch(old,dSdt_new,0,cur_time,Source_Type,0,NUM_STATE);
    
#ifdef REACTIONS
    {
	MultiFab& React_new = get_new_data(Reactions_Type);
	int ncomp = React_new.nComp();

        FillPatch(old,React_new,0,cur_time,Reactions_Type,0,ncomp);
    }

#endif

#ifdef LEVELSET
    MultiFab& LS_new = get_new_data(LS_State_Type);
    int nGrowRegrid = 0;
    
    FillPatch(old,LS_new,nGrowRegrid,cur_time,LS_State_Type,0,1);
    
    // FIXME: Assumes that interpolated coarse data should rather just be setvald
    LStype.setVal(3); // This means we don't care about these points
    LStype.copy(oldlev->LStype);
    
    // Reinitialize to allow narrowband to push into new cell area
    // ...need to build narrowband structure prior to using it though
    for (MFIter mfi(LStype); mfi.isValid(); ++mfi)
      {
	IntFab& type = LStype[mfi];
	IntFab& nband = LSnband[mfi];
	IntFab& mine = LSmine[mfi];
	const Box& box = mfi.validbox();
	int nbandsize = nband.box().numPts();
	int minesize = mine.box().numPts();
	
	// Set nband data based on type
	BL_FORT_PROC_CALL(LS_NARROWBAND,ls_narrowband)
	     (BL_TO_FORTRAN(type),
	      nband.dataPtr(), &nbandsize,
	      mine.dataPtr(), &minesize,
	      box.loVect(), box.hiVect());
      }
    reinit_phi(cur_time);
#endif
}

//
// This version inits the data on a new level that did not
// exist before regridding.
//
void
Castro::init ()
{
    BL_PROFILE("Castro::init()");

    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[State_Type].curTime();
    Real prev_time = getLevel(level-1).state[State_Type].prevTime();

    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time,dt_old,dt);
    MultiFab& S_new = get_new_data(State_Type);
    FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, NUM_STATE);

    // Set E in terms of e + kinetic energy
    // enforce_consistent_e(S_new);

#ifdef RADIATION
    if (do_radiation) {
      MultiFab& Er_new = get_new_data(Rad_Type);
      int ncomp = Er_new.nComp();
      FillCoarsePatch(Er_new, 0, cur_time, Rad_Type, 0, ncomp);
    }
#endif

#ifdef GRAVITY
    if (do_grav) {
	MultiFab& phi_new = get_new_data(PhiGrav_Type);
	FillCoarsePatch(phi_new, 0, cur_time, PhiGrav_Type, 0, 1);
    }
#endif

    MultiFab& dSdt_new = get_new_data(Source_Type);
    FillCoarsePatch(dSdt_new, 0, cur_time, Source_Type, 0, NUM_STATE);
    
#ifdef ROTATION
    if (do_rotation) {
      MultiFab& phirot_new = get_new_data(PhiRot_Type);
      FillCoarsePatch(phirot_new, 0, cur_time, PhiRot_Type, 0, 1);
    }
#endif

#ifdef LEVELSET
    FillCoarsePatch(get_new_data(LS_State_Type),0,cur_time,LS_State_Type,0,1);
#endif

}

Real
Castro::initialTimeStep ()
{
    Real dummy_dt = 0.0;
    Real init_dt  = 0.0;

    if (initial_dt > 0.0) 
    {
       init_dt = initial_dt;
    } 
    else 
    {
       init_dt = init_shrink*estTimeStep(dummy_dt);
    }

    return init_dt;
}

Real
Castro::estTimeStep (Real dt_old)
{
    BL_PROFILE("Castro::estTimeStep()");

    if (fixed_dt > 0.0)
        return fixed_dt;

    // This is just a dummy value to start with 
    Real estdt  = 1.0e+200;

    const MultiFab& stateMF = get_new_data(State_Type);

    const Real* dx = geom.CellSize();    
    
#ifdef DIFFUSION
    if (do_hydro or diffuse_temp) 
#else
    if (do_hydro)
#endif
    {

#ifdef RADIATION
      if (Radiation::rad_hydro_combined) {

	  // Compute radiation + hydro limited timestep.
	
#ifdef _OPENMP
#pragma omp parallel
#endif
	  {
	      Real dt = 1.e200;

	      const MultiFab& radMF = get_new_data(Rad_Type);
	      FArrayBox gPr;

	      for (MFIter mfi(stateMF, true); mfi.isValid(); ++mfi) 
	      {
	          const Box& tbox = mfi.tilebox();
	          const Box& vbox = mfi.validbox();
		  
		  gPr.resize(tbox);
		  radiation->estimate_gamrPr(stateMF[mfi], radMF[mfi], gPr, dx, vbox);
	  
		  BL_FORT_PROC_CALL(CA_ESTDT_RAD, ca_estdt_rad)
		      (BL_TO_FORTRAN(stateMF[mfi]),
	               BL_TO_FORTRAN(gPr),
	               tbox.loVect(),tbox.hiVect(),dx,&dt);
              }
#ifdef _OPENMP
#pragma omp critical (castro_estdt_rad)	      
#endif
	      {
	          estdt = std::min(estdt,dt);
              }
          }
      }
      else 
      {
#endif   

	  // Compute hydro-limited timestep.
	if (do_hydro)
	  {

#ifdef _OPENMP
#pragma omp parallel
#endif
	    {
	      Real dt = 1.e200;
	      
	      for (MFIter mfi(stateMF,true); mfi.isValid(); ++mfi)
		{
		  const Box& box = mfi.tilebox();
		  
		  BL_FORT_PROC_CALL(CA_ESTDT,ca_estdt)
		      (ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
		       BL_TO_FORTRAN_3D(stateMF[mfi]),
		       ZFILL(dx),&dt);
		}
#ifdef _OPENMP
#pragma omp critical (castro_estdt)	      
#endif
	      {
		estdt = std::min(estdt,dt);
	      }
	    }
	  }
	    
#ifdef DIFFUSION
	// Diffusion-limited timestep
	if (diffuse_temp)
	  {
#ifdef _OPENMP
#pragma omp parallel
#endif
	    {
	      Real dt = 1.e200;
	      
	      for (MFIter mfi(stateMF,true); mfi.isValid(); ++mfi)
		{
		  const Box& box = mfi.tilebox();

		  BL_FORT_PROC_CALL(CA_ESTDT_DIFFUSION,ca_estdt_diffusion)
		      (ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
		       BL_TO_FORTRAN_3D(stateMF[mfi]),
		       ZFILL(dx),&dt);
		}
#ifdef _OPENMP
#pragma omp critical (castro_estdt)	      
#endif
	      {
		estdt = std::min(estdt,dt);
	      }
	    }
	  }
#endif  // diffusion

#ifdef RADIATION
      }
#endif

       ParallelDescriptor::ReduceRealMin(estdt);
       estdt *= cfl;
       if (verbose && ParallelDescriptor::IOProcessor()) 
           std::cout << "...estimated hydro-limited timestep at level " << level << ": " << estdt << std::endl;
    }

#ifdef REACTIONS
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& reactions_new = get_new_data(Reactions_Type);

    if (do_react) {
    
        // Compute burning-limited timestep.
    
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            Real dt = 1.e200;
    
	    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
	    {
	        const Box& box = mfi.validbox();
		BL_FORT_PROC_CALL(CA_ESTDT_BURNING,ca_estdt_burning)
                    (BL_TO_FORTRAN_3D(S_new[mfi]),
		     BL_TO_FORTRAN_3D(reactions_new[mfi]),
		     ARLIM_3D(box.loVect()),ARLIM_3D(box.hiVect()),ZFILL(dx),&dt);

	    }
#ifdef _OPENMP
#pragma omp critical (castro_estdt_burning)
#endif
	    {
	        estdt = std::min(estdt,dt);
	    }
	      
        }
    
	ParallelDescriptor::ReduceRealMin(estdt);

	if (verbose && ParallelDescriptor::IOProcessor()) 
	  std::cout << "...estimated burning-limited and hydro-limited timestep at level " << level << ": " << estdt << std::endl;
    }
#endif

#ifdef RADIATION
    if (do_radiation) radiation->EstTimeStep(estdt, level);
#endif

#ifdef PARTICLES
#ifdef GRAVITY
    ParticleEstTimeStep(estdt);
#endif
#endif

    if (verbose && ParallelDescriptor::IOProcessor())
        cout << "Castro::estTimeStep at level " << level << ":  estdt = " << estdt << '\n';

    return estdt;
}

void
Castro::computeNewDt (int                   finest_level,
                      int                   sub_cycle,
                      Array<int>&           n_cycle,
                      const Array<IntVect>& ref_ratio,
                      Array<Real>&          dt_min,
                      Array<Real>&          dt_level,
                      Real                  stop_time,
                      int                   post_regrid_flag)
{
    BL_PROFILE("Castro::computeNewDt()");

    //
    // We are at the start of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0)
        return;

    int i;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        Castro& adv_level = getLevel(i);
        dt_min[i] = adv_level.estTimeStep(dt_level[i]);
    }

    if (fixed_dt <= 0.0)
    {
       if (post_regrid_flag == 1) 
       {
          //
          // Limit dt's by pre-regrid dt
          //
          for (i = 0; i <= finest_level; i++)
          {
              dt_min[i] = std::min(dt_min[i],dt_level[i]);
          }
       } 
       else 
       {
          //
          // Limit dt's by change_max * old dt
          //
          for (i = 0; i <= finest_level; i++)
          {
             if (verbose && ParallelDescriptor::IOProcessor())
                 if (dt_min[i] > change_max*dt_level[i])
                 {
                        cout << "Castro::compute_new_dt : limiting dt at level "
                             << i << '\n';
                        cout << " ... new dt computed: " << dt_min[i]
                             << '\n';
                        cout << " ... but limiting to: "
                             << change_max * dt_level[i] << " = " << change_max
                             << " * " << dt_level[i] << '\n';
                 }
              dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
          }
       } 
    }

    //
    // Find the minimum over all levels
    //
    for (i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_min[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

void
Castro::computeInitialDt (int                   finest_level,
                          int                   sub_cycle,
                          Array<int>&           n_cycle,
                          const Array<IntVect>& ref_ratio,
                          Array<Real>&          dt_level,
                          Real                  stop_time)
{
    BL_PROFILE("Castro::computeInitialDt()");

    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
        return;

    int i;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    ///TODO/DEBUG: This will need to change for optimal subcycling.
    for (i = 0; i <= finest_level; i++)
    {
        dt_level[i] = getLevel(i).initialTimeStep();
        n_factor   *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_level[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

void
Castro::post_timestep (int iteration)
{
    BL_PROFILE("Castro::post_timestep()");

    //
    // Integration cycle on fine level grids is complete
    // do post_timestep stuff here.
    //
    int finest_level = parent->finestLevel();

#ifdef RADIATION
    if (do_radiation && (level < finest_level)) {
        // computeTemp is not needed before or after this call because
        // setup for deferred sync does not touch state, only flux registers.
        radiation->deferred_sync_setup(level);

	if (do_reflux) {
	    radiation->reflux(level);
	    // Since radiation->reflux does not touch the fluid state,
	    // we do need to recompute Temp here.
	}
    }
#endif

#ifdef GRAVITY
    // Check the whole hierarchy before the syncs
    if (do_grav && level == 0 && gravity->test_results_of_solves() == 1 &&
        gravity->get_gravity_type() == "PoissonGrav")
    {
       if (verbose && ParallelDescriptor::IOProcessor())
          std::cout << "before gravity_sync " << std::endl;
       gravity->test_composite_phi(level);
    }
#endif

    if (do_reflux && level < finest_level) {

        MultiFab& S_new_crse = get_new_data(State_Type);

#ifdef GRAVITY
        MultiFab drho_and_drhoU;
        if (do_grav && gravity->get_gravity_type() == "PoissonGrav")  {
           // Define the update to rho and rhoU due to refluxing.
           drho_and_drhoU.define(grids,3+1,0,Fab_allocate);
           MultiFab::Copy(drho_and_drhoU,S_new_crse,Density,0,3+1,0);
           drho_and_drhoU.mult(-1.0);
        }
#endif

        reflux();

        // We need to do this before anything else because refluxing changes the values of coarse cells
        //    underneath fine grids with the assumption they'll be over-written by averaging down
        if (level < finest_level)
           avgDown();

        // This needs to be done after any changes to the state from refluxing.
        enforce_nonnegative_species(S_new_crse);

#ifdef GRAVITY
        if (do_grav && gravity->get_gravity_type() == "PoissonGrav" && gravity->NoSync() == 0)  {

	    MultiFab::Add(drho_and_drhoU, S_new_crse, Density, 0, 3+1, 0);

            MultiFab dphi(grids,1,0,Fab_allocate);
            dphi.setVal(0.);

            gravity->reflux_phi(level,dphi);
                
            // Compute (cross-level) gravity sync based on drho, dphi
            PArray<MultiFab> grad_delta_phi_cc(finest_level-level+1,PArrayManage); 
            for (int lev = level; lev <= finest_level; lev++) {
               grad_delta_phi_cc.set(lev-level,
                                     new MultiFab(getLevel(lev).boxArray(),3,0,Fab_allocate));
               grad_delta_phi_cc[lev-level].setVal(0.);
            }

	    int ncycle = parent->nCycle(level);
	    gravity->gravity_sync(level,finest_level,iteration,ncycle,drho_and_drhoU,dphi,grad_delta_phi_cc);

            for (int lev = level; lev <= finest_level; lev++)  
            {
              Real dt_lev = parent->dtLevel(lev);
              MultiFab&  S_new_lev = getLevel(lev).get_new_data(State_Type);
              Real cur_time = state[State_Type].curTime();

              const BoxArray& ba = getLevel(lev).boxArray();
              MultiFab grad_phi_cc(ba,3,0,Fab_allocate);
              gravity->get_new_grav_vector(lev,grad_phi_cc,cur_time);

#ifdef _OPENMP
#pragma omp parallel	      
#endif
	      {
		  FArrayBox sync_src;
		  FArrayBox dstate;

		  for (MFIter mfi(S_new_lev,true); mfi.isValid(); ++mfi)
		  {
		      const Box& bx = mfi.tilebox();
		      dstate.resize(bx,3+1);
		      if (lev == level) {
			  dstate.copy(drho_and_drhoU[mfi],bx);
		      } else {
			  dstate.setVal(0.); 
		      }
		      
		      // Compute sync source
		      sync_src.resize(bx,3+1);
		      BL_FORT_PROC_CALL(CA_SYNCGSRC,ca_syncgsrc)
			  (ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
			   BL_TO_FORTRAN_3D(grad_phi_cc[mfi]),
			   BL_TO_FORTRAN_3D(grad_delta_phi_cc[lev-level][mfi]),
			   BL_TO_FORTRAN_3D(S_new_lev[mfi]),
			   BL_TO_FORTRAN_3D(dstate),
			   BL_TO_FORTRAN_3D(sync_src),
			   dt_lev);

		      sync_src.mult(0.5*dt_lev);
		      S_new_lev[mfi].plus(sync_src,bx,0,Xmom,3);
		      S_new_lev[mfi].plus(sync_src,bx,0,Eden,1);
		  }
	      }
            }

            // Check the whole hierarchy after the syncs
            if (level == 0 && gravity->test_results_of_solves() == 1 && 
                              gravity->get_gravity_type() == "PoissonGrav")
            {
               if (verbose && ParallelDescriptor::IOProcessor())
                  std::cout << "after gravity_sync " << std::endl;
               gravity->test_composite_phi(level);
            }
        }
#endif
    }

    if (level < finest_level)
        avgDown();

    if (level == 0)
    {
        int nstep = parent->levelSteps(0);

        if ((sum_interval > 0) && (nstep%sum_interval == 0) )
            sum_integrated_quantities();

#ifdef GRAVITY
        if (moving_center) write_center();
#endif
    }

#ifdef RADIATION
    if (level == 0) {
      if (do_radiation) {
	for (int lev = finest_level; lev >= 0; lev--) {
	  radiation->analytic_solution(lev);
	}
      }
    }

    // diagnostic stuff
    
    if (level == 0) 
      do_energy_diagnostics();
#endif

    // Re-compute temperature after all the other updates.
    MultiFab& S_new = getLevel(level).get_new_data(State_Type);
    computeTemp(S_new);
}

void
Castro::post_restart ()
{
    BL_PROFILE("Castro::post_restart()");

   Real cur_time = state[State_Type].curTime();

#ifdef PARTICLES
    ParticlePostRestart(parent->theRestartFile());
#endif

#ifdef GRAVITY
    if (do_grav)
    {
        if (level == 0)
        {
            // Passing numpts_1d at level 0 
            int numpts_1d = get_numpts ();
            gravity->set_numpts_in_gravity(numpts_1d);

            for (int lev = 0; lev <= parent->finestLevel(); lev++)
            {
                AmrLevel& this_level  = getLevel(lev);
                Castro& cs_level = getLevel(lev);
                gravity->install_level(lev,&this_level,
                                       cs_level.Volume(),cs_level.Area());
            }

            if (moving_center == 1)
            {
               MultiFab&  S_new = get_new_data(State_Type);
               define_new_center(S_new,cur_time);
            }

            gravity->set_mass_offset(cur_time);

            if ( gravity->get_gravity_type() == "PoissonGrav")
            {
                if (gravity->NoComposite() != 1)
                {
 		   int use_previous_phi = 1;
		   gravity->multilevel_solve_for_new_phi(0,parent->finestLevel(),use_previous_phi);
                   if (gravity->test_results_of_solves() == 1)
                       gravity->test_composite_phi(level);
                }
#ifdef PARTICLES
                if (do_dm_particles)
                {
                    // Do solve if we haven't already done it above
                    if (gravity->NoComposite() == 1)
                       gravity->multilevel_solve_for_new_phi(0,parent->finestLevel());

                    for (int k = 0; k <= parent->finestLevel(); k++)
                    {
                        const BoxArray& ba = getLevel(k).boxArray();
                        MultiFab grav_vec_new(ba,3,0,Fab_allocate);
                        gravity->get_new_grav_vector(k,grav_vec_new,cur_time);
                    }
                }
#endif

            }

            if (grown_factor > 1)
                post_grown_restart();
        }
    }
#endif

#ifdef ROTATION
    MultiFab& phirot_new = get_new_data(PhiRot_Type);
    MultiFab& rot_new = get_new_data(Rotation_Type);
    MultiFab& S_new = get_new_data(State_Type);
    if (do_rotation)
      fill_rotation_field(phirot_new, rot_new, S_new, cur_time);
    else {
      phirot_new.setVal(0.0);
      rot_new.setVal(0.0);
    }
#endif    
    
#ifdef DIFFUSION
      // diffusion is a static object, only alloc if not already there
      if (diffusion == 0)
        diffusion = new Diffusion(parent,&phys_bc);

      if (level == 0) 
         for (int lev = 0; lev <= parent->finestLevel(); lev++) {
            AmrLevel& this_level = getLevel(lev);
                Castro& cs_level = getLevel(lev);
            diffusion->install_level(lev,&this_level,
                                     cs_level.Volume(),cs_level.Area());
         }
#endif

#ifdef REACTIONS
    MultiFab &React_new = get_new_data(Reactions_Type);
    React_new.setVal(0.0);
#endif

      set_special_tagging_flag(cur_time);
}

void
Castro::postCoarseTimeStep (Real cumtime)
{
    //
    // postCoarseTimeStep() is only called by level 0.
    //
#ifdef GRAVITY
    if (do_grav)
        gravity->set_mass_offset(cumtime, 0);
#endif

#ifdef PARTICLES
    if (Castro::theDMPC() && particle_move_type == "Random")
        ParticleMoveRandom();
#endif
}

void
Castro::post_regrid (int lbase,
                     int new_finest)
{
    delete fine_mask;
    fine_mask = 0;

#ifdef PARTICLES
    if (level == lbase) ParticleRedistribute();
#endif

#ifdef GRAVITY
    if (do_grav)
    {
       const Real cur_time = state[State_Type].curTime();
       if ( (level == lbase) && cur_time > 0.)  
       {
	   if ( gravity->get_gravity_type() == "PoissonGrav" && (gravity->NoComposite() != 1) ) {
	       int use_previous_phi = 1;
	       gravity->multilevel_solve_for_new_phi(level,new_finest,use_previous_phi);
	   }
       }
    }
#endif
}

void
Castro::post_init (Real stop_time)
{
    BL_PROFILE("Castro::post_init()");

    if (level > 0)
        return;

    // Send refinement data to Fortran

    int max_level = parent->maxLevel();
    int nlevs = max_level + 1;
    
    Real dx_level[3*nlevs];

    const Real* dx_coarse = geom.CellSize();    

    for (int dir = 0; dir < 3; dir++)
      dx_level[dir] = (ZFILL(dx_coarse))[dir];
    
    for (int lev = 1; lev <= max_level; lev++) {
      IntVect ref_ratio = parent->refRatio(lev-1);
      
      for (int dir = 0; dir < 3; dir++)
	if (dir < BL_SPACEDIM)
	  dx_level[3 * lev + dir] = dx_level[3 * (lev - 1) + dir] / ref_ratio[dir];
	else
	  dx_level[3 * lev + dir] = 0.0;
    }

    BL_FORT_PROC_CALL(SET_REFINEMENT_PARAMS,set_refinement_params)
      (max_level, dx_level);
    
    //
    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    //
    int finest_level = parent->finestLevel();
    for (int k = finest_level-1; k>= 0; k--)
        getLevel(k).avgDown();

    Real cur_time = state[State_Type].curTime();
    
#ifdef GRAVITY
    if (do_grav) {
       if (gravity->get_gravity_type() == "PoissonGrav") {

          // Calculate offset before first multilevel solve.
          gravity->set_mass_offset(cur_time);

          if (gravity->NoComposite() != 1)  {
             gravity->multilevel_solve_for_new_phi(level,finest_level);
             if (gravity->test_results_of_solves() == 1)
                gravity->test_composite_phi(level);
          }
       }
 
       // Make this call just to fill the initial state data.
       for (int k = 0; k <= parent->finestLevel(); k++)
       {
          BoxArray ba = getLevel(k).boxArray();
          MultiFab grav_vec_new(ba,3,0,Fab_allocate);
          gravity->get_new_grav_vector(k,grav_vec_new,cur_time);
       }
    }
#endif

#ifdef ROTATION
    MultiFab& phirot_new = get_new_data(PhiRot_Type);
    MultiFab& rot_new = get_new_data(Rotation_Type);
    MultiFab& S_new = get_new_data(State_Type);
    if (do_rotation)
      fill_rotation_field(phirot_new, rot_new, S_new, cur_time);
    else {
      phirot_new.setVal(0.0);
      rot_new.setVal(0.0);
    }
#endif
    
#ifdef RADIATION
    if (do_radiation) {
      // The option of whether to do a multilevel initialization is
      // controlled within the radiation class.

      radiation->post_init(level);

      for (int k = finest_level-1; k>= 0; k--) {
        getLevel(k).avgDown(Rad_Type);
      }

      do_energy_diagnostics();

      // re-estimate the initial timestep using the initialized
      // fine-level data and radiation field.
      //if (level == 0)
      //post_init_estDT();
    }
#endif

// Allow the user to define their own post_init functions.

#ifdef do_problem_post_init

    problem_post_init();

#endif

    if ( (sum_interval > 0) && (parent->levelSteps(0)%sum_interval == 0) )
        sum_integrated_quantities();

#ifdef GRAVITY
    if (level == 0 && moving_center == 1)
       write_center();
#endif
}

void
Castro::post_grown_restart ()
{
    if (level > 0)
        return;

    int finest_level = parent->finestLevel();
    Real cur_time = state[State_Type].curTime();
    
#ifdef GRAVITY
    if (do_grav) {
       if (gravity->get_gravity_type() == "PoissonGrav") {

          // Calculate offset before first multilevel solve.
          gravity->set_mass_offset(cur_time);

          if (gravity->NoComposite() != 1)  {
             gravity->multilevel_solve_for_new_phi(level,finest_level);
             if (gravity->test_results_of_solves() == 1)
                gravity->test_composite_phi(level);
          }
       }
 
       // Make this call just to fill the initial state data.
       for (int k = 0; k <= parent->finestLevel(); k++)
       {
          BoxArray ba = getLevel(k).boxArray();
          MultiFab grav_vec_new(ba,3,0,Fab_allocate);
          gravity->get_new_grav_vector(k,grav_vec_new,cur_time);
       }
    }
#endif

#ifdef ROTATION
    MultiFab& phirot_new = get_new_data(PhiRot_Type);
    MultiFab& rot_new = get_new_data(Rotation_Type);
    MultiFab& S_new = get_new_data(State_Type);
    if (do_rotation)
      fill_rotation_field(phirot_new, rot_new, S_new, cur_time);
    else {
      phirot_new.setVal(0.0);
      rot_new.setVal(0.0);
    }
#endif    
    
#ifdef RADIATION
    if (do_radiation) {
      // The option of whether to do a multilevel initialization is
      // controlled within the radiation class.

      radiation->post_init(level);

      for (int k = finest_level-1; k>= 0; k--) {
        getLevel(k).avgDown(Rad_Type);
      }

      do_energy_diagnostics();

      // re-estimate the initial timestep using the initialized
      // fine-level data and radiation field.
      //if (level == 0)
      //post_init_estDT();
    }
#endif
}

int
Castro::okToContinue ()
{
    if (level > 0)
        return 1;

    int test = 1;
    if (parent->dtLevel(0) < dt_cutoff) test = 0;

    return test; 
}

#ifdef AUX_UPDATE
void
Castro::advance_aux(Real time, Real dt)
{
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... special update for auxiliary variables \n";

    MultiFab&  S_old = get_old_data(State_Type);
    MultiFab&  S_new = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_old,true); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.tilebox();
        FArrayBox& old_fab = S_old[mfi];
        FArrayBox& new_fab = S_new[mfi];
	BL_FORT_PROC_CALL(CA_AUXUPDATE,ca_auxupdate)
	     (BL_TO_FORTRAN(old_fab),
	      BL_TO_FORTRAN(new_fab),
              box.loVect(), box.hiVect(),
              &dt);
    }
}
#endif

#ifdef LEVELSET

void
Castro::advance_levelset(Real time, Real dt)
{
    BL_PROFILE("Castro::advance_levelset()");

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... update levelset\n";

    MultiFab&  LS_new = get_new_data(LS_State_Type);
    LS_new.copy(get_old_data(LS_State_Type));

    for (MFIter mfi(LStype); mfi.isValid(); ++mfi)
    {
        IntFab& type = LStype[mfi];
        IntFab& nband = LSnband[mfi];
        IntFab& mine = LSmine[mfi];
        const Box& box = mfi.validbox();
        int nbandsize = nband.box().numPts();
        int minesize = mine.box().numPts();
        
        // Set nband data based on type
	BL_FORT_PROC_CALL(LS_NARROWBAND,ls_narrowband)
	     (BL_TO_FORTRAN(type),
              nband.dataPtr(), &nbandsize,
              mine.dataPtr(), &minesize,
              box.loVect(), box.hiVect());
    }
    
    Real phidt = dt;
    Real phit = dt;
    
    while(phit > 0) {
        
        int nGrowLS = 1;
        int nCompLS = 1;
        const Real* dx = geom.CellSize();
        
        phidt = phit;
        
	BoxLib::fill_boundary(LStype, 0, 1, geom);
        
        for (FillPatchIterator fpi(*this,LS_new,nGrowLS,state[LS_State_Type].curTime(),LS_State_Type,0,nCompLS);
             fpi.isValid(); ++fpi)
        {
            
            const FArrayBox& lsfab = fpi();  // Note that this is a temporary fab with nCompLS components filled by FPI
            const Box& box = fpi.validbox();
            IntFab& type = LStype[fpi];
            
            IntFab& nband = LSnband[fpi];
            IntFab& mine = LSmine[fpi];
            int nbandsize = nband.box().numPts();
            int minesize = mine.box().numPts();
            
            const FArrayBox& uadv = u_gdnv[0][fpi];
            const FArrayBox& vadv = u_gdnv[1][fpi];
#if (BL_SPACEDIM == 3)
            const FArrayBox& wadv = u_gdnv[2][fpi];
#endif

            Real lscfl;
            BL_FORT_PROC_CALL(LS_CFL,ls_cfl)
                  (&lscfl, 
                   BL_TO_FORTRAN(lsfab), BL_TO_FORTRAN(uadv),  BL_TO_FORTRAN(vadv),
#if (BL_SPACEDIM == 3)
                   BL_TO_FORTRAN(wadv), 
#endif
                   nband.dataPtr(), &nbandsize, mine.dataPtr(), &minesize,
                   box.loVect(), box.hiVect(), &phit, dx,
                   BL_TO_FORTRAN(type));
            phidt = std::min(phidt,lscfl);
        }
        ParallelDescriptor::ReduceRealMin(phidt);
        
        phit = phit - phidt;
	if (verbose && ParallelDescriptor::IOProcessor())
	  {
	    std::cout<<"phit"<<phit<<"\n";
	    std::cout<<"phidt"<<phidt<<"\n";
	  }
        bool reinit = false;
        
        for (FillPatchIterator fpi(*this,LS_new,nGrowLS,state[LS_State_Type].curTime(),LS_State_Type,0,nCompLS);
             fpi.isValid(); ++fpi)
        {
            const FArrayBox& lsfab = fpi();  // Note that this is a temporary fab with nCompLS components filled by FPI
            FArrayBox& phinew = LS_new[fpi];  // This is a fab of the actual state, with phi in the 0 component
            const Box& box = fpi.validbox();
            const FArrayBox& uadv = u_gdnv[0][fpi];
            const FArrayBox& vadv = u_gdnv[1][fpi];
#if (BL_SPACEDIM == 3)
            const FArrayBox& wadv = u_gdnv[2][fpi];
#endif

            IntFab& type = LStype[fpi];
            
            IntFab& nband = LSnband[fpi];
            IntFab& mine = LSmine[fpi];
            int nbandsize = nband.box().numPts();
            int minesize = mine.box().numPts();
            
            // Advance level set for this fab, result directly into state
            int reinit_flag;
	    BL_FORT_PROC_CALL(LS_PHIUPD,ls_phiupd)
                  (&reinit_flag,BL_TO_FORTRAN(lsfab), BL_TO_FORTRAN(phinew), 
                   BL_TO_FORTRAN(uadv), BL_TO_FORTRAN(vadv), 
#if (BL_SPACEDIM == 3)
                   BL_TO_FORTRAN(wadv), 
#endif
                   nband.dataPtr(), &nbandsize, mine.dataPtr(), &minesize,
                   box.loVect(), box.hiVect(), &phidt, dx,
                   BL_TO_FORTRAN(type)); 
            reinit |= reinit_flag;
        }
        
        ParallelDescriptor::ReduceBoolOr(reinit);
        
        if (reinit)
        {
            reinit_phi(time);
        }
    }
    
    delete [] u_gdnv;

}

void
Castro::reinit_phi(Real time)
{
    BL_PROFILE("Castro::reinit_phi()");

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "...... reinitialzing levelset\n";

    MultiFab& LS_new = get_new_data(LS_State_Type);
    const Real* dx   = geom.CellSize();
    int nGrowLS = 2;
    int nCompLS = 1;

    BoxLib::fill_boundary(LStype, 0, 1, geom);
        
    // Load valid region of phi
    Array<int> intfacep, intfacen, heap, heaploc;
    for (FillPatchIterator fpi(*this,LS_new,nGrowLS,state[LS_State_Type].curTime(),LS_State_Type,0,nCompLS);
         fpi.isValid(); ++fpi)
    {
        const FArrayBox& lsfab = fpi();  // Note that this is a temporary fab with nCompLS components filled by FPI
        FArrayBox& phinew = LS_new[fpi];  // This is a fab of the actual state, with phi in the 0 component
        IntFab& type = LStype[fpi];
        const Box& box = fpi.validbox();
        
        int intfacesize = box.numPts();
        int intfacenump = 0;
        int intfacenumn = 0;
        intfacep.resize(BL_SPACEDIM*intfacesize);
        intfacen.resize(BL_SPACEDIM*intfacesize);
        
        IntFab& nband = LSnband[fpi];
        IntFab& mine = LSmine[fpi];
        int nbandsize = nband.box().numPts();

        int nbandnum = 0;
        heap.resize(BL_SPACEDIM*nbandsize);
        int  typesize = type.box().numPts();
        heaploc.resize(typesize);

        // Set type for all cells in band to "outside band"
	BL_FORT_PROC_CALL(LS_RETYPIFY,ls_retypify)
            (BL_TO_FORTRAN(type), nband.dataPtr(), &nbandsize);

        // Set list of positive and negative narrowband points to be filled by FASTMARCH
	BL_FORT_PROC_CALL(LS_FINDINTERFACE,ls_findinterface)
            (BL_TO_FORTRAN(lsfab),  BL_TO_FORTRAN(phinew),  
             BL_TO_FORTRAN(type),  
             box.loVect(), box.hiVect(), dx, &intfacenump, &intfacenumn,
             intfacep.dataPtr(), intfacen.dataPtr(),
             nband.dataPtr(), &nbandsize, &intfacesize);

        // Fill in narrow band (+ve/-ve sides), set mines
        if(intfacenump > 0)
        {
            int positive = 1;
	    BL_FORT_PROC_CALL(LS_FASTMARCH,ls_fastmarch)
                 (BL_TO_FORTRAN(phinew),  BL_TO_FORTRAN(type),  
                  box.loVect(), box.hiVect(), dx, &intfacenump, intfacep.dataPtr(),
                  nband.dataPtr(), &nbandsize, &nbandnum,
                  mine.dataPtr(), &positive,&intfacesize,
                  heap.dataPtr(), heaploc.dataPtr());
        }
        
        if(intfacenumn > 0)
        {
            int negative = -1;
	    BL_FORT_PROC_CALL(LS_FASTMARCH,ls_fastmarch)
                 (BL_TO_FORTRAN(phinew),  BL_TO_FORTRAN(type),  
                  box.loVect(), box.hiVect(), dx,  &intfacenumn,  intfacen.dataPtr(),
            	  nband.dataPtr(), &nbandsize, &nbandnum,
            	  mine.dataPtr(), &negative,&intfacesize,
            	  heap.dataPtr(), heaploc.dataPtr());
        }
    }

    bool notdone = true;
    bool redo = false;
    
    // Check grow region and see if anything changes due to neighboring grids
    while (notdone)
    {
	BoxLib::fill_boundary(LS_new, geom);
	BoxLib::fill_boundary(LStype, 0, 1, geom);
        
        for (MFIter mfi(LS_new); mfi.isValid(); ++mfi)
        {
            FArrayBox& phinew = LS_new[mfi];  // This is a fab of the actual state, with phi in the 0 component
            IntFab& type = LStype[mfi];
            const Box& box = mfi.validbox();
            
            IntFab& nband = LSnband[mfi];
            int nbandsize = nband.box().numPts();
            
            int nbandnum = 0;
            
            int  numtemptype = type.box().numPts();            
            heaploc.resize(numtemptype);
            
	    BL_FORT_PROC_CALL(LS_NBANDNUMIFY,ls_nbandnumify)
               (nband.dataPtr(), &nbandsize,&nbandnum);
            
            int positive = 1;
            int redo_flag;
	    BL_FORT_PROC_CALL(LS_FASTMARCH2,ls_fastmarch2)
                 (&redo_flag, BL_TO_FORTRAN(phinew),BL_TO_FORTRAN(type),
                  box.loVect(), box.hiVect(), dx,
                  nband.dataPtr(), &nbandsize, &nbandnum, &positive, heaploc.dataPtr());
            redo |= redo_flag;
            
            int negative = -1;
	    BL_FORT_PROC_CALL(LS_FASTMARCH2,ls_fastmarch2)
                 (&redo_flag, BL_TO_FORTRAN(phinew),BL_TO_FORTRAN(type),
                  box.loVect(), box.hiVect(), dx,
                  nband.dataPtr(), &nbandsize, &nbandnum, &negative, heaploc.dataPtr());
            redo |= redo_flag;
        }
        ParallelDescriptor::ReduceBoolOr(redo);
        notdone = redo;
        redo = false;
    }

    for (MFIter mfi(LStype); mfi.isValid(); ++mfi)
    {
        IntFab& type = LStype[mfi];
        IntFab& nband = LSnband[mfi];
        IntFab& mine = LSmine[mfi];
        const Box& box = mfi.validbox();
        
        int nbandsize = nband.box().numPts();
        int minesize = mine.box().numPts();
        
	BL_FORT_PROC_CALL(LS_MINE,ls_mine)
             (BL_TO_FORTRAN(type),
              nband.dataPtr(), &nbandsize, mine.dataPtr(), &minesize,
              box.loVect(), box.hiVect());
    }
}

#endif

void
Castro::time_center_source_terms(MultiFab& S_new, MultiFab& ext_src_old, MultiFab &ext_src_new, Real dt)
{
    BL_PROFILE("Castro::time_center_source_terms()");

    // Subtract off half of the old source term, and add half of the new.

    ext_src_old.mult(-0.5*dt);
    ext_src_new.mult( 0.5*dt);
    MultiFab::Add(S_new,ext_src_old,0,0,S_new.nComp(),0);
    MultiFab::Add(S_new,ext_src_new,0,0,S_new.nComp(),0);
}

#ifdef SGS
void
Castro::getSource (Real time, Real dt, MultiFab& state, MultiFab& ext_src, MultiFab& sgs_state, MultiFab* sgs_fluxes)
{
   const Real* dx = geom.CellSize();

   ext_src.setVal(0.0,ext_src.nGrow());

   sgs_state.setVal(0.0);

   // Set these to zero so we always add them up correctly
   for (int dir = 0; dir < BL_SPACEDIM; dir++) 
       sgs_fluxes[dir].setVal(0.0);

   // We need to use temporary FABs to hold the fluxes for the old source
   //   because sgs_fluxes has no ghost cells but the fluxes array in the ext_src routine
   //   needs to have a layer of ghost cells for fluxes
   FArrayBox fluxx, fluxy, fluxz;

   for (MFIter mfi(ext_src); mfi.isValid(); ++mfi)
   {
        const Box& bx = grids[mfi.index()];
	
        Box bxx = bx; bxx.surroundingNodes(0); bxx.grow(1);
        fluxx.resize(bxx,NUM_STATE);

        Box bxy = bx; bxy.surroundingNodes(1); bxy.grow(1);
        fluxy.resize(bxy,NUM_STATE);

        Box bxz = bx; bxz.surroundingNodes(2); bxz.grow(1);
        fluxz.resize(bxz,NUM_STATE);

        BL_FORT_PROC_CALL(CA_EXT_SRC,ca_ext_src)
            (bx.loVect(), bx.hiVect(),
             BL_TO_FORTRAN(state[mfi]),
             BL_TO_FORTRAN(fluxx),
             BL_TO_FORTRAN(fluxy),
             BL_TO_FORTRAN(fluxz),
             BL_TO_FORTRAN(ext_src[mfi]),
             BL_TO_FORTRAN_N(sgs_mf[mfi],0),
             BL_TO_FORTRAN_N(sgs_mf[mfi],1),
             BL_TO_FORTRAN_N(sgs_mf[mfi],2),
             dx,&time,&dt);
  
        sgs_fluxes[0][mfi].copy(fluxx,0,0,NUM_STATE);
        sgs_fluxes[1][mfi].copy(fluxy,0,0,NUM_STATE);
        sgs_fluxes[2][mfi].copy(fluxz,0,0,NUM_STATE);
   }
   geom.FillPeriodicBoundary(ext_src,0,NUM_STATE);
}

#else

void
Castro::getSource (Real time, Real dt, MultiFab& state_old, MultiFab& state_new, MultiFab& ext_src, int ng)
{
   const Real* dx = geom.CellSize();
   const Real* prob_lo = geom.ProbLo();

   ext_src.setVal(0.0);

#ifdef _OPENMP
#pragma omp parallel
#endif    
   for (MFIter mfi(ext_src,true); mfi.isValid(); ++mfi)
     {
       const Box& bx = mfi.growntilebox(ng);
#ifdef DIMENSION_AGNOSTIC	   
       BL_FORT_PROC_CALL(CA_EXT_SRC,ca_ext_src)
	 (ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
	  BL_TO_FORTRAN_3D(state_old[mfi]),
	  BL_TO_FORTRAN_3D(state_new[mfi]),
	  BL_TO_FORTRAN_3D(ext_src[mfi]),
	  ZFILL(prob_lo),ZFILL(dx),&time,&dt);
#else	   
       BL_FORT_PROC_CALL(CA_EXT_SRC,ca_ext_src)
	 (bx.loVect(), bx.hiVect(),
	  BL_TO_FORTRAN(state_old[mfi]),
	  BL_TO_FORTRAN(state_new[mfi]),
	  BL_TO_FORTRAN(ext_src[mfi]),
	  prob_lo,dx,&time,&dt);
#endif	   
     }
}

#endif

#ifdef DIFFUSION
#ifdef TAU
void
Castro::define_tau (MultiFab& tau_diff, MultiFab& grav_vector, Real time)
{
   MultiFab& S_old = get_old_data(State_Type);
   const int ncomp = S_old.nComp();

   const Geometry& fine_geom = parent->Geom(parent->finestLevel());
   const Real*       dx_fine = fine_geom.CellSize();

   for (FillPatchIterator fpi(*this,S_old,NUM_GROW,time,State_Type,Density,ncomp);
                          fpi.isValid();++fpi)
   {
        Box bx(fpi.validbox());
        int i = fpi.index();
        BL_FORT_PROC_CALL(CA_DEFINE_TAU,ca_define_tau)
                 (bx.loVect(), bx.hiVect(),
                  BL_TO_FORTRAN(tau_diff[fpi]),
                  BL_TO_FORTRAN(fpi()),
                  BL_TO_FORTRAN(grav_vector[fpi]),
                  dx_fine);
   }
}
#endif

void
Castro::getTempDiffusionTerm (Real time, MultiFab& TempDiffTerm, MultiFab* tau)
{
    BL_PROFILE("Castro::getTempDiffusionTerm()");

   MultiFab& S_old = get_old_data(State_Type);
   if (verbose && ParallelDescriptor::IOProcessor()) 
      std::cout << "Calculating diffusion term at time " << time << std::endl;

#ifdef TAU
   if (tau == 0) 
     std::cerr << "ERROR:tau must be defined if USE_TAU = TRUE " << std::endl;
#else
   if (tau != 0) 
     std::cerr << "ERROR:tau must == 0 if USE_TAU = FALSE " << std::endl;
#endif

   // Fill coefficients at this level.
   PArray<MultiFab> coeffs(BL_SPACEDIM,PArrayManage);
   PArray<MultiFab> coeffs_temporary(3,PArrayManage); // This is what we pass to the dimension-agnostic Fortran
   for (int dir = 0; dir < 3; dir++) {
       if (dir < BL_SPACEDIM) {
	 coeffs.set(dir,new MultiFab(getEdgeBoxArray(dir), 1, 0, Fab_allocate));
	 coeffs_temporary.set(dir,new MultiFab(getEdgeBoxArray(dir), 1, 0, Fab_allocate));
       } else {
	 coeffs_temporary.set(dir,new MultiFab(grids, 1, 0, Fab_allocate));
       }
   }

   const Geometry& fine_geom = parent->Geom(parent->finestLevel());
   const Real*       dx_fine = fine_geom.CellSize();

   // Fill temperature at this level.
   MultiFab Temperature(grids,1,1,Fab_allocate);

   {
       FillPatchIterator fpi(*this,S_old,1,time,State_Type,0,NUM_STATE);
       MultiFab& state_old = fpi.get_mf();
       
       MultiFab::Copy(Temperature, state_old, Temp, 0, 1, 1);

       for (MFIter mfi(state_old); mfi.isValid(); ++mfi)
       {
	   const Box& bx = grids[mfi.index()];

	   BL_FORT_PROC_CALL(CA_FILL_TEMP_COND,ca_fill_temp_cond)
  	       (ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		BL_TO_FORTRAN_3D(state_old[mfi]),
#ifdef TAU
		BL_TO_FORTRAN_3D((*tau)[mfi]),
#endif
		BL_TO_FORTRAN_3D(coeffs_temporary[0][mfi]),
		BL_TO_FORTRAN_3D(coeffs_temporary[1][mfi]),
		BL_TO_FORTRAN_3D(coeffs_temporary[2][mfi]),
  	        ZFILL(dx_fine));
       }
   }

   // Now copy the temporary array results back to the
   // correctly dimensioned coeffs array.

   for (int dir = 0; dir < BL_SPACEDIM; dir++)
     MultiFab::Copy(coeffs[dir], coeffs_temporary[dir], 0, 0, 1, 0);
   
   if (Geometry::isAnyPeriodic())
     for (int d = 0; d < BL_SPACEDIM; d++)
       geom.FillPeriodicBoundary(coeffs[d]);

   MultiFab CrseTemp;
   if (level > 0) {
       // Fill temperature at next coarser level, if it exists.
       const BoxArray& crse_grids = getLevel(level-1).boxArray();
       CrseTemp.define(crse_grids,1,1,Fab_allocate);
       FillPatch(getLevel(level-1),CrseTemp,1,time,State_Type,Temp,1);
   }

   diffusion->applyop(level,Temperature,CrseTemp,TempDiffTerm,coeffs);

   // Extrapolate to ghost cells
   if (TempDiffTerm.nGrow() > 0) {
       for (MFIter mfi(TempDiffTerm); mfi.isValid(); ++mfi)
       {
	   const Box& bx = mfi.validbox();
	   BL_FORT_PROC_CALL(CA_TEMPDIFFEXTRAP,ca_tempdiffextrap)
               (ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		BL_TO_FORTRAN_3D(TempDiffTerm[mfi]));
       }
   }
}

void
Castro::getSpecDiffusionTerm (Real time, MultiFab& SpecDiffTerm)
{
    BL_PROFILE("Castro::getSpecDiffusionTerm()");

   MultiFab& S_old = get_old_data(State_Type);
   if (verbose && ParallelDescriptor::IOProcessor()) 
      std::cout << "Calculating species diffusion term at time " << time << std::endl;

   // Fill coefficients at this level.
   PArray<MultiFab> coeffs(BL_SPACEDIM,PArrayManage);
   PArray<MultiFab> coeffs_temporary(3,PArrayManage); // This is what we pass to the dimension-agnostic Fortran
   for (int dir = 0; dir < 3; dir++) {
       if (dir < BL_SPACEDIM) {
	 coeffs.set(dir,new MultiFab(getEdgeBoxArray(dir), 1, 0, Fab_allocate));
	 coeffs_temporary.set(dir,new MultiFab(getEdgeBoxArray(dir), 1, 0, Fab_allocate));
       } else {
	 coeffs_temporary.set(dir,new MultiFab(grids, 1, 0, Fab_allocate));
       }
   }

   const Geometry& fine_geom = parent->Geom(parent->finestLevel());
   const Real*       dx_fine = fine_geom.CellSize();

   FillPatchIterator fpi(*this,S_old,1,time,State_Type,0,NUM_STATE);
   MultiFab& state_old = fpi.get_mf();

   for (MFIter mfi(state_old); mfi.isValid(); ++mfi)
   {
       const Box& bx = grids[mfi.index()];

       BL_FORT_PROC_CALL(CA_FILL_SPEC_COEFF,ca_fill_spec_coeff)
  	       (ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		BL_TO_FORTRAN_3D(state_old[mfi]),
		BL_TO_FORTRAN_3D(coeffs_temporary[0][mfi]),
		BL_TO_FORTRAN_3D(coeffs_temporary[1][mfi]),
		BL_TO_FORTRAN_3D(coeffs_temporary[2][mfi]),
  	        ZFILL(dx_fine));
   }

   // Now copy the temporary array results back to the
   // correctly dimensioned coeffs array.
   for (int dir = 0; dir < BL_SPACEDIM; dir++)
     MultiFab::Copy(coeffs[dir], coeffs_temporary[dir], 0, 0, 1, 0);
   
   if (Geometry::isAnyPeriodic())
     for (int d = 0; d < BL_SPACEDIM; d++)
       geom.FillPeriodicBoundary(coeffs[d]);

   // Create MultiFabs that only hold the data for one species at a time.
   MultiFab Species(grids,1,1,Fab_allocate);
   MultiFab     SDT(grids,SpecDiffTerm.nComp(),SpecDiffTerm.nGrow(),Fab_allocate);
   MultiFab CrseSpec, CrseDen;
   if (level > 0) {
       const BoxArray& crse_grids = getLevel(level-1).boxArray();
       CrseSpec.define(crse_grids,1,1,Fab_allocate);
   }

   // Fill one species at a time at this level.
   for (int ispec = 0; ispec < NumSpec; ispec++)
   {
       MultiFab::Copy  (Species, state_old, FirstSpec+ispec, 0, 1, 1);
       MultiFab::Divide(Species, state_old, Density        , 0, 1, 1);

       // Fill temperature at next coarser level, if it exists.
       if (level > 0) 
       {
           FillPatch(getLevel(level-1),CrseSpec,1,time,State_Type,FirstSpec+ispec,1);
           FillPatch(getLevel(level-1),CrseDen ,1,time,State_Type,Density        ,1);
           MultiFab::Divide(CrseSpec, CrseDen, 0, 0, 1, 1);
       }

       diffusion->applyop(level,Species,CrseSpec,SDT,coeffs);

       // Extrapolate to ghost cells
       if (SDT.nGrow() > 0) {
           for (MFIter mfi(SDT); mfi.isValid(); ++mfi)
           {
    	       const Box& bx = mfi.validbox();
    	       BL_FORT_PROC_CALL(CA_TEMPDIFFEXTRAP,ca_tempdiffextrap)
                   (ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		    BL_TO_FORTRAN_3D(SDT[mfi]));
           }
       }
       // Copy back into SpecDiffTerm from the temporary SDT
       MultiFab::Copy(SpecDiffTerm, SDT, 0, ispec, 1, 1);
   }
}

#if (BL_SPACEDIM == 1)
// **********************************************
// Note: this currently just gets the term that looks like div(2 mu grad(u)) which is 
//       only part of the viscous term.  We assume that the coefficient that is filled is "2 mu"
// **********************************************
void
Castro::getViscousTerm (Real time, MultiFab& ViscousTermforMomentum, MultiFab& ViscousTermforEnergy)
{
    BL_PROFILE("Castro::getViscousTerm()");

   if (verbose && ParallelDescriptor::IOProcessor()) 
      std::cout << "Calculating viscous term at time " << time << std::endl;

   getFirstViscousTerm(time,ViscousTermforMomentum);

   MultiFab SecndTerm(grids,ViscousTermforMomentum.nComp(),ViscousTermforMomentum.nGrow(),Fab_allocate);
   getSecndViscousTerm(time,SecndTerm);
   MultiFab::Add(ViscousTermforMomentum, SecndTerm, 0, 0, ViscousTermforMomentum.nComp(), 0);

   getViscousTermForEnergy(time,ViscousTermforEnergy);
}

void
Castro::getFirstViscousTerm (Real time, MultiFab& ViscousTerm)
{
   MultiFab& S_old = get_old_data(State_Type);

   // Fill coefficients at this level.
   PArray<MultiFab> coeffs(BL_SPACEDIM,PArrayManage);
   PArray<MultiFab> coeffs_temporary(3,PArrayManage); // This is what we pass to the dimension-agnostic Fortran
   for (int dir = 0; dir < 3; dir++) {
       if (dir < BL_SPACEDIM) {
	 coeffs.set(dir,new MultiFab(getEdgeBoxArray(dir), 1, 0, Fab_allocate));
	 coeffs_temporary.set(dir,new MultiFab(getEdgeBoxArray(dir), 1, 0, Fab_allocate));
       } else {
	 coeffs_temporary.set(dir,new MultiFab(grids, 1, 0, Fab_allocate));
       }
   }

   const Geometry& fine_geom = parent->Geom(parent->finestLevel());
   const Real*       dx_fine = fine_geom.CellSize();

   // Fill velocity at this level.
   MultiFab Vel(grids,1,1,Fab_allocate);

   FillPatchIterator fpi(*this,S_old,1,time,State_Type,0,NUM_STATE);
   MultiFab& state_old = fpi.get_mf();

   // Remember this is just 1-d 
   MultiFab::Copy  (Vel, state_old, Xmom   , 0, 1, 1);
   MultiFab::Divide(Vel, state_old, Density, 0, 1, 1);

   for (MFIter mfi(state_old); mfi.isValid(); ++mfi)
   {
       const Box& bx = grids[mfi.index()];

       BL_FORT_PROC_CALL(CA_FILL_FIRST_VISC_COEFF,ca_fill_first_visc_coeff)
  	       (ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		BL_TO_FORTRAN_3D(state_old[mfi]),
		BL_TO_FORTRAN_3D(coeffs_temporary[0][mfi]),
		BL_TO_FORTRAN_3D(coeffs_temporary[1][mfi]),
		BL_TO_FORTRAN_3D(coeffs_temporary[2][mfi]),
  	        ZFILL(dx_fine));
   }

   // Now copy the temporary array results back to the
   // correctly dimensioned coeffs array.
   for (int dir = 0; dir < BL_SPACEDIM; dir++)
     MultiFab::Copy(coeffs[dir], coeffs_temporary[dir], 0, 0, 1, 0);
   
   if (Geometry::isAnyPeriodic())
     for (int d = 0; d < BL_SPACEDIM; d++)
       geom.FillPeriodicBoundary(coeffs[d]);

   MultiFab CrseVel, CrseDen;
   if (level > 0) {
       // Fill temperature at next coarser level, if it exists.
       const BoxArray& crse_grids = getLevel(level-1).boxArray();
       CrseVel.define(crse_grids,1,1,Fab_allocate);
       FillPatch(getLevel(level-1),CrseVel ,1,time,State_Type,Xmom   ,1);
       FillPatch(getLevel(level-1),CrseDen ,1,time,State_Type,Density,1);
       MultiFab::Divide(CrseVel, CrseDen, 0, 0, 1, 1);
   }
   diffusion->applyop(level,Vel,CrseVel,ViscousTerm,coeffs);

   // Extrapolate to ghost cells
   if (ViscousTerm.nGrow() > 0) {
       for (MFIter mfi(ViscousTerm); mfi.isValid(); ++mfi)
       {
	   const Box& bx = mfi.validbox();
	   BL_FORT_PROC_CALL(CA_TEMPDIFFEXTRAP,ca_tempdiffextrap)
               (ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		BL_TO_FORTRAN_3D(ViscousTerm[mfi]));
       }
   }
}

void
Castro::getSecndViscousTerm (Real time, MultiFab& ViscousTerm)
{
   MultiFab& S_old = get_old_data(State_Type);

   // Fill coefficients at this level.
   PArray<MultiFab> coeffs(BL_SPACEDIM,PArrayManage);
   PArray<MultiFab> coeffs_temporary(3,PArrayManage); // This is what we pass to the dimension-agnostic Fortran
   for (int dir = 0; dir < 3; dir++) {
       if (dir < BL_SPACEDIM) {
	 coeffs.set(dir,new MultiFab(getEdgeBoxArray(dir), 1, 0, Fab_allocate));
	 coeffs_temporary.set(dir,new MultiFab(getEdgeBoxArray(dir), 1, 0, Fab_allocate));
       } else {
	 coeffs_temporary.set(dir,new MultiFab(grids, 1, 0, Fab_allocate));
       }
   }

   const Geometry& fine_geom = parent->Geom(parent->finestLevel());
   const Real*       dx_fine = fine_geom.CellSize();

   // Fill velocity at this level.
   MultiFab Vel(grids,1,1,Fab_allocate);

   FillPatchIterator fpi(*this,S_old,1,time,State_Type,0,NUM_STATE);
   MultiFab& state_old = fpi.get_mf();

   // Remember this is just 1-d 
   MultiFab::Copy  (Vel, state_old, Xmom   , 0, 1, 1);
   MultiFab::Divide(Vel, state_old, Density, 0, 1, 1);

   for (MFIter mfi(state_old); mfi.isValid(); ++mfi)
   {
       const Box& bx = grids[mfi.index()];

       BL_FORT_PROC_CALL(CA_FILL_SECND_VISC_COEFF,ca_fill_secnd_visc_coeff)
  	       (ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		BL_TO_FORTRAN_3D(state_old[mfi]),
		BL_TO_FORTRAN_3D(coeffs_temporary[0][mfi]),
		BL_TO_FORTRAN_3D(coeffs_temporary[1][mfi]),
		BL_TO_FORTRAN_3D(coeffs_temporary[2][mfi]),
  	        ZFILL(dx_fine));
   }

   // Now copy the temporary array results back to the
   // correctly dimensioned coeffs array.
   for (int dir = 0; dir < BL_SPACEDIM; dir++)
     MultiFab::Copy(coeffs[dir], coeffs_temporary[dir], 0, 0, 1, 0);
   
   if (Geometry::isAnyPeriodic())
     for (int d = 0; d < BL_SPACEDIM; d++)
       geom.FillPeriodicBoundary(coeffs[d]);

   MultiFab CrseVel, CrseDen;
   if (level > 0) {
       // Fill temperature at next coarser level, if it exists.
       const BoxArray& crse_grids = getLevel(level-1).boxArray();
       CrseVel.define(crse_grids,1,1,Fab_allocate);
       FillPatch(getLevel(level-1),CrseVel ,1,time,State_Type,Xmom   ,1);
       FillPatch(getLevel(level-1),CrseDen ,1,time,State_Type,Density,1);
       MultiFab::Divide(CrseVel, CrseDen, 0, 0, 1, 1);
   }
   diffusion->applyViscOp(level,Vel,CrseVel,ViscousTerm,coeffs);

   // Extrapolate to ghost cells
   if (ViscousTerm.nGrow() > 0) {
       for (MFIter mfi(ViscousTerm); mfi.isValid(); ++mfi)
       {
	   const Box& bx = mfi.validbox();
	   BL_FORT_PROC_CALL(CA_TEMPDIFFEXTRAP,ca_tempdiffextrap)
               (ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		BL_TO_FORTRAN_3D(ViscousTerm[mfi]));
       }
   }
}

void
Castro::getViscousTermForEnergy (Real time, MultiFab& ViscousTerm)
{
   MultiFab& S_old = get_old_data(State_Type);

   // Fill coefficients at this level.
   PArray<MultiFab> coeffs(BL_SPACEDIM,PArrayManage);
   PArray<MultiFab> coeffs_temporary(3,PArrayManage); // This is what we pass to the dimension-agnostic Fortran
   for (int dir = 0; dir < 3; dir++) {
       if (dir < BL_SPACEDIM) {
	 coeffs.set(dir,new MultiFab(getEdgeBoxArray(dir), 1, 0, Fab_allocate));
	 coeffs_temporary.set(dir,new MultiFab(getEdgeBoxArray(dir), 1, 0, Fab_allocate));
       } else {
	 coeffs_temporary.set(dir,new MultiFab(grids, 1, 0, Fab_allocate));
       }
   }

   const Geometry& fine_geom = parent->Geom(parent->finestLevel());
   const Real*       dx_fine = fine_geom.CellSize();

   FillPatchIterator fpi(*this,S_old,2,time,State_Type,0,NUM_STATE);
   MultiFab& state_old = fpi.get_mf();

   // Remember this is just 1-d 
   int coord_type = Geometry::Coord();
   for (MFIter mfi(state_old); mfi.isValid(); ++mfi)
   {
       const Box& bx = grids[mfi.index()];

       BL_FORT_PROC_CALL(CA_COMPUTE_DIV_TAU_U,ca_compute_div_tau_u)
  	       (ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		BL_TO_FORTRAN_3D(ViscousTerm[mfi]),
		BL_TO_FORTRAN_3D(state_old[mfi]),
  	        ZFILL(dx_fine),&coord_type);
   }

   // Extrapolate to ghost cells
   if (ViscousTerm.nGrow() > 0) {
       for (MFIter mfi(ViscousTerm); mfi.isValid(); ++mfi)
       {
	   const Box& bx = mfi.validbox();
	   BL_FORT_PROC_CALL(CA_TEMPDIFFEXTRAP,ca_tempdiffextrap)
               (ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		BL_TO_FORTRAN_3D(ViscousTerm[mfi]));
       }
   }
}
#endif
#endif

void
Castro::reflux ()
{
    BL_PROFILE("Castro::reflux()");

    BL_ASSERT(level<parent->finestLevel());

    const Real strt = ParallelDescriptor::second();

    getFluxReg(level+1).Reflux(get_new_data(State_Type),volume,1.0,0,0,NUM_STATE,geom);

#ifdef SGS
    getSGSFluxReg(level+1).Reflux(get_new_data(State_Type),volume,1.0,0,0,NUM_STATE,geom);
#endif

#ifdef RADIATION
    if (Radiation::rad_hydro_combined) {
      getRADFluxReg(level+1).Reflux(get_new_data(Rad_Type),volume,1.0,0,0,Radiation::nGroups,geom);
    }
#endif

    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "Castro::reflux() at level " << level << " : time = " << end << std::endl;
#ifdef BL_LAZY
	});
#endif
    }
}

void
Castro::avgDown ()
{
    BL_PROFILE("Castro::avgDown()");

  if (level == parent->finestLevel()) return;

  avgDown(State_Type);

#ifdef GRAVITY
  avgDown(Gravity_Type);
  avgDown(PhiGrav_Type);
#endif

#ifdef ROTATION
  avgDown(Rotation_Type);
  avgDown(PhiRot_Type);
#endif

  avgDown(Source_Type);
  
#ifdef REACTIONS
  avgDown(Reactions_Type);
#endif

#ifdef LEVELSET
  avgDown(LS_State_Type);
#endif

#ifdef RADIATION
  if (do_radiation) {
    avgDown(Rad_Type);
  }
#endif

}

void
Castro::enforce_nonnegative_species (MultiFab& S_new)
{
#ifdef _OPENMP
#pragma omp parallel
#endif    
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.tilebox();
       BL_FORT_PROC_CALL(CA_ENFORCE_NONNEGATIVE_SPECIES,ca_enforce_nonnegative_species)
           (BL_TO_FORTRAN(S_new[mfi]),bx.loVect(),bx.hiVect());
    }
}

void
Castro::enforce_consistent_e (MultiFab& S)
{
#ifdef _OPENMP
#pragma omp parallel
#endif    
    for (MFIter mfi(S,true); mfi.isValid(); ++mfi)
    {
        const Box& box     = mfi.tilebox();
        const int* lo      = box.loVect();
        const int* hi      = box.hiVect();
        BL_FORT_PROC_CALL(CA_ENFORCE_CONSISTENT_E,ca_enforce_consistent_e)
          (ARLIM_3D(lo), ARLIM_3D(hi), BL_TO_FORTRAN_3D(S[mfi]));
    }
}

void
Castro::avgDown (int state_indx)
{
    BL_PROFILE("Castro::avgDown(state_indx)");

    if (level == parent->finestLevel()) return;

    Castro& fine_lev = getLevel(level+1);

    const Geometry& fgeom = fine_lev.geom;
    const Geometry& cgeom =          geom;

    MultiFab&  S_crse   = get_new_data(state_indx);
    MultiFab&  S_fine   = fine_lev.get_new_data(state_indx);

    BoxLib::average_down(S_fine, S_crse,
			 fgeom, cgeom,
			 0, S_fine.nComp(), fine_ratio);
}

void
Castro::allocOldData ()
{
    for (int k = 0; k < NUM_STATE_TYPE; k++)
        state[k].allocOldData();
}

void
Castro::removeOldData()
{
    AmrLevel::removeOldData();
}

void
Castro::errorEst (TagBoxArray& tags,
                  int          clearval,
                  int          tagval,
                  Real         time,
                  int          n_error_buf,
                  int          ngrow)
{
    BL_PROFILE("Castro::errorEst()");

    const int*  domain_lo = geom.Domain().loVect();
    const int*  domain_hi = geom.Domain().hiVect();
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();

    for (int j = 0; j < err_list.size(); j++)
    {
        MultiFab* mf = derive(err_list[j].name(), time, err_list[j].nGrow());

        BL_ASSERT(!(mf == 0));

#ifdef _OPENMP
#pragma omp parallel
#endif
	{
	    Array<int>  itags;

	    for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
	    {
		// FABs
		FArrayBox&  datfab  = (*mf)[mfi];
		TagBox&     tagfab  = tags[mfi];

		// tile box
		const Box&  tilebx  = mfi.tilebox();

		// physical tile box
		const RealBox& pbx  = RealBox(tilebx,geom.CellSize(),geom.ProbLo());

		//fab box
		const Box&  datbox  = datfab.box();

		// We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
		// So we are going to get a temporary integer array.
		tagfab.get_itags(itags, tilebx);

		// data pointer and index space
		int*        tptr    = itags.dataPtr();
		const int*  tlo     = tilebx.loVect();
		const int*  thi     = tilebx.hiVect();
		//
		const int*  lo      = tlo;
		const int*  hi      = thi;
		//
		const Real* xlo     = pbx.lo();
		//
		Real*       dat     = datfab.dataPtr();
		const int*  dlo     = datbox.loVect();
		const int*  dhi     = datbox.hiVect();
		const int   ncomp   = datfab.nComp();
		
		err_list[j].errFunc()(tptr, tlo, thi, &tagval,
				      &clearval, dat, dlo, dhi,
				      lo,hi, &ncomp, domain_lo, domain_hi,
				      dx, xlo, prob_lo, &time, &level);
		//
		// Now update the tags in the TagBox.
		//
                tagfab.tags_and_untags(itags, tilebx);
	    }
	}

        delete mf;
    }

    // Now we'll tag any user-specified zones using the full state array.

    MultiFab& S_new = get_new_data(State_Type);


#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Array<int>  itags;

	for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
	{
	    // tile box
	    const Box&  tilebx  = mfi.tilebox();

            TagBox&     tagfab  = tags[mfi];

	    // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
	    // So we are going to get a temporary integer array.
	    tagfab.get_itags(itags, tilebx);

            // data pointer and index space
	    int*        tptr    = itags.dataPtr();
	    const int*  tlo     = tilebx.loVect();
	    const int*  thi     = tilebx.hiVect();

#ifdef DIMENSION_AGNOSTIC
	    BL_FORT_PROC_CALL(SET_PROBLEM_TAGS, set_problem_tags)
	                     (tptr,  ARLIM_3D(tlo), ARLIM_3D(thi),
			      BL_TO_FORTRAN_3D(S_new[mfi]),
			      &tagval, &clearval, 
			      ARLIM_3D(tilebx.loVect()), ARLIM_3D(tilebx.hiVect()), 
			      ZFILL(dx), ZFILL(prob_lo), &time, &level);
#else	    
	    BL_FORT_PROC_CALL(SET_PROBLEM_TAGS, set_problem_tags)
	                     (tptr,  ARLIM(tlo), ARLIM(thi),
			      BL_TO_FORTRAN(S_new[mfi]),
			      &tagval, &clearval, 
			      tilebx.loVect(), tilebx.hiVect(), 
			      dx, prob_lo, &time, &level);
#endif
	    
	    //
	    // Now update the tags in the TagBox.
	    //
            tagfab.tags_and_untags(itags, tilebx);
	}
    }
}

MultiFab*
Castro::derive (const std::string& name,
                Real           time,
                int            ngrow)
{
#ifdef NEUTRINO
  if (name.substr(0,4) == "Neut") {
    // Extract neutrino energy group number from name string and
    // pass to fortran so that derive will have access to it:
    int is = atoi(name.c_str() + name.find('s') + 1);
    int ig = atoi(name.c_str() + name.find('g') + 1);

    BL_ASSERT(is < Radiation::nNeutrinoSpecies);
    for (int n = 0; n < is; n++) {
      ig += Radiation::nNeutrinoGroups[n];
    }

    BL_FORT_PROC_CALL(CA_SETGROUP,ca_setgroup)
      (ig);
  }
#endif

#ifdef PARTICLES
  return ParticleDerive(name,time,ngrow);
#else
   return AmrLevel::derive(name,time,ngrow);
#endif
}

void
Castro::derive (const std::string& name,
                Real           time,
                MultiFab&      mf,
                int            dcomp)
{
#ifdef NEUTRINO
  if (name.substr(0,4) == "Neut") {
    // Extract neutrino energy group number from name string and
    // pass to fortran so that derive will have access to it:
    int is = atoi(name.c_str() + name.find('s') + 1);
    int ig = atoi(name.c_str() + name.find('g') + 1);

    BL_ASSERT(is < Radiation::nNeutrinoSpecies);
    for (int n = 0; n < is; n++) {
      ig += Radiation::nNeutrinoGroups[n];
    }

    BL_FORT_PROC_CALL(CA_SETGROUP,ca_setgroup)
      (ig);
  }
#endif

    AmrLevel::derive(name,time,mf,dcomp);
}

void
Castro::network_init ()
{
   BL_FORT_PROC_CALL(CA_NETWORK_INIT,ca_network_init) ();
}

void
Castro::extern_init ()
{
  // initialize the external runtime parameters -- these will
  // live in the probin

  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "reading extern runtime parameters ..." << std::endl;
  }

  int probin_file_length = probin_file.length();
  Array<int> probin_file_name(probin_file_length);

  for (int i = 0; i < probin_file_length; i++)
    probin_file_name[i] = probin_file[i];

  BL_FORT_PROC_CALL(CA_EXTERN_INIT,ca_extern_init) 
    (probin_file_name.dataPtr(),
     &probin_file_length);
}

#ifdef SGS
void
Castro::reset_old_sgs(Real dt)
{
   MultiFab&   S_old = get_old_data(State_Type);
   MultiFab& sgs_old = get_old_data(SGS_Type);

   int is_old = 1; 

   for (MFIter mfi(S_old); mfi.isValid(); ++mfi)
   {
       const Box& bx = mfi.validbox();

       BL_FORT_PROC_CALL(CA_RESET_SGS,ca_reset_sgs)
           (BL_TO_FORTRAN(S_old[mfi]),
            BL_TO_FORTRAN_N(sgs_old[mfi],0),
            BL_TO_FORTRAN_N(sgs_old[mfi],0),
            bx.loVect(),bx.hiVect(),verbose,is_old,dt);
   }
}
void
Castro::reset_new_sgs(Real dt)
{
   MultiFab&   S_new = get_new_data(State_Type);
   MultiFab& sgs_old = get_old_data(SGS_Type);
   MultiFab& sgs_new = get_new_data(SGS_Type);

   int is_old = 0; 

   for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
   {
       const Box& bx = mfi.validbox();

       BL_FORT_PROC_CALL(CA_RESET_SGS,ca_reset_sgs)
           (BL_TO_FORTRAN(S_new[mfi]),
            BL_TO_FORTRAN_N(sgs_old[mfi],0),
            BL_TO_FORTRAN_N(sgs_new[mfi],0),
            bx.loVect(),bx.hiVect(),verbose,is_old,dt);
   }
}
#endif

void
Castro::reset_internal_energy(MultiFab& S_new)
{
    Real sum  = 0.;
    Real sum0 = 0.;

    if (parent->finestLevel() == 0 && print_energy_diagnostics)
    {
        // Pass in the multifab and the component
        sum0 = volWgtSumMF(&S_new,Eden,true);
    }

    // Synchronize (rho e) and (rho E) so they are consistent with each other
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        BL_FORT_PROC_CALL(RESET_INTERNAL_E,reset_internal_e)
	    (ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()), 
	     BL_TO_FORTRAN_3D(S_new[mfi]),
	     print_fortran_warnings);
    }

    if (parent->finestLevel() == 0 && print_energy_diagnostics)
    {
        // Pass in the multifab and the component
        sum = volWgtSumMF(&S_new,Eden,true);
#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
	ParallelDescriptor::ReduceRealSum(sum0);
	ParallelDescriptor::ReduceRealSum(sum);
        if (ParallelDescriptor::IOProcessor() && std::abs(sum-sum0) > 0)
            std::cout << "(rho E) added from reset terms                 : " << sum-sum0 << " out of " << sum0 << std::endl;
#ifdef BL_LAZY
	});
#endif
    }
}

void
Castro::computeTemp(MultiFab& State)
{
#ifdef RADIATION

  int resetEint = 1;
  radiation->computeTemp(State, resetEint);

#else
    reset_internal_energy(State);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(State,true); mfi.isValid(); ++mfi)
    { 
      const Box& bx = mfi.tilebox();
      BL_FORT_PROC_CALL(COMPUTE_TEMP,compute_temp)
	(ARLIM_3D(bx.loVect()),ARLIM_3D(bx.hiVect()),BL_TO_FORTRAN_3D(State[mfi]));
    }

#endif
}

void
Castro::set_special_tagging_flag(Real time)
{
   if (!do_special_tagging) return;

   MultiFab& S_new = get_new_data(State_Type);
   Real max_den = S_new.norm0(Density);

   int flag_was_changed = 0;
   BL_FORT_PROC_CALL(CA_SET_SPECIAL_TAGGING_FLAG,
                     ca_set_special_tagging_flag)(max_den,&flag_was_changed);
   if (ParallelDescriptor::IOProcessor()) {
      if (flag_was_changed == 1) {
        std::ofstream os("Bounce_time",std::ios::out);
        os << "T_Bounce " << time << std::endl;
        os.close();
      } 
   } 
}

int
Castro::get_numpts ()
{
     int numpts_1d;

     Box bx(geom.Domain());
     long nx = bx.size()[0];

#if (BL_SPACEDIM == 1)
     numpts_1d = nx;
#elif (BL_SPACEDIM == 2)
     long ny = bx.size()[1];
     Real ndiagsq = Real(nx*nx + ny*ny);
     numpts_1d = int(sqrt(ndiagsq))+2*NUM_GROW;  
#elif (BL_SPACEDIM == 3)
     long ny = bx.size()[1];
     long nz = bx.size()[2];
     Real ndiagsq = Real(nx*nx + ny*ny + nz*nz);
     numpts_1d = int(sqrt(ndiagsq))+2*NUM_GROW;
#endif 

     if (verbose && ParallelDescriptor::IOProcessor())
         std::cout << "Castro::numpts_1d at level  " << level << " is " << numpts_1d << std::endl;

     return numpts_1d;
}

void
Castro::make_radial_data(int is_new)
{
#if (BL_SPACEDIM > 1)

 // We only call this for level = 0
   BL_ASSERT(level == 0);
   
   int numpts_1d = get_numpts();

   Array<Real> radial_vol(numpts_1d,0);

   const Real* dx = geom.CellSize();
   Real  dr = dx[0];

   if (is_new == 1) {
      MultiFab& S = get_new_data(State_Type);
      int nc = S.nComp();
      Array<Real> radial_state(numpts_1d*nc,0);
      for (MFIter mfi(S); mfi.isValid(); ++mfi)
      {
         Box bx(mfi.validbox());
         BL_FORT_PROC_CALL(CA_COMPUTE_AVGSTATE,ca_compute_avgstate)
             (bx.loVect(), bx.hiVect(),dx,&dr,&nc,
              BL_TO_FORTRAN(     S[mfi]),radial_state.dataPtr(),
              BL_TO_FORTRAN(volume[mfi]),radial_vol.dataPtr(),
              geom.ProbLo(),&numpts_1d);
      }

      ParallelDescriptor::ReduceRealSum(radial_vol.dataPtr(),numpts_1d);
      ParallelDescriptor::ReduceRealSum(radial_state.dataPtr(),numpts_1d*nc);

      int first = 0; 
      int np_max = 0; 
      for (int i = 0; i < numpts_1d; i++) {
         if (radial_vol[i] > 0.) 
         {
            for (int j = 0; j < nc; j++) {
              radial_state[nc*i+j] /= radial_vol[i];
            }
         } else if (first == 0) {
            np_max = i; 
            first  = 1; 
         }
      }

      Array<Real> radial_state_short(np_max*nc,0);

      for (int i = 0; i < np_max; i++) {
         for (int j = 0; j < nc; j++) {
           radial_state_short[nc*i+j] = radial_state[nc*i+j];
         }
      }

      Real new_time = state[State_Type].curTime();
      BL_FORT_PROC_CALL(SET_NEW_OUTFLOW_DATA, set_new_outflow_data)
        (radial_state_short.dataPtr(),&new_time,&np_max,&nc);
   }
   else
   {
      MultiFab& S = get_old_data(State_Type);
      int nc = S.nComp();
      Array<Real> radial_state(numpts_1d*nc,0);
      for (MFIter mfi(S); mfi.isValid(); ++mfi)
      {
         Box bx(mfi.validbox());
         BL_FORT_PROC_CALL(CA_COMPUTE_AVGSTATE,ca_compute_avgstate)
             (bx.loVect(), bx.hiVect(),dx,&dr,&nc,
              BL_TO_FORTRAN(     S[mfi]),radial_state.dataPtr(),
              BL_TO_FORTRAN(volume[mfi]),radial_vol.dataPtr(),
              geom.ProbLo(),&numpts_1d);
      }

      ParallelDescriptor::ReduceRealSum(radial_vol.dataPtr(),numpts_1d);
      ParallelDescriptor::ReduceRealSum(radial_state.dataPtr(),numpts_1d*nc);

      int first = 0; 
      int np_max = 0; 
      for (int i = 0; i < numpts_1d; i++) {
         if (radial_vol[i] > 0.) 
         {
            for (int j = 0; j < nc; j++) {
              radial_state[nc*i+j] /= radial_vol[i];
            }
         } else if (first == 0) {
            np_max = i; 
            first  = 1; 
         }
      }

      Array<Real> radial_state_short(np_max*nc,0);

      for (int i = 0; i < np_max; i++) {
         for (int j = 0; j < nc; j++) {
           radial_state_short[nc*i+j] = radial_state[nc*i+j];
         }
      }

      Real old_time = state[State_Type].prevTime();
      BL_FORT_PROC_CALL(SET_OLD_OUTFLOW_DATA, set_old_outflow_data)
        (radial_state_short.dataPtr(),&old_time,&np_max,&nc);
   }

#endif
}

#ifdef GRAVITY
void
Castro::define_new_center(MultiFab& S, Real time)
{
    Real center[BL_SPACEDIM];
    const Real* dx = geom.CellSize();

    IntVect max_index = S.maxIndex(Density,0);
    Box bx(max_index,max_index);
    bx.grow(1);
    BoxArray ba(bx);
    MultiFab mf(ba,1,0);
    int owner = mf.DistributionMap()[0];

    // Define a cube 3-on-a-side around the point with the maximum density
    FillPatch(*this,mf,0,time,State_Type,Density,1);

    int mi[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) mi[i] = max_index[i];

    // Find the position of the "center" by interpolating from data at cell centers
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) 
    {
        BL_FORT_PROC_CALL(FIND_CENTER,find_center)(mf[mfi].dataPtr(),&center[0],mi,dx,geom.ProbLo());
    }
    // Now broadcast to everyone else.
    ParallelDescriptor::Bcast(&center[0], BL_SPACEDIM, owner);

    // Make sure if R-Z that center stays exactly on axis
    if ( Geometry::IsRZ() ) center[0] = 0;  

    BL_FORT_PROC_CALL(SET_CENTER,set_center)(&center[0]);
}

void
Castro::write_center ()
{
    int ndatalogs = parent->NumDataLogs();
 
    if ( (moving_center==1) && (ndatalogs > 0) && ParallelDescriptor::IOProcessor())
    {
       std::ostream& data_logc = parent->DataLog(0);
 
       int nstep = parent->levelSteps(0);
       Real time = state[State_Type].curTime();

       Real center[BL_SPACEDIM];
       BL_FORT_PROC_CALL(GET_CENTER,get_center)(center);
 
       if (time == 0.0) {
           data_logc << std::setw( 8) <<  "   nstep";
           data_logc << std::setw(14) <<  "         time  ";
           data_logc << std::setw(14) <<  "         center" << std::endl;;
       }
 
           data_logc << std::setw( 8) <<  nstep;
           data_logc << std::setw(14) <<  std::setprecision(6) <<  time;
           data_logc << std::setw(14) <<  std::setprecision(6) << center[0];
#if (BL_SPACEDIM >= 2)
           data_logc << std::setw(14) <<  std::setprecision(6) << center[1];
#endif
#if (BL_SPACEDIM == 3)
           data_logc << std::setw(14) <<  std::setprecision(6) << center[2];
#endif
           data_logc << std::endl;
    }
}

#endif


Real
Castro::getCPUTime()
{
  
  int numCores = ParallelDescriptor::NProcs();
#ifdef _OPENMP
  numCores = numCores*omp_get_max_threads();
#endif    

  Real T = numCores*(ParallelDescriptor::second() - startCPUTime) + 
    previousCPUTimeUsed;

  return T;
}


MultiFab*
Castro::build_fine_mask()
{
    BL_ASSERT(level > 0); // because we are building a mask for the coarser level

    if (fine_mask != 0) return fine_mask;

    BoxArray baf = parent->boxArray(level);
    baf.coarsen(crse_ratio);

    const BoxArray& bac = parent->boxArray(level-1);
    fine_mask = new MultiFab(bac,1,0);
    fine_mask->setVal(1.0);
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*fine_mask); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*fine_mask)[mfi];

	std::vector< std::pair<int,Box> > isects = baf.intersections(fab.box());

	for (int ii = 0; ii < isects.size(); ii++)
	{
	    fab.setVal(0.0,isects[ii].second,0);
	}
    }

    return fine_mask;
}
