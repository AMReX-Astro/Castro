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

#include <AMReX_Utility.H>
#include <AMReX_CONSTANTS.H>
#include <Castro.H>
#include <Castro_F.H>
#include <Derive_F.H>
#include <Castro_error_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_ParmParse.H>

#ifdef RADIATION
#include "Radiation.H"
#include "RAD_F.H"
#endif

#ifdef AMREX_PARTICLES
#include <AMReX_Particles_F.H>
#endif

#ifdef SELF_GRAVITY
#include "Gravity.H"
#endif

#ifdef DIFFUSION
#include "Diffusion.H"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace amrex;

bool         Castro::signalStopJob = false;

ErrorList    Castro::err_list;
int          Castro::num_err_list_default = 0;
int          Castro::radius_grow   = 1;
BCRec        Castro::phys_bc;
int          Castro::NUM_STATE     = -1;
int          Castro::NUM_GROW      = -1;

int          Castro::lastDtPlotLimited = 0;
Real         Castro::lastDtBeforePlotLimiting = 0.0;

Real         Castro::frac_change   = 1.e200;

int          Castro::Density       = -1;
int          Castro::Eden          = -1;
int          Castro::Eint          = -1;
int          Castro::Temp          = -1;
int          Castro::Xmom          = -1;
int          Castro::Ymom          = -1;
int          Castro::Zmom          = -1;
#ifdef HYBRID_MOMENTUM
int          Castro::Rmom          = -1;
int          Castro::Lmom          = -1;
int          Castro::Pmom          = -1;
#endif

int          Castro::QRHO = -1;
int          Castro::QU = -1;
int          Castro::QV = -1;
int          Castro::QW = -1;
int          Castro::QGAME = -1;
int          Castro::QPRES = -1;
int          Castro::QREINT = -1;
int          Castro::QTEMP = -1;
int          Castro::QFA = -1;
int          Castro::QFS = -1;
int          Castro::QFX = -1;
#ifdef MHD
int          Castro::QMAGX = -1;
int          Castro::QMAGY = -1;
int          Castro::QMAGZ = -1;
#endif
#ifdef RADIATION
int          Castro::QPTOT = -1;
int          Castro::QREITOT = -1;
int          Castro::QRAD = -1;
#endif

int          Castro::GDRHO = -1;
int          Castro::GDU = -1;
int          Castro::GDV = -1;
int          Castro::GDW = -1;
int          Castro::GDPRES = -1;
int          Castro::GDGAME = -1;
#ifdef RADIATION
int          Castro::GDLAMS = -1;
int          Castro::GDERADS = -1;
#endif

int          Castro::NumSpec       = 0;
int          Castro::FirstSpec     = -1;

int          Castro::NumAux        = 0;
int          Castro::FirstAux      = -1;

int          Castro::NumAdv        = 0;
int          Castro::FirstAdv      = -1;

#ifdef SHOCK_VAR
int          Castro::Shock         = -1;
#endif

int          Castro::NQSRC         = -1;
int          Castro::NQAUX         = -1;
int          Castro::NQ            = -1;

int          Castro::NGDNV         = -1;

Real         Castro::num_zones_advanced = 0.0;

Vector<std::string> Castro::source_names;

int          Castro::MOL_STAGES;
Vector< Vector<Real> > Castro::a_mol;
Vector<Real> Castro::b_mol;
Vector<Real> Castro::c_mol;


#include <castro_defaults.H>

#ifdef SELF_GRAVITY
// the gravity object
Gravity*     Castro::gravity  = 0;
#endif

#ifdef DIFFUSION
// the diffusion object
Diffusion*    Castro::diffusion  = 0;
#endif

#ifdef RADIATION
int          Castro::do_radiation = -1;

// the radiation object
Radiation*   Castro::radiation = 0;
#endif


std::string  Castro::probin_file = "probin";


#if BL_SPACEDIM == 1
#ifndef AMREX_USE_CUDA
IntVect      Castro::hydro_tile_size(1024);
#else
IntVect      Castro::hydro_tile_size(1024);
#endif
IntVect      Castro::no_tile_size(1024);
#elif BL_SPACEDIM == 2
#ifndef AMREX_USE_CUDA
IntVect      Castro::hydro_tile_size(1024,16);
#else
IntVect      Castro::hydro_tile_size(1024,64);
#endif
IntVect      Castro::no_tile_size(1024,1024);
#else
#ifndef AMREX_USE_CUDA
IntVect      Castro::hydro_tile_size(1024,16,16);
#else
IntVect      Castro::hydro_tile_size(1024,64,64);
#endif
IntVect      Castro::no_tile_size(1024,1024,1024);
#endif

// this will be reset upon restart
Real         Castro::previousCPUTimeUsed = 0.0;

Real         Castro::startCPUTime = 0.0;

int          Castro::Knapsack_Weight_Type = -1;
int          Castro::num_state_type = 0;

// Castro::variableSetUp is in Castro_setup.cpp
// variableCleanUp is called once at the end of a simulation
void
Castro::variableCleanUp ()
{
#ifdef SELF_GRAVITY
  if (gravity != 0) {
    if (verbose > 1 && ParallelDescriptor::IOProcessor()) {
      std::cout << "Deleting gravity in variableCleanUp..." << '\n';
    }
    delete gravity;
    gravity = 0;
  }
#endif

#ifdef DIFFUSION
  if (diffusion != 0) {
    if (verbose > 1 && ParallelDescriptor::IOProcessor()) {
      std::cout << "Deleting diffusion in variableCleanUp..." << '\n';
    }
    delete diffusion;
    diffusion = 0;
  }
#endif

#ifdef RADIATION
  if (radiation != 0) { int report = (verbose || radiation->verbose);
    if (report && ParallelDescriptor::IOProcessor()) {
      std::cout << "Deleting radiation in variableCleanUp..." << '\n';
    }
    delete radiation;
    radiation = 0;
    if (report && ParallelDescriptor::IOProcessor()) {
      std::cout << "                                        done" << std::endl;
    }
  }
#endif

#ifdef AMREX_PARTICLES
  delete TracerPC;
  TracerPC = 0;
#endif

    desc_lst.clear();

    ca_finalize_meth_params();

    network_finalize();
    eos_finalize();
#ifdef SPONGE
    sponge_finalize();
#endif
    amrinfo_finalize();
}

void
Castro::read_params ()
{
    static bool done = false;

    if (done) return;

    done = true;

    ParmParse pp("castro");

#include <castro_queries.H>

    // Get boundary conditions
    Vector<int> lo_bc(BL_SPACEDIM), hi_bc(BL_SPACEDIM);
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
                    amrex::Error();
                }
                if (hi_bc[dir] != Interior)
                {
                    std::cerr << "Castro::read_params:periodic in direction "
                              << dir
                              << " but high BC is not Interior\n";
                    amrex::Error();
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
                amrex::Error();
            }
            if (hi_bc[dir] == Interior)
            {
                std::cerr << "Castro::read_params:interior bc in direction "
                          << dir
                          << " but not periodic\n";
                amrex::Error();
            }
        }
    }

    if ( Geometry::IsRZ() && (lo_bc[0] != Symmetry) ) {
        std::cerr << "ERROR:Castro::read_params: must set r=0 boundary condition to Symmetry for r-z\n";
        amrex::Error();
    }

#if (BL_SPACEDIM == 1)
    if ( Geometry::IsSPHERICAL() )
    {
      if ( (lo_bc[0] != Symmetry) && (Geometry::ProbLo(0) == 0.0) )
      {
        std::cerr << "ERROR:Castro::read_params: must set r=0 boundary condition to Symmetry for spherical\n";
        amrex::Error();
      }
    }
#elif (BL_SPACEDIM == 2)
    if ( Geometry::IsSPHERICAL() )
      {
	amrex::Abort("We don't support spherical coordinate systems in 2D");
      }
#elif (BL_SPACEDIM == 3)
    if ( Geometry::IsRZ() )
      {
	amrex::Abort("We don't support cylindrical coordinate systems in 3D");
      }
    else if ( Geometry::IsSPHERICAL() )
      {
	amrex::Abort("We don't support spherical coordinate systems in 3D");
      }
#endif

#ifdef HYBRID_MOMENTUM
    // We do not support hybrid advection when using the HLLC solver.

    if (riemann_solver == 2) {
        amrex::Abort("HLLC Riemann solver unsupported when using hybrid momentum.");
    }
#endif


    // sanity checks

    if (grown_factor < 1)
       amrex::Error("grown_factor must be integer >= 1");

    if (cfl <= 0.0 || cfl > 1.0)
      amrex::Error("Invalid CFL factor; must be between zero and one.");

    // The timestep retry mechanism is currently incompatible with MOL.

    if (time_integration_method != CornerTransportUpwind && use_retry)
        amrex::Error("Method of lines integration is incompatible with the timestep retry mechanism.");

    // fourth order implies do_ctu=0
    if (fourth_order == 1 && time_integration_method == CornerTransportUpwind)
      {
	if (ParallelDescriptor::IOProcessor())
	    amrex::Error("WARNING: fourth_order requires a different time_integration_method");
      }

    // The CUDA MOL implementation is only supported in 3D right now.
#if defined(AMREX_USE_CUDA) && (AMREX_SPACEDIM < 3)
    if (time_integration_method != CornerTransportUpwind) {
        amrex::Error("Only the CTU advance is supported for 1D/2D when using CUDA.");
    }
#endif

    if (rebalanced_splitting && time_integration_method != CornerTransportUpwind) {
        amrex::Error("Rebalanced splitting is only supported for the CTU advance.");
    }

    if (hybrid_riemann == 1 && BL_SPACEDIM == 1)
      {
        std::cerr << "hybrid_riemann only implemented in 2- and 3-d\n";
        amrex::Error();
      }

    if (hybrid_riemann == 1 && (Geometry::IsSPHERICAL() || Geometry::IsRZ() ))
      {
        std::cerr << "hybrid_riemann should only be used for Cartesian coordinates\n";
        amrex::Error();
      }



    // Make sure not to call refluxing if we're not actually doing any hydro.
    if (do_hydro == 0) do_reflux = 0;

    if (max_dt < fixed_dt)
      {
	std::cerr << "cannot have max_dt < fixed_dt\n";
	amrex::Error();
      }

#ifdef AMREX_PARTICLES
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

    // The CUDA MOL implementation doesn't currently do radiation.
#ifdef AMREX_USE_CUDA
    if (do_radiation && time_integration_method != CornerTransportUpwind) {
        amrex::Error("Radiation is currently unsupported for MOL when using CUDA.");
    }
#endif
#endif

#ifdef ROTATION
    if (do_rotation) {
      if (rotational_period <= 0.0) {
	std::cerr << "Error:Castro::Rotation enabled but rotation period less than zero\n";
	amrex::Error();
      }
    }
    if (Geometry::IsRZ())
      rot_axis = 2;
#if (BL_SPACEDIM == 1)
      if (do_rotation) {
	std::cerr << "ERROR:Castro::Rotation not implemented in 1d\n";
	amrex::Error();
      }
#endif
#endif

   StateDescriptor::setBndryFuncThreadSafety(bndry_func_thread_safe);

   ParmParse ppa("amr");
   ppa.query("probin_file",probin_file);

    Vector<int> tilesize(BL_SPACEDIM);
    if (pp.queryarr("hydro_tile_size", tilesize, 0, BL_SPACEDIM))
    {
	for (int i=0; i<BL_SPACEDIM; i++) hydro_tile_size[i] = tilesize[i];
    }

}

Castro::Castro ()
    :
    prev_state(num_state_type)
{
}

Castro::Castro (Amr&            papa,
                int             lev,
                const Geometry& level_geom,
                const BoxArray& bl,
		const DistributionMapping& dm,
                Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,dm,time),
    prev_state(num_state_type)
{
    buildMetrics();

    initMFs();

    for (int i = 0; i < n_lost; i++) {
      material_lost_through_boundary_cumulative[i] = 0.0;
      material_lost_through_boundary_temp[i] = 0.0;
    }

#ifdef SELF_GRAVITY

   // Initialize to zero here in case we run with do_grav = false.
   MultiFab& new_grav_mf = get_new_data(Gravity_Type);
   new_grav_mf.setVal(0.0);

   if (do_grav) {
      // gravity is a static object, only alloc if not already there
      if (gravity == 0)
	gravity = new Gravity(parent,parent->finestLevel(),&phys_bc,Density);

      // Passing numpts_1d at level 0
      if (!Geometry::isAllPeriodic() && gravity != 0)
      {
         int numpts_1d = get_numpts();

	 // For 1D, we need to add ghost cells to the numpts
	 // given to us by Castro.

#if (BL_SPACEDIM == 1)
	 numpts_1d += 2 * NUM_GROW;
#endif

         gravity->set_numpts_in_gravity(numpts_1d);
      }

      gravity->install_level(level,this,volume,area);

      if (verbose && level == 0 &&  ParallelDescriptor::IOProcessor())
         std::cout << "Setting the gravity type to " << gravity->get_gravity_type() << std::endl;

#ifdef SELF_GRAVITY
      if (gravity->get_gravity_type() == "PoissonGrav" && gravity->NoComposite() != 0 && gravity->NoSync() == 0)
      {
	  std::cerr << "Error: not meaningful to have gravity.no_sync == 0 without having gravity.no_composite == 0.";
	  amrex::Error();
      }
#endif

       // We need to initialize this to zero since certain bc types don't overwrite the potential NaNs
       // ghost cells because they are only multiplying them by a zero coefficient.
       MultiFab& phi_new = get_new_data(PhiGrav_Type);
       phi_new.setVal(0.0,phi_new.nGrow());

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

   MultiFab& sources_new = get_new_data(Source_Type);
   sources_new.setVal(0.0, sources_new.nGrow());

#ifdef REACTIONS

   // Initialize reaction data to zero.

   MultiFab& reactions_new = get_new_data(Reactions_Type);
   reactions_new.setVal(0.0);

   // Initialize balance vector to zero.

   if (rebalanced_splitting) {
       MultiFab& balance_new = get_new_data(Balance_Type);
       balance_new.setVal(0.0);
   }

#endif

#ifdef SDC
#ifdef REACTIONS
   // Initialize reactions source term to zero.

   MultiFab& react_src_new = get_new_data(SDC_React_Type);
   react_src_new.setVal(0.0, NUM_GROW);
#endif
#endif

   if (Knapsack_Weight_Type > 0) {
    get_new_data(Knapsack_Weight_Type).setVal(1.0);
   }

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
      radiation->regrid(level, grids, dmap);
    }
#endif

}

Castro::~Castro ()
{
#ifdef RADIATION
    if (radiation != 0) {
      //radiation->cleanup(level);
      radiation->close(level);
    }
#endif
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

    volume.clear();
    volume.define(grids,dmap,1,NUM_GROW);
    geom.GetVolume(volume);

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        area[dir].clear();
	area[dir].define(getEdgeBoxArray(dir),dmap,1,NUM_GROW);
        geom.GetFaceArea(area[dir],dir);
    }
    for (int dir = BL_SPACEDIM; dir < 3; dir++)
    {
        area[dir].clear();
        area[dir].define(grids, dmap, 1, 0);
        area[dir].setVal(0.0);
    }

    dLogArea[0].clear();
#if (BL_SPACEDIM <= 2)
    geom.GetDLogA(dLogArea[0],grids,dmap,0,NUM_GROW);
#endif

    if (level == 0) setGridInfo();

    wall_time_start = 0.0;
}

// Initialize the MultiFabs and flux registers that live as class members.

void
Castro::initMFs()
{
    fluxes.resize(3);

    for (int dir = 0; dir < BL_SPACEDIM; ++dir)
	fluxes[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, NUM_STATE, 0));

    for (int dir = BL_SPACEDIM; dir < 3; ++dir)
	fluxes[dir].reset(new MultiFab(get_new_data(State_Type).boxArray(), dmap, NUM_STATE, 0));

    mass_fluxes.resize(3);

    for (int dir = 0; dir < BL_SPACEDIM; ++dir)
	mass_fluxes[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, 1, 0));

    for (int dir = BL_SPACEDIM; dir < 3; ++dir)
	mass_fluxes[dir].reset(new MultiFab(get_new_data(State_Type).boxArray(), dmap, 1, 0));

#if (BL_SPACEDIM <= 2)
    if (!Geometry::IsCartesian())
	P_radial.define(getEdgeBoxArray(0), dmap, 1, 0);
#endif

    // Keep track of which components of the momentum flux have pressure
    if (AMREX_SPACEDIM == 1 || (AMREX_SPACEDIM == 2 && Geometry::IsRZ())) {
        mom_flux_has_p[0][0] = false;
    }
    else {
        mom_flux_has_p[0][0] = true;
    }

    mom_flux_has_p[0][1] = false;
    mom_flux_has_p[0][2] = false;

    mom_flux_has_p[1][0] = false;
    mom_flux_has_p[1][1] = true;
    mom_flux_has_p[1][2] = false;

    mom_flux_has_p[2][0] = false;
    mom_flux_has_p[2][1] = false;
    mom_flux_has_p[2][2] = true;


#ifdef RADIATION
    if (Radiation::rad_hydro_combined) {
	rad_fluxes.resize(BL_SPACEDIM);
        for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
	    rad_fluxes[dir].reset(new MultiFab(getEdgeBoxArray(dir), dmap, Radiation::nGroups, 0));
	}
    }
#endif

    if (do_reflux && level > 0) {

	flux_reg.define(grids, dmap, crse_ratio, level, NUM_STATE);
	flux_reg.setVal(0.0);

#if (BL_SPACEDIM < 3)
	if (!Geometry::IsCartesian()) {
	    pres_reg.define(grids, dmap, crse_ratio, level, 1);
	    pres_reg.setVal(0.0);
	}
#endif

#ifdef RADIATION
	if (Radiation::rad_hydro_combined) {
	    rad_flux_reg.define(grids, dmap, crse_ratio, level, Radiation::nGroups);
	    rad_flux_reg.setVal(0.0);
	}
#endif

#ifdef SELF_GRAVITY
	if (do_grav && gravity->get_gravity_type() == "PoissonGrav" && gravity->NoSync() == 0) {
	    phi_reg.define(grids, dmap, crse_ratio, level, 1);
	    phi_reg.setVal(0.0);
	}
#endif

    }

    // Set the flux register scalings.

    if (do_reflux) {

	flux_crse_scale = -1.0;
	flux_fine_scale = 1.0;

	// The fine pressure scaling depends on dimensionality,
	// as the dimensionality determines the number of
	// adjacent zones. In 1D the face is a point so
	// there's only one fine neighbor for a given coarse
	// face; in 2D there's crse_ratio[1] faces adjacent
	// to a face perpendicular to the radial dimension;
	// and in 3D there would be crse_ratio**2, though
	// we do not separate the pressure out in 3D. Note
	// that the scaling by dt has already been handled
	// in the construction of the P_radial array.

	// The coarse pressure scaling is the same as for the
	// fluxes, we want the total refluxing contribution
	// over the full set of fine timesteps to equal P_radial.

#if (BL_SPACEDIM == 1)
	pres_crse_scale = -1.0;
	pres_fine_scale = 1.0;
#elif (BL_SPACEDIM == 2)
	pres_crse_scale = -1.0;
	pres_fine_scale = 1.0 / crse_ratio[1];
#endif

    }

    post_step_regrid = 0;

    lastDtRetryLimited = false;
    lastDtFromRetry = 1.e200;

}

void
Castro::setTimeLevel (Real time,
                      Real dt_old,
                      Real dt_new)
{
    AmrLevel::setTimeLevel(time,dt_old,dt_new);
}

void
Castro::setGridInfo ()
{

    // Send refinement data to Fortran. We do it here
    // because now the grids have been initialized and
    // we need this data for setting up the problem.
    // Note that this routine will always get called
    // on level 0, even if we are doing a restart,
    // so it is safe to put this here.

    if (level == 0) {

      const int max_level = parent->maxLevel();
      const int nlevs = max_level + 1;

      Real dx_level[3*nlevs];
      int domlo_level[3*nlevs];
      int domhi_level[3*nlevs];
      int ref_ratio_to_f[3*nlevs];
      int n_error_buf_to_f[nlevs];
      int blocking_factor_to_f[nlevs];

      const Real* dx_coarse = geom.CellSize();

      const int* domlo_coarse = geom.Domain().loVect();
      const int* domhi_coarse = geom.Domain().hiVect();

      for (int dir = 0; dir < 3; dir++) {
	dx_level[dir] = (ZFILL(dx_coarse))[dir];

	domlo_level[dir] = (ARLIM_3D(domlo_coarse))[dir];
	domhi_level[dir] = (ARLIM_3D(domhi_coarse))[dir];

	// Refinement ratio and error buffer on finest level are meaningless,
	// and we want them to be zero on the finest level because some
	// of the algorithms depend on this feature.

	ref_ratio_to_f[dir + 3 * (nlevs - 1)] = 0;
	n_error_buf_to_f[nlevs-1] = 0;
      }

      for (int lev = 0; lev <= max_level; lev++)
	blocking_factor_to_f[lev] = parent->blockingFactor(lev)[0];

      for (int lev = 1; lev <= max_level; lev++) {
	IntVect ref_ratio = parent->refRatio(lev-1);

	// Note that we are explicitly calculating here what the
	// data would be on refined levels rather than getting the
	// data directly from those levels, because some potential
	// refined levels may not exist at the beginning of the simulation.

	for (int dir = 0; dir < 3; dir++)
	  if (dir < BL_SPACEDIM) {
	    dx_level[3 * lev + dir] = dx_level[3 * (lev - 1) + dir] / ref_ratio[dir];
	    int ncell = (domhi_level[3 * (lev - 1) + dir] - domlo_level[3 * (lev - 1) + dir] + 1) * ref_ratio[dir];
	    domlo_level[3 * lev + dir] = domlo_level[dir];
	    domhi_level[3 * lev + dir] = domlo_level[3 * lev + dir] + ncell - 1;
	    ref_ratio_to_f[3 * (lev - 1) + dir] = ref_ratio[dir];
	  } else {
	    dx_level[3 * lev + dir] = 0.0;
	    domlo_level[3 * lev + dir] = 0;
	    domhi_level[3 * lev + dir] = 0;
	    ref_ratio_to_f[3 * (lev - 1) + dir] = 0;
	  }

	n_error_buf_to_f[lev - 1] = parent->nErrorBuf(lev - 1);
      }

      ca_set_grid_info(max_level, dx_level, domlo_level, domhi_level,
		       ref_ratio_to_f, n_error_buf_to_f, blocking_factor_to_f);

    }

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
	amrex::Abort("We don't support dx != dy");
      }
#elif (BL_SPACEDIM == 3)
    const Real SMALL = 1.e-13;
    if ( (fabs(dx[0] - dx[1]) > SMALL*dx[0]) || (fabs(dx[0] - dx[2]) > SMALL*dx[0]) )
      {
	amrex::Abort("We don't support dx != dy != dz");
      }
#endif

    ca_set_amr_info(level, -1, -1, -1.0, -1.0);

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

#ifdef SDC
#ifdef REACTIONS
   MultiFab& react_src_new = get_new_data(SDC_React_Type);
   react_src_new.setVal(0.0, NUM_GROW);
#endif
#endif

   if (Knapsack_Weight_Type > 0) {
       get_new_data(Knapsack_Weight_Type).setVal(1.0);
   }

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

#ifdef AMREX_DIMENSION_AGNOSTIC
          BL_FORT_PROC_CALL(CA_INITDATA,ca_initdata)
          (level, cur_time, ARLIM_3D(lo), ARLIM_3D(hi), ns,
  	   BL_TO_FORTRAN_ANYD(S_new[mfi]), ZFILL(dx),
  	   ZFILL(gridloc.lo()), ZFILL(gridloc.hi()));
#else
          BL_FORT_PROC_CALL(CA_INITDATA,ca_initdata)
  	  (level, cur_time, lo, hi, ns,
  	   BL_TO_FORTRAN(S_new[mfi]), dx,
  	   gridloc.lo(), gridloc.hi());
#endif

	  // Generate the initial hybrid momenta based on this user data.

#ifdef HYBRID_MOMENTUM
	  ca_init_hybrid_momentum(lo, hi, BL_TO_FORTRAN_ANYD(S_new[mfi]));
#endif

          // Verify that the sum of (rho X)_i = rho at every cell

          ca_check_initial_species(AMREX_ARLIM_3D(lo), AMREX_ARLIM_3D(hi),
				   BL_TO_FORTRAN_ANYD(S_new[mfi]));
       }
       enforce_consistent_e(S_new);

       // thus far, we assume that all initialization has worked on cell-centers
       // (to second-order, these are cell-averages, so we're done in that case).
       // For fourth-order, we need to convert to cell-averages now.
#ifndef AMREX_USE_CUDA
       if (fourth_order) {
         Sborder.define(grids, dmap, NUM_STATE, NUM_GROW);
         AmrLevel::FillPatch(*this, Sborder, NUM_GROW, cur_time, State_Type, 0, NUM_STATE);

         // note: this cannot be tiled

         for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
           {
             const Box& box     = mfi.validbox();

             ca_make_fourth_in_place(BL_TO_FORTRAN_BOX(box),
                                     BL_TO_FORTRAN_FAB(Sborder[mfi]));
           }

         // now copy back the averages
         MultiFab::Copy(S_new, Sborder, 0, 0, NUM_STATE, 0);
         Sborder.clear();
       }
#endif

       // Do a FillPatch so that we can get the ghost zones filled.

       int ng = S_new.nGrow();

       if (ng > 0)
	   AmrLevel::FillPatch(*this, S_new, ng, cur_time, State_Type, 0, S_new.nComp());
    }

    int is_new = 1;
    clean_state(is_new, S_new.nGrow());

#ifdef RADIATION
    if (do_radiation) {
      for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
          int i = mfi.index();

	  if (radiation->verbose > 2) {
	    std::cout << "Calling RADINIT at level " << level << ", grid "
		 << i << std::endl;
	  }

          RealBox    gridloc(grids[mfi.index()],
                             geom.CellSize(), geom.ProbLo());
          const Box& box = mfi.validbox();
          const int* lo  = box.loVect();
          const int* hi  = box.hiVect();

	  Rad_new[mfi].setVal(0.0);

#ifdef AMREX_DIMENSION_AGNOSTIC
	  BL_FORT_PROC_CALL(CA_INITRAD,ca_initrad)
	      (level, cur_time, ARLIM_3D(lo), ARLIM_3D(hi), Radiation::nGroups,
	       BL_TO_FORTRAN_ANYD(Rad_new[mfi]), ZFILL(dx),
	       ZFILL(gridloc.lo()), ZFILL(gridloc.hi()));
#else
	  BL_FORT_PROC_CALL(CA_INITRAD,ca_initrad)
	      (level, cur_time, lo, hi, Radiation::nGroups,
	       BL_TO_FORTRAN(Rad_new[mfi]),dx,
	       gridloc.lo(),gridloc.hi());
#endif

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


#ifdef SELF_GRAVITY
#if (BL_SPACEDIM > 1)
    if ( (level == 0) && (spherical_star == 1) ) {
       const int nc = S_new.nComp();
       const int n1d = get_numpts();
       allocate_outflow_data(&n1d,&nc);
       int is_new = 1;
       make_radial_data(is_new);
    }
#endif

    MultiFab& G_new = get_new_data(Gravity_Type);
    G_new.setVal(0.);

    MultiFab& phi_new = get_new_data(PhiGrav_Type);
    phi_new.setVal(0.);
#endif

    MultiFab& source_new = get_new_data(Source_Type);
    source_new.setVal(0., source_new.nGrow());

#ifdef ROTATION
    MultiFab& rot_new = get_new_data(Rotation_Type);
    rot_new.setVal(0.);

    MultiFab& phirot_new = get_new_data(PhiRot_Type);
    phirot_new.setVal(0.);
#endif

#ifdef AMREX_PARTICLES
    if (level == 0)
	init_particles();
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

    for (int s = 0; s < num_state_type; ++s) {
	MultiFab& state_MF = get_new_data(s);
	FillPatch(old, state_MF, state_MF.nGrow(), cur_time, s, 0, state_MF.nComp());
    }

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

    Real time = cur_time;

    // If we just triggered a regrid, we need to account for the fact that
    // the data on the coarse level has already been advanced.

    if (getLevel(level-1).post_step_regrid)
	time = prev_time;

    setTimeLevel(time,dt_old,dt);

    for (int s = 0; s < num_state_type; ++s) {
	MultiFab& state_MF = get_new_data(s);
	FillCoarsePatch(state_MF, 0, time, s, 0, state_MF.nComp());
    }
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

    ca_set_amr_info(level, -1, -1, -1.0, -1.0);

    Real estdt = max_dt;

    const MultiFab& stateMF = get_new_data(State_Type);

    const Real* dx = geom.CellSize();

    std::string limiter = "castro.max_dt";

    // Start the hydro with the max_dt value, but divide by CFL
    // to account for the fact that we multiply by it at the end.
    // This ensures that if max_dt is more restrictive than the hydro
    // criterion, we will get exactly max_dt for a timestep.

    Real estdt_hydro = max_dt / cfl;

#ifdef DIFFUSION
    if (do_hydro or diffuse_temp or diffuse_enth)
#else
    if (do_hydro)
#endif
    {

#ifdef RADIATION
      if (Radiation::rad_hydro_combined) {

	  // Compute radiation + hydro limited timestep.

#ifdef _OPENMP
#pragma omp parallel reduction(min:estdt_hydro)
#endif
        {
          Real dt = max_dt / cfl;

          const MultiFab& radMF = get_new_data(Rad_Type);
          FArrayBox gPr;

          for (MFIter mfi(stateMF, true); mfi.isValid(); ++mfi)
            {
              const Box& tbox = mfi.tilebox();
              const Box& vbox = mfi.validbox();

              gPr.resize(tbox);
              radiation->estimate_gamrPr(stateMF[mfi], radMF[mfi], gPr, dx, vbox);

              ca_estdt_rad(tbox.loVect(),tbox.hiVect(),
                           BL_TO_FORTRAN(stateMF[mfi]),
                           BL_TO_FORTRAN(gPr),
                           dx,&dt);
            }
          estdt_hydro = std::min(estdt_hydro, dt);
        }

      }
      else
      {
#endif

	  // Compute hydro-limited timestep.
	if (do_hydro)
	  {

#ifdef _OPENMP
#pragma omp parallel reduction(min:estdt_hydro)
#endif
	    {
	      Real dt = max_dt / cfl;

	      for (MFIter mfi(stateMF,true); mfi.isValid(); ++mfi)
		{
		  const Box& box = mfi.tilebox();

#pragma gpu
		  ca_estdt(AMREX_INT_ANYD(box.loVect()), AMREX_INT_ANYD(box.hiVect()),
			   BL_TO_FORTRAN_ANYD(stateMF[mfi]),
			   AMREX_REAL_ANYD(dx),
                           AMREX_MFITER_REDUCE_MIN(&dt));
		}
              estdt_hydro = std::min(estdt_hydro, dt);
            }
	  }

#ifdef DIFFUSION
	// Diffusion-limited timestep
	// Note that the diffusion uses the same CFL safety factor
	// as the main hydrodynamics timestep limiter.
	if (diffuse_temp)
	{
#ifdef _OPENMP
#pragma omp parallel reduction(min:estdt_hydro)
#endif
          {
            Real dt = max_dt / cfl;

            for (MFIter mfi(stateMF,true); mfi.isValid(); ++mfi)
              {
                const Box& box = mfi.tilebox();
                ca_estdt_temp_diffusion(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
                                        BL_TO_FORTRAN_ANYD(stateMF[mfi]),
                                        ZFILL(dx),&dt);
              }
            estdt_hydro = std::min(estdt_hydro, dt);
          }
	}
	if (diffuse_enth)
	{
#ifdef _OPENMP
#pragma omp parallel reduction(min:estdt_hydro)
#endif
          {
            Real dt = max_dt / cfl;

            for (MFIter mfi(stateMF,true); mfi.isValid(); ++mfi)
              {
                const Box& box = mfi.tilebox();
                ca_estdt_enth_diffusion(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
                                        BL_TO_FORTRAN_ANYD(stateMF[mfi]),
                                        ZFILL(dx),&dt);
              }
            estdt_hydro = std::min(estdt_hydro, dt);
          }
	}
#endif  // diffusion

#ifdef RADIATION
      }
#endif

       ParallelDescriptor::ReduceRealMin(estdt_hydro);
       estdt_hydro *= cfl;
       if (verbose && ParallelDescriptor::IOProcessor())
           std::cout << "...estimated hydro-limited timestep at level " << level << ": " << estdt_hydro << std::endl;

       // Determine if this is more restrictive than the maximum timestep limiting

       if (estdt_hydro < estdt) {
	 limiter = "hydro";
	 estdt = estdt_hydro;
       }
    }

#ifdef REACTIONS
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& R_new = get_new_data(Reactions_Type);

    // Dummy value to start with
    Real estdt_burn = max_dt;

    if (do_react) {

        // Compute burning-limited timestep.

#ifdef _OPENMP
#pragma omp parallel reduction(min:estdt_burn)
#endif
      {
        Real dt = max_dt;

        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
          {
            const Box& box = mfi.validbox();

            if (state[State_Type].hasOldData() && state[Reactions_Type].hasOldData()) {

              MultiFab& S_old = get_old_data(State_Type);
              MultiFab& R_old = get_old_data(Reactions_Type);

              ca_estdt_burning(ARLIM_3D(box.loVect()),ARLIM_3D(box.hiVect()),
                               BL_TO_FORTRAN_ANYD(S_old[mfi]),
                               BL_TO_FORTRAN_ANYD(S_new[mfi]),
                               BL_TO_FORTRAN_ANYD(R_old[mfi]),
                               BL_TO_FORTRAN_ANYD(R_new[mfi]),
                               ZFILL(dx),&dt_old,&dt);

            } else {

              ca_estdt_burning(ARLIM_3D(box.loVect()),ARLIM_3D(box.hiVect()),
                               BL_TO_FORTRAN_ANYD(S_new[mfi]),
                               BL_TO_FORTRAN_ANYD(S_new[mfi]),
                               BL_TO_FORTRAN_ANYD(R_new[mfi]),
                               BL_TO_FORTRAN_ANYD(R_new[mfi]),
                               ZFILL(dx),&dt_old,&dt);

            }

          }
        estdt_burn = std::min(estdt_burn,dt);
      }

      ParallelDescriptor::ReduceRealMin(estdt_burn);

      if (verbose && ParallelDescriptor::IOProcessor() && estdt_burn < max_dt)
        std::cout << "...estimated burning-limited timestep at level " << level << ": " << estdt_burn << std::endl;

      // Determine if this is more restrictive than the hydro limiting

      if (estdt_burn < estdt) {
        limiter = "burning";
        estdt = estdt_burn;
      }
    }
#endif

#ifdef RADIATION
    if (do_radiation) radiation->EstTimeStep(estdt, level);
#endif

    if (verbose && ParallelDescriptor::IOProcessor())
      std::cout << "Castro::estTimeStep (" << limiter << "-limited) at level " << level << ":  estdt = " << estdt << '\n';

    return estdt;
}

void
Castro::computeNewDt (int                   finest_level,
                      int                   sub_cycle,
                      Vector<int>&           n_cycle,
                      const Vector<IntVect>& ref_ratio,
                      Vector<Real>&          dt_min,
                      Vector<Real>&          dt_level,
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

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
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
          for (int i = 0; i <= finest_level; i++)
          {
              dt_min[i] = std::min(dt_min[i],dt_level[i]);
          }
       }
       else
       {
          //
          // Limit dt's by change_max * old dt,
          // if we didn't limit the last timestep
          // to hit a plotfile interval.
          //
          if (lastDtPlotLimited) {

              dt_min[0] = std::min(dt_min[0], lastDtBeforePlotLimiting);

              lastDtPlotLimited = 0;
              lastDtBeforePlotLimiting = 0.0;

          }
          else {

              for (int i = 0; i <= finest_level; i++)
              {
                  if (verbose && ParallelDescriptor::IOProcessor())
                      if (dt_min[i] > change_max*dt_level[i])
                      {
                          std::cout << "Castro::compute_new_dt : limiting dt at level "
                                    << i << '\n';
                          std::cout << " ... new dt computed: " << dt_min[i]
                                    << '\n';
                          std::cout << " ... but limiting to: "
                                    << change_max * dt_level[i] << " = " << change_max
                                    << " * " << dt_level[i] << '\n';
                      }
                  dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
              }

          }
       }
    }

    //
    // If we limited the last step by a retry,
    // apply that here if the retry-recommended
    // timestep is smaller than what we calculated.
    //
    for (int i = 0; i <= finest_level; ++i) {
        if (getLevel(i).lastDtRetryLimited == 1) {
            if (getLevel(i).lastDtFromRetry < dt_min[i]) {
                if (verbose && ParallelDescriptor::IOProcessor()) {
                    std::cout << " ... limiting dt at level " << i << " to: "
                              << getLevel(i).lastDtFromRetry << " = retry-limited timestep\n";
                }
                dt_min[i] = getLevel(i).lastDtFromRetry;
            }
        }
    }

    //
    // Find the minimum over all levels
    //
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_min[i]);
    }

    //
    // Optionally, limit dt's by the value of
    // plot_per or small_plot_per.
    //
    if (plot_per_is_exact) {

        const Real plot_per = parent->plotPer();

        if (plot_per > 0.0) {

            const Real cur_time = state[State_Type].curTime();

            // Calculate the new dt by comparing to the dt needed to get
            // to the next multiple of plot_per.

            const Real dtMod = std::fmod(cur_time, plot_per);

            Real newPlotDt;

            // Note that if we are just about exactly on a multiple of plot_per,
            // then we need to be careful to avoid floating point issues.

            if (std::abs(dtMod - plot_per) <= std::numeric_limits<Real>::epsilon()) {
                newPlotDt = plot_per + (plot_per - dtMod);
            }
            else {
                newPlotDt = plot_per - dtMod;
            }

            if (newPlotDt < dt_0) {
                lastDtPlotLimited = 1;
                lastDtBeforePlotLimiting = dt_0;
                dt_0 = newPlotDt;

                // Avoid taking timesteps that are so small that
                // they may cause problems in the hydrodynamics.

                const Real epsDt = 1.e-4 * lastDtBeforePlotLimiting;
                dt_0 = std::max(dt_0, epsDt);

                if (verbose)
                    amrex::Print() << " ... limiting dt to " << dt_0 << " to hit the next plot interval.\n";
            }

        }

    }

    if (small_plot_per_is_exact) {

        const Real small_plot_per = parent->smallplotPer();

        if (small_plot_per > 0.0) {

            const Real cur_time = state[State_Type].curTime();

            // Same logic as for plot_per_is_exact.

            const Real dtMod = std::fmod(cur_time, small_plot_per);

            Real newSmallPlotDt;

            if (std::abs(dtMod - small_plot_per) <= std::numeric_limits<Real>::epsilon()) {
                newSmallPlotDt = small_plot_per + (small_plot_per - dtMod);
            }
            else {
                newSmallPlotDt = small_plot_per - dtMod;
            }

            if (newSmallPlotDt < dt_0) {
                lastDtPlotLimited = 1;
                lastDtBeforePlotLimiting = dt_0;
                dt_0 = newSmallPlotDt;

                const Real epsDt = 1.e-4 * lastDtBeforePlotLimiting;
                dt_0 = std::max(dt_0, epsDt);

                if (verbose)
                    amrex::Print() << " ... limiting dt to " << dt_0 << " to hit the next smallplot interval.\n";
            }

        }

    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = std::numeric_limits<Real>::epsilon();
    Real cur_time = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) >= (stop_time - eps)) {
            dt_0 = stop_time - cur_time;
            if (verbose)
                amrex::Print() << " ... limiting dt to " << dt_0 << " to hit the stop_time.\n";
        }
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

void
Castro::computeInitialDt (int                   finest_level,
                          int                   sub_cycle,
                          Vector<int>&           n_cycle,
                          const Vector<IntVect>& ref_ratio,
                          Vector<Real>&          dt_level,
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
    // TODO/DEBUG: This will need to change for optimal subcycling.
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

    // Pass some information about the state of the simulation to a Fortran module.

    ca_set_amr_info(level, iteration, -1, -1.0, -1.0);

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

    // Now do the refluxing. If we're using gravity it
    // will also do the sync solve associated with the reflux.

    if (do_reflux && level < parent->finestLevel())
	reflux(level, level+1);

    // Ensure consistency with finer grids.

    if (level < finest_level)
	avgDown();


    // Clean up any aberrant state data generated by the reflux and average-down,
    // and then update quantities like temperature to be consistent.
    int is_new=1;
    MultiFab& S_new = get_new_data(State_Type);
    clean_state(is_new, S_new.nGrow());

    // Flush Fortran output

    if (verbose)
	flush_output();

#ifdef DO_PROBLEM_POST_TIMESTEP

    // Provide a hook for the user to do things after all of
    // the normal updates have been applied. The user is
    // responsible for any actions after this point, like
    // doing a computeTemp call if they change the state data.

    problem_post_timestep();

#endif

    if (level == 0)
    {
        int nstep = parent->levelSteps(0);
	Real dtlev = parent->dtLevel(0);
	Real cumtime = parent->cumTime() + dtlev;

	bool sum_int_test = false;

	if (sum_interval > 0) {

	  if (nstep%sum_interval == 0)
	    sum_int_test = true;

	}

	bool sum_per_test = false;

	if (sum_per > 0.0) {

	  const int num_per_old = floor((cumtime - dtlev) / sum_per);
	  const int num_per_new = floor((cumtime        ) / sum_per);

	  if (num_per_old != num_per_new)
	    sum_per_test = true;

	}

        if (sum_int_test || sum_per_test)
	  sum_integrated_quantities();

#ifdef SELF_GRAVITY
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

#ifdef AMREX_PARTICLES
    if (TracerPC)
    {
	const int ncycle = parent->nCycle(level);
	//
	// Don't redistribute/timestamp on the final subiteration except on the coarsest grid.
	//
	if (iteration < ncycle || level == 0)
	{
	    int ngrow = (level == 0) ? 0 : iteration;

	    TracerPC->Redistribute(level, parent->finestLevel(), ngrow);

	    TimestampParticles(ngrow+1);
	}
    }
#endif
}

void
Castro::post_restart ()
{
   BL_PROFILE("Castro::post_restart()");

   Real cur_time = state[State_Type].curTime();

#ifdef AMREX_PARTICLES
   ParticlePostRestart(parent->theRestartFile());
#endif

#ifdef SELF_GRAVITY
    if (do_grav)
    {
        if (level == 0)
        {
            // Passing numpts_1d at level 0
            int numpts_1d = get_numpts ();

#if (BL_SPACEDIM == 1)
	    numpts_1d += 2 * NUM_GROW;
#endif

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

		   // Update the maximum density, used in setting the solver tolerance.

		   gravity->update_max_rhs();

		   gravity->multilevel_solve_for_new_phi(0,parent->finestLevel(),use_previous_phi);
                   if (gravity->test_results_of_solves() == 1)
                       gravity->test_composite_phi(level);
                }
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

#ifdef DO_PROBLEM_POST_RESTART
    problem_post_restart();
#endif

}

void
Castro::postCoarseTimeStep (Real cumtime)
{
    // postCoarseTimeStep() is only called by level 0.
    BL_ASSERT(level == 0);
    AmrLevel::postCoarseTimeStep(cumtime);
#ifdef SELF_GRAVITY
    if (do_grav)
        gravity->set_mass_offset(cumtime, 0);
#endif
}

void
Castro::check_for_post_regrid (Real time)
{

    BL_PROFILE("Castro::check_for_post_regrid()");

    // Check whether we have any zones at this time signifying that they
    // need to be tagged that do not have corresponding zones on the
    // fine level.

    if (level < parent->maxLevel()) {

	TagBoxArray tags(grids, dmap);

	tags.setVal(TagBox::CLEAR);

	for (int i = 0; i < err_list.size(); ++i)
            apply_tagging_func(tags, TagBox::CLEAR, TagBox::SET, time, i);

        apply_problem_tags(tags, TagBox::CLEAR, TagBox::SET, time);

	// Globally collate the tags.

	Vector<IntVect> tvec;

	tags.collate(tvec);

	// If we requested any tags at all, we have a potential trigger for a regrid.

	int num_tags = tvec.size();

	bool missing_on_fine_level = false;

	if (num_tags > 0) {

	    if (level == parent->finestLevel()) {

		// If there is no level above us at all, we know a regrid is needed.

		missing_on_fine_level = true;

	    } else {

		// If there is a level above us, we need to check whether
		// every tagged zone has a corresponding entry in the fine
		// grid. If not, we need to trigger a regrid so we can get
		// those other zones refined too.

		const BoxArray& fgrids = getLevel(level+1).boxArray();

		for (int i = 0; i < tvec.size(); ++i) {

		    Box c_bx(tvec[i], tvec[i]);
		    Box f_bx = c_bx.refine(parent->refRatio(level));

		    if (!fgrids.contains(f_bx)) {
			missing_on_fine_level = true;
			break;
		    }

		}

	    }

	}

	if (missing_on_fine_level) {
	    post_step_regrid = 1;

	    if (amrex::ParallelDescriptor::IOProcessor()) {
		std::cout << "\n"
			  << "Current refinement pattern insufficient at level " << level << ".\n"
			  << "Performing a regrid to obtain more refinement.\n";

	    }
	}

    }

}

void
Castro::post_regrid (int lbase,
                     int new_finest)
{
    fine_mask.clear();

#ifdef AMREX_PARTICLES
    if (TracerPC && level == lbase) {
	TracerPC->Redistribute(lbase);
    }
#endif

#ifdef SELF_GRAVITY
    if (do_grav)
    {

	if (use_post_step_regrid && getLevel(lbase).post_step_regrid && gravity->get_gravity_type() == "PoissonGrav") {

	   if (level > lbase) {

	       // In the case where we're coming here during a regrid that occurs
	       // after a timestep, we only want to interpolate the gravitational
	       // field from the old time. The state data will already have been
	       // filled, so all we need to do is interpolate the grad_phi data.

	       // Instantiate a bare physical BC function for grad_phi. It doesn't do anything
	       // since the fine levels for Poisson gravity do not touch the physical boundary.

	       GradPhiPhysBCFunct gp_phys_bc;

	       // We need to use a nodal interpolater.

	       Interpolater* gp_interp = &node_bilinear_interp;

	       Vector<MultiFab*> grad_phi_coarse = amrex::GetVecOfPtrs(gravity->get_grad_phi_prev(level-1));
	       Vector<MultiFab*> grad_phi_fine = amrex::GetVecOfPtrs(gravity->get_grad_phi_curr(level));

	       Real time = getLevel(lbase).get_state_data(Gravity_Type).prevTime();

	       // For the BCs, we will use the Gravity_Type BCs for convenience, but these will
	       // not do anything because we do not fill on physical boundaries.

	       const Vector<BCRec>& gp_bcs = getLevel(level).get_desc_lst()[Gravity_Type].getBCs();

	       for (int n = 0; n < BL_SPACEDIM; ++n) {
		   amrex::InterpFromCoarseLevel(*grad_phi_fine[n], time, *grad_phi_coarse[n],
						0, 0, 1,
						parent->Geom(level-1), parent->Geom(level),
						gp_phys_bc, 0, gp_phys_bc, 0, parent->refRatio(level-1),
						gp_interp, gp_bcs, 0);
	       }

	   }

       } else {

	    const Real cur_time = state[State_Type].curTime();
	    if ( (level == lbase) && cur_time > 0.)
	    {
		if ( gravity->get_gravity_type() == "PoissonGrav" && (gravity->NoComposite() != 1) ) {
		    int use_previous_phi = 1;

		    // Update the maximum density, used in setting the solver tolerance.

		    if (level == 0)
			gravity->update_max_rhs();

		    gravity->multilevel_solve_for_new_phi(level,new_finest,use_previous_phi);

		}

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

    //
    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    //
    int finest_level = parent->finestLevel();
    for (int k = finest_level-1; k>= 0; k--)
        getLevel(k).avgDown();

#ifdef SELF_GRAVITY

    if (do_grav) {

       Real cur_time = state[State_Type].curTime();

       if (gravity->get_gravity_type() == "PoissonGrav") {

	  // Update the maximum density, used in setting the solver tolerance.

	  gravity->update_max_rhs();

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
          MultiFab& grav_new = getLevel(k).get_new_data(Gravity_Type);
          gravity->get_new_grav_vector(k,grav_new,cur_time);
       }
    }
#endif

#ifdef ROTATION
    MultiFab& phirot_new = get_new_data(PhiRot_Type);
    MultiFab& rot_new = get_new_data(Rotation_Type);
    MultiFab& S_new = get_new_data(State_Type);
    if (do_rotation) {
      Real cur_time = state[State_Type].curTime();
      fill_rotation_field(phirot_new, rot_new, S_new, cur_time);
    }
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

#ifdef DO_PROBLEM_POST_INIT

    problem_post_init();

#endif

        int nstep = parent->levelSteps(0);
	Real dtlev = parent->dtLevel(0);
	Real cumtime = parent->cumTime();
	if (cumtime != 0.0) cumtime += dtlev;

	bool sum_int_test = false;

	if (sum_interval > 0) {

	  if (nstep%sum_interval == 0)
	    sum_int_test = true;

	}

	bool sum_per_test = false;

	if (sum_per > 0.0) {

	  const int num_per_old = floor((cumtime - dtlev) / sum_per);
	  const int num_per_new = floor((cumtime        ) / sum_per);

	  if (num_per_old != num_per_new)
	    sum_per_test = true;

	}

        if (sum_int_test || sum_per_test)
	  sum_integrated_quantities();

#ifdef SELF_GRAVITY
    if (level == 0 && moving_center == 1)
       write_center();
#endif
}

void
Castro::post_grown_restart ()
{
    if (level > 0)
        return;

#ifdef SELF_GRAVITY
    if (do_grav) {
	int finest_level = parent->finestLevel();
	Real cur_time = state[State_Type].curTime();

	if (gravity->get_gravity_type() == "PoissonGrav") {

	  // Update the maximum density, used in setting the solver tolerance.

	  gravity->update_max_rhs();

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
          MultiFab& grav_new = getLevel(k).get_new_data(Gravity_Type);
          gravity->get_new_grav_vector(k,grav_new,cur_time);
       }
    }
#endif

#ifdef ROTATION
    MultiFab& phirot_new = get_new_data(PhiRot_Type);
    MultiFab& rot_new = get_new_data(Rotation_Type);
    MultiFab& S_new = get_new_data(State_Type);
    if (do_rotation) {
      Real cur_time = state[State_Type].curTime();
      fill_rotation_field(phirot_new, rot_new, S_new, cur_time);
    }
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

      int finest_level = parent->finestLevel();

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

    if (signalStopJob) {
      test = 0;
      if (ParallelDescriptor::IOProcessor())
	std::cout << " Signalling a stop of the run due to signalStopJob = true." << std::endl;
    }
    else if (parent->dtLevel(0) < dt_cutoff) {
      test = 0;
      if (ParallelDescriptor::IOProcessor())
	std::cout << " Signalling a stop of the run because dt < dt_cutoff." << std::endl;
    }

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
	void ca_auxupdate(BL_TO_FORTRAN(old_fab),
			  BL_TO_FORTRAN(new_fab),
			  box.loVect(), box.hiVect(),
			  &dt);
    }
}
#endif


void
Castro::FluxRegCrseInit() {

    if (level == parent->finestLevel()) return;

    Castro& fine_level = getLevel(level+1);

    for (int i = 0; i < BL_SPACEDIM; ++i)
	fine_level.flux_reg.CrseInit(*fluxes[i], i, 0, 0, NUM_STATE, flux_crse_scale);

#if (BL_SPACEDIM <= 2)
    if (!Geometry::IsCartesian())
	fine_level.pres_reg.CrseInit(P_radial, 0, 0, 0, 1, pres_crse_scale);
#endif

#ifdef RADIATION
    if (Radiation::rad_hydro_combined)
	for (int i = 0; i < BL_SPACEDIM; ++i)
	    fine_level.rad_flux_reg.CrseInit(*rad_fluxes[i], i, 0, 0, Radiation::nGroups, flux_crse_scale);
#endif

}


void
Castro::FluxRegFineAdd() {

    if (level == 0) return;

    for (int i = 0; i < BL_SPACEDIM; ++i)
	flux_reg.FineAdd(*fluxes[i], i, 0, 0, NUM_STATE, flux_fine_scale);

#if (BL_SPACEDIM <= 2)
    if (!Geometry::IsCartesian())
	getLevel(level).pres_reg.FineAdd(P_radial, 0, 0, 0, 1, pres_fine_scale);
#endif

#ifdef RADIATION
    if (Radiation::rad_hydro_combined)
	for (int i = 0; i < BL_SPACEDIM; ++i)
	    getLevel(level).rad_flux_reg.FineAdd(*rad_fluxes[i], i, 0, 0, Radiation::nGroups, flux_fine_scale);
#endif

}


void
Castro::reflux(int crse_level, int fine_level)
{
    BL_PROFILE("Castro::reflux()");

    BL_ASSERT(fine_level > crse_level);

    const Real strt = ParallelDescriptor::second();

#ifdef SELF_GRAVITY
    int nlevs = fine_level - crse_level + 1;

    Vector<std::unique_ptr<MultiFab> > drho(nlevs);
    Vector<std::unique_ptr<MultiFab> > dphi(nlevs);

    if (do_grav && gravity->get_gravity_type() == "PoissonGrav" && gravity->NoSync() == 0)  {

	for (int lev = crse_level; lev <= fine_level; ++lev) {

	    const auto& amrlevel = getLevel(lev);
	    const auto& ba = amrlevel.boxArray();
	    const auto& dm = amrlevel.DistributionMap();

	    drho[lev - crse_level].reset(new MultiFab(ba, dm, 1, 0));
	    dphi[lev - crse_level].reset(new MultiFab(ba, dm, 1, 0));

	    drho[lev - crse_level]->setVal(0.0);
	    dphi[lev - crse_level]->setVal(0.0);

	}

    }
#endif

    FluxRegister* reg;

    for (int lev = fine_level; lev > crse_level; --lev) {

	reg = &getLevel(lev).flux_reg;

	Castro& crse_lev = getLevel(lev-1);
	Castro& fine_lev = getLevel(lev);

	MultiFab& crse_state = crse_lev.get_new_data(State_Type);

	// Clear out the data that's not on coarse-fine boundaries so that this register only
	// modifies the fluxes on coarse-fine interfaces.

	reg->ClearInternalBorders(crse_lev.geom);

	// Trigger the actual reflux on the coarse level now.

	reg->Reflux(crse_state, crse_lev.volume, 1.0, 0, 0, NUM_STATE, crse_lev.geom);

	// Store the density change, for the gravity sync.

#ifdef SELF_GRAVITY
	int ilev = lev - crse_level - 1;

	if (do_grav && gravity->get_gravity_type() == "PoissonGrav" && gravity->NoSync() == 0) {
	    reg->Reflux(*drho[ilev], crse_lev.volume, 1.0, 0, Density, 1, crse_lev.geom);
	    amrex::average_down(*drho[ilev + 1], *drho[ilev], 0, 1, getLevel(lev).crse_ratio);
	}
#endif

	// Also update the coarse fluxes MultiFabs using the reflux data. This should only make
	// a difference if we re-evaluate the source terms later.

	Vector<std::unique_ptr<MultiFab> > temp_fluxes(3);

	if (update_sources_after_reflux) {

	    for (int i = 0; i < BL_SPACEDIM; ++i) {
		temp_fluxes[i].reset(new MultiFab(crse_lev.fluxes[i]->boxArray(),
						  crse_lev.fluxes[i]->DistributionMap(),
						  crse_lev.fluxes[i]->nComp(), crse_lev.fluxes[i]->nGrow()));
		temp_fluxes[i]->setVal(0.0);
	    }
	    for (OrientationIter fi; fi; ++fi) {
		const FabSet& fs = (*reg)[fi()];
		int idir = fi().coordDir();
		fs.copyTo(*temp_fluxes[idir], 0, 0, 0, temp_fluxes[idir]->nComp());
	    }
	    for (int i = 0; i < BL_SPACEDIM; ++i) {
		MultiFab::Add(*crse_lev.fluxes[i], *temp_fluxes[i], 0, 0, crse_lev.fluxes[i]->nComp(), 0);
                MultiFab::Add(*crse_lev.mass_fluxes[i], *temp_fluxes[i], Density, 0, 1, 0);
		temp_fluxes[i].reset();
	    }

	}

	// We no longer need the flux register data, so clear it out.

	reg->setVal(0.0);

#if (BL_SPACEDIM <= 2)
	if (!Geometry::IsCartesian()) {

	    reg = &getLevel(lev).pres_reg;

	    MultiFab dr(crse_lev.grids, crse_lev.dmap, 1, 0);
	    dr.setVal(crse_lev.geom.CellSize(0));

	    reg->ClearInternalBorders(crse_lev.geom);

	    reg->Reflux(crse_state, dr, 1.0, 0, Xmom, 1, crse_lev.geom);

	    if (update_sources_after_reflux) {

		temp_fluxes[0].reset(new MultiFab(crse_lev.P_radial.boxArray(),
						  crse_lev.P_radial.DistributionMap(),
						  crse_lev.P_radial.nComp(), crse_lev.P_radial.nGrow()));
		temp_fluxes[0]->setVal(0.0);

                for (OrientationIter fi; fi; ++fi)
		{
		    const FabSet& fs = (*reg)[fi()];
		    int idir = fi().coordDir();
		    if (idir == 0) {
			fs.copyTo(*temp_fluxes[idir], 0, 0, 0, temp_fluxes[idir]->nComp());
		    }
                }

		MultiFab::Add(crse_lev.P_radial, *temp_fluxes[0], 0, 0, crse_lev.P_radial.nComp(), 0);
		temp_fluxes[0].reset();

	    }

	    reg->setVal(0.0);

	}
#endif

#ifdef RADIATION

	// This follows the same logic as the pure hydro fluxes; see above for details.

	if (Radiation::rad_hydro_combined) {

	    reg = &getLevel(lev).rad_flux_reg;

	    reg->ClearInternalBorders(crse_lev.geom);

	    reg->Reflux(crse_lev.get_new_data(Rad_Type), crse_lev.volume, 1.0, 0, 0, Radiation::nGroups, crse_lev.geom);

	    if (update_sources_after_reflux) {

		for (int i = 0; i < BL_SPACEDIM; ++i) {
		    temp_fluxes[i].reset(new MultiFab(crse_lev.rad_fluxes[i]->boxArray(),
						      crse_lev.rad_fluxes[i]->DistributionMap(),
						      crse_lev.rad_fluxes[i]->nComp(), crse_lev.rad_fluxes[i]->nGrow()));
		    temp_fluxes[i]->setVal(0.0);
		}
		for (OrientationIter fi; fi; ++fi) {
		    const FabSet& fs = (*reg)[fi()];
		    int idir = fi().coordDir();
		    fs.copyTo(*temp_fluxes[idir], 0, 0, 0, temp_fluxes[idir]->nComp());
		}
		for (int i = 0; i < BL_SPACEDIM; ++i) {
		    MultiFab::Add(*crse_lev.rad_fluxes[i], *temp_fluxes[i], 0, 0, crse_lev.rad_fluxes[i]->nComp(), 0);
		    temp_fluxes[i].reset();
		}

	    }

	    reg->setVal(0.0);

	}

#endif

#ifdef SELF_GRAVITY
	if (do_grav && gravity->get_gravity_type() == "PoissonGrav" && gravity->NoSync() == 0)  {

	    reg = &getLevel(lev).phi_reg;

	    // Note that the scaling by the area here is corrected for by dividing by the
	    // cell volume in the reflux. In this way we get a discrete divergence that
	    // is analogous to the divergence of the flux in the hydrodynamics. See Equation
	    // 37 in the Castro I paper. The dimensions of dphi are therefore actually
	    // phi / cm**2, which makes it correct for the RHS of the Poisson equation.

	    for (int i = 0; i < BL_SPACEDIM; ++i) {
		reg->CrseInit(*(gravity->get_grad_phi_curr(lev-1)[i]), crse_lev.area[i], i, 0, 0, 1, -1.0);
		reg->FineAdd(*(gravity->get_grad_phi_curr(lev)[i]), fine_lev.area[i], i, 0, 0, 1, 1.0);
	    }

	    reg->Reflux(*dphi[ilev], crse_lev.volume, 1.0, 0, 0, 1, crse_lev.geom);

	    amrex::average_down(*dphi[ilev + 1], *dphi[ilev], 0, 1, getLevel(lev).crse_ratio);

	    reg->setVal(0.0);

	}
#endif

    }

    // Do the sync solve across all levels.

#ifdef SELF_GRAVITY
    if (do_grav && gravity->get_gravity_type() == "PoissonGrav" && gravity->NoSync() == 0)
	gravity->gravity_sync(crse_level, fine_level, amrex::GetVecOfPtrs(drho), amrex::GetVecOfPtrs(dphi));
#endif

    // Now subtract the new-time updates to the state data,
    // recompute it, and add it back. This corrects for the fact
    // that the new-time data was calculated the first time around
    // using a state that hadn't yet been refluxed. Note that this
    // needs to come after the gravity sync solve because the gravity
    // source depends on an up-to-date value of phi. We'll do this
    // on the fine level in addition to the coarser levels, because
    // global sources like gravity or source terms that rely on
    // ghost zone fills like diffusion depend on the data in the
    // coarser levels.

    if (update_sources_after_reflux) {

	for (int lev = fine_level; lev >= crse_level; --lev) {

            MultiFab& S_old = getLevel(lev).get_old_data(State_Type);
	    MultiFab& S_new = getLevel(lev).get_new_data(State_Type);
            MultiFab& source = getLevel(lev).get_new_data(Source_Type);
	    Real time = getLevel(lev).state[State_Type].curTime();
	    Real dt_advance = getLevel(lev).dt_advance; // Note that this may be shorter than the full timestep due to subcycling.
            Real dt_amr = parent->dtLevel(lev); // The full timestep expected by the Amr class.

            if (getLevel(lev).apply_sources()) {

                int is_new=1;
                getLevel(lev).apply_source_to_state(is_new, S_new, source, -dt_advance);
            }

            // Temporarily restore the last iteration's old data for the purposes of recalculating the corrector.
            // This is only necessary if we've done subcycles on that level.

            if (use_retry && dt_advance < dt_amr && getLevel(lev).keep_prev_state) {

                for (int k = 0; k < num_state_type; k++) {

                    if (getLevel(lev).prev_state[k]->hasOldData()) {

                        // Use the new-time data as a temporary buffer. Ideally this would be done
                        // as a pointer swap, but we cannot assume that the distribution mapping
                        // is the same between the current state and the original state, due to
                        // possible regrids that have occurred in between.

                        MultiFab& old = getLevel(lev).get_old_data(k);
                        MultiFab::Copy(getLevel(lev).prev_state[k]->newData(), old, 0, 0, old.nComp(), old.nGrow());
                        MultiFab::Copy(old, getLevel(lev).prev_state[k]->oldData(), 0, 0, old.nComp(), old.nGrow());

                        getLevel(lev).state[k].setTimeLevel(time, dt_advance, 0.0);
                        getLevel(lev).prev_state[k]->setTimeLevel(time, dt_amr, 0.0);

                    }

                }

            }

            if (getLevel(lev).apply_sources()) {

                getLevel(lev).do_new_sources(source, S_old, S_new, time, dt_advance);

                int is_new=1;
                getLevel(lev).apply_source_to_state(is_new, S_new, source, dt_advance);

            }

            if (use_retry && dt_advance < dt_amr && getLevel(lev).keep_prev_state) {

                for (int k = 0; k < num_state_type; k++) {

                    if (getLevel(lev).prev_state[k]->hasOldData()) {

                        // Now retrieve the original old time data.

                        MultiFab& old = getLevel(lev).get_old_data(k);
                        MultiFab::Copy(old, getLevel(lev).prev_state[k]->newData(), 0, 0, old.nComp(), old.nGrow());

                        getLevel(lev).state[k].setTimeLevel(time, dt_amr, 0.0);
                        getLevel(lev).prev_state[k]->setTimeLevel(time, dt_advance, 0.0);

                    }

                }

                // Now deallocate the old data, it is no longer needed.

                if (lev == 0 || lev > level) {
                    amrex::FillNull(getLevel(lev).prev_state);
                    getLevel(lev).keep_prev_state = false;
                }

            }

	}

    }

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

#ifdef SELF_GRAVITY
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

#ifdef SDC
#ifdef REACTIONS
  avgDown(SDC_React_Type);
#endif
#endif

#ifdef RADIATION
  if (do_radiation) {
    avgDown(Rad_Type);
  }
#endif

  if (rebalanced_splitting) {
    avgDown(Balance_Type);
  }

}

void
Castro::normalize_species (MultiFab& S_new, int ng)
{

    BL_PROFILE("Castro::normalize_species()");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.growntilebox(ng);

#pragma gpu
       ca_normalize_species(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                            BL_TO_FORTRAN_ANYD(S_new[mfi]));
    }
}

void
Castro::enforce_consistent_e (MultiFab& S)
{

    BL_PROFILE("Castro::enforce_consistent_e()");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S,true); mfi.isValid(); ++mfi)
    {
        const Box& box     = mfi.tilebox();
        const int* lo      = box.loVect();
        const int* hi      = box.hiVect();

        ca_enforce_consistent_e(ARLIM_3D(lo), ARLIM_3D(hi), BL_TO_FORTRAN_ANYD(S[mfi]));
    }
}

Real
Castro::enforce_min_density (MultiFab& S_old, MultiFab& S_new, int ng)
{

    BL_PROFILE("Castro::enforce_min_density()");

    // This routine sets the density in S_new to be larger than the density floor.
    // Note that it will operate everywhere on S_new, including ghost zones.
    // S_old is present so that, after the hydro call, we know what the old density
    // was so that we have a reference for comparison. If you are calling it elsewhere
    // and there's no meaningful reference state, just pass in the same MultiFab twice.

    // The return value is the the negative fractional change in the state that has the
    // largest magnitude. If there is no reference state, this is meaningless.

    Real dens_change = 1.e0;

    MultiFab reset_source;

    if (print_update_diagnostics)
    {

	// Before we do anything, make a copy of the state.

	reset_source.define(S_new.boxArray(), S_new.DistributionMap(), S_new.nComp(), 0);

	MultiFab::Copy(reset_source, S_new, 0, 0, S_new.nComp(), 0);

    }

#ifdef _OPENMP
#pragma omp parallel reduction(min:dens_change)
#endif
    for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi) {

	const Box& bx = mfi.growntilebox(ng);

	const FArrayBox& stateold = S_old[mfi];
	FArrayBox& statenew = S_new[mfi];
	const FArrayBox& vol      = volume[mfi];

	ca_enforce_minimum_density
            (AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
             BL_TO_FORTRAN_ANYD(stateold),
             BL_TO_FORTRAN_ANYD(statenew),
             BL_TO_FORTRAN_ANYD(vol),
             &dens_change,
             verbose);

    }

    if (print_update_diagnostics)
    {

	// Evaluate what the effective reset source was.

	MultiFab::Subtract(reset_source, S_new, 0, 0, S_old.nComp(), 0);

	bool local = true;
	Vector<Real> reset_update = evaluate_source_change(reset_source, 1.0, local);

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
	    ParallelDescriptor::ReduceRealSum(reset_update.dataPtr(), reset_update.size(), ParallelDescriptor::IOProcessorNumber());

	    if (ParallelDescriptor::IOProcessor()) {
		if (std::abs(reset_update[0]) != 0.0) {
		    std::cout << std::endl << "  Contributions to the state from negative density resets:" << std::endl;

		    print_source_change(reset_update);
		}
	    }

#ifdef BL_LAZY
        });
#endif

    }

    return dens_change;

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

    amrex::average_down(S_fine, S_crse,
			 fgeom, cgeom,
			 0, S_fine.nComp(), fine_ratio);
}

void
Castro::allocOldData ()
{
    for (int k = 0; k < num_state_type; k++)
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

    ca_set_amr_info(level, -1, -1, -1.0, -1.0);

    Real t = time;

    // If we are forcing a post-timestep regrid,
    // note that we need to use the new time here,
    // not the old time.

    if (post_step_regrid)
	t = get_state_data(State_Type).curTime();

    // Apply each of the specified tagging functions.

    for (int j = 0; j < num_err_list_default; j++)
	apply_tagging_func(tags, clearval, tagval, t, j);

    // Now apply the user-specified tagging functions.
    // Include problem-specific hooks before and after.

    problem_pre_tagging_hook(tags, clearval, tagval, t);

    for (int j = num_err_list_default; j < err_list.size(); j++)
        apply_tagging_func(tags, clearval, tagval, t, j);

    // Now we'll tag any user-specified zones using the full state array.

    apply_problem_tags(tags, clearval, tagval, time);

    problem_post_tagging_hook(tags, clearval, tagval, t);
}



void
Castro::apply_problem_tags (TagBoxArray& tags,
                            int          clearval,
                            int          tagval,
                            Real         time)
{

    BL_PROFILE("Castro::apply_problem_tags()");

    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();

    MultiFab& S_new = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Vector<int>  itags;

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

#ifdef AMREX_DIMENSION_AGNOSTIC
	    set_problem_tags(ARLIM_3D(tilebx.loVect()), ARLIM_3D(tilebx.hiVect()),
                             tptr, ARLIM_3D(tlo), ARLIM_3D(thi),
			     BL_TO_FORTRAN_ANYD(S_new[mfi]),
			     &tagval, &clearval,
			     ZFILL(dx), ZFILL(prob_lo), &time, &level);
#else
	    set_problem_tags(tilebx.loVect(), tilebx.hiVect(),
                             tptr, ARLIM(tlo), ARLIM(thi),
			     BL_TO_FORTRAN(S_new[mfi]),
			     &tagval, &clearval,
		             dx, prob_lo, &time, &level);
#endif

	    //
	    // Now update the tags in the TagBox.
	    //
            tagfab.tags_and_untags(itags, tilebx);
	}
    }

}



void
Castro::apply_tagging_func(TagBoxArray& tags, int clearval, int tagval, Real time, int j)
{

    BL_PROFILE("Castro::apply_tagging_func()");

    const int*  domain_lo = geom.Domain().loVect();
    const int*  domain_hi = geom.Domain().hiVect();
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();

    auto mf = derive(err_list[j].name(), time, err_list[j].nGrow());

    BL_ASSERT(mf);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Vector<int>  itags;

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

}



std::unique_ptr<MultiFab>
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

    ca_setgroup(ig);
  }
#endif

#ifdef AMREX_PARTICLES
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

    ca_setgroup(ig);
  }
#endif

    AmrLevel::derive(name,time,mf,dcomp);
}

void
Castro::amrinfo_init ()
{
   ca_amrinfo_init();
}

void
Castro::amrinfo_finalize()
{
   ca_amrinfo_finalize();
}

void
Castro::network_init ()
{
   ca_network_init();
}

void
Castro::network_finalize ()
{
   ca_network_finalize();
}

void
Castro::eos_finalize ()
{
   ca_eos_finalize();
}

void
Castro::extern_init ()
{
  // initialize the external runtime parameters -- these will
  // live in the probin

  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "reading extern runtime parameters ..." << std::endl;
  }

  const int probin_file_length = probin_file.length();
  Vector<int> probin_file_name(probin_file_length);

  for (int i = 0; i < probin_file_length; i++)
    probin_file_name[i] = probin_file[i];

  ca_extern_init(probin_file_name.dataPtr(),&probin_file_length);
}

void
Castro::reset_internal_energy(MultiFab& S_new)
{

    BL_PROFILE("Castro::reset_internal_energy()");

    MultiFab old_state;

    // Make a copy of the state so we can evaluate how much changed.

    if (print_update_diagnostics)
    {
	old_state.define(S_new.boxArray(), S_new.DistributionMap(), S_new.nComp(), 0);
        MultiFab::Copy(old_state, S_new, 0, 0, S_new.nComp(), 0);
    }

    int ng = S_new.nGrow();

    // Ensure (rho e) isn't too small or negative
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(ng);

#pragma gpu
        ca_reset_internal_e(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
			    BL_TO_FORTRAN_ANYD(S_new[mfi]),
			    print_fortran_warnings);
    }

    // Flush Fortran output

    if (verbose)
      flush_output();

    if (print_update_diagnostics)
    {
	// Evaluate what the effective reset source was.

	MultiFab reset_source(S_new.boxArray(), S_new.DistributionMap(), S_new.nComp(), 0);

	MultiFab::Copy(reset_source, S_new, 0, 0, S_new.nComp(), 0);

	MultiFab::Subtract(reset_source, old_state, 0, 0, old_state.nComp(), 0);

	bool local = true;
	Vector<Real> reset_update = evaluate_source_change(reset_source, 1.0, local);

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
	    ParallelDescriptor::ReduceRealSum(reset_update.dataPtr(), reset_update.size(), ParallelDescriptor::IOProcessorNumber());

	    if (ParallelDescriptor::IOProcessor()) {
		if (std::abs(reset_update[Eint]) != 0.0) {
		    std::cout << std::endl << "  Contributions to the state from negative energy resets:" << std::endl;

		    print_source_change(reset_update);
		}
	    }

#ifdef BL_LAZY
	});
#endif
    }
}

void
Castro::computeTemp(MultiFab& State, int ng)
{

  BL_PROFILE("Castro::computeTemp()");

  reset_internal_energy(State);

#ifdef RADIATION
  FArrayBox temp;
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(State,true); mfi.isValid(); ++mfi)
    {
      const Box& bx = mfi.growntilebox(ng);

#ifdef RADIATION
      if (Radiation::do_real_eos == 0) {
	temp.resize(bx);
	temp.copy(State[mfi],bx,Eint,bx,0,1);

	ca_compute_temp_given_cv
	  (bx.loVect(), bx.hiVect(),
	   BL_TO_FORTRAN(temp),
	   BL_TO_FORTRAN(State[mfi]),
	   &Radiation::const_c_v, &Radiation::c_v_exp_m, &Radiation::c_v_exp_n);

	State[mfi].copy(temp,bx,0,bx,Temp,1);
      } else {
#endif
#pragma gpu
	ca_compute_temp(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
			BL_TO_FORTRAN_ANYD(State[mfi]));
#ifdef RADIATION
      }
#endif
    }
}



void
Castro::apply_source_term_predictor()
{

    BL_PROFILE("Castro::apply_source_term_predictor()");

    // Optionally predict the source terms to t + dt/2,
    // which is the time-level n+1/2 value, To do this we use a
    // lagged predictor estimate: dS/dt_n = (S_n - S_{n-1}) / dt, so
    // S_{n+1/2} = S_n + (dt / 2) * dS/dt_n. We'll add the S_n
    // terms later; now we add the second term. We defer the
    // multiplication by dt / 2 to initialize_do_advance since
    // for a retry we may not yet know at this point what the
    // advance timestep is.

    // Note that since for dS/dt we want (S^{n+1} - S^{n}) / dt,
    // we only need to take twice the new-time source term, since in
    // the predictor-corrector approach, the new-time source term is
    // 1/2 * S^{n+1} - 1/2 * S^{n}. This is untrue in general for the
    // non-momentum sources, so for safety we'll only apply it to the
    // momentum sources.

    // Note that if the old data doesn't exist yet (e.g. it is
    // the first step of the simulation) FillPatch will just
    // return the new data, so this is a valid operation and
    // the result will be zero, so there is no source term
    // prediction in the first step.

    const Real old_time = get_state_data(Source_Type).prevTime();
    const Real new_time = get_state_data(Source_Type).curTime();

    const Real dt_old = new_time - old_time;

    AmrLevel::FillPatchAdd(*this, sources_for_hydro, NUM_GROW, new_time, Source_Type, Xmom, 3, Xmom);

    sources_for_hydro.mult(2.0 / dt_old, NUM_GROW);

}



void
Castro::swap_state_time_levels(const Real dt)
{

    for (int k = 0; k < num_state_type; k++) {

	// The following is a hack to make sure that we only
	// ever have new data for certain state types that only
	// ever need new time data; by doing a swap now, we'll
	// guarantee that allocOldData() does nothing. We do
	// this because we never need the old data, so we
	// don't want to allocate memory for it.

#ifdef SDC
#ifdef REACTIONS
        if (k == SDC_React_Type)
            state[k].swapTimeLevels(0.0);
#endif
#endif

        if (rebalanced_splitting) {
            if (k == Balance_Type) {
                state[k].swapTimeLevels(0.0);
            }
        }

        state[k].allocOldData();

        state[k].swapTimeLevels(dt);

    }

}



#ifdef SELF_GRAVITY
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

   BL_PROFILE("Castro::make_radial_data()");

 // We only call this for level = 0
   BL_ASSERT(level == 0);

   int numpts_1d = get_numpts();

   Vector<Real> radial_vol(numpts_1d,0);

   const Real* dx = geom.CellSize();
   Real  dr = dx[0];

   if (is_new == 1) {
      MultiFab& S = get_new_data(State_Type);
      const int nc = S.nComp();
      Vector<Real> radial_state(numpts_1d*nc,0);
      for (MFIter mfi(S); mfi.isValid(); ++mfi)
      {
         Box bx(mfi.validbox());
         ca_compute_avgstate(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),ZFILL(dx),&dr,&nc,
			     BL_TO_FORTRAN_ANYD(     S[mfi]),radial_state.dataPtr(),
			     BL_TO_FORTRAN_ANYD(volume[mfi]),radial_vol.dataPtr(),
			     ZFILL(geom.ProbLo()),&numpts_1d);
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

      Vector<Real> radial_state_short(np_max*nc,0);

      for (int i = 0; i < np_max; i++) {
         for (int j = 0; j < nc; j++) {
           radial_state_short[nc*i+j] = radial_state[nc*i+j];
         }
      }

      const Real new_time = state[State_Type].curTime();
      set_new_outflow_data(radial_state_short.dataPtr(),&new_time,&np_max,&nc);
   }
   else
   {
      MultiFab& S = get_old_data(State_Type);
      const int nc = S.nComp();
      Vector<Real> radial_state(numpts_1d*nc,0);
      for (MFIter mfi(S); mfi.isValid(); ++mfi)
      {
         Box bx(mfi.validbox());
         ca_compute_avgstate(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),ZFILL(dx),&dr,&nc,
			     BL_TO_FORTRAN_ANYD(     S[mfi]),radial_state.dataPtr(),
			     BL_TO_FORTRAN_ANYD(volume[mfi]),radial_vol.dataPtr(),
			     ZFILL(geom.ProbLo()),&numpts_1d);
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

      Vector<Real> radial_state_short(np_max*nc,0);

      for (int i = 0; i < np_max; i++) {
         for (int j = 0; j < nc; j++) {
           radial_state_short[nc*i+j] = radial_state[nc*i+j];
         }
      }

      const Real old_time = state[State_Type].prevTime();
      set_old_outflow_data(radial_state_short.dataPtr(),&old_time,&np_max,&nc);
   }

#endif
}

void
Castro::define_new_center(MultiFab& S, Real time)
{
    BL_PROFILE("Castro::define_new_center()");

    Real center[3];
    const Real* dx = geom.CellSize();

    IntVect max_index = S.maxIndex(Density,0);
    Box bx(max_index,max_index);
    bx.grow(1);
    BoxArray ba(bx);
    int owner = ParallelDescriptor::IOProcessorNumber();
    DistributionMapping dm { Vector<int>(1,owner) };
    MultiFab mf(ba,dm,1,0);

    // Define a cube 3-on-a-side around the point with the maximum density
    FillPatch(*this,mf,0,time,State_Type,Density,1);

    int mi[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) mi[i] = max_index[i];

    // Find the position of the "center" by interpolating from data at cell centers
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        ca_find_center(mf[mfi].dataPtr(),&center[0],ARLIM_3D(mi),ZFILL(dx),ZFILL(geom.ProbLo()));
    }
    // Now broadcast to everyone else.
    ParallelDescriptor::Bcast(&center[0], BL_SPACEDIM, owner);

    // Make sure if R-Z that center stays exactly on axis
    if ( Geometry::IsRZ() ) center[0] = 0;

    ca_set_center(ZFILL(center));
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

       Real center[3];
       ca_get_center(center);

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


MultiFab&
Castro::build_fine_mask()
{
    BL_PROFILE("Castro::build_fine_mask()");

    BL_ASSERT(level > 0); // because we are building a mask for the coarser level

    if (!fine_mask.empty()) return fine_mask;

    BoxArray baf = parent->boxArray(level);
    baf.coarsen(crse_ratio);

    const BoxArray& bac = parent->boxArray(level-1);
    const DistributionMapping& dmc = parent->DistributionMap(level-1);
    fine_mask.define(bac,dmc,1,0);
    fine_mask.setVal(1.0);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(fine_mask); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = fine_mask[mfi];

	const std::vector< std::pair<int,Box> >& isects = baf.intersections(fab.box());

	for (int ii = 0; ii < isects.size(); ++ii)
	{
	    fab.setVal(0.0,isects[ii].second,0);
	}
    }

    return fine_mask;
}

iMultiFab&
Castro::build_interior_boundary_mask (int ng)
{
    BL_PROFILE("Castro::build_interior_boundary_mask()");

    for (int i = 0; i < ib_mask.size(); ++i)
    {
	if (ib_mask[i]->nGrow() == ng) {
	    return *ib_mask[i];
	}
    }

    //  If we got here, we need to build a new one
    ib_mask.push_back(std::unique_ptr<iMultiFab>(new iMultiFab(grids, dmap, 1, ng)));

    iMultiFab& imf = *ib_mask.back();

    int ghost_covered_by_valid = 0;
    int other_cells = 1; // uncovered ghost, valid, and outside domain cells are set to 1

    imf.BuildMask(geom.Domain(), geom.periodicity(),
		  ghost_covered_by_valid, other_cells, other_cells, other_cells);

    return imf;
}

// Fill a version of the state with ng ghost zones from the state data.

void
Castro::expand_state(MultiFab& S, Real time, int iclean, int ng)
{
  BL_PROFILE("Castro::expand_state()");

  // S is the multifab we are filling with State_Type StateData,
  // including a ghost cell fill at the end.  Before we do the fill,
  // we do a clean_state on the State_Type data.  iclean = 0 means
  // clean the "old" time data, iclean = 1 means clean the "new time
  // data", and iclean = 2 means clean both the old and new time data
  // (in case we are interpolating in time).  For any other value, no
  // cleaning is done

  // note: we don't clean ghost cells here, since we are doing a fill
  // right afterwards

  if (iclean == 0) {
    clean_state(iclean, 0);
  } else if (iclean == 1) {
    clean_state(iclean, 0);
  } else if (iclean == 2) {
    clean_state(0, 0);
    clean_state(1, 0);
  }

  BL_ASSERT(S.nGrow() >= ng);

  AmrLevel::FillPatch(*this, S, ng, time, State_Type, 0, NUM_STATE);
}


void
Castro::check_for_nan(MultiFab& state, int check_ghost)
{
  BL_PROFILE("Castro::check_for_nan()");

  int ng = 0;
  if (check_ghost == 1) {
    ng = state.nComp();
  }

  if (state.contains_nan(Density,state.nComp(),ng,true))
    {
      for (int i = 0; i < state.nComp(); i++)
        {
	  if (state.contains_nan(Density + i, 1, ng, true))
            {
	      std::string abort_string = std::string("State has NaNs in the ") + desc_lst[State_Type].name(i) + std::string(" component::check_for_nan()");
	      amrex::Abort(abort_string.c_str());
            }
        }
    }
}

// Convert a MultiFab with conservative state data u to a primitive MultiFab q.
#ifdef SDC
void
Castro::cons_to_prim(MultiFab& u, MultiFab& q, MultiFab& qaux)
{

    BL_PROFILE("Castro::cons_to_prim()");

    BL_ASSERT(u.nComp() == NUM_STATE);
    BL_ASSERT(q.nComp() == NQ);
    BL_ASSERT(u.nGrow() >= q.nGrow());

    int ng = q.nGrow();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(u, true); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.growntilebox(ng);

	ca_ctoprim(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
		   u[mfi].dataPtr(), ARLIM_3D(u[mfi].loVect()), ARLIM_3D(u[mfi].hiVect()),
		   q[mfi].dataPtr(), ARLIM_3D(q[mfi].loVect()), ARLIM_3D(q[mfi].hiVect()),
		   qaux[mfi].dataPtr(), ARLIM_3D(qaux[mfi].loVect()), ARLIM_3D(qaux[mfi].hiVect()));

    }

}
#endif


// Given State_Type state data, perform a number of cleaning steps to make
// sure the data is sensible. The return value is the same as the return
// value of enforce_min_density.

Real
Castro::clean_state(MultiFab& state) {

    BL_PROFILE("Castro::clean_state(state)");

    // Enforce a minimum density.

    MultiFab temp_state(state.boxArray(), state.DistributionMap(), state.nComp(), state.nGrow());

    MultiFab::Copy(temp_state, state, 0, 0, state.nComp(), state.nGrow());


#ifndef AMREX_USE_CUDA
    Real frac_change = enforce_min_density(temp_state, state, state.nGrow());
#else
    Real frac_change = 1.e200;
#endif

    // Ensure all species are normalized.

    normalize_species(state, state.nGrow());

    // Sync the linear and hybrid momenta.

#ifdef HYBRID_MOMENTUM
    hybrid_sync(state);
#endif

    // Compute the temperature (note that this will also reset
    // the internal energy for consistency with the total energy).

    computeTemp(state, state.nGrow());

    return frac_change;

}


Real
Castro::clean_state(int is_new, int ng) {

  BL_PROFILE("Castro::clean_state(is_new, ng)");

  // this is the "preferred" clean_state interface -- it will work
  // directly on the StateData.  is_new=0 means the old data is used,
  // is_new=1 means the new data is used.


  MultiFab& state = is_new == 1 ? get_new_data(State_Type) : get_old_data(State_Type);

  MultiFab temp_state(state.boxArray(), state.DistributionMap(), state.nComp(), ng);

  MultiFab::Copy(temp_state, state, 0, 0, state.nComp(), ng);

  Real frac_change = clean_state(is_new, temp_state, ng);

  return frac_change;

}


Real
Castro::clean_state(int is_new, MultiFab& state_old, int ng) {

    BL_PROFILE("Castro::clean_state(is_new, state_old, ng)");

  MultiFab& state = is_new == 1 ? get_new_data(State_Type) : get_old_data(State_Type);

  // Enforce a minimum density.
#ifndef AMREX_USE_CUDA
    Real frac_change = enforce_min_density(state_old, state, ng);
#else
    Real frac_change = 1.e200;
#endif


  // Ensure all species are normalized.
  normalize_species(state, ng);

  // Sync the linear and hybrid momenta.
#ifdef HYBRID_MOMENTUM
  hybrid_sync(state);
#endif

  // Compute the temperature (note that this will also reset
  // the internal energy for consistency with the total energy).
  computeTemp(state, ng);

  return frac_change;

}
