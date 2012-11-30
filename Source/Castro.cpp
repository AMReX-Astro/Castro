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

#include "buildInfo.H"

static int  sum_interval = -1;
static Real fixed_dt     = -1.0;
static Real initial_dt   = -1.0;
static Real dt_cutoff    = 0.0;

bool         Castro::dump_old      = false;

int          Castro::verbose       = 0;
Real         Castro::cfl           = 0.8;
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
int          Castro::LastSpec      = -1;

int          Castro::NumAux        = 0;
int          Castro::FirstAux      = -1;
int          Castro::LastAux       = -1;

int          Castro::NumAdv        = 0;
int          Castro::FirstAdv      = -1;
int          Castro::LastAdv       = -1;

Real         Castro::difmag        = 0.1;
Real         Castro::small_dens    = -1.e200;
Real         Castro::small_temp    = -1.e200;
Real         Castro::small_pres    = -1.e200;
Real         Castro::gamma         =  0.0;

int          Castro::do_hydro = -1;
int          Castro::do_react = -1;
int          Castro::add_ext_src = 0;

#ifdef POINTMASS
Real         Castro::point_mass    = 0.0;
#endif

#ifdef GRAVITY
Gravity*     Castro::gravity  = 0;
int          Castro::plot_phiGrav = 0;
int          Castro::do_grav  = -1;
#endif

#ifdef DIFFUSION
Diffusion*    Castro::diffusion  = 0;
int           Castro::diffuse_temp = 0;
#endif

#ifdef RADIATION
int          Castro::do_radiation = -1;
Radiation*   Castro::radiation = 0;
#endif

#ifdef ROTATION
int          Castro::do_rotation = -1;
Real         Castro::rotational_frequency = 0.0;
#endif

int          Castro::grown_factor = 1;
int          Castro::star_at_center = -1;
int          Castro::moving_center = 0;
int          Castro::allow_untagging = 0;
int          Castro::normalize_species = 0;
int          Castro::fix_mass_flux = 0;
int          Castro::allow_negative_energy = 1;
int          Castro::do_special_tagging = 0;
int          Castro::ppm_type = 1;
int          Castro::use_colglaz = 0;
int          Castro::spherical_star = 0;
int          Castro::do_sponge  = 0;

int          Castro::print_fortran_warnings  = 0;
int          Castro::print_energy_diagnostics  = 0;

int          Castro::show_center_of_mass = 0;

#ifdef SGS
Real         Castro::sum_turb_src = 0.0;
#endif

std::string  Castro::job_name = "";


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
    pp.get("init_shrink",init_shrink);
    pp.get("cfl",cfl);
    pp.query("change_max",change_max);
    pp.query("fixed_dt",fixed_dt);
    pp.query("initial_dt",initial_dt);
    pp.query("sum_interval",sum_interval);
    pp.query("do_reflux",do_reflux);
    do_reflux = (do_reflux ? 1 : 0);
    pp.get("dt_cutoff",dt_cutoff);

    pp.query("dump_old",dump_old);

    pp.query("difmag",difmag);
    pp.query("small_dens",small_dens);
    pp.query("small_temp",small_temp);
    pp.query("small_pres",small_pres);
    pp.query("gamma",gamma);

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
    pp.get("do_react",do_react);
    pp.query("add_ext_src",add_ext_src);

#ifdef DIFFUSION
    pp.query("diffuse_temp",diffuse_temp);
#endif

    pp.query("grown_factor",grown_factor);
    if (grown_factor < 1) 
       BoxLib::Error("grown_factor must be integer >= 1");

    pp.query("star_at_center",star_at_center);

    pp.query("moving_center",moving_center);

    pp.query("allow_untagging",allow_untagging);
    pp.query("normalize_species",normalize_species);
    pp.query("fix_mass_flux",fix_mass_flux);
    pp.query("allow_negative_energy",allow_negative_energy);
    pp.query("do_special_tagging",do_special_tagging);
    pp.query("ppm_type", ppm_type);
    pp.query("use_colglaz",use_colglaz);
    pp.query("spherical_star",spherical_star);
    pp.query("do_sponge",do_sponge);

    pp.query("show_center_of_mass",show_center_of_mass);
    pp.query("print_energy_diagnostics",print_energy_diagnostics);
    pp.query("print_fortran_warnings",print_fortran_warnings);

    if (use_colglaz == 1 && BL_SPACEDIM > 1) 
    {
        std::cerr << "use_colglaz = 1 only implemented for 1-d \n";
        BoxLib::Error();
    }

    // Make sure not to call refluxing if we're not actually doing any hydro.
    if (do_hydro == 0) do_reflux = 0;

#ifdef GRAVITY
    pp.get("do_grav",do_grav);
    pp.query("plot_phiGrav",plot_phiGrav);

#if (BL_SPACEDIM == 1)
    if (do_grav && !Geometry::IsSPHERICAL()) {
        std::cerr << "ERROR:Castro::Gravity in 1D assumes that the coordinate system is spherical\n";
        BoxLib::Error();
    }
#endif
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
    if (do_rotation) pp.get("rotational_frequency",rotational_frequency);
#if (BL_SPACEDIM == 1)
      if (do_rotation) {
	std::cerr << "ERROR:Castro::Rotation not implemented in 1d\n";
	BoxLib::Error();
      }
#endif
#endif

   pp.query("job_name",job_name);  
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
        flux_reg = new FluxRegister(grids,crse_ratio,level,NUM_STATE);

#ifdef SGS
    sgs_flux_reg = 0;
    if (level > 0 && do_reflux)
        sgs_flux_reg = new FluxRegister(grids,crse_ratio,level,NUM_STATE);
#endif

#ifdef RADIATION    
    rad_flux_reg = 0;
    if (Radiation::rad_hydro_combined && level > 0 && do_reflux) {
      rad_flux_reg = new FluxRegister(grids,crse_ratio,level,Radiation::nGroups);
    }
#endif

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
         int numpts_1d = get_numpts(level, geom.Domain());
         gravity->set_numpts_in_gravity(numpts_1d);
      }
#endif

      gravity->install_level(level,this,volume,area);

      if (verbose && level == 0 &&  ParallelDescriptor::IOProcessor()) 
         std::cout << "Setting the gravity type to " << gravity->get_gravity_type() << std::endl;
   }
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
}

void
Castro::restart (Amr&     papa,
                 istream& is,
                 bool     bReadSpecial)
{
    AmrLevel::restart(papa,is,bReadSpecial);

    buildMetrics();

    BL_ASSERT(flux_reg == 0);
    if (level > 0 && do_reflux)
        flux_reg = new FluxRegister(grids,crse_ratio,level,NUM_STATE);

#ifdef SGS
    BL_ASSERT(sgs_flux_reg == 0);
    if (level > 0 && do_reflux)
        sgs_flux_reg = new FluxRegister(grids,crse_ratio,level,NUM_STATE);
#endif

#ifdef RADIATION
    BL_ASSERT(rad_flux_reg == 0);
    if (Radiation::rad_hydro_combined && level > 0 && do_reflux)
      rad_flux_reg = new FluxRegister(grids,crse_ratio,level,Radiation::nGroups);
#endif

    const Real* dx  = geom.CellSize();

    if ( (grown_factor > 1) && (parent->maxLevel() < 1) )
    {
       std::cout << "grown_factor is " << grown_factor << std::endl;
       std::cout << "max_level is " << parent->maxLevel() << std::endl;
       BoxLib::Error("Must have max_level > 0 if doing special restart with grown_factor");
    }

    if (grown_factor > 1 && level == 0)
    {
       if (verbose && ParallelDescriptor::IOProcessor())
          std::cout << "Doing special restart with grown_factor " << grown_factor << std::endl;

       MultiFab& S_new = get_new_data(State_Type);
       Real cur_time   = state[State_Type].curTime();

       Box orig_domain;
       if (star_at_center == 0) {
          orig_domain = BoxLib::coarsen(geom.Domain(),grown_factor);
       } else if (star_at_center == 1) {

          Box domain(geom.Domain());
          int d,lo=0,hi=0;
          if (Geometry::IsRZ()) {
             if (grown_factor != 2) 
                BoxLib::Abort("Must have grown_factor = 2");

             d = 0;
             int dlen =  domain.size()[d];
             lo = 0;
             hi = dlen/2;
             orig_domain.setSmall(d,lo);
             orig_domain.setBig(d,hi);

             d = 1;
             dlen =  domain.size()[d];
             lo =   dlen/4    ;
             hi = 3*dlen/4 - 1;
             orig_domain.setSmall(d,lo);
             orig_domain.setBig(d,hi);

          } else {
             for (int d = 0; d < BL_SPACEDIM; d++) 
             {
                int dlen =  domain.size()[d];
                if (grown_factor == 2) {
                   lo =   dlen/4    ;
                   hi = 3*dlen/4 - 1;
                } else if (grown_factor == 3) {
                   lo =   (dlen)/3    ;
                   hi = 2*(dlen)/3 - 1;
                } else { 
                   BoxLib::Abort("Must have grown_factor = 2 or 3");
                }
                orig_domain.setSmall(d,lo);
                orig_domain.setBig(d,hi);
             }
          }
       } else {
          if (ParallelDescriptor::IOProcessor())
             std::cout << "... invalid value of star_at_center: " << star_at_center << std::endl;
          BoxLib::Abort();
       }

       int ns = NUM_STATE;

       for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
       {
           RealBox    gridloc = RealBox(grids[mfi.index()],geom.CellSize(),geom.ProbLo()); 
           const Box& bx      = mfi.validbox();
           const int* lo      = bx.loVect();
           const int* hi      = bx.hiVect();

           if (! orig_domain.contains(bx)) {
              BL_FORT_PROC_CALL(CA_INITDATA,ca_initdata)
                (level, cur_time, lo, hi, ns,
                 BL_TO_FORTRAN(S_new[mfi]), dx,
                 gridloc.lo(), gridloc.hi());
           }
       }
    }

    if (grown_factor > 1 && level == 1)
        getLevel(0).avgDown();

#if (BL_SPACEDIM > 1)
    if ( (level == 0) && (spherical_star == 1) ) {
       MultiFab& S_new = get_new_data(State_Type);
       int nc = S_new.nComp();
       int n1d = get_numpts(level,geom.Domain());
       BL_FORT_PROC_CALL(ALLOCATE_OUTFLOW_DATA,allocate_outflow_data)(&n1d,&nc);
       int is_new = 1;
       make_radial_data(is_new);
    }
#endif

#ifdef GRAVITY
    if (do_grav && level == 0) {
       BL_ASSERT(gravity == 0);
       gravity = new Gravity(parent,parent->finestLevel(),&phys_bc,Density);
    }
#endif

#ifdef DIFFUSION
    if (level == 0) {
       BL_ASSERT(diffusion == 0);
       diffusion = new Diffusion(parent,&phys_bc);
    }
#endif

#ifdef RADIATION
    if (do_radiation) {
      if (radiation == 0) {
        // radiation is a static object, only alloc if not already there
        int rad_restart = 1; // disables quasi-steady initialization
        radiation = new Radiation(parent, this, rad_restart);
      }
      radiation->regrid(level, grids);
      radiation->restart(level, parent->theRestartFile(), is);
    }
#endif
}

void
Castro::checkPoint(const std::string& dir,
                   std::ostream&  os,
                   VisMF::How     how,
                   bool dump_old_default)
{
  AmrLevel::checkPoint(dir, os, how, dump_old);

#ifdef RADIATION
  if (do_radiation) {
    radiation->checkPoint(level, dir, os, how);
  }
#endif

#ifdef PARTICLES
  ParticleCheckPoint(dir);
#endif
}

std::string
Castro::thePlotFileType () const
{
    //
    // Increment this whenever the writePlotFile() format changes.
    //
    static const std::string the_plot_file_type("HyperCLaw-V1.1");

    return the_plot_file_type;
}

void
Castro::setPlotVariables ()
{
  AmrLevel::setPlotVariables();

#ifdef RADIATION
  if (Radiation::nNeutrinoSpecies > 0 &&
      Radiation::plot_neutrino_group_energies_total == 0) {
    char rad_name[10];
    for (int j = 0; j < Radiation::nNeutrinoSpecies; j++) {
      for (int i = 0; i < Radiation::nNeutrinoGroups[j]; i++) {
        sprintf(rad_name, "rads%dg%d", j, i);
        parent->deleteStatePlotVar(rad_name);
      }
    }
  }
#endif

  ParmParse pp("castro");

  bool plot_X;

  if (pp.query("plot_X",plot_X))
  {
      if (plot_X)
      {
          //
	  // Get the number of species from the network model.
          //
	  BL_FORT_PROC_CALL(GET_NUM_SPEC, get_num_spec)(&NumSpec);
          //
	  // Get the species names from the network model.
          //
	  for (int i = 0; i < NumSpec; i++)
          {
              int len = 20;
              Array<int> int_spec_names(len);
              //
              // This call return the actual length of each string in "len"
              //
              BL_FORT_PROC_CALL(GET_SPEC_NAMES, get_spec_names)
                  (int_spec_names.dataPtr(),&i,&len);
              char* spec_name = new char[len+1];
              for (int j = 0; j < len; j++) 
                  spec_name[j] = int_spec_names[j];
              spec_name[len] = '\0';
	      string spec_string = "X(";
              spec_string += spec_name;
              spec_string += ')';
	      parent->addDerivePlotVar(spec_string);
              delete [] spec_name;
	  }
      }
  }
}

void
Castro::writePlotFile (const std::string& dir,
                       ostream&       os,
                       VisMF::How     how)
{
    int i, n;
    //
    // The list of indices of State to write to plotfile.
    // first component of pair is state_type,
    // second component of pair is component # within the state_type
    //
    std::vector<std::pair<int,int> > plot_var_map;
    for (int typ = 0; typ < desc_lst.size(); typ++)
        for (int comp = 0; comp < desc_lst[typ].nComp();comp++)
            if (parent->isStatePlotVar(desc_lst[typ].name(comp)) &&
                desc_lst[typ].getType() == IndexType::TheCellType())
                plot_var_map.push_back(std::pair<int,int>(typ,comp));

    int num_derive = 0;
    std::list<std::string> derive_names;
    const std::list<DeriveRec>& dlist = derive_lst.dlist();

    for (std::list<DeriveRec>::const_iterator it = dlist.begin();
	 it != dlist.end();
	 ++it)
    {
        if (parent->isDerivePlotVar(it->name()))
	{
            derive_names.push_back(it->name());
            num_derive++;
	}
    }

    int n_data_items = plot_var_map.size() + num_derive;
#ifdef GRAVITY
    if (do_grav) 
      if (gravity->get_gravity_type() == "PoissonGrav" && plot_phiGrav) n_data_items++;
#endif

    Real cur_time = state[State_Type].curTime();

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        //
        // The first thing we write out is the plotfile type.
        //
        os << thePlotFileType() << '\n';

        if (n_data_items == 0)
            BoxLib::Error("Must specify at least one valid data item to plot");

        os << n_data_items << '\n';

	//
	// Names of variables -- first state, then derived
	//
	for (i =0; i < plot_var_map.size(); i++)
        {
	    int typ = plot_var_map[i].first;
	    int comp = plot_var_map[i].second;
	    os << desc_lst[typ].name(comp) << '\n';
        }

	for ( std::list<std::string>::iterator it = derive_names.begin();
	      it != derive_names.end(); ++it)
        {
	    const DeriveRec* rec = derive_lst.get(*it);
            os << rec->variableName(0) << '\n';
        }
#ifdef GRAVITY
        if (do_grav && plot_phiGrav && gravity->get_gravity_type() == "PoissonGrav")
            os << "phiGrav" << '\n';
#endif

        os << BL_SPACEDIM << '\n';
        os << parent->cumTime() << '\n';
        int f_lev = parent->finestLevel();
        os << f_lev << '\n';
        for (i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbLo(i) << ' ';
        os << '\n';
        for (i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbHi(i) << ' ';
        os << '\n';
        for (i = 0; i < f_lev; i++)
            os << parent->refRatio(i)[0] << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->Geom(i).Domain() << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->levelSteps(i) << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
        {
            for (int k = 0; k < BL_SPACEDIM; k++)
                os << parent->Geom(i).CellSize()[k] << ' ';
            os << '\n';
        }
        os << (int) Geometry::Coord() << '\n';
        os << "0\n"; // Write bndry data.

#ifdef RADIATION
	if (do_radiation && Radiation::do_multigroup) {
	  std::ofstream groupfile;
	  std::string FullPathGroupFile = dir;
	  FullPathGroupFile += "/RadiationGroups";
	  groupfile.open(FullPathGroupFile.c_str(), std::ios::out);

	  radiation->write_groups(groupfile);

	  groupfile.close();
	}
#endif

        // job_info file with details about the run
	std::ofstream jobInfoFile;
	std::string FullPathJobInfoFile = dir;
	std::string PrettyLine = "===============================================================================\n";

	FullPathJobInfoFile += "/job_info";
	jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);	

	// job information
	jobInfoFile << PrettyLine;
	jobInfoFile << " Job Information\n";
	jobInfoFile << PrettyLine;
	
	jobInfoFile << "job name: " << job_name << "\n\n";

	jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
	jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif
	jobInfoFile << "\n\n";

        // plotfile information
	jobInfoFile << PrettyLine;
	jobInfoFile << " Plotfile Information\n";
	jobInfoFile << PrettyLine;

	time_t now = time(0);

	// Convert now to tm struct for local timezone
	tm* localtm = localtime(&now);
	jobInfoFile   << "output data / time: " << asctime(localtm);

	char currentDir[FILENAME_MAX];
	if (getcwd(currentDir, FILENAME_MAX)) {
	  jobInfoFile << "output dir:         " << currentDir << "\n";
	}

	jobInfoFile << "\n\n";


        // build information
	jobInfoFile << PrettyLine;
	jobInfoFile << " Build Information\n";
	jobInfoFile << PrettyLine;

	jobInfoFile << "build date:    " << buildInfoGetBuildDate() << "\n";
	jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << "\n";
	jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << "\n";
	jobInfoFile << "BoxLib dir:    " << buildInfoGetBoxlibDir() << "\n";

	jobInfoFile << "\n";
	
	jobInfoFile << "COMP:  " << buildInfoGetComp() << "\n";
	jobInfoFile << "FCOMP: " << buildInfoGetFcomp() << "\n";

	jobInfoFile << "\n";

	jobInfoFile << "EOS:     " << buildInfoGetAux(1) << "\n";
	jobInfoFile << "network: " << buildInfoGetAux(2) << "\n";

	jobInfoFile << "\n";

	const char* githash1 = buildInfoGetGitHash(1);
	const char* githash2 = buildInfoGetGitHash(2);
	const char* githash3 = buildInfoGetGitHash(3);
	if (strlen(githash1) > 0) {
	  jobInfoFile << "Castro   git hash: " << githash1 << "\n";
	}
	if (strlen(githash2) > 0) {
	  jobInfoFile << "BoxLib   git hash: " << githash2 << "\n";
	}
	if (strlen(githash3) > 0) {	
	  jobInfoFile << "AstroDev git hash: " << githash3 << "\n";
	}

	jobInfoFile << "\n\n";


	// runtime parameters
	jobInfoFile << PrettyLine;
	jobInfoFile << " Inputs File Parameters\n";
	jobInfoFile << PrettyLine;
	
	ParmParse::dumpTable(jobInfoFile, true);

	jobInfoFile.close();
	

    }
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";
    char buf[64];
    sprintf(buf, "Level_%d", level);
    std::string Level = buf;
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if (!FullPath.empty() && FullPath[FullPath.size()-1] != '/')
        FullPath += '/';
    FullPath += Level;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(FullPath, 0755))
            BoxLib::CreateDirectoryFailed(FullPath);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
        os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
        os << parent->levelSteps(level) << '\n';

        for (i = 0; i < grids.size(); ++i)
        {
            RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
            for (n = 0; n < BL_SPACEDIM; n++)
                os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
        }
        //
        // The full relative pathname of the MultiFabs at this level.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        if (n_data_items > 0)
        {
            std::string PathNameInHeader = Level;
            PathNameInHeader += BaseName;
            os << PathNameInHeader << '\n';
        }
    }
    //
    // We combine all of the multifabs -- state, derived, etc -- into one
    // multifab -- plotMF.
    // NOTE: we are assuming that each state variable has one component,
    // but a derived variable is allowed to have multiple components.
    int       cnt   = 0;
    const int nGrow = 0;
    MultiFab  plotMF(grids,n_data_items,nGrow);
    MultiFab* this_dat = 0;
    //
    // Cull data from state variables -- use no ghost cells.
    //
    for (i = 0; i < plot_var_map.size(); i++)
    {
	int typ  = plot_var_map[i].first;
	int comp = plot_var_map[i].second;
	this_dat = &state[typ].newData();
	MultiFab::Copy(plotMF,*this_dat,comp,cnt,1,nGrow);
	cnt++;
    }
    //
    // Cull data from derived variables.
    // 
    if (derive_names.size() > 0)
    {
	for (std::list<std::string>::iterator it = derive_names.begin();
	     it != derive_names.end(); ++it) 
	{
	    MultiFab* derive_dat = derive(*it,cur_time,nGrow);
	    MultiFab::Copy(plotMF,*derive_dat,0,cnt,1,nGrow);
	    delete derive_dat;
	    cnt++;
	}
    }

#ifdef GRAVITY
    //
    // Copy phi into plotMF.
    //
    if (do_grav && plot_phiGrav && gravity->get_gravity_type() == "PoissonGrav")
        MultiFab::Copy(plotMF,*gravity->get_phi_curr(level),0,cnt++,1,nGrow);
#endif

    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;
    VisMF::Write(plotMF,TheFullPath,how,true);

#ifdef PARTICLES
    //
    // Write the particles in a plotfile directory 
    // Particles are only written if particles.write_in_plotfile = 1 in inputs file.
    //
    ParticlePlotFile(dir);
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
    //
    // Loop over grids, call FORTRAN function to init with data.
    //
    int ns          = NUM_STATE;
    const Real* dx  = geom.CellSize();
    MultiFab& S_new = get_new_data(State_Type);
    Real cur_time   = state[State_Type].curTime();

    // make sure dx = dy = dz -- that's all we guarantee to support
    const Real SMALL = 1.e-13;
#if (BL_SPACEDIM == 2)
    if (fabs(dx[0] - dx[1]) > SMALL)
      {
	BoxLib::Abort("We don't support dx != dy");
      }
#elif (BL_SPACEDIM == 3)
    if ( (fabs(dx[0] - dx[1]) > SMALL) || (fabs(dx[0] - dx[2]) > SMALL) )
      {
	BoxLib::Abort("We don't support dx != dy != dz");
      }
#endif


    if (verbose && ParallelDescriptor::IOProcessor())
       std::cout << "Initializing the data at level " << level << std::endl;

#ifdef RADIATION
    // rad quantities are in the state even if (do_radiation == 0)
    MultiFab &Rad_new = get_new_data(Rad_Type);
    // extra state quantity for graphic debugging:
    MultiFab &Test_new = get_new_data(Test_Type);
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
          RealBox    gridloc = RealBox(grids[mfi.index()],geom.CellSize(),geom.ProbLo());
          const Box& box     = mfi.validbox();
          const int* lo      = box.loVect();
          const int* hi      = box.hiVect();
  
          // Temp unused for GammaLaw,
          //   set it here so that pltfiles have defined numbers
          S_new[mfi].setVal(0.,Temp);
  
          BL_FORT_PROC_CALL(CA_INITDATA,ca_initdata)
  	  (level, cur_time, lo, hi, ns,
  	   BL_TO_FORTRAN(S_new[mfi]), dx,
  	   gridloc.lo(), gridloc.hi());

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

	  BL_FORT_PROC_CALL(CA_INITRAD,ca_initrad)
	      (level, cur_time, lo, hi, Radiation::nGroups,
	       BL_TO_FORTRAN(Rad_new[i]),dx,
	       gridloc.lo(),gridloc.hi());

	  if (Radiation::nNeutrinoSpecies > 0 && Radiation::nNeutrinoGroups[0] == 0) {
	      // Hack: running photon radiation through neutrino solver
            Rad_new[i].mult(Radiation::Etorad, 
                            0, Radiation::nGroups);
	  }

          if (Rad_new.nComp() > Radiation::nGroups) {
            // Initialize flux components to 0
            Rad_new[i].setVal(0.0, box, Radiation::nGroups,
                              Rad_new.nComp() - Radiation::nGroups);
          }

          if (Test_new.nComp() > 0) {
            Test_new[i].setVal(0.0);
          }
      }
    }
#endif // RADIATION

#endif // MAESTRO_INIT

    set_special_tagging_flag(cur_time);

#if (BL_SPACEDIM > 1)
    if ( (level == 0) && (spherical_star == 1) ) {
       int nc = S_new.nComp();
       int n1d = get_numpts(level,geom.Domain());
       BL_FORT_PROC_CALL(ALLOCATE_OUTFLOW_DATA,allocate_outflow_data)(&n1d,&nc);
       int is_new = 1;
       make_radial_data(is_new);
    }
#endif

#ifdef GRAVITY
    // Set these to zero so they're defined for the plotfile.
    if (!do_grav) {
       MultiFab& G_new = get_new_data(Gravity_Type);
       G_new.setVal(0.);
    }
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
	   BL_TO_FORTRAN(LS_new[mfi.index()]),
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
    
    for (FillPatchIterator fpi(old,S_new,0,cur_time,State_Type,0,NUM_STATE);
          fpi.isValid();
          ++fpi)
    {
        S_new[fpi].copy(fpi());
    }

    // Set E in terms of e + kinetic energy
    // enforce_consistent_e(S_new);

#ifdef RADIATION
    if (do_radiation) {
      MultiFab& Er_new = get_new_data(Rad_Type);
      int ncomp = Er_new.nComp();
      
      for (FillPatchIterator fpi(old,Er_new,0,cur_time,Rad_Type,0,ncomp);
          fpi.isValid();
          ++fpi)
      {
        Er_new[fpi].copy(fpi());
      }
    }
#endif

#ifdef REACTIONS
    MultiFab& React_new = get_new_data(Reactions_Type);
    int ncomp = React_new.nComp();
      
      for (FillPatchIterator fpi(old,React_new,0,cur_time,Reactions_Type,0,ncomp);
          fpi.isValid();
          ++fpi)
      {
        React_new[fpi].copy(fpi());
      }

#endif

#ifdef LEVELSET
    MultiFab& LS_new = get_new_data(LS_State_Type);
    int nGrowRegrid = 0;
    
    for (FillPatchIterator fpi(old,LS_new,nGrowRegrid,cur_time,LS_State_Type,0,1);
	 fpi.isValid(); ++fpi)
      {
	LS_new[fpi.index()].copy(fpi());
      }
    
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
    if (fixed_dt > 0.0)
        return fixed_dt;

    // This is just a dummy value to start with 
    Real estdt  = 1.0e+200;

    const MultiFab& stateMF = get_new_data(State_Type);

    if (do_hydro) 
    {
      const Real* dx    = geom.CellSize();

#ifdef RADIATION
      if (Radiation::rad_hydro_combined) {
	const MultiFab& radMF = get_new_data(Rad_Type);
	for (MFIter mfi(stateMF); mfi.isValid(); ++mfi) {
	  int i = mfi.index();

	  const Box& box = mfi.validbox();
	  Real dt = estdt;
	  
	  FArrayBox gPr(box);
	  radiation->estimate_gamrPr(stateMF[i], radMF[i], gPr, dx, box);

	  BL_FORT_PROC_CALL(CA_ESTDT_RAD, ca_estdt_rad)
	    (BL_TO_FORTRAN(stateMF[i]),
	     BL_TO_FORTRAN(gPr),
	     box.loVect(),box.hiVect(),dx,&dt);
	  
	  estdt = std::min(estdt,dt);
	}
      }
      else 
      {
#endif
       for (MFIter mfi(stateMF); mfi.isValid(); ++mfi)
       {
           const Box& box = mfi.validbox();
           Real dt = estdt;

   	      BL_FORT_PROC_CALL(CA_ESTDT,ca_estdt)
                  (BL_TO_FORTRAN(stateMF[mfi]),
                   box.loVect(),box.hiVect(),
                   dx,&dt);
   
           estdt = std::min(estdt,dt);
       }

#ifdef RADIATION
      }
#endif

       ParallelDescriptor::ReduceRealMin(estdt);
       estdt *= cfl;
       if (verbose && ParallelDescriptor::IOProcessor()) 
           std::cout << "...estdt from     hydro at level " << level << ": " << estdt << std::endl;
    }

#ifdef REACTIONS
#if 0
    const MultiFab& reactMF = get_new_data(Reactions_Type);

    // If these quantities aren't defined yet then we haven't taken a step at this level yet
    int initial_step = 0;

    Real cur_time  = state[State_Type].curTime();
    if (cur_time == 0.0) initial_step = 1;

    Real estdt_burning = 1.e200;
    Real dt            = 1.e200;
    for (MFIter mfi(stateMF); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
        BL_FORT_PROC_CALL(CA_ESTDT_BURNING,ca_estdt_burning)
            (BL_TO_FORTRAN(stateMF[mfi]),
             BL_TO_FORTRAN(reactMF[mfi]),
             box.loVect(),box.hiVect(),&dt_old,&dt,
             &initial_step);

        estdt_burning = std::min(estdt_burning,dt);
    }

    ParallelDescriptor::ReduceRealMin(estdt_burning);

    if (verbose && ParallelDescriptor::IOProcessor())
        cout << "hydro timestep at level " << level << ":  estdt = " << estdt << '\n';

    if (verbose && ParallelDescriptor::IOProcessor())
        cout << "burning timestep at level " << level << ":  estdt_burning = " << estdt_burning << '\n';

    estdt = std::min(estdt, estdt_burning);

#endif
#endif

#ifdef RADIATION
    radiation->EstTimeStep(estdt, level);
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
    //
    // Integration cycle on fine level grids is complete
    // do post_timestep stuff here.
    //
    int finest_level = parent->finestLevel();

#ifdef RADIATION
    if (do_radiation && (level < finest_level)) {
      if (Radiation::do_deferred_sync == 1) {
        // computeTemp is not needed before or after this call because
        // setup for deferred sync does not touch state, only flux registers.
        radiation->deferred_sync_setup(level);

	if (do_reflux) {
	  radiation->reflux(level);

          // Since radiation->reflux does not touch the fluid state,
          // we do need to recompute Temp here.
	}
      }
      else {

        radiation->sync_solve(level);

	if (do_reflux) {
	  radiation->reflux(level);
	}

        // Rad sync touches these levels, so recompute temp on all of them.
        for (int lev = level; lev <= finest_level; lev++) {
          MultiFab& S_new = getLevel(lev).get_new_data(State_Type);
          computeTemp(S_new);
        }

      }
    }
#endif

#ifdef GRAVITY
    // Check the whole hierarchy before the syncs
    if (level == 0 && gravity->test_results_of_solves() == 1 &&
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
           drho_and_drhoU.define(grids,BL_SPACEDIM+1,0,Fab_allocate);
           MultiFab::Copy(drho_and_drhoU,S_new_crse,Density,0,BL_SPACEDIM+1,0);
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

           for (MFIter mfi(drho_and_drhoU); mfi.isValid(); ++mfi)
              drho_and_drhoU[mfi].plus(S_new_crse[mfi],Density,0,BL_SPACEDIM+1);

            MultiFab dphi(grids,1,0,Fab_allocate);
            dphi.setVal(0.);

            gravity->reflux_phi(level,dphi);
                
            // Compute (cross-level) gravity sync based on drho, dphi
            PArray<MultiFab> grad_delta_phi_cc(finest_level-level+1,PArrayManage); 
            for (int lev = level; lev <= finest_level; lev++) {
               grad_delta_phi_cc.set(lev-level,
                                     new MultiFab(getLevel(lev).boxArray(),BL_SPACEDIM,0,Fab_allocate));
               grad_delta_phi_cc[lev-level].setVal(0.);
            }
            gravity->gravity_sync(level,finest_level,drho_and_drhoU,dphi,grad_delta_phi_cc);

            // Create sync src terms at coarser level (i.e. "level")
            FArrayBox grad_phi_cc, sync_src;
            FArrayBox dstate;

            for (int lev = level; lev <= finest_level; lev++)  
            {
              Real dt_lev = parent->dtLevel(lev);
              MultiFab&  S_new_lev = getLevel(lev).get_new_data(State_Type);
              Real cur_time = state[State_Type].curTime();

              const BoxArray& ba = getLevel(lev).boxArray();
              MultiFab grad_phi_cc(ba,BL_SPACEDIM,0,Fab_allocate);
              gravity->get_new_grav_vector(lev,grad_phi_cc,cur_time);

              for (MFIter mfi(S_new_lev); mfi.isValid(); ++mfi)
              {
                const Box bx = mfi.validbox();
                dstate.resize(bx,BL_SPACEDIM+1);
                if (lev == level) {
                   dstate.copy(drho_and_drhoU[mfi]);
                } else {
                   dstate.setVal(0.); 
                }

                // Compute sync source
                sync_src.resize(bx,BL_SPACEDIM+1);
                int i = mfi.index();
                BL_FORT_PROC_CALL(CA_SYNCGSRC,ca_syncgsrc)
                    (bx.loVect(), bx.hiVect(),
                     BL_TO_FORTRAN(grad_phi_cc[i]),
                     BL_TO_FORTRAN(grad_delta_phi_cc[lev-level][i]),
                     BL_TO_FORTRAN(S_new_lev[i]),
                     BL_TO_FORTRAN(dstate),
                     BL_TO_FORTRAN(sync_src),
                     dt_lev);

                sync_src.mult(0.5*dt_lev);
                S_new_lev[mfi].plus(sync_src,0,Xmom,BL_SPACEDIM);
                S_new_lev[mfi].plus(sync_src,0,Eden,1);
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
            int numpts_1d = get_numpts (level, geom.Domain());
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
                   gravity->multilevel_solve_for_phi(0,parent->finestLevel());
                   if (gravity->test_results_of_solves() == 1)
                       gravity->test_composite_phi(level);
                }
#ifdef PARTICLES
                if (do_dm_particles)
                {
                    // Do solve if we haven't already done it above
                    if (gravity->NoComposite() == 1)
                       gravity->multilevel_solve_for_phi(0,parent->finestLevel());

                    for (int k = 0; k <= parent->finestLevel(); k++)
                    {
                        const BoxArray& ba = getLevel(k).boxArray();
                        MultiFab grav_vec_new(ba,BL_SPACEDIM,0,Fab_allocate);
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
#ifdef PARTICLES
    ParticleMoveRandom();
#endif
}

void
Castro::post_regrid (int lbase,
                     int new_finest)
{

#ifdef PARTICLES
    if (level == lbase) ParticleRedistribute();
#endif

#ifdef GRAVITY
    if (do_grav)
    {
       const Real cur_time = state[State_Type].curTime();
       if ( (level == lbase) && cur_time > 0.)  
       {
          if ( gravity->get_gravity_type() == "PoissonGrav" && (gravity->NoComposite() != 1) )  
              gravity->multilevel_solve_for_phi(level,new_finest);
       }
    }
#endif

}

void
Castro::post_init (Real stop_time)
{
    if (level > 0)
        return;
    //
    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    //
    int finest_level = parent->finestLevel();
    for (int k = finest_level-1; k>= 0; k--)
        getLevel(k).avgDown();

#ifdef GRAVITY
    Real cur_time = state[State_Type].curTime();
    if (do_grav) {
       if (gravity->get_gravity_type() == "PoissonGrav") {

          // Calculate offset before first multilevel solve.
          gravity->set_mass_offset(cur_time);

          if (gravity->NoComposite() != 1)  {
             gravity->multilevel_solve_for_phi(level,finest_level);
             if (gravity->test_results_of_solves() == 1)
                gravity->test_composite_phi(level);
          }
       }
 
       // Make this call just to fill the initial state data.
       for (int k = 0; k <= parent->finestLevel(); k++)
       {
          BoxArray ba = getLevel(k).boxArray();
          MultiFab grav_vec_new(ba,BL_SPACEDIM,0,Fab_allocate);
          gravity->get_new_grav_vector(k,grav_vec_new,cur_time);
       }
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
#ifdef GRAVITY
    Real cur_time = state[State_Type].curTime();
    if (do_grav) {
       if (gravity->get_gravity_type() == "PoissonGrav") {

          // Calculate offset before first multilevel solve.
          gravity->set_mass_offset(cur_time);

          if (gravity->NoComposite() != 1)  {
             gravity->multilevel_solve_for_phi(level,finest_level);
             if (gravity->test_results_of_solves() == 1)
                gravity->test_composite_phi(level);
          }
       }
 
       // Make this call just to fill the initial state data.
       for (int k = 0; k <= parent->finestLevel(); k++)
       {
          BoxArray ba = getLevel(k).boxArray();
          MultiFab grav_vec_new(ba,BL_SPACEDIM,0,Fab_allocate);
          gravity->get_new_grav_vector(k,grav_vec_new,cur_time);
       }
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

    for (MFIter mfi(S_old); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();
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
        
        LStype.FillBoundary();
        BoxLib::FillPeriodicBoundary(geom,LStype,0,1);
        
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
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "...... reinitialzing levelset\n";

    MultiFab& LS_new = get_new_data(LS_State_Type);
    const Real* dx   = geom.CellSize();
    int nGrowLS = 2;
    int nCompLS = 1;
    
    LStype.FillBoundary();
    BoxLib::FillPeriodicBoundary(geom,LStype,0,1);
    
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
        LS_new.FillBoundary();
        geom.FillPeriodicBoundary(LS_new);
        LStype.FillBoundary();
        BoxLib::FillPeriodicBoundary(geom,LStype,0,1);
        
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
    // Subtract off half of the old source term, and add half of the new.

    ext_src_old.mult(-0.5*dt);
    ext_src_new.mult( 0.5*dt);
    MultiFab::Add(S_new,ext_src_old,0,0,S_new.nComp(),0);
    MultiFab::Add(S_new,ext_src_new,0,0,S_new.nComp(),0);
}

#ifdef SGS
void
Castro::getOldSource (Real old_time, Real dt, MultiFab&  ext_src, MultiFab* sgs_fluxes)
{
   const Real* dx = geom.CellSize();

   MultiFab& S_old = get_old_data(State_Type);
   const int ncomp = S_old.nComp();

   ext_src.setVal(0.0,ext_src.nGrow());

   MultiFab& sgs_mf = get_old_data(SGS_Type);
   sgs_mf.setVal(0.0);

   // Set these to zero so we always add them up correctly
   for (int dir = 0; dir < BL_SPACEDIM; dir++) 
       sgs_fluxes[dir].setVal(0.0);

   // We need to use temporary FABs to hold the fluxes for the old source
   //   because sgs_fluxes has no ghost cells but the fluxes array in the ext_src routine
   //   needs to have a layer of ghost cells for fluxes
   FArrayBox fluxx, fluxy, fluxz;

   for (FillPatchIterator Old_fpi(*this,S_old,4,old_time,State_Type,Density,ncomp);
                          Old_fpi.isValid();++Old_fpi)
   {
        const Box& bx = grids[Old_fpi.index()];
 
        Box bxx = bx; bxx.surroundingNodes(0); bxx.grow(1);
        fluxx.resize(bxx,NUM_STATE);

        Box bxy = bx; bxy.surroundingNodes(1); bxy.grow(1);
        fluxy.resize(bxy,NUM_STATE);

        Box bxz = bx; bxz.surroundingNodes(2); bxz.grow(1);
        fluxz.resize(bxz,NUM_STATE);

        BL_FORT_PROC_CALL(CA_EXT_SRC,ca_ext_src)
            (bx.loVect(), bx.hiVect(),
             BL_TO_FORTRAN(  Old_fpi()),
             BL_TO_FORTRAN(fluxx),
             BL_TO_FORTRAN(fluxy),
             BL_TO_FORTRAN(fluxz),
             BL_TO_FORTRAN(ext_src[Old_fpi]),
             BL_TO_FORTRAN_N(sgs_mf[Old_fpi],0),
             BL_TO_FORTRAN_N(sgs_mf[Old_fpi],1),
             BL_TO_FORTRAN_N(sgs_mf[Old_fpi],2),
             dx,&old_time,&dt);

        sgs_fluxes[0][Old_fpi.index()].copy(fluxx,0,0,NUM_STATE);
        sgs_fluxes[1][Old_fpi.index()].copy(fluxy,0,0,NUM_STATE);
        sgs_fluxes[2][Old_fpi.index()].copy(fluxz,0,0,NUM_STATE);
   }
   geom.FillPeriodicBoundary(ext_src,0,NUM_STATE);
}

void
Castro::getNewSource (Real new_time, Real dt, MultiFab& ext_src, MultiFab* sgs_fluxes)
{
   const Real* dx = geom.CellSize();

   MultiFab& S_new = get_new_data(State_Type);
   const int ncomp = S_new.nComp();

   ext_src.setVal(0.0,ext_src.nGrow());

   MultiFab& sgs_mf = get_new_data(SGS_Type);
   sgs_mf.setVal(0.0);

   // Set these to zero so we always add them up correctly
   for (int dir = 0; dir < BL_SPACEDIM; dir++) 
       sgs_fluxes[dir].setVal(0.0);

   for (FillPatchIterator New_fpi(*this,S_new,4,new_time,State_Type,Density,ncomp);
                          New_fpi.isValid();++New_fpi)
   {
        const Box& bx = grids[New_fpi.index()];

        BL_FORT_PROC_CALL(CA_EXT_SRC,ca_ext_src)
            (bx.loVect(), bx.hiVect(),
             BL_TO_FORTRAN(  New_fpi()),
             BL_TO_FORTRAN(sgs_fluxes[0][New_fpi]),
             BL_TO_FORTRAN(sgs_fluxes[1][New_fpi]),
#if (BL_SPACEDIM == 3)
             BL_TO_FORTRAN(sgs_fluxes[2][New_fpi]),
#endif
             BL_TO_FORTRAN(ext_src[New_fpi]),
             BL_TO_FORTRAN_N(sgs_mf[New_fpi],0),
             BL_TO_FORTRAN_N(sgs_mf[New_fpi],1),
             BL_TO_FORTRAN_N(sgs_mf[New_fpi],2),
             dx,&new_time,&dt);
   }
   geom.FillPeriodicBoundary(ext_src,0,NUM_STATE);

}

#else

void
Castro::getOldSource (Real old_time, Real dt, MultiFab&  ext_src)
{
   const Real* dx = geom.CellSize();
   const Real* prob_lo   = geom.ProbLo();

   MultiFab& S_old = get_old_data(State_Type);
   const int ncomp = S_old.nComp();

   ext_src.setVal(0.0,ext_src.nGrow());

   MultiFab levelArea[BL_SPACEDIM];
   for (int i = 0; i < BL_SPACEDIM ; i++)
       geom.GetFaceArea(levelArea[i],grids,i,NUM_GROW);

   for (FillPatchIterator Old_fpi(*this,S_old,4,old_time,State_Type,Density,ncomp);
                          Old_fpi.isValid();++Old_fpi)
   {
        const Box& bx = grids[Old_fpi.index()];
        BL_FORT_PROC_CALL(CA_EXT_SRC,ca_ext_src)
            (bx.loVect(), bx.hiVect(),
             BL_TO_FORTRAN(  Old_fpi()),
             BL_TO_FORTRAN(  Old_fpi()),
             BL_TO_FORTRAN(ext_src[Old_fpi]),
             prob_lo,dx,&old_time,&dt);
   }
   geom.FillPeriodicBoundary(ext_src,0,NUM_STATE);
}

void
Castro::getNewSource (Real old_time, Real new_time, Real dt, MultiFab& ext_src)
{
   const Real* dx = geom.CellSize();
   const Real* prob_lo   = geom.ProbLo();

   MultiFab& S_old = get_old_data(State_Type);
   const int ncomp = S_old.nComp();

   ext_src.setVal(0.0,ext_src.nGrow());

   for (FillPatchIterator Old_fpi(*this,S_old,4,old_time,State_Type,Density,ncomp),
                          New_fpi(*this,S_old,4,new_time,State_Type,Density,ncomp);
                          Old_fpi.isValid() && New_fpi.isValid();++Old_fpi,++New_fpi)
   {
        const Box& bx = grids[Old_fpi.index()];
        BL_FORT_PROC_CALL(CA_EXT_SRC,ca_ext_src)
            (bx.loVect(), bx.hiVect(),
             BL_TO_FORTRAN(  Old_fpi()),
             BL_TO_FORTRAN(  New_fpi()),
             BL_TO_FORTRAN(ext_src[Old_fpi]),
             prob_lo,dx,&new_time,&dt);
   }
   geom.FillPeriodicBoundary(ext_src,0,NUM_STATE);
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
                  BL_TO_FORTRAN(tau_diff[i]),
                  BL_TO_FORTRAN(fpi()),
                  BL_TO_FORTRAN(grav_vector[i]),
                  dx_fine);
   }
}
#endif

void
Castro::getTempDiffusionTerm (Real time, MultiFab& TempDiffTerm, MultiFab* tau)
{
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

   // Fill temperature at this level.
   MultiFab Temperature(grids,1,1,Fab_allocate);
   for (FillPatchIterator fpi(*this,S_old,1,time,State_Type,Temp,1);
                          fpi.isValid();++fpi)
   {
       Temperature[fpi].copy(fpi(),fpi().box());
   }

   // Fill coefficients at this level.
   PArray<MultiFab> coeffs(BL_SPACEDIM,PArrayManage);
   for (int dir = 0; dir < BL_SPACEDIM ; dir++) {
      coeffs.set(dir,new MultiFab);
      BoxArray edge_boxes(grids);
      edge_boxes.surroundingNodes(dir);
      coeffs[dir].define(edge_boxes,1,0,Fab_allocate);
   }

   const Geometry& fine_geom = parent->Geom(parent->finestLevel());
   const Real*       dx_fine = fine_geom.CellSize();

   for (FillPatchIterator fpi(*this,S_old,1,time,State_Type,0,NUM_STATE);
                          fpi.isValid();++fpi)
   {
       Box bx(fpi.validbox());
       int i = fpi.index();
       BL_FORT_PROC_CALL(CA_FILL_TEMP_COND,ca_fill_temp_cond)
                (bx.loVect(), bx.hiVect(),
                 BL_TO_FORTRAN(fpi()),
#ifdef TAU
                 BL_TO_FORTRAN((*tau)[i]),
#endif
                 D_DECL(BL_TO_FORTRAN(coeffs[0][i]),
                        BL_TO_FORTRAN(coeffs[1][i]),
                        BL_TO_FORTRAN(coeffs[2][i])),
                 dx_fine);
   }

   if (Geometry::isAnyPeriodic())
     for (int d = 0; d < BL_SPACEDIM; d++)
       geom.FillPeriodicBoundary(coeffs[d]);

   if (level == 0) {
      diffusion->applyop(level,Temperature,TempDiffTerm,coeffs);
   } else if (level > 0) {
      // Fill temperature at next coarser level, if it exists.
      const BoxArray& crse_grids = getLevel(level-1).boxArray();
      MultiFab CrseTemp(crse_grids,1,1,Fab_allocate);
      for (FillPatchIterator fpi(getLevel(level-1),CrseTemp,1,time,State_Type,Temp,1);
          fpi.isValid();
          ++fpi)
      {
        CrseTemp[fpi].copy(fpi());
      }
      diffusion->applyop(level,Temperature,CrseTemp,TempDiffTerm,coeffs);
   }
}
#endif

void
Castro::reflux ()
{
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

        ParallelDescriptor::ReduceRealMax(end,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Castro::reflux() at level " << level << " : time = " << end << std::endl;
    }
}

void
Castro::avgDown ()
{
  if (level == parent->finestLevel()) return;

  avgDown(State_Type);

#ifdef GRAVITY
  avgDown(Gravity_Type);
#endif

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
    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
       const Box bx = mfi.validbox();
       BL_FORT_PROC_CALL(CA_ENFORCE_NONNEGATIVE_SPECIES,ca_enforce_nonnegative_species)
           (BL_TO_FORTRAN(S_new[mfi]),bx.loVect(),bx.hiVect());
    }
}

void
Castro::enforce_consistent_e (MultiFab& S)
{
    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        const Box& box     = mfi.validbox();
        const int* lo      = box.loVect();
        const int* hi      = box.hiVect();
        BL_FORT_PROC_CALL(CA_ENFORCE_CONSISTENT_E,ca_enforce_consistent_e)
          (lo, hi, BL_TO_FORTRAN(S[mfi]));
    }
}

void
Castro::avgDown (int state_indx)
{
    if (level == parent->finestLevel()) return;

    Castro& fine_lev = getLevel(level+1);
    MultiFab&  S_crse   = get_new_data(state_indx);
    MultiFab&  S_fine   = fine_lev.get_new_data(state_indx);
    MultiFab&  fvolume  = fine_lev.volume;
    const int  ncomp    = S_fine.nComp();

    BL_ASSERT(S_crse.boxArray() == volume.boxArray());
    BL_ASSERT(fvolume.boxArray() == S_fine.boxArray());
    //
    // Coarsen() the fine stuff on processors owning the fine data.
    //
    BoxArray crse_S_fine_BA(S_fine.boxArray().size());

    for (int i = 0; i < S_fine.boxArray().size(); ++i)
    {
        crse_S_fine_BA.set(i,BoxLib::coarsen(S_fine.boxArray()[i],fine_ratio));
    }

    MultiFab crse_S_fine(crse_S_fine_BA,ncomp,0);
    MultiFab crse_fvolume(crse_S_fine_BA,1,0);

    crse_fvolume.copy(volume);

    for (MFIter mfi(S_fine); mfi.isValid(); ++mfi)
    {
        const int        i        = mfi.index();
        const Box&       ovlp     = crse_S_fine_BA[i];
        FArrayBox&       crse_fab = crse_S_fine[i];
        const FArrayBox& crse_vol = crse_fvolume[i];
        const FArrayBox& fine_fab = S_fine[i];
        const FArrayBox& fine_vol = fvolume[i];

	BL_FORT_PROC_CALL(CA_AVGDOWN,ca_avgdown)
            (BL_TO_FORTRAN(crse_fab), ncomp,
             BL_TO_FORTRAN(crse_vol),
             BL_TO_FORTRAN(fine_fab),
             BL_TO_FORTRAN(fine_vol),
             ovlp.loVect(),ovlp.hiVect(),fine_ratio.getVect());
    }

    S_crse.copy(crse_S_fine);
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
    const int*  domain_lo = geom.Domain().loVect();
    const int*  domain_hi = geom.Domain().hiVect();
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();
    Array<int>  itags;

    for (int j = 0; j < err_list.size(); j++)
    {
        MultiFab* mf = derive(err_list[j].name(), time, err_list[j].nGrow());

        BL_ASSERT(!(mf == 0));

        for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
        {
            int         idx     = mfi.index();
            RealBox     gridloc = RealBox(grids[idx],geom.CellSize(),geom.ProbLo());
            itags               = tags[idx].tags();
            int*        tptr    = itags.dataPtr();
            const int*  tlo     = tags[idx].box().loVect();
            const int*  thi     = tags[idx].box().hiVect();
            const int*  lo      = mfi.validbox().loVect();
            const int*  hi      = mfi.validbox().hiVect();
            const Real* xlo     = gridloc.lo();
            Real*       dat     = (*mf)[mfi].dataPtr();
            const int*  dlo     = (*mf)[mfi].box().loVect();
            const int*  dhi     = (*mf)[mfi].box().hiVect();
            const int   ncomp   = (*mf)[mfi].nComp();

            err_list[j].errFunc()(tptr, ARLIM(tlo), ARLIM(thi), &tagval,
				  &clearval, dat, ARLIM(dlo), ARLIM(dhi),
				  lo,hi, &ncomp, domain_lo, domain_hi,
				  dx, xlo, prob_lo, &time, &level);
            //
            // Don't forget to set the tags in the TagBox.
            //
            if (allow_untagging == 1) 
            {
               tags[idx].tags_and_untags(itags);
            } else {
               tags[idx].tags(itags);
            }
        }

        delete mf;
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

Real
Castro::sumDerive (const std::string& name,
                   Real           time)
{
    Real sum     = 0.0;
    MultiFab* mf = derive(name, time, 0);

    BL_ASSERT(!(mf == 0));

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        if (level < parent->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                fab.setVal(0,isects[ii].second,0);
            }
        }

        sum += fab.sum(0);
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

//
// Helper function for Castro::SyncInterp().
//

static
void
set_bc_new (int*            bc_new,
            int             n,
            int             src_comp,
            const int*      clo,
            const int*      chi,
            const int*      cdomlo,
            const int*      cdomhi,
            const BoxArray& cgrids,
            int**           bc_orig_qty)
            
{
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        int bc_index = (n+src_comp)*(2*BL_SPACEDIM) + dir;
        bc_new[bc_index]             = INT_DIR;
        bc_new[bc_index+BL_SPACEDIM] = INT_DIR;
 
        if (clo[dir] < cdomlo[dir] || chi[dir] > cdomhi[dir])
        {
            for (int crse = 0; crse < cgrids.size(); crse++)
            {
                const int* c_lo = cgrids[crse].loVect();
                const int* c_hi = cgrids[crse].hiVect();

                if (clo[dir] < cdomlo[dir] && c_lo[dir] == cdomlo[dir])
                    bc_new[bc_index] = bc_orig_qty[crse][bc_index];
                if (chi[dir] > cdomhi[dir] && c_hi[dir] == cdomhi[dir])
                    bc_new[bc_index+BL_SPACEDIM] = bc_orig_qty[crse][bc_index+BL_SPACEDIM]; 
            }
        }
    }
}

//
// Interpolate a cell-centered Sync correction from a
// coarse level (c_lev) to a fine level (f_lev).
//
// This routine interpolates the num_comp components of CrseSync
// (starting at src_comp) and either increments or puts the result into
// the num_comp components of FineSync (starting at dest_comp)
// The components of bc_orig_qty corespond to the quantities of CrseSync.
//

void
Castro::SyncInterp (MultiFab&      CrseSync,
                    int            c_lev,
                    MultiFab&      FineSync,
                    int            f_lev,
                    IntVect&       ratio,
                    int            src_comp,
                    int            dest_comp,
                    int            num_comp,
                    int            increment,
                    Real           dt_clev, 
                    int**          bc_orig_qty,
                    SyncInterpType which_interp,
                    int            state_comp)
{
    BL_ASSERT(which_interp >= 0 && which_interp <= 5);

    Interpolater* interpolater = 0;

    switch (which_interp)
    {
    case PC_T:           interpolater = &pc_interp;           break;
    case CellCons_T:     interpolater = &cell_cons_interp;    break;
    case CellConsLin_T:  interpolater = &lincc_interp;        break;
    case CellConsProt_T: interpolater = &protected_interp;    break;
    default:
        BoxLib::Abort("Castro::SyncInterp(): how did this happen");
    }

    Castro&   fine_level = getLevel(f_lev);
    const BoxArray& fgrids     = fine_level.boxArray();
    const Geometry& fgeom      = parent->Geom(f_lev);
    const BoxArray& cgrids     = getLevel(c_lev).boxArray();
    const Geometry& cgeom      = parent->Geom(c_lev);
    const Real*     dx_crse    = cgeom.CellSize();
    Box             cdomain    = BoxLib::coarsen(fgeom.Domain(),ratio);
    const int*      cdomlo     = cdomain.loVect();
    const int*      cdomhi     = cdomain.hiVect();
    int*            bc_new     = new int[2*BL_SPACEDIM*(src_comp+num_comp)];

    BoxArray cdataBA(fgrids.size());

    for (int i = 0; i < fgrids.size(); i++)
        cdataBA.set(i,interpolater->CoarseBox(fgrids[i],ratio));
    //
    // Note: The boxes in cdataBA may NOT be disjoint !!!
    //
    MultiFab cdataMF(cdataBA,num_comp,0);

    cdataMF.setVal(0);

    cdataMF.copy(CrseSync, src_comp, 0, num_comp);
    //
    // Set physical boundary conditions in cdataMF.
    //
    // HACK HACK HACK -- for now to get it to compile
#if 1
    for (MFIter mfi(cdataMF); mfi.isValid(); ++mfi)
    {
        int         i       = mfi.index();
        RealBox     gridloc = RealBox(fine_level.boxArray()[i],
                                      fine_level.Geom().CellSize(),
                                      fine_level.Geom().ProbLo());
        FArrayBox&  cdata   = cdataMF[mfi];
        const int*  clo     = cdata.loVect();
        const int*  chi     = cdata.hiVect();
        const Real* xlo     = gridloc.lo();

        for (int n = 0; n < num_comp; n++)
        {
            set_bc_new(bc_new,n,src_comp,clo,chi,cdomlo,cdomhi,cgrids,bc_orig_qty);

            BL_FORT_PROC_CALL(FILCC,filcc)
                (BL_TO_FORTRAN(cdata),
                 cdomlo, cdomhi, dx_crse, xlo,
                 &(bc_new[2*BL_SPACEDIM*(n+src_comp)]));
        }
    }
#endif
    cgeom.FillPeriodicBoundary(cdataMF, 0, num_comp);
    //
    // Interpolate from cdataMF to fdata and update FineSync.
    // Note that FineSync and cdataMF will have the same distribution
    // since the length of their BoxArrays are equal.
    //
    FArrayBox    fdata;
    Array<BCRec> bc_interp(num_comp);

    MultiFab* fine_stateMF = 0;
    if (interpolater == &protected_interp)
    {
        fine_stateMF = &(getLevel(f_lev).get_new_data(State_Type));
    }

    for (MFIter mfi(cdataMF); mfi.isValid(); ++mfi)
    {
        int        i     = mfi.index();
        FArrayBox& cdata = cdataMF[mfi];
        const int* clo   = cdata.loVect();
        const int* chi   = cdata.hiVect();

        fdata.resize(fgrids[i], num_comp);
        //
        // Set the boundary condition array for interpolation.
        //
        for (int n = 0; n < num_comp; n++)
        {
            set_bc_new(bc_new,n,src_comp,clo,chi,cdomlo,cdomhi,cgrids,bc_orig_qty);
        }

        for (int n = 0; n < num_comp; n++)
        {
            for (int dir = 0; dir < BL_SPACEDIM; dir++)
            {
                int bc_index = (n+src_comp)*(2*BL_SPACEDIM) + dir;
                bc_interp[n].setLo(dir,bc_new[bc_index]);
                bc_interp[n].setHi(dir,bc_new[bc_index+BL_SPACEDIM]);
            }
        }

        interpolater->interp(cdata,0,fdata,0,num_comp,fgrids[i],ratio,
                             cgeom,fgeom,bc_interp,src_comp,State_Type);

        if (increment)
        {
            fdata.mult(dt_clev);

            if (interpolater == &protected_interp)
            {
              cdata.mult(dt_clev);
              FArrayBox& fine_state = (*fine_stateMF)[i];
              interpolater->protect(cdata,0,fdata,0,fine_state,state_comp,
                                    num_comp,fgrids[i],ratio,
                                    cgeom,fgeom,bc_interp);
              Real dt_clev_inv = 1./dt_clev;
              cdata.mult(dt_clev_inv);
            }
            
            FineSync[i].plus(fdata,0,dest_comp,num_comp);
        }
        else
        {
            FineSync[i].copy(fdata,0,dest_comp,num_comp);
        }
    }

    delete [] bc_new;
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
  std::string probin_file = "probin";

  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "reading extern runtime parameters ..." << std::endl;
  }

  ParmParse pp("amr");
  if (pp.contains("probin_file"))
  {
    pp.get("probin_file",probin_file);
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
       const Box bx = mfi.validbox();

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
       const Box bx = mfi.validbox();

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
        sum0 = volWgtSumMF(&S_new,Eden);
    }

    // Synchronize (rho e) and (rho E) so they are consistent with each other
    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        const Box bx = mfi.validbox();

        BL_FORT_PROC_CALL(RESET_INTERNAL_E,reset_internal_e)
            (BL_TO_FORTRAN(S_new[mfi]),
             bx.loVect(), bx.hiVect(), print_fortran_warnings);
    }

    if (parent->finestLevel() == 0 && print_energy_diagnostics)
    {
        // Pass in the multifab and the component
        sum = volWgtSumMF(&S_new,Eden);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "(rho E) added from reset terms      : " << sum-sum0 << " out of " << sum0 << std::endl;
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

    for (MFIter mfi(State); mfi.isValid(); ++mfi)
    { const Box bx = mfi.validbox();
	BL_FORT_PROC_CALL(COMPUTE_TEMP,compute_temp)
	  (bx.loVect(),bx.hiVect(),BL_TO_FORTRAN(State[mfi]));
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
Castro::get_numpts (int level, Box bx)
{
     int numpts_1d;

     int nx = bx.size()[0];

#if (BL_SPACEDIM == 1)
     numpts_1d = nx;
#elif (BL_SPACEDIM == 2)
     int ny = bx.size()[1];

     if (level == 0)
     {
        Real ndiagsq = Real(nx*nx + ny*ny);
        numpts_1d = int(sqrt(ndiagsq))+6;
     }
     else
     {
        numpts_1d = std::min(nx,ny);
     }
#elif (BL_SPACEDIM == 3)
     int ny = bx.size()[1];
     int nz = bx.size()[2];
     if (level == 0)
     {
        Real ndiagsq = Real(nx*nx + ny*ny + nz*nz);
        numpts_1d = int(sqrt(ndiagsq))+6;
     }
     else
     {
        numpts_1d = std::min(nx,ny);
        numpts_1d = std::min(numpts_1d,nz);
     }
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
   
   int numpts_1d = get_numpts(level, geom.Domain());

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
    for (FillPatchIterator fpi(*this,mf,0,time,State_Type,Density,1); fpi.isValid(); ++fpi)
    {
        mf[fpi].copy(fpi());
    }

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
