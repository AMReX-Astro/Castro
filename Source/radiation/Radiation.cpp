#include <AMReX_LO_BCTYPES.H>
#include <AMReX_ParmParse.H>
#include <Radiation.H>
#include <RadSolve.H>
#include <rad_util.H>
#include <filt_prim.H>

#include <opacity.H>

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <sstream>

using namespace amrex;

Radiation::Solver_Type Radiation::SolverType = Radiation::InvalidSolver;

Real Radiation::radtoE = 0.;
//Real Radiation::radtoJ = 0.;
Real Radiation::Etorad = 0.;
Real Radiation::radfluxtoF = 0.;

int Radiation::do_multigroup = 0;
int Radiation::nGroups = NGROUPS;
int Radiation::accelerate = 1;
int Radiation::rad_hydro_combined = 0;
int Radiation::Er_Lorentz_term = 1;

int Radiation::icomp_lambda  = -1;
int Radiation::icomp_kp      = -1;
int Radiation::icomp_kr      = -1;
int Radiation::icomp_lab_Er  = -1;
int Radiation::icomp_lab_Fr  = -1;
int Radiation::icomp_com_Fr  = -1;
int Radiation::nplotvar      = 0;
Vector<std::string>  Radiation::plotvar_names;
int Radiation::filter_lambda_T = 0;
int Radiation::filter_lambda_S = 0;
int Radiation::filter_prim_int = 0;
int Radiation::filter_prim_T = 4;
int Radiation::filter_prim_S = 0;

// These physical constants get their values in the Radiation constructor:
Real Radiation::convert_MeV_erg = 0.0;
Real Radiation::clight          = 0.0;
Real Radiation::hPlanck         = 0.0;
Real Radiation::kBoltz          = 0.0;
Real Radiation::Avogadro        = 0.0;

Real Radiation::c               = 0.0;
Real Radiation::sigma           = 0.0;
Real Radiation::aRad            = 0.0;

int Radiation::current_group_number = -1;
std::string Radiation::current_group_name = "Radiation";

int Radiation::pure_hydro = 0;

// static initialization, must be called before Castro::variableSetUp

void Radiation::read_static_params()
{
  ParmParse pp("radiation");

  {
    int solver_type = Radiation::SolverType;
    pp.get("SolverType", solver_type);
    SolverType = static_cast<Solver_Type>(solver_type);
  }

  if (Radiation::SolverType == Radiation::MGFLDSolver) {
    do_multigroup = 1;
  }
  else {
    do_multigroup = 0;
  }

  // check that we're not using a single group solver if NGROUPS > 1
  if (!do_multigroup && Radiation::nGroups > 1) {
      amrex::Error("Radiation::nGroups > 1 but using single group solver");
  }

  if (Radiation::SolverType == Radiation::MGFLDSolver) {
    accelerate = 2;
  }
  pp.query("accelerate", accelerate);

  if (Radiation::SolverType == Radiation::SGFLDSolver ||
      Radiation::SolverType == Radiation::MGFLDSolver ) {
    Radiation::rad_hydro_combined = 1;
    pp.query("rad_hydro_combined", rad_hydro_combined);
  }

  if (SolverType == SGFLDSolver) {
    radiation::fspace_advection_type = 1;
    if (radiation::comoving) {
      Er_Lorentz_term = 0;
    }
    else {
      Er_Lorentz_term = 1;
      pp.query("Er_Lorentz_term", Er_Lorentz_term);
    }
  }

  pp.query("filter_lambda_T", filter_lambda_T);
  filter_lambda_S = filter_lambda_T - 1;
  pp.query("filter_lambda_S", filter_lambda_S);
  if (filter_lambda_T > 4) {
    amrex::Error("filter_lambda_T > 4");
  }
  if (filter_lambda_T < 0) {
    amrex::Error("filter_lambda_T < 0");
  }
  if (filter_lambda_T > 0) {
    if (filter_lambda_S >= filter_lambda_T) {
      amrex::Error("Invalid filter_lambda_S; S must be less than T when T > 0.");
    }
  }

  pp.query("filter_prim_int", filter_prim_int);
  pp.query("filter_prim_T", filter_prim_T);
  filter_prim_S = filter_prim_T - 1;
  pp.query("filter_prim_S", filter_prim_S);
  if (filter_prim_T > 4) {
    amrex::Error("filter_prim_T > 4");
  }
  if (filter_prim_T < 0) {
    amrex::Error("filter_prim_T < 0");
  }
  if (filter_prim_T > 0) {
    if (filter_prim_S >= filter_prim_T) {
      amrex::Error("Invalid filter_prim_S; S must be less than T when T > 0.");
    }
  }


  if (Radiation::SolverType == Radiation::MGFLDSolver) {

    Radiation::nGroups = NGROUPS;

    // sanity check -- in the old style, we allowed the number of groups to be
    // set at compile time, so ensure that, if the user set this, it matches what
    // we expect
    int test_groups = -1;
    pp.query("nGroups", test_groups);

    if (test_groups > 0 && test_groups != Radiation::nGroups) {
      amrex::Error("you set the number of groups at runtime, but this does not match the compiled value");
    }

  }
  else if (Radiation::SolverType != Radiation::SingleGroupSolver &&
           Radiation::SolverType != Radiation::SGFLDSolver) {
      amrex::Error("Unknown Radiation::SolverType");
  }


  // set up the extra plot variables
  {
      if (radiation::plot_lambda) {
          icomp_lambda = plotvar_names.size();

          if (!do_multigroup || radiation::limiter == 0) {
              plotvar_names.push_back("lambda");
          } else {
              for (int g=0; g<nGroups; ++g) {
                  std::ostringstream ss;
                  ss << "lambda" << g;
                  plotvar_names.push_back(ss.str());
              }
          }
      }
      if (radiation::plot_kappa_p) {
          icomp_kp = plotvar_names.size();
          if (!do_multigroup) {
              plotvar_names.push_back("kappa_P");
          } else {
              for (int g=0; g<nGroups; ++g) {
                  std::ostringstream ss;
                  ss << "kappa_P" << g;
                  plotvar_names.push_back(ss.str());
              }
          }
      }
      if (radiation::plot_kappa_r) {
          icomp_kr = plotvar_names.size();
          if (!do_multigroup) {
              plotvar_names.push_back("kappa_R");
          } else {
              for (int g=0; g<nGroups; ++g) {
                  std::ostringstream ss;
                  ss << "kappa_R" << g;
                  plotvar_names.push_back(ss.str());
              }
          }
      }
      if (radiation::plot_lab_Er) {
          icomp_lab_Er = plotvar_names.size();
          if (!do_multigroup) {
              plotvar_names.push_back("Erlab");
          } else {
              for (int g=0; g<nGroups; ++g) {
                  std::ostringstream ss;
                  ss << "Erlab" << g;
                  plotvar_names.push_back(ss.str());
              }
          }
      }
      if (radiation::plot_lab_flux) {
          icomp_lab_Fr = plotvar_names.size();
          std::string frame = "lab";
          Vector<std::string> dimname;
          dimname.push_back("x");
          dimname.push_back("y");
          dimname.push_back("z");
          if (!do_multigroup) {
              for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
                  std::ostringstream ss;
                  ss << "Fr" << frame << dimname[idim];
                  plotvar_names.push_back(ss.str());
              }
          } else {
              for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
                  for (int g=0; g<nGroups; ++g) {
                      std::ostringstream ss;
                      ss << "Fr" << frame << g << dimname[idim];
                      plotvar_names.push_back(ss.str());
                  }
              }
          }
      }
      if (radiation::plot_com_flux) {
          icomp_com_Fr = plotvar_names.size();
          std::string frame = "com";
          Vector<std::string> dimname;
          dimname.push_back("x");
          dimname.push_back("y");
          dimname.push_back("z");
          if (!do_multigroup) {
              for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
                  std::ostringstream ss;
                  ss << "Fr" << frame << dimname[idim];
                  plotvar_names.push_back(ss.str());
              }
          } else {
              for (int idim=0; idim<AMREX_SPACEDIM; ++idim) {
                  for (int g=0; g<nGroups; ++g) {
                      std::ostringstream ss;
                      ss << "Fr" << frame << g << dimname[idim];
                      plotvar_names.push_back(ss.str());
                  }
              }
          }
      }
      nplotvar = plotvar_names.size();
  }

  amrex::Print() << "radiation initialized, nGroups = " << Radiation::nGroups << std::endl;

}

Radiation::Radiation(Amr* Parent, Castro* castro, int restart)
  : parent(Parent)
{
  // castro is passed in, rather than obtained from parent, because this
  // routine will be called in some cases before any AmrLevels have
  // been installed into the parent's array of levels.

  ParmParse pp("radiation");

  do_sync = 1; pp.query("do_sync", do_sync);

  {

    clight = C::c_light;
    hPlanck = C::hplanck;
    kBoltz = C::k_B;
    Avogadro = C::n_A;
    convert_MeV_erg = 1.e6_rt * C::ev2erg;

    aRad = 4.*C::sigma_SB / C::c_light;

    c        = clight;
    sigma    = C::sigma_SB;

    if (!do_multigroup) {
        // In single group and abstract test problems we can play with
        // c and sigma independent of physical reality, but messing with
        // them in a multigroup problem is likely to be bad.
        pp.query("c", c);
        pp.query("sigma", sigma);
    }

    // Set Hypre flux factors here. Since this only occurs once,
    // every instance of Hypre must use the same factor (or
    // be responsible for changing it internally).

    HypreABec::fluxFactor() = c;
    HypreMultiABec::fluxFactor() = c;

  }

  radtoE = 1.0;
  //    radtoJ = c/(4.*M_PI);
  Etorad = 1.0;
  radfluxtoF = 1.0;

  reltol   = 1.e-6;          pp.query("reltol", reltol);
  if (SolverType == SGFLDSolver || SolverType == MGFLDSolver) {
    abstol = 0.0;
  }
  else {
    abstol   = 1.e-6;
  }
  pp.query("abstol", abstol);
  maxiter  = 50;             pp.query("maxiter", maxiter);
  miniter  =  1;             pp.query("miniter", miniter);
  convergence_check_type = 0;
  pp.query("convergence_check_type", convergence_check_type);
  if (SolverType == SGFLDSolver && radiation::limiter == 1) {
    amrex::Abort("SGFLDSolver does not support limiter = 1");
  }
  if (SolverType == MGFLDSolver && radiation::limiter == 1) {
    amrex::Abort("MGFLDSolver does not support limiter = 1");
  }

  inner_update_limiter = 0;
  pp.query("inner_update_limiter", inner_update_limiter);

  update_opacity    = 1000;

  if (SolverType == SGFLDSolver || SolverType == MGFLDSolver) {
    update_planck     = 1000;
    update_rosseland  = 1000;
    update_limiter    = 1000;
  }
  else {
    update_planck     = 50;
    update_rosseland  = 50;
    update_limiter    = 4;
  }
  pp.query("update_planck", update_planck);
  pp.query("update_rosseland", update_rosseland);
  pp.query("update_opacity", update_opacity);
  pp.query("update_limiter", update_limiter);

  dT  = 1.0;                 pp.query("delta_temp", dT);

  // for inner iterations of neutrino J equation
  relInTol = 1.e-4;          pp.query("relInTol", relInTol);
  if (SolverType == SGFLDSolver || SolverType == MGFLDSolver) {
    absInTol = 0.0;
  }
  else {
    absInTol = 1.e-4;
  }
  pp.query("absInTol", absInTol);
  maxInIter = 30;            pp.query("maxInIter", maxInIter);
  minInIter =  1;            pp.query("minInIter", minInIter);

  skipAccelAllowed = 0;
  pp.query("skipAccelAllowed", skipAccelAllowed);

  matter_update_type = 0;
  pp.query("matter_update_type", matter_update_type);

  n_bisect = 1000;
  pp.query("n_bisect", n_bisect);
  dedT_fac = 1.0;
  pp.query("dedT_fac", dedT_fac);

  inner_convergence_check = 2;
  pp.query("inner_convergence_check", inner_convergence_check);

  delta_e_rat_dt_tol = 100.0;
  pp.query("delta_e_rat_dt_tol", delta_e_rat_dt_tol);
  delta_T_rat_dt_tol = 100.0;
  pp.query("delta_T_rat_dt_tol", delta_T_rat_dt_tol);

  underfac = 1.0;    pp.query("underfac", underfac);

  use_WiensLaw = 0;
  pp.query("use_WiensLaw", use_WiensLaw);
  Tf_Wien = -1.0;
  pp.query("Tf_Wien", Tf_Wien);

  verbose  = 0;      pp.query("v", verbose);  pp.query("verbose", verbose);


  do_kappa_stm_emission = 0;
  pp.query("do_kappa_stm_emission", do_kappa_stm_emission);

  use_dkdT = 0;
  pp.query("use_dkdT", use_dkdT);

  if (verbose > 2) {
    Vector<int> temp;
    if (pp.queryarr("spot",temp,0,AMREX_SPACEDIM)) {
      IntVect tempi(temp);
      spot = tempi;
    }
    if (ParallelDescriptor::IOProcessor()) std::cout << "Spot: " << spot << std::endl;
  }

  if (verbose > 0 && ParallelDescriptor::IOProcessor()) {
    std::cout << "Creating Radiation object" << std::endl;
  }
  if (verbose >= 1 && ParallelDescriptor::IOProcessor()) {
    std::cout << "processors = " << ParallelDescriptor::NProcs() << std::endl;
    std::cout << "do_sync            = " << do_sync << std::endl;

    std::cout << "c        = " << c << std::endl;
    std::cout << "sigma    = " << sigma << std::endl;
    std::cout << "reltol   = " << reltol << std::endl;
    std::cout << "abstol   = " << abstol << std::endl;
    std::cout << "maxiter  = " << maxiter << std::endl;
    std::cout << "relInTol = " << relInTol << std::endl;
    std::cout << "absInTol = " << absInTol << std::endl;
    std::cout << "maxInIter = " << maxInIter << std::endl;
    std::cout << "delta_e_rat_dt_tol = " << delta_e_rat_dt_tol << std::endl;
    std::cout << "delta_T_rat_dt_tol = " << delta_T_rat_dt_tol << std::endl;
    std::cout << "limiter  = " << radiation::limiter << std::endl;
    std::cout << "closure  = " << radiation::closure << std::endl;
    std::cout << "update_limiter   = " << update_limiter << std::endl;
    std::cout << "update_planck    = " << update_planck << std::endl;
    std::cout << "update_rosseland = " << update_rosseland << std::endl;
    std::cout << "delta_temp = " << dT << std::endl;
    std::cout << "underfac = " << underfac << std::endl;
    std::cout << "do_multigroup = " << do_multigroup << std::endl;
    std::cout << "accelerate = " << accelerate << std::endl;
    std::cout << "verbose  = " << verbose << std::endl;
    if (SolverType == SingleGroupSolver) {
      std::cout << "SolverType = 0: SingleGroupSolver " << std::endl;
    }
    else if (SolverType == SGFLDSolver) {
      std::cout << "SolverType = 5: SGFLDSolver " << std::endl;
    }
    else if (SolverType == MGFLDSolver) {
      std::cout << "SolverType = 6: MGFLDSolver " << std::endl;
    }
    if (SolverType == MGFLDSolver || SolverType == SGFLDSolver) {
      std::cout << "rad_hydro_combined = " << rad_hydro_combined << std::endl;
      std::cout << "comoving = " << radiation::comoving << std::endl;
    }
    if (SolverType == MGFLDSolver) {
      std::cout << "fspace_advection_type = " << radiation::fspace_advection_type << std::endl;
    }
    if (SolverType == SGFLDSolver && radiation::comoving == 0) {
      std::cout << "Er_Lorentz_term = " << Er_Lorentz_term << std::endl;
    }
  }

  if (do_multigroup) {
    get_groups(verbose);
  }
  else {
    // xnu is a dummy for single group
    xnu.resize(2, 1.0);
    nugroup.resize(1, 1.0);
  }

  // current implementation of the Radiation boundary condition reads
  // incoming flux information in the RadBndry constructor.  we just
  // set the boundary condition type here:

  Vector<int> lo_bc(AMREX_SPACEDIM), hi_bc(AMREX_SPACEDIM);
  pp.getarr("lo_bc",lo_bc,0,AMREX_SPACEDIM);
  pp.getarr("hi_bc",hi_bc,0,AMREX_SPACEDIM);
  for (int i = 0; i < AMREX_SPACEDIM; i++) {
    rad_bc.setLo(i,lo_bc[i]);
    rad_bc.setHi(i,hi_bc[i]);
    if (verbose > 1 && ParallelDescriptor::IOProcessor()) {
      std::cout << "dimension " << i << " rad boundary conditions = "
            << lo_bc[i] << ", " << hi_bc[i] << std::endl;
    }
  }

  // size flux register arrays and persistent MultiFabs:

  int levels = parent->maxLevel() + 1; // maximum allowable number of levels

  flux_cons.resize(levels);
  flux_cons_old.resize(levels);
  flux_trial.resize(levels);

  dflux.resize(levels);

  plotvar.resize(levels);

  delta_t_old.resize(levels, 0.0);

  delta_e_rat_level.resize(levels, 0.0);
  delta_T_rat_level.resize(levels, 0.0);

  pp.query("pure_hydro", pure_hydro);

}

void Radiation::regrid(int level, const BoxArray& grids, const DistributionMapping& dmap)
{
  BL_PROFILE("Radiation::Regrid");
  if (verbose > 0 && ParallelDescriptor::IOProcessor()) {
    std::cout << "Regridding radiation object at level " << level
         << "..." << std::endl;
  }
  if (level > 0) {
    IntVect crse_ratio = parent->refRatio(level-1);

    flux_cons[level].reset(new FluxRegister(grids, dmap, crse_ratio, level, nGroups));
    flux_cons[level]->setVal(0.0);

    // For deferred sync, flux_cons_old does not need to be defined here.
    // It will be set in the deferred_sync_setup routine.

    flux_trial[level].reset(new FluxRegister(grids, dmap, crse_ratio, level, nGroups));
    flux_trial[level]->setVal(0.0);

  }

  dflux[level].reset(new MultiFab(grids, dmap, 1, 0));

  if (nplotvar > 0) {
      plotvar[level].reset(new MultiFab(grids, dmap, nplotvar, 0));
      plotvar[level]->setVal(0.0);
  }

  // This array will not be used on the finest level.  I create it here,
  // though, in case a finer level is created before this level is next
  // regridded:

  if (verbose > 1 && ParallelDescriptor::IOProcessor()) {
    std::cout << "                                         done" << std::endl;
  }
}

void Radiation::close(int level)
{
  // Only appropriate when a level disappears, otherwise see regrid:
  if (level > parent->finestLevel()) {
    if (verbose > 0 && ParallelDescriptor::IOProcessor()) {
      std::cout << "Clearing radiation object at level " << level
            << "..." << std::endl;
    }
    flux_cons[level].reset();
    flux_trial[level].reset();

    // flux_cons_old is not deleted here because if it exists it still
    // has energy in it.  It will be deleted once it is finally used.

    dflux[level].reset();

    plotvar[level].reset();

    if (verbose > 1 && ParallelDescriptor::IOProcessor()) {
      std::cout << "                                       done" << std::endl;
    }

    // When a level is closed, then we have to make sure there is
    // consistent flux information in level-1 to start the next
    // time step.  We need to do this if we are in the middle of
    // a level-2 timestep, or if there is no level-2.  In either
    // case the operations in init_flux are appropriate.

    BL_ASSERT(level > 0);
    int clev = level - 1;

    // Check in case more than one level was removed:
    if (clev == parent->finestLevel()) {
      int ncycle    = parent->nCycle(clev);
      int iteration = parent->levelSteps(clev);
      iteration = (clev > 0) ? iteration % ncycle : iteration;
      if (iteration > 0) {
        init_flux(clev, ncycle);
      }
    }
  }
}

void Radiation::restart(int level, const BoxArray& grids,
                        const DistributionMapping& dmap,
                        const std::string& dir, std::istream& is)
{
  //
  // With the deferred sync option, we have to restart the rad flux register
  //

  std::string Path, aString;

  do {
    is >> aString;
    if (aString.find("delta_e_rat") == 0) {
      is >> delta_e_rat_level[level];
    }
    else if (aString.find("delta_T_rat") == 0) {
      is >> delta_T_rat_level[level];
    }
    else {
      Path = aString;
    }
  } while (Path.empty());

  //
  // Read flux register only if present in the chkfile.
  //
  std::string Flag;
  is >> Flag;

  if (Flag == "Present") {
    BL_ASSERT(level > 0);
    //
    // Read delta_t associated with this flux information.
    //
    is >> delta_t_old[level-1];
    //
    // Prepend the name of the chkfile directory.
    //
    std::string FullPathName = dir;
    if (!dir.empty() && dir[dir.length()-1] != '/')
        FullPathName += '/';
    FullPathName += Path;
    //
    // Input conservation flux register.
    //
    const IntVect& crse_ratio = parent->refRatio(level-1);
    flux_cons_old[level].reset(new FluxRegister(grids, dmap, crse_ratio, level, nGroups));
    flux_cons_old[level]->read(FullPathName, is);
  }
}

void Radiation::checkPoint(int level,
                           const std::string& dir,
                           std::ostream&  os,
                           VisMF::How     how)
{
  //
  // With the deferred sync option, we have to restart the rad flux register
  //

  char buf[64];

  //
  // Write deltas to header for timestep control.
  //
  if (ParallelDescriptor::IOProcessor()) {
    int oldprec = os.precision(20);
    auto DeltaString = std::format("delta_e_rat_level[{}]= ", level);
    os << DeltaString << delta_e_rat_level[level] << '\n';
    DeltaString = std::format("delta_T_rat_level[{}]= ", level);
    os << DeltaString << delta_T_rat_level[level] << '\n';
    os.precision(oldprec);
  }

  // Path name construction stolen from AmrLevel::checkPoint

  std::string Level = std::format("Level_{}", level);

  //
  // Write name of conservation flux register to header.
  //

  std::string PathNameInHeader = Level;
  PathNameInHeader += "/RadFlux";
  if (ParallelDescriptor::IOProcessor()) {
    os << PathNameInHeader;
  }

  if (flux_cons_old[level]) {
    //
    // Conservation flux register exists.
    //
    if (ParallelDescriptor::IOProcessor()) {
      BL_ASSERT(level > 0);
      int oldprec = os.precision(20);
      os << " Present  " << delta_t_old[level-1] << '\n';
      os.precision(oldprec);
    }
    //
    // This is the full pathname for the written FluxRegister.
    //
    std::string FullPathName = dir;
    if (!FullPathName.empty() &&
        FullPathName[FullPathName.length()-1] != '/') {
      FullPathName += '/';
    }
    FullPathName += Level;
    FullPathName += "/RadFlux";
    //
    // Output conservation flux register.
    //
    flux_cons_old[level]->write(FullPathName, os /* , how */ );
  }
  else {
    //
    // Conservation flux register does not exist.
    //
    if (ParallelDescriptor::IOProcessor()) {
      os << " Absent\n";
    }
  }
}

void Radiation::post_init(int level)
{
    return;
}

void Radiation::pre_timestep(int level)
{
  BL_PROFILE("Radiation::pre_timestep");
  int fine_level = parent->finestLevel();
  int ncycle     = parent->nCycle(level);

  static int done = 0;
  if (level < fine_level) {
      // For deferred sync, we may have moved flux_cons into flux_cons_old
      // and not rebuilt flux_cons, so check for that here.

      int flevel = level + 1;
      if (!flux_cons[flevel]) {
          const BoxArray& grids = parent->getLevel(flevel).boxArray();
          const DistributionMapping& dmap = parent->getLevel(flevel).DistributionMap();
          const IntVect& crse_ratio = parent->refRatio(level);
          flux_cons[flevel].reset(new FluxRegister(grids, dmap, crse_ratio, flevel, nGroups));
          flux_cons[flevel]->setVal(0.0);
      }
  }

  // If we aren't doing a multilevel solve, we still need to initialize
  // dflux and load the flux registers at each level.  For a single-level
  // calculation this only needs to be done once per run, whether
  // at initialization or at restart.  We can't trust iteration to
  // tell us, since at restart iteration is not 0.
  int iteration = parent->levelSteps(level);
  iteration = (level > 0) ? iteration % ncycle : done;
  if (level < fine_level || iteration == 0) {
      init_flux(level, ncycle);
      done = 1;
  }
}

void Radiation::init_flux(int level, int ncycle)
{
  BL_PROFILE("Radiation::init_flux");
  if (verbose > 0 && ParallelDescriptor::IOProcessor()) {
    std::cout << "Radiation flux initialization at level " << level
         << "..." << std::endl;
  }

  int fine_level = parent->finestLevel();

  dflux[level]->setVal(0.0);

  if (level < fine_level) {
      flux_cons[level+1]->setVal(0.0);
  }

  if (verbose > 1 && ParallelDescriptor::IOProcessor()) {
    std::cout << "                                           done" << std::endl;
  }
}

// Overwrites temperature with exchange term, exch = temp on input:

// This version used by single group, multigroup

void Radiation::compute_exchange(MultiFab& exch,
                                 MultiFab& Er,
                                 MultiFab& fkp, int igroup)
{
    BL_PROFILE("Radiation::compute_exchange");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(exch, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        auto exch_arr = exch[mfi].array();
        auto Er_arr = Er[mfi].array();
        auto fkp_arr = fkp[mfi].array();

        Real lsigma = sigma;
        Real lc = c;

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            exch_arr(i,j,k) = fkp_arr(i,j,k) * (4.e0_rt * lsigma * std::pow(exch_arr(i,j,k), 4) - lc * Er_arr(i,j,k));
        });
    }
}

void Radiation::compute_eta(MultiFab& eta, MultiFab& etainv,
                            MultiFab& state, MultiFab& temp,
                            MultiFab& fkp, MultiFab& Er,
                            Real delta_t, Real c,
                            Real underrel, int lag_planck, int igroup)
{
    BL_PROFILE("Radiation::compute_eta");

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox c_v;

        for (MFIter mfi(eta, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.tilebox();

            if (lag_planck) {

                Array4<Real> const eta_arr = eta.array(mfi);
                Array4<Real> const fkp_arr = fkp.array(mfi);
                AMREX_PARALLEL_FOR_3D(bx, i, j, k,
                                      { eta_arr(i,j,k) = fkp_arr(i,j,k); });

            }
            else {

                // This is the only case where we need a direct call for
                // Planck mean as a function of temperature.

                auto eta_arr = eta[mfi].array();
                auto state_arr = state[mfi].array();

                const Real nu = nugroup[igroup];
                const Real dT_loc = dT;

                amrex::ParallelFor(bx,
                [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
                {
                    Real rho = state_arr(i,j,k,URHO);
                    Real temp = state_arr(i,j,k,UTEMP) + dT_loc;
                    Real Ye;
                    if (NumAux > 0)  {
                        Ye = state_arr(i,j,k,UFX);
                    } else {
                        Ye = 0.e0_rt;
                    }

                    Real kp, kr;
                    bool comp_kp = true;
                    bool comp_kr = false;
                    opacity(kp, kr, rho, temp, Ye, nu, comp_kp, comp_kr);

                    eta_arr(i,j,k,igroup) = kp;
                });

            }

            c_v.resize(bx);
            Elixir c_v_elix = c_v.elixir();

            get_c_v(c_v, temp[mfi], state[mfi], bx);

            auto eta_arr = eta[mfi].array();
            auto etainv_arr = etainv[mfi].array();
            auto frho_arr = state[mfi].array(URHO);
            auto temp_arr = temp[mfi].array();
            auto c_v_arr = c_v.array();
            auto fkp_arr = fkp[mfi].array();
            auto Er_arr = Er[mfi].array(igroup);

            Real dT_loc = dT;

            const Real fac1 = 16.e0_rt * sigma * delta_t;
            const Real fac0 = 0.25e0_rt * fac1 / dT;
            const Real fac2 = delta_t * c / dT;

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
            {
                Real d;

                if (lag_planck != 0)
                {
                    // assume eta and fkp are the same
                    d = fac1 * fkp_arr(i,j,k) * std::pow(temp_arr(i,j,k), 3);
                }
                else
                {
                    d = fac0 * (eta_arr(i,j,k) * std::pow(temp_arr(i,j,k) + dT_loc, 4) -
                                fkp_arr(i,j,k) * std::pow(temp_arr(i,j,k), 4)) -
                        fac2 * (eta_arr(i,j,k) - fkp_arr(i,j,k)) * Er_arr(i,j,k);
                    // alternate form, sometimes worse, sometimes better:
                    //   d = fac1 * fkp_arr(i,j,k) * std::pow(temp_arr(i,j,k), 3) +
                    //       fac0 * (eta_arr(i,j,k) - fkp_arr(i,j,k)) * std::pow(temp(i,j,k), 4) -
                    //       fac2 * (eta_arr(i,j,k) - fkp_arr(i,j,k)) * Er_arr(i,j,k);
                    // another alternate form (much worse):
                    //   d = fac1 * fkp_arr(i,j,k) * std::pow(temp_arr(i,j,k) + dtTloc, 3) +
                    //       fac0 * (eta_arr(i,j,k) - fkp_arr(i,j,k)) * std::pow(temp(i,j,k) + dT_loc, 4) -
                    //       fac2 * (eta_arr(i,j,k) - fkp_arr(i,j,k)) * Er_arr(i,j,k);
                }

                Real frc = frho_arr(i,j,k) * c_v_arr(i,j,k) + 1.0e-50_rt;
                eta_arr(i,j,k) = d / (d + frc);
                etainv_arr(i,j,k) = underrel * frc / (d + frc);
                eta_arr(i,j,k) = 1.e0_rt - etainv_arr(i,j,k);
                // eta_arr(i,j,k) = 1.e0_rt - underrel * (1.e0_rt - eta_arr(i,j,k));
            });
        }
    }
}

void Radiation::internal_energy_update(Real& relative, Real& absolute,
                                       MultiFab& frhoes,
                                       MultiFab& frhoem,
                                       MultiFab& eta,
                                       MultiFab& etainv,
                                       MultiFab& dflux_old,
                                       MultiFab& dflux_new,
                                       MultiFab& exch,
                                       Real delta_t)
{
  BL_PROFILE("Radiation::internal_energy_update");

  ReduceOps<ReduceOpMax, ReduceOpMax> reduce_op;
  ReduceData<Real, Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  Real theta = 1.0;

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(eta, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

      const auto eta_arr = eta[mfi].array();
      const auto etainv_arr = etainv[mfi].array();
      const auto frhoem_arr = frhoem[mfi].array();
      const auto exch_arr = exch[mfi].array();
      const auto dfo = dflux_old[mfi].array();
      const auto dfn = dflux_new[mfi].array();
      auto frhoes_arr = frhoes[mfi].array();

      reduce_op.eval(bx, reduce_data,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
      {
          Real chg = 0.e0_rt;
          Real tot = 0.e0_rt;

          Real tmp = eta_arr(i,j,k) * frhoes_arr(i,j,k) +
                     etainv_arr(i,j,k) *
                     (frhoem_arr(i,j,k) -
                      delta_t * ((1.e0_rt - theta) *
                                 (dfo(i,j,k) - dfn(i,j,k)) +
                                 exch_arr(i,j,k)));

          chg = std::abs(tmp - frhoes_arr(i,j,k));
          tot = std::abs(frhoes_arr(i,j,k));

          frhoes_arr(i,j,k) = tmp;

          Real absres = chg;
          Real relres = chg / (tot + 1.e-50_rt);

          return {relres, absres};
      });
  }

  ReduceTuple hv = reduce_data.value();

  relative = amrex::get<0>(hv);
  absolute = amrex::get<1>(hv);

  ParallelDescriptor::ReduceRealMax(relative);
  ParallelDescriptor::ReduceRealMax(absolute);
}

void Radiation::internal_energy_update(Real& relative, Real& absolute,
                                       MultiFab& frhoes,
                                       MultiFab& frhoem,
                                       MultiFab& eta,
                                       MultiFab& etainv,
                                       MultiFab& dflux_old,
                                       MultiFab& dflux_new,
                                       MultiFab& exch,
                                       MultiFab& Dterm,
                                       Real delta_t)
{
  BL_PROFILE("Radiation::internal_energy_update_d");

  relative = 0.0;
  absolute = 0.0;

  ReduceOps<ReduceOpMax, ReduceOpMax> reduce_op;
  ReduceData<Real, Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  const Real theta = 1.0;
  const Real tiny = 1.e-50_rt;

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(eta,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box &reg = mfi.tilebox();

      const auto eta_arr = eta[mfi].array();
      const auto etainv_arr = etainv[mfi].array();
      const auto frhoem_arr = frhoem[mfi].array();
      const auto exch_arr = exch[mfi].array();
      const auto dfo = dflux_old[mfi].array();
      const auto dfn = dflux_new[mfi].array();
      const auto dterm_arr = Dterm[mfi].array();
      auto frhoes_arr = frhoes[mfi].array();

      reduce_op.eval(reg, reduce_data,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
      {

          Real chg = 0.e0_rt;
          Real tot = 0.e0_rt;

          Real tmp = eta_arr(i,j,k) * frhoes_arr(i,j,k) +
              etainv_arr(i,j,k) * (frhoem_arr(i,j,k) -
                                   delta_t * ((1.e0_rt - theta) *
                                              (dfo(i,j,k) - dfn(i,j,k)) +
                                              exch_arr(i,j,k))) +
              delta_t * dterm_arr(i,j,k);

          chg = std::abs(tmp - frhoes_arr(i,j,k));
          tot = std::abs(frhoes_arr(i,j,k));

          frhoes_arr(i,j,k) = tmp;

          Real absres = chg;
          Real relres = chg / (tot + tiny);

          return {relres, absres};

      });

  }

  ReduceTuple hv = reduce_data.value();

  relative = amrex::get<0>(hv);
  absolute = amrex::get<1>(hv);

  ParallelDescriptor::ReduceRealMax(relative);
  ParallelDescriptor::ReduceRealMax(absolute);

}

void Radiation::nonconservative_energy_update(Real& relative, Real& absolute,
                                              MultiFab& frhoes,
                                              MultiFab& frhoem,
                                              MultiFab& eta,
                                              MultiFab& etainv,
                                              MultiFab& Er_new,
                                              MultiFab& dflux_old,
                                              MultiFab& dflux_new,
                                              MultiFab& temp,
                                              MultiFab& fkp,
                                              MultiFab& state,
                                              Real delta_t)
{
  BL_PROFILE("Radiation::nonconservative_energy_update");

  ReduceOps<ReduceOpMax, ReduceOpMax> reduce_op;
  ReduceData<Real, Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  const Real theta = 1.0;
  const Real tiny = 1.e-50_rt;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
      FArrayBox c_v;

      for (MFIter mfi(eta,true); mfi.isValid(); ++mfi) {
          const Box &reg = mfi.tilebox();

          c_v.resize(reg);
          Elixir elix_c_v = c_v.elixir();

          get_c_v(c_v, temp[mfi], state[mfi], reg);

          const auto state_arr = state[mfi].array();
          const auto temp_arr = temp[mfi].array();
          const auto fkp_arr = fkp[mfi].array();
          const auto er_arr = Er_new[mfi].array();
          const auto eta_arr = eta[mfi].array();
          const auto etainv_arr = etainv[mfi].array();
          const auto frhoem_arr = frhoem[mfi].array();
          const auto frhoes_arr = frhoes[mfi].array();
          const auto dfo_arr = dflux_old[mfi].array();
          const auto dfn_arr = dflux_new[mfi].array();
          const auto c_v_arr = c_v.array();

          reduce_op.eval(reg, reduce_data,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
          {

              Real chg = 0.e0_rt;
              Real tot = 0.e0_rt;

              Real frhocv = state_arr(i,j,k,URHO) * c_v_arr(i,j,k);

              Real dbdt = 16.0_rt * C::sigma_SB * std::pow(temp_arr(i,j,k), 3);
              Real b = 4.0_rt * C::sigma_SB * std::pow(temp_arr(i,j,k), 4);
              Real exch = fkp_arr(i,j,k) * (b - C::c_light * er_arr(i,j,k));
              Real tmp = eta_arr(i,j,k) * frhoes_arr(i,j,k) +
                  etainv_arr(i,j,k) * (frhoem_arr(i,j,k) -
                                       delta_t * ((1.0_rt - theta) * (dfo_arr(i,j,k) - dfn_arr(i,j,k)) + exch));

              if (frhocv > tiny && tmp > frhoes_arr(i,j,k)) {
                  Real db = (tmp - frhoes_arr(i,j,k)) * dbdt / frhocv;
                  tmp = std::pow((b + db) / (4.0_rt * C::sigma_SB), 0.25_rt);
                  tmp = frhoes_arr(i,j,k) + frhocv * (tmp - temp_arr(i,j,k));
              }

              chg = std::abs(tmp - frhoes_arr(i,j,k));
              tot = std::abs(frhoes_arr(i,j,k));
              frhoes_arr(i,j,k) = tmp;

              Real absres = chg;
              Real relres = chg / (tot + tiny);

              return {absres, relres};
          });

      }

  } // OpenMP

  ReduceTuple hv = reduce_data.value();

  absolute = amrex::get<0>(hv);
  relative = amrex::get<1>(hv);

  ParallelDescriptor::ReduceRealMax(relative);
  ParallelDescriptor::ReduceRealMax(absolute);
}

void Radiation::state_update(MultiFab& state, MultiFab& frhoes)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        auto state_arr = state[mfi].array();
        auto frhoes_arr = frhoes[mfi].array();

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            Real kin = state_arr(i,j,k,UEDEN) - state_arr(i,j,k,UEINT);
            state_arr(i,j,k,UEINT) = frhoes_arr(i,j,k);
            state_arr(i,j,k,UEDEN) = frhoes_arr(i,j,k) + kin;

            // frhoes will be overwritten with temperature here

            if (state_arr(i,j,k,UEINT) <= 0.e0_rt)
            {
                frhoes_arr(i,j,k) = small_temp;
            }
            else
            {
                Real rhoInv = 1.e0_rt / state_arr(i,j,k,URHO);

                eos_re_t eos_state;
                eos_state.rho = state_arr(i,j,k,URHO);
                eos_state.T   = state_arr(i,j,k,UTEMP);
                eos_state.e   = state_arr(i,j,k,UEINT) * rhoInv;
                for (int n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] = state_arr(i,j,k,UFS+n) * rhoInv;
                }
#if NAUX_NET > 0
                for (int n = 0; n < NumAux; ++n) {
                    eos_state.aux[n] = state_arr(i,j,k,UFX+n) * rhoInv;
                }
#endif

                eos(eos_input_re, eos_state);

                frhoes_arr(i,j,k) = eos_state.T;

                state_arr(i,j,k,UTEMP) = frhoes_arr(i,j,k);
            }
        });
    }
}


void Radiation::extrapolateBorders(MultiFab& f, int indx)
{
    BL_PROFILE("Radiation::extrapolateBorders");

    BL_ASSERT(f.nGrow() >= 1);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(f); mfi.isValid(); ++mfi) {
        // Note no tiling in the current implementation.
        const Box& bx = mfi.validbox();
        const Box& grownbx = amrex::grow(bx, 1);

        Array4<Real> const f_arr = f[mfi].array(indx);

        amrex::LoopOnCpu(grownbx, [=] (int i, int j, int k) noexcept
        {
            // Note that the results on the corners will be the same
            // regardless of which order we do the loop in.

            if (i == bx.smallEnd(0)) {
                f_arr(i-1,j,k) = 2.0_rt * f_arr(i,j,k) - f_arr(i+1,j,k);
            }

            if (i == bx.bigEnd(0)) {
                f_arr(i+1,j,k) = 2.0_rt * f_arr(i,j,k) - f_arr(i-1,j,k);
            }

#if AMREX_SPACEDIM >= 2
            if (j == bx.smallEnd(1)) {
                f_arr(i,j-1,k) = 2.0_rt * f_arr(i,j,k) - f_arr(i,j+1,k);
            }

            if (j == bx.bigEnd(1)) {
                f_arr(i,j+1,k) = 2.0_rt * f_arr(i,j,k) - f_arr(i,j-1,k);
            }
#endif

#if AMREX_SPACEDIM == 3
            if (k == bx.smallEnd(2)) {
                f_arr(i,j,k-1) = 2.0_rt * f_arr(i,j,k) - f_arr(i,j,k+1);
            }

            if (k == bx.bigEnd(2)) {
                f_arr(i,j,k+1) = 2.0_rt * f_arr(i,j,k) - f_arr(i,j,k-1);
            }
#endif
        });
    }
}


void Radiation::getBndryData(RadBndry& bd, MultiFab& Er,
                             Real time, int level)
{
  BL_PROFILE("Radiation::getBndryData");
  Castro      *castro = (Castro*)&parent->getLevel(level);
  const BoxArray& grids = castro->boxArray();
  const DistributionMapping& dmap = castro->DistributionMap();

  if (level == 0) {
    bd.setBndryValues(Er, Rad, 0, 1, rad_bc); // Rad=0
  }
  else {
    BoxArray cgrids(grids);
    IntVect crse_ratio = parent->refRatio(level-1);
    cgrids.coarsen(crse_ratio);
    BndryRegister crse_br(cgrids, dmap, 0, 1, 1, 1);
    crse_br.setVal(1.0e30);
    filBndry(crse_br, level-1, time);

    bd.setBndryValues(crse_br, 0, Er, Rad, 0, 1, crse_ratio, rad_bc);
  }

  // We do this last, in case Er has ghost cells which get written into
  // the boundary values:

  bd.setTime(time);

  bd.setBndryFluxConds(rad_bc);
}

void Radiation::getBndryDataMG(MGRadBndry& mgbd, MultiFab& Er,
                               Real time, int level)
{
  BL_PROFILE("Radiation::getBndryDataMG");
  Castro *castro = dynamic_cast<Castro*>(&parent->getLevel(level));
  const BoxArray& grids = castro->boxArray();
  const DistributionMapping& dmap = castro->DistributionMap();

  if(level == 0) {
    mgbd.setBndryValues(Er, 0, 0, Radiation::nGroups, rad_bc);
  }
  else {
    BoxArray cgrids(grids);
    IntVect crse_ratio = parent->refRatio(level-1);
    cgrids.coarsen(crse_ratio);
    BndryRegister crse_br(cgrids, dmap, 0, 1, 1, Radiation::nGroups);
    crse_br.setVal(1.0e30);
    filBndry(crse_br, level-1, time);

    mgbd.setBndryValues(crse_br, 0, Er, 0,
                        0, Radiation::nGroups, crse_ratio, rad_bc);
  }

  mgbd.setTime(time);
  mgbd.setBndryFluxConds(rad_bc);
}

void Radiation::getBndryDataMG_ga(MGRadBndry& mgbd, MultiFab& Er, int level)
{
  BL_PROFILE("Radiation::getBndryDataMG_ga");
  Castro *castro = dynamic_cast<Castro*>(&parent->getLevel(level));
  const BoxArray& grids = castro->boxArray();
  const DistributionMapping& dmap = castro->DistributionMap();

  if(level == 0) {
    mgbd.setBndryValues(Er, 0, 0, 1, rad_bc);
  }
  else {
    BoxArray cgrids(grids);
    IntVect crse_ratio = parent->refRatio(level-1);
    cgrids.coarsen(crse_ratio);
    BndryRegister crse_br(cgrids, dmap, 0, 1, 1, 1);
    crse_br.setVal(0.0);

    mgbd.setBndryValues(crse_br, 0, Er, 0,
                        0, 1, crse_ratio, rad_bc);
  }

  mgbd.setBndryFluxConds(rad_bc);
}

void Radiation::filBndry(BndryRegister& bdry, int level, Real time)
{
  BL_PROFILE("Radiation::filBndry");
  // in this routine "level" is the coarse level

  Castro      *castro = (Castro*)&parent->getLevel(level);
  const BoxArray& grids = castro->boxArray();
  const DistributionMapping& dmap = castro->DistributionMap();
  const Geometry& geom  = parent->Geom(level);

  Real old_time = castro->get_state_data(Rad_Type).prevTime();
  Real new_time = castro->get_state_data(Rad_Type).curTime();
  Real eps = (new_time > old_time) ? 0.001*(new_time - old_time) : 1.0;

  BL_ASSERT( (time > old_time-eps) && (time < new_time + eps));

  MultiFab& S_new = castro->get_new_data(Rad_Type);
  // the next line is OK even if S_old is not defined yet
  MultiFab& S_old = castro->get_old_data(Rad_Type);

  if (!geom.isAnyPeriodic()) {
    int n_ghost = 0;

    if (time > new_time - eps) {
      //bdry.copyFrom(S_new, n_ghost, 0, 0, 1);
      bdry.copyFrom(S_new, n_ghost, 0, 0, Radiation::nGroups);
    }
    else if (time < old_time + eps) {
      //bdry.copyFrom(S_old, n_ghost, 0, 0, 1);
      bdry.copyFrom(S_old, n_ghost, 0, 0, Radiation::nGroups);
    }
    else {
      Real a = (new_time - time) / (new_time - old_time);
      Real b = (time - old_time) / (new_time - old_time);
      //bdry.linComb(a, S_old, 0, b, S_new, 0, 0, 1, n_ghost);
      bdry.linComb(a, S_old, 0, b, S_new, 0, 0, Radiation::nGroups, n_ghost);
    }
  }
  else {

    // older version, ghost cell seems unnecessary (perhaps a holdover from
    // when we put boundary conditions in ghost cells):

    // (later) resurrected because ghost cell needed for periodic case:

    int need_old_data = (time <= new_time - eps);
    int need_new_data = (time >= old_time + eps);

    int n_grow = 1;

    MultiFab sold_tmp;
    MultiFab snew_tmp;

    if (need_old_data) {
      sold_tmp.define(grids, dmap, 1, n_grow);
      sold_tmp.setVal(0.0); // need legal numbers for linComb below
      MultiFab::Copy(sold_tmp, S_old, Rad, 0, 1, 0);
      sold_tmp.FillBoundary(geom.periodicity());
    }

    if (need_new_data) {
      snew_tmp.define(grids, dmap, 1, n_grow);
      snew_tmp.setVal(0.0); // need legal numbers for linComb below
      MultiFab::Copy(snew_tmp, S_new, Rad, 0, 1, 0);
      snew_tmp.FillBoundary(geom.periodicity());
    }

    int n_ghost = 1;

    if (time > new_time - eps) {
      bdry.copyFrom(snew_tmp, n_ghost, 0, 0, 1);
    }
    else if (time < old_time + eps) {
      bdry.copyFrom(sold_tmp, n_ghost, 0, 0, 1);
    }
    else {
      Real a = (new_time - time) / (new_time - old_time);
      Real b = (time - old_time) / (new_time - old_time);
      bdry.linComb(a, sold_tmp, 0, b, snew_tmp, 0, 0, 1, n_ghost);
    }
  }
}

void Radiation::get_c_v(FArrayBox& c_v, FArrayBox& temp, FArrayBox& state,
                        const Box& reg)
{
    BL_PROFILE("Radiation::get_c_v");

    auto c_v_arr = c_v.array();
    auto temp_arr = temp.array();
    auto state_arr = state.array();

    amrex::ParallelFor(reg,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        Real rhoInv = 1.e0_rt / state_arr(i,j,k,URHO);

        eos_re_t eos_state;
        eos_state.rho = state_arr(i,j,k,URHO);
        eos_state.T   = temp_arr(i,j,k);
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = state_arr(i,j,k,UFS+n) * rhoInv;
        }
#if NAUX_NET > 0
        for (int n = 0; n < NumAux; ++n) {
            eos_state.aux[n] = state_arr(i,j,k,UFX+n) * rhoInv;
        }
#endif

        eos(eos_input_rt, eos_state);

        c_v_arr(i,j,k) = eos_state.cv;
    });
    Gpu::synchronize();
}

void Radiation::get_frhoe(MultiFab& frhoe,
                          MultiFab& state)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter si(state,TilingIfNotGPU()); si.isValid(); ++si) {
        const Box& reg = si.tilebox();

        auto frhoe_arr = frhoe[si].array();
        auto state_arr = state[si].array();

        amrex::ParallelFor(reg,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            frhoe_arr(i,j,k) = state_arr(i,j,k,UEINT);
        });
    }
}

// temp contains frhoe on input:

void Radiation::get_planck_and_temp(MultiFab& fkp,
                                    MultiFab& temp,
                                    MultiFab& state,
                                    int igroup, Real delta_t)
{
    BL_PROFILE("Radiation::get_planck_and_temp");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        auto fkp_arr = fkp[mfi].array(igroup);
        auto temp_arr = temp[mfi].array();
        auto state_arr = state[mfi].array();

        // Get T from rhoe; overwrite temp with T

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            if (temp_arr(i,j,k) <= 0.e0_rt)
            {
                temp_arr(i,j,k) = small_temp;
            }
            else
            {
                Real rhoInv = 1.e0_rt / state_arr(i,j,k,URHO);

                eos_re_t eos_state;
                eos_state.rho = state_arr(i,j,k,URHO);
                eos_state.T   = state_arr(i,j,k,UTEMP);
                eos_state.e   = temp_arr(i,j,k) * rhoInv;
                for (int n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] = state_arr(i,j,k,UFS+n) * rhoInv;
                }
#if NAUX_NET > 0
                for (int n = 0; n < NumAux; ++n) {
                    eos_state.aux[n] = state_arr(i,j,k,UFX+n) * rhoInv;
                }
#endif

                eos(eos_input_re, eos_state);

                temp_arr(i,j,k) = eos_state.T;
            }
        });

        const Real nu = nugroup[igroup];

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            Real rho = state_arr(i,j,k,URHO);
            Real temp = state_arr(i,j,k,UTEMP);
            Real Ye;
            if (NumAux > 0)  {
                Ye = state_arr(i,j,k,UFX);
            } else {
                Ye = 0.e0_rt;
            }

            Real kp, kr;
            bool comp_kp = true;
            bool comp_kr = false;
            opacity(kp, kr, rho, temp, Ye, nu, comp_kp, comp_kr);

            fkp_arr(i,j,k) = kp;
        });

        int ncomp = temp[mfi].nComp();

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            const Real temp_floor = 1.e-10_rt;

            for (int n = 0; n < ncomp; ++n) {
                if (temp_arr(i,j,k,n) < temp_floor) {
                    temp_arr(i,j,k,n) = temp_floor;
                }
            }
        });
    }
}

// Uses filPatch to fill state data in a ghost cell around each grid
// so that kappa_r can be constructed everywhere.  Values across
// physical boundaries will not be used, however.

void Radiation::get_rosseland(MultiFab& kappa_r,
                              AmrLevel* castro,
                              int igroup)
{
  BL_PROFILE("Radiation::get_rosseland");

  BL_ASSERT(kappa_r.nGrow() == 1);

  if(Radiation::nGroups == 0) {
      kappa_r.setVal(0.);
      return;
  }

  Real time = castro->get_state_data(State_Type).curTime();
  MultiFab& S_new = castro->get_new_data(State_Type);
  int nstate = S_new.nComp();

  FillPatchIterator fpi(*castro, S_new, 1, time, State_Type, 0, nstate);
  MultiFab& state = fpi.get_mf();

  const Real nu = nugroup[igroup];

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
      FArrayBox frhoe;
      for(MFIter mfi(kappa_r, TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
          const Box& bx = mfi.growntilebox();
          frhoe.resize(bx, 1);
          Elixir frhoe_elix = frhoe.elixir();

          auto frhoe_arr = frhoe.array();
          auto state_arr = state[mfi].array();
          auto kpr = kappa_r[mfi].array();

          amrex::ParallelFor(bx,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
          {
              // frhoe will be overwritten with temperature here

              if (state_arr(i,j,k,UEINT) <= 0.e0_rt)
              {
                  frhoe_arr(i,j,k) = small_temp;
              }
              else
              {
                  Real rhoInv = 1.e0_rt / state_arr(i,j,k,URHO);

                  eos_re_t eos_state;
                  eos_state.rho = state_arr(i,j,k,URHO);
                  eos_state.T   = state_arr(i,j,k,UTEMP);
                  eos_state.e   = state_arr(i,j,k,UEINT) * rhoInv;
                  for (int n = 0; n < NumSpec; ++n) {
                      eos_state.xn[n] = state_arr(i,j,k,UFS+n) * rhoInv;
                  }
#if NAUX_NET > 0
                  for (int n = 0; n < NumAux; ++n) {
                      eos_state.aux[n] = state_arr(i,j,k,UFX+n) * rhoInv;
                  }
#endif

                  eos(eos_input_re, eos_state);

                  frhoe_arr(i,j,k) = eos_state.T;

                  state_arr(i,j,k,UTEMP) = frhoe_arr(i,j,k);
              }

              Real rho = state_arr(i,j,k,URHO);
              Real temp = state_arr(i,j,k,UTEMP);
              Real Ye;
              if (NumAux > 0) {
                  Ye = state_arr(i,j,k,UFX);
              } else {
                  Ye = 0.e0_rt;
              }

              Real kp, kr;
              bool comp_kp = false;
              bool comp_kr = true;

              opacity(kp, kr, rho, temp, Ye, nu, comp_kp, comp_kr);

              kpr(i,j,k,igroup) = kr;
          });
      }
  }
}

// Only updates interior cells of current level, leaves kappa_r unchanged
// in ghost cells bordering coarse grids or physical boundaries.

void Radiation::update_rosseland_from_temp(MultiFab& kappa_r,
                                           MultiFab& temp,
                                           MultiFab& state,
                                           const Geometry& geom,
                                           int igroup)
{
  BL_PROFILE("update_rosseland_from_temp");

  BL_ASSERT(kappa_r.nGrow() == 1);
  BL_ASSERT(temp.nGrow()    == 0);
  BL_ASSERT(kappa_r.nComp() == Radiation::nGroups);

  const Real nu = nugroup[igroup];

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

      auto state_arr = state[mfi].array();
      auto temp_arr = temp[mfi].array();
      auto kpr = kappa_r[mfi].array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
      {
          state_arr(i,j,k,UTEMP) = temp_arr(i,j,k);

          Real rho = state_arr(i,j,k,URHO);
          Real temp = state_arr(i,j,k,UTEMP);
          Real Ye;
          if (NumAux > 0) {
              Ye = state_arr(i,j,k,UFX);
          } else {
              Ye = 0.e0_rt;
          }

          Real kp, kr;
          bool comp_kp = false;
          bool comp_kr = true;

          opacity(kp, kr, rho, temp, Ye, nu, comp_kp, comp_kr);

          kpr(i,j,k,igroup) = kr;
      });
  }

  kappa_r.FillBoundary(geom.periodicity());
}

void Radiation::SGFLD_compute_rosseland(MultiFab& kappa_r, const MultiFab& state)
{
  BL_PROFILE("Radiation::SGFLD_compute_rosseland (MultiFab)");

  const int igroup = 0;
  const Real nu = nugroup[igroup];

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(kappa_r, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& kbox = mfi.growntilebox();

      auto state_arr = state[mfi].array();
      auto kpr = kappa_r[mfi].array();

      amrex::ParallelFor(kbox,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
      {
          Real rho = state_arr(i,j,k,URHO);
          Real temp = state_arr(i,j,k,UTEMP);
          Real Ye;
          if (NumAux > 0) {
              Ye = state_arr(i,j,k,UFX);
          } else {
              Ye = 0.e0_rt;
          }

          Real kp, kr;
          bool comp_kp = false;
          bool comp_kr = true;

          opacity(kp, kr, rho, temp, Ye, nu, comp_kp, comp_kr);

          kpr(i,j,k,igroup) = kr;
      });
  }
}

void Radiation::SGFLD_compute_rosseland(FArrayBox& kappa_r, const FArrayBox& state)
{
  BL_PROFILE("Radiation::SGFLD_compute_rosseland (FArrayBox)");

  const Box& kbox = kappa_r.box();

  const int igroup = 0;
  const Real nu = nugroup[igroup];

  auto state_arr = state.array();
  auto kpr = kappa_r.array();

  amrex::ParallelFor(kbox,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
  {
      Real rho = state_arr(i,j,k,URHO);
      Real temp = state_arr(i,j,k,UTEMP);
      Real Ye;
      if (NumAux > 0) {
          Ye = state_arr(i,j,k,UFX);
      } else {
          Ye = 0.e0_rt;
      }

      Real kp, kr;
      bool comp_kp = false;
      bool comp_kr = true;

      opacity(kp, kr, rho, temp, Ye, nu, comp_kp, comp_kr);

      kpr(i,j,k,igroup) = kr;
  });
  Gpu::synchronize();
}

void Radiation::deferred_sync_setup(int crse_level)
{
  BL_PROFILE("Radiation::deferred_sync_setup");
  if (verbose && ParallelDescriptor::IOProcessor()) {
    std::cout << "Radiation deferred sync setup for coarse level " << crse_level << "..." << std::endl;
  }

  // sync is only called if this is not the finest level

  int level = crse_level + 1;

  const Orientation ori(0, Orientation::low);

  int do_swap = (flux_cons_old[level] &&
                 (flux_cons_old[level]->coarsenedBoxes() ==
                  flux_cons[level]->coarsenedBoxes()) &&
                 ((*flux_cons_old[level])[ori].DistributionMap() ==
                  (*flux_cons[level])[ori].DistributionMap()));

  if (do_swap) {
      // flux_cons_old exists and is based on same grids as flux_cons
      std::swap(flux_cons_old[level], flux_cons[level]);
  }
  else {
      flux_cons_old[level] = std::move(flux_cons[level]);
      // leave flux_cons undefined because we may be about to regrid anyway
  }

  delta_t_old[crse_level] = parent->dtLevel(crse_level);

  if (verbose && ParallelDescriptor::IOProcessor()) {
    if (do_sync) {
      std::cout << "                                                   done"
           << std::endl;
    }
    else {
      std::cout << "                                                   zeroed out"
           << std::endl;
    }
  }
}

void Radiation::deferred_sync(int level, MultiFab& rhs, int indx)
{
  int fine_level = parent->finestLevel();
  Castro *castro = dynamic_cast<Castro*>(&parent->getLevel(level));
  const BoxArray& grids = castro->boxArray();
  const DistributionMapping& dmap = castro->DistributionMap();
  Real delta_t = parent->dtLevel(level);

  if (level < parent->maxLevel() &&
      flux_cons_old[level+1]) {

    FluxRegister& sync_flux = *flux_cons_old[level+1];

    if (indx == 0) {
      // clean up fine-fine interfaces (does all groups at once)
        sync_flux.ClearInternalBorders(castro->Geom());
    }

    Real scale = delta_t_old[level] / delta_t;
    if (!do_sync) scale = 0.0;

    sync_flux.Reflux(rhs, scale, indx, 0, 1, castro->Geom());

    // If there is a sync source from a still finer level,
    // test for an intersection that will not be caught by
    // an intermediate level and print a message to standard
    // output if one is found.

    if (level+1 < parent->maxLevel() &&
        flux_cons_old[level+2]) {

      BoxArray fgrids;

      if (level+1 <= fine_level) {

        // We know here that level+1 still exists.
        // (If it does not exist, fgrids remains empty. The operations
        // used below are well-defined for an empty BoxArray.)

        fgrids = parent->getLevel(level+1).boxArray();
      }

      const IntVect& rat = sync_flux.refRatio();
      BoxArray old_boxes = sync_flux.coarsenedBoxes();
      old_boxes.refine(rat); // now refined to level+1

      // If level+1 grids have not changed, any higher level
      // flux registers must be properly nested within them
      // and can be ignored.

      if (fgrids != old_boxes) {

        // level+1 grids have changed.  Check whether level+1
        // still contains the region any higher registers
        // will reflux into:

        IntVect ref_rat = IntVect::TheUnitVector();

        for (int flev = level+1; flev < parent->maxLevel(); flev++) {

          // ref_rat is the ratio between flev and level

          ref_rat *= parent->refRatio(flev-1);

          if (flux_cons_old[flev+1]) {

            if (flev > level+1) {
              fgrids.refine(parent->refRatio(flev-1));
              // fgrids is now at refinement level flev
            }

            FluxRegister& ff_sync = *flux_cons_old[flev+1];
            BoxArray ffgr = ff_sync.coarsenedBoxes();
            ffgr.grow(1); // these are the cells to reflux into

            // we don't reflux into cells outside the domain
            const Box& domain = parent->Geom(flev).Domain();

            // this may be inefficient, so do it box by box instead
            //ffgr = amrex::intersect(ffgr, domain);
            for (int i = 0; i < ffgr.size(); i++) {
              ffgr.set(i, (ffgr[i] & domain));
            }

            // For large parallel jobs, the following
            // BoxArray::contains() function may be
            // prohibitively expensive.
            // We could optimize by only doing it once and
            // reusing the result for all the groups.

            if (fgrids.contains(ffgr)) {
              // The higher level flux register is still contained
              // within level+1, so we can ignore it and any
              // finer levels here.

              break;
            }
            else {
              // Either level+1 does not exist, or it does not
              // contain the flux register from level flev+1.

              // I don't think we can use FluxRegister to reflux
              // from flev+1 to level, not even if we constructed
              // an intermediate FluxRegister at some other resolution
              // and initialized it from ff_sync (as we do below in
              // the coarse-to-fine direction).  The problem is that
              // for an extreme difference in levels, the faces in
              // ff_sync may not lie on cell faces in level, and the
              // FluxRegister mechanism assumes that they do (and
              // can't be fooled into doing the right thing anyway).

              // Fortunately this is a rare case---it apparently can
              // only happen when there is a change in the refinement
              // criteria.  What I'm going to do is construct a
              // MultiFab at level flev, reflux ff_sync into it, and
              // then average the result down to level.  This is much
              // simpler than it would be to
              // re-implement reflux from the FluxRegister class for
              // this special circumstance.

              // ClearInternalBorders may not have been done yet.
              // We'll do it again when we advance finer levels up to
              // flev, so there is some duplication of effort, but
              // this is an unusual situation.

              if (indx == 0) {
                  // clean up fine-fine interfaces (does all groups at once)
                  ff_sync.ClearInternalBorders(parent->Geom(flev));
              }

              BoxArray coarsened_grids(ffgr);

              // coarsened_grids is now at flev resolution, and
              // large enough to contain all of the cells that
              // ff_sync would coarsen into.  (It's tempting to just
              // use a refined copy of grids at level, but that
              // could be too large to refine.)

              coarsened_grids.coarsen(ref_rat);

              // The way the coarsen method is defined, if the fine
              // BoxArray does not coarsen evenly by the specified ratio,
              // the coarsened result will be just a little bit larger
              // so that it still fully contains the fine region.

              BoxArray refined_grids(coarsened_grids);
              refined_grids.refine(ref_rat);

              // refined_grids is now back to flev resolution, but is
              // enlarged so that it evenly covers whole cells at the
              // resolution of level.

              const DistributionMapping& dm_ff_sync = ff_sync.DistributionMap();

              MultiFab flev_data(refined_grids, dm_ff_sync, 1, 0);
              flev_data.setVal(0.0);

              // Reflux into this expanded (flev resolution) MultiFab.

              // The units of the data in ff_sync are
              // <conserved quantity> per (flev) timestep.
              // The averaging process will preserve the volume
              // integral of <conserved quantity>.

              // So the remaining scale factor is determined only
              // by the ratio of the (flev) timestep for which
              // this flux data was originally prepared, to the
              // current timestep.

              Real scale = delta_t_old[flev] / delta_t;
              if (!do_sync) scale = 0.0;

              ff_sync.Reflux(flev_data, scale, indx, 0, 1,
                             parent->Geom(flev));

              // Coarsen to level resolution.

              MultiFab rhs_tmp(grids, dmap, 1, 0);
              rhs_tmp.setVal(0.0);     // clear garbage

              // The data we are coarsening already has the metric factor built in.
              amrex::average_down(flev_data, rhs_tmp, 0, 1, ref_rat);

              // Add coarsened result into rhs.
              rhs.plus(rhs_tmp, 0, 1, 0);
            }
          }
        } // end loop over finer levels
      }

    } // end refluxing from levels higher than level+1

  } // end refluxing from all finer levels

  // If there is a sync source from this level or a coarser
  // level that now intersects this level, reflux it too:

  IntVect ref_rat = IntVect::TheUnitVector();

  for (int flev = level; flev > 0; flev--) {

    // delta_t_old is used here as a hack indicating that the
    // data in this flux register is "fresh".  See the end of
    // the update routine where delta_t_old is zeroed out after
    // the data becomes stale.

    if (flux_cons_old[flev] &&
        delta_t_old[flev-1] > 0.0) {

      FluxRegister& crse_sync_flux = *flux_cons_old[flev];

      const IntVect& rat = crse_sync_flux.refRatio();
      BoxArray old_boxes = crse_sync_flux.coarsenedBoxes();
      old_boxes.refine(rat); // now refined to level flev

      const BoxArray& fgrids = parent->getLevel(flev).boxArray();
      if (fgrids == old_boxes) {

        // If the grids at level flev have not changed, then the
        // current level must still be properly nested within it
        // and we don't need to go any further.

        break;
      }
      else {

        // ref_rat is the ratio between level and flev

        old_boxes.refine(ref_rat); // now refined to current level

        // flev grids have changed, so check if the current level
        // is still properly nested within them.  If so, no
        // refluxing is necessary from this or coarser levels.

        // For large parallel jobs, the following
        // BoxArray::contains() function may be
        // prohibitively expensive.
        // We could optimize by only doing it once and
        // reusing the result for all the groups.

        if (old_boxes.contains(grids)) {
          break;
        }

        // We now know we have to reflux from crse_sync_flux

        const DistributionMapping& dm_crse_sync_flux = crse_sync_flux.DistributionMap();

        FluxRegister ref_sync_flux(old_boxes, dm_crse_sync_flux,
                                   IntVect::TheUnitVector(),
                                   level+1, // this is not used
                                   1);
        ref_rat *= rat;

        // ref_rat is now the ratio between level and flev-1

        // Fill in each Fab in ref_sync_flux from
        // the corresponding Fab in crse_sync_flux, refined
        // with piecewise-constant (area-adjusted)
        // interpolation with ratio given by ref_rat.

        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
          const Orientation lo_face(dir,Orientation::low);
          const Orientation hi_face(dir,Orientation::high);

          for (FabSetIter fsi(ref_sync_flux[lo_face]);
               fsi.isValid(); ++fsi) {

            rfface(ref_sync_flux[lo_face][fsi].array(),
                   crse_sync_flux[lo_face][fsi].array(indx),
                   dir, ref_rat);
          }

          for (FabSetIter fsi(ref_sync_flux[hi_face]);
               fsi.isValid(); ++fsi) {

            rfface(ref_sync_flux[hi_face][fsi].array(),
                   crse_sync_flux[hi_face][fsi].array(indx),
                   dir, ref_rat);
          }
        }

        // ClearInternalBorders has already been done---it
        // was done when refluxing this data to level flev-1.

        // The units of the data in crse_sync_flux are
        // <conserved quantity> per (coarse) timestep.  The
        // RFFACE routine divides each value by the number of
        // times it is being duplicated, leaving the same
        // total sum as before.

        // So the remaining scale factor is determined only
        // by the ratio of the (coarse) timestep for which
        // this flux data was originally prepared, to the
        // current timestep.

        Real scale = delta_t_old[flev-1] / delta_t;
        if (!do_sync) scale = 0.0;

        ref_sync_flux.Reflux(rhs, scale, 0, 0, 1, castro->Geom());

      } // end if "this register may overlap current level"

    } // end if "this register exists and has data"

  } // end loop over coarser flux register levels

  // All refluxing into the rhs is now done.
}

void Radiation::reflux(int level)
{
  return;
}

// Computes the scaled gradient for use in flux limiters

void Radiation::scaledGradient(int level,
                               Array<MultiFab, AMREX_SPACEDIM>& R,
                               MultiFab& kappa_r, int kcomp,
                               MultiFab& Er, int igroup,
                               int nGrow_Er, int Rcomp)
{
  BL_PROFILE("Radiation::scaledGradient");
  BL_ASSERT(kappa_r.nGrow() == 1);

  MultiFab Erbtmp;
  if (nGrow_Er == 0) { // default value
    if (radiation::limiter > 0) {
      const BoxArray& grids = parent->boxArray(level);
      const DistributionMapping& dmap = parent->DistributionMap(level);
      Erbtmp.define(grids,dmap,1,1);
      Erbtmp.setVal(-1.0);
      MultiFab::Copy(Erbtmp, Er, igroup, 0, 1, 0);

      // Values in ghost cells are set to -1, indicating that one-sided
      // differences should be used in computing the gradient term for
      // the flux limiter.  In order to make the solution independent
      // of the grid layout, we now go back and overwrite values in
      // those cells bordering grids at the same level:

      Erbtmp.FillBoundary(parent->Geom(level).periodicity());
    }
  }

  MultiFab& Erborder = (nGrow_Er==0) ? Erbtmp : Er;

  auto dx = parent->Geom(level).CellSizeArray();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {

      for (MFIter mfi(R[idim], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

          const Box& nbx = mfi.tilebox();  // note that R is edge based

          auto R_arr = R[idim][mfi].array(Rcomp);

          if (radiation::limiter == 0) {

              R[idim][mfi].setVal<RunOn::Device>(0.0, Rcomp);

          }
          else {

              int include_cross_terms = 0;

              if (radiation::limiter == 1) {
                  include_cross_terms = 0;
              } else if (radiation::limiter == 2) {
                  include_cross_terms = 1;
              } else {
                  amrex::Abort("Unknown limiter");
              }

              auto kap_arr = kappa_r[mfi].array(kcomp);
              auto Er_arr = Erborder[mfi].array();

              amrex::ParallelFor(nbx,
              [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
              {
                  Real dxInv[3] = {0.0};

                  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                      dxInv[d] = 1.e0_rt / dx[d];
                  }

                  Real dal, dar, dbl, dbr;

                  const Real tiny = 1.e-50_rt;

                  if (idim == 0)
                  {
                      if (include_cross_terms == 1)
                      {
#if (AMREX_SPACEDIM >= 2)
                          dal = Er_arr(i-1,j+1,k) - Er_arr(i-1,j-1,k);
                          dar = Er_arr(i  ,j+1,k) - Er_arr(i  ,j-1,k);

                          if      (Er_arr(i-1,j-1,k) == -1.e0_rt)
                          {
                              dal = 2.e0_rt * (Er_arr(i-1,j+1,k) - Er_arr(i-1,j  ,k));
                          }
                          else if (Er_arr(i-1,j+1,k) == -1.e0_rt)
                          {
                              dal = 2.e0_rt * (Er_arr(i-1,j  ,k) - Er_arr(i-1,j-1,k));
                          }

                          if      (Er_arr(i  ,j-1,k) == -1.e0_rt)
                          {
                              dar = 2.e0_rt * (Er_arr(i  ,j+1,k) - Er_arr(i  ,j  ,k));
                          }
                          else if (Er_arr(i  ,j+1,k) == -1.e0_rt)
                          {
                              dar = 2.e0_rt * (Er_arr(i  ,j  ,k) - Er_arr(i  ,j-1,k));
                          }
#else
                          dal = 0.e0_rt;
                          dar = 0.e0_rt;
#endif

#if (AMREX_SPACEDIM == 3)
                          dbl = Er_arr(i-1,j,k+1) - Er_arr(i-1,j,k-1);
                          dbr = Er_arr(i  ,j,k+1) - Er_arr(i  ,j,k-1);

                          if      (Er_arr(i-1,j,k-1) == -1.e0_rt)
                          {
                              dbl = 2.e0_rt * (Er_arr(i-1,j,k+1) - Er_arr(i-1,j,k  ));
                          }
                          else if (Er_arr(i-1,j,k+1) == -1.e0_rt)
                          {
                              dbl = 2.e0_rt * (Er_arr(i-1,j,k  ) - Er_arr(i-1,j,k-1));
                          }

                          if      (Er_arr(i  ,j,k-1) == -1.e0_rt)
                          {
                              dbr = 2.e0_rt * (Er_arr(i  ,j,k+1) - Er_arr(i  ,j,k  ));
                          }
                          else if (Er_arr(i  ,j,k+1) == -1.e0_rt)
                          {
                              dbr = 2.e0_rt * (Er_arr(i  ,j,k  ) - Er_arr(i  ,j,k-1));
                          }
#else
                          dbl = 0.e0_rt;
                          dbr = 0.e0_rt;
#endif

                      }
                      else
                      {
                          dal = 0.e0_rt;
                          dar = 0.e0_rt;
                          dbl = 0.e0_rt;
                          dbr = 0.e0_rt;
                      }

                      Real rg;

                      if (Er_arr(i-1,j,k) == -1.e0_rt)
                      {
                          rg = std::pow((Er_arr(i+1,j,k) - Er_arr(i,j,k)) * dxInv[0], 2) +
                               std::pow(0.5_rt * dar * dxInv[1], 2) +
                               std::pow(0.5_rt * dbr * dxInv[2], 2);
                      }
                      else if (Er_arr(i,j,k) == -1.e0_rt)
                      {
                          rg = std::pow((Er_arr(i-1,j,k) - Er_arr(i-2,j,k)) * dxInv[0], 2) +
                               std::pow(0.5_rt * dal * dxInv[1], 2) +
                               std::pow(0.5_rt * dbl * dxInv[2], 2);
                      }
                      else
                      {
                          rg = std::pow((Er_arr(i,j,k) - Er_arr(i-1,j,k)) * dxInv[0], 2) +
                               std::pow((1.0_rt / 4.0_rt) * (dal + dar) * dxInv[1], 2) +
                               std::pow((1.0_rt / 4.0_rt) * (dbl + dbr) * dxInv[2], 2);
                      }

                      Real kap = kavg(kap_arr(i-1,j,k), kap_arr(i,j,k), dx[0], -1);
                      R_arr(i,j,k) = std::sqrt(rg) / (kap * amrex::max(Er_arr(i-1,j,k), Er_arr(i,j,k), tiny));

                  }
                  else if (idim == 1)
                  {
                      if (include_cross_terms == 1)
                      {
                          dal = Er_arr(i+1,j-1,k  ) - Er_arr(i-1,j-1,k  );
                          dar = Er_arr(i+1,j  ,k  ) - Er_arr(i-1,j  ,k  );

                          if      (Er_arr(i-1,j-1,k  ) == -1.e0_rt)
                          {
                              dal = 2.e0_rt * (Er_arr(i+1,j-1,k  ) - Er_arr(i  ,j-1,k  ));
                          }
                          else if (Er_arr(i+1,j-1,k  ) == -1.e0_rt)
                          {
                              dal = 2.e0_rt * (Er_arr(i  ,j-1,k  ) - Er_arr(i-1,j-1,k  ));
                          }

                          if      (Er_arr(i-1,j  ,k  ) == -1.e0_rt)
                          {
                              dar = 2.e0_rt * (Er_arr(i+1,j  ,k  ) - Er_arr(i  ,j  ,k  ));
                          }
                          else if (Er_arr(i+1,j  ,k  ) == -1.e0_rt)
                          {
                              dar = 2.e0_rt * (Er_arr(i  ,j  ,k  ) - Er_arr(i-1,j  ,k  ));
                          }

#if (AMREX_SPACEDIM == 3)
                          dbl = Er_arr(i  ,j-1,k+1) - Er_arr(i  ,j-1,k-1);
                          dbr = Er_arr(i  ,j  ,k+1) - Er_arr(i  ,j  ,k-1);

                          if      (Er_arr(i  ,j-1,k-1) == -1.e0_rt)
                          {
                              dbl = 2.e0_rt * (Er_arr(i  ,j-1,k+1) - Er_arr(i  ,j-1,k  ));
                          }
                          else if (Er_arr(i  ,j-1,k+1) == -1.e0_rt)
                          {
                              dbl = 2.e0_rt * (Er_arr(i  ,j-1,k  ) - Er_arr(i  ,j-1,k-1));
                          }

                          if      (Er_arr(i  ,j  ,k-1) == -1.e0_rt)
                          {
                              dbr = 2.e0_rt * (Er_arr(i  ,j  ,k+1) - Er_arr(i  ,j,  k  ));
                          }
                          else if (Er_arr(i  ,j  ,k+1) == -1.e0_rt)
                          {
                              dbr = 2.e0_rt * (Er_arr(i  ,j  ,k  ) - Er_arr(i  ,j,  k-1));
                          }
#else
                          dbl = 0.e0_rt;
                          dbr = 0.e0_rt;
#endif
                      }
                      else
                      {
                          dal = 0.e0_rt;
                          dar = 0.e0_rt;
                          dbl = 0.e0_rt;
                          dbr = 0.e0_rt;
                      }

                      Real rg;

                      if (Er_arr(i,j-1,k) == -1.e0_rt)
                      {
                          rg = std::pow((Er_arr(i,j+1,k) - Er_arr(i,j,k)) * dxInv[1], 2) +
                               std::pow(0.5_rt * dar * dxInv[0], 2) +
                               std::pow(0.5_rt * dbr * dxInv[2], 2);
                      }
                      else if (Er_arr(i,j,k) == -1.e0_rt)
                      {
                          rg = std::pow((Er_arr(i,j-1,k) - Er_arr(i,j-2,k)) * dxInv[1], 2) +
                               std::pow(0.5_rt * dal * dxInv[0], 2) +
                               std::pow(0.5_rt * dbl * dxInv[2], 2);
                      }
                      else
                      {
                          rg = std::pow((Er_arr(i,j,k) - Er_arr(i,j-1,k)) * dxInv[1], 2) +
                               std::pow((1.0_rt / 4.0_rt) * (dal + dar) * dxInv[0], 2) +
                               std::pow((1.0_rt / 4.0_rt) * (dbl + dbr) * dxInv[2], 2);
                      }

                      Real kap = kavg(kap_arr(i,j-1,k), kap_arr(i,j,k), dx[1], -1);
                      R_arr(i,j,k) = std::sqrt(rg) / (kap * amrex::max(Er_arr(i,j-1,k), Er_arr(i,j,k), tiny));
                  }
                  else
                  {
                      if (include_cross_terms == 1)
                      {
                          dal = Er_arr(i+1,j  ,k-1) - Er_arr(i-1,j  ,k-1);
                          dar = Er_arr(i+1,j  ,k  ) - Er_arr(i-1,j  ,k  );

                          if      (Er_arr(i-1,j  ,k-1) == -1.e0_rt)
                          {
                              dal = 2.e0_rt * (Er_arr(i+1,j  ,k-1) - Er_arr(i  ,j  ,k-1));
                          }
                          else if (Er_arr(i+1,j  ,k-1) == -1.e0_rt)
                          {
                              dal = 2.e0_rt * (Er_arr(i  ,j  ,k-1) - Er_arr(i-1,j  ,k-1));
                          }

                          if      (Er_arr(i-1,j  ,k  ) == -1.e0_rt)
                          {
                              dar = 2.e0_rt * (Er_arr(i+1,j  ,k  ) - Er_arr(i  ,j  ,k  ));
                          }
                          else if (Er_arr(i+1,j  ,k  ) == -1.e0_rt)
                          {
                              dar = 2.e0_rt * (Er_arr(i  ,j  ,k  ) - Er_arr(i-1,j  ,k  ));
                          }

                          dbl = Er_arr(i  ,j+1,k-1) - Er_arr(i  ,j-1,k-1);
                          dbr = Er_arr(i  ,j+1,k  ) - Er_arr(i  ,j-1,k  );

                          if      (Er_arr(i  ,j-1,k-1) == -1.e0_rt)
                          {
                              dbl = 2.e0_rt * (Er_arr(i  ,j+1,k-1) - Er_arr(i  ,j  ,k-1));
                          }
                          else if (Er_arr(i  ,j+1,k-1) == -1.e0_rt)
                          {
                              dbl = 2.e0_rt * (Er_arr(i  ,j  ,k-1) - Er_arr(i  ,j-1,k-1));
                          }

                          if      (Er_arr(i  ,j-1,k  ) == -1.e0_rt)
                          {
                              dbr = 2.e0_rt * (Er_arr(i  ,j+1,k  ) - Er_arr(i  ,j,  k  ));
                          }
                          else if (Er_arr(i  ,j+1,k  ) == -1.e0_rt)
                          {
                              dbr = 2.e0_rt * (Er_arr(i  ,j  ,k  ) - Er_arr(i  ,j-1,k  ));
                          }
                      }
                      else
                      {
                          dal = 0.e0_rt;
                          dar = 0.e0_rt;
                          dbl = 0.e0_rt;
                          dbr = 0.e0_rt;
                      }

                      Real rg;

                      if (Er_arr(i,j,k-1) == -1.e0_rt)
                      {
                          rg = std::pow((Er_arr(i,j,k+1) - Er_arr(i,j,k)) * dxInv[2], 2) +
                               std::pow(0.5_rt * dar * dxInv[0], 2) +
                               std::pow(0.5_rt * dbr * dxInv[1], 2);
                      }
                      else if (Er_arr(i,j,k) == -1.e0_rt)
                      {
                          rg = std::pow((Er_arr(i,j,k-1) - Er_arr(i,j,k-2)) * dxInv[2], 2) +
                               std::pow(0.5_rt * dal * dxInv[0], 2) +
                               std::pow(0.5_rt * dbl * dxInv[1], 2);
                      }
                      else
                      {
                          rg = std::pow((Er_arr(i,j,k) - Er_arr(i,j,k-1)) * dxInv[2], 2) +
                               std::pow((1.0_rt / 4.0_rt) * (dal + dar) * dxInv[0], 2) +
                               std::pow((1.0_rt / 4.0_rt) * (dbl + dbr) * dxInv[1], 2);
                      }

                      Real kap = kavg(kap_arr(i,j,k-1), kap_arr(i,j,k), dx[2], -1);
                      R_arr(i,j,k) = std::sqrt(rg) / (kap * amrex::max(Er_arr(i,j,k-1), Er_arr(i,j,k), tiny));
                  }
              });
          }
      }
  }
}

// On input, lambda should contain scaled gradient.
// On output this will be overwritten with the flux limiter.

void Radiation::fluxLimiter(int level,
                            Array<MultiFab, AMREX_SPACEDIM>& lambda,
                            int lamcomp)
{
    BL_PROFILE("Radiation:fluxLimiter");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        for (MFIter mfi(lambda[idim], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();

            auto lambda_arr = lambda[idim][mfi].array(lamcomp);

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
            {
                lambda_arr(i,j,k) = FLDlambda(lambda_arr(i,j,k));
            });
        }
    }
}

void Radiation::get_rosseland_v_dcf(MultiFab& kappa_r, MultiFab& v, MultiFab& dcf,
                                    Real delta_t, Real c,
                                    AmrLevel* castro, int igroup)
{
    BL_ASSERT(kappa_r.nGrow() == 1);
    BL_ASSERT(      v.nGrow() == 1);
    BL_ASSERT(    dcf.nGrow() == 1);

    int nstate = castro->get_new_data(State_Type).nComp();
    Real time = castro->get_state_data(State_Type).curTime();

    FillPatchIterator fpi_r(*castro, kappa_r, 1, time, Rad_Type, 0, 1);
    MultiFab& Er = fpi_r.get_mf();

    FillPatchIterator fpi_s(*castro, kappa_r, 1, time, State_Type, 0, nstate);
    MultiFab& S = fpi_s.get_mf();

    Real sigma_loc = sigma;
    Real c_loc = c;
    Real dT_loc = dT;

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox temp, c_v, kp, kp2;
        for (MFIter mfi(kappa_r,true); mfi.isValid(); ++mfi)
        {
            const Box& reg = mfi.growntilebox();

            temp.resize(reg);
            Elixir temp_elix = temp.elixir();

            kp.resize(reg);
            Elixir kp_elix = kp.elixir();

            kp2.resize(reg);
            Elixir kp2_elix = kp2.elixir();

            c_v.resize(reg);
            Elixir c_v_elix = c_v.elixir();

            auto temp_arr = temp.array();
            auto S_arr = S[mfi].array();
            auto v_arr = v[mfi].array();
            auto er_arr = Er[mfi].array();
            auto kr_arr = kappa_r[mfi].array();
            auto kp_arr = kp.array();
            auto kp2_arr = kp2.array();
            auto c_v_arr = c_v.array();
            auto dcf_arr = dcf[mfi].array();

            amrex::ParallelFor(reg,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
            {
                // Get T from rhoe

                if (S_arr(i,j,k,UEINT) <= 0.e0_rt)
                {
                    temp_arr(i,j,k) = small_temp;
                }
                else
                {
                    Real rhoInv = 1.e0_rt / S_arr(i,j,k,URHO);

                    eos_re_t eos_state;
                    eos_state.rho = S_arr(i,j,k,URHO);
                    eos_state.T   = S_arr(i,j,k,UTEMP);
                    eos_state.e   = S_arr(i,j,k,UEINT) * rhoInv;
                    for (int n = 0; n < NumSpec; ++n) {
                        eos_state.xn[n] = S_arr(i,j,k,UFS+n) * rhoInv;
                    }
#if NAUX_NET > 0
                    for (int n = 0; n < NumAux; ++n) {
                        eos_state.aux[n] = S_arr(i,j,k,UFX+n) * rhoInv;
                    }
#endif

                    eos(eos_input_re, eos_state);

                    temp_arr(i,j,k) = eos_state.T;
                }
            });
            Gpu::synchronize();

            get_c_v(c_v, temp, S[mfi], reg);

            S[mfi].copy<RunOn::Device>(temp,reg,0,reg,UTEMP,1);

            const Real nu = nugroup[igroup];

            amrex::ParallelFor(reg,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
            {
                Real rho = S_arr(i,j,k,URHO);
                Real temp = S_arr(i,j,k,UTEMP);
                Real Ye;
                if (NumAux > 0)  {
                    Ye = S_arr(i,j,k,UFX);
                } else {
                    Ye = 0.e0_rt;
                }

                Real kp, kr;
                bool comp_kp = true;
                bool comp_kr = true;
                opacity(kp, kr, rho, temp, Ye, nu, comp_kp, comp_kr);

                kp_arr(i,j,k) = kp;
                kr_arr(i,j,k) = kr;

                temp += dT_loc;

                opacity(kp, kr, rho, temp, Ye, nu, comp_kp, comp_kr);

                kp2_arr(i,j,k) = kp;
            });

            amrex::ParallelFor(reg,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
            {
                const Real fac0 = 4.e0_rt * sigma_loc * delta_t / dT_loc;
                const Real fac2 = c_loc * delta_t / dT_loc;

                v_arr(i,j,k,0) = S_arr(i,j,k,UMX) / S_arr(i,j,k,URHO);
#if AMREX_SPACEDIM >= 2
                v_arr(i,j,k,1) = S_arr(i,j,k,UMY) / S_arr(i,j,k,URHO);
#endif
#if AMREX_SPACEDIM == 3
                v_arr(i,j,k,2) = S_arr(i,j,k,UMZ) / S_arr(i,j,k,URHO);
#endif

                Real alpha = fac0 * (kp2_arr(i,j,k) * std::pow(temp_arr(i,j,k) + dT_loc, 4) -
                                     kp_arr(i,j,k) * std::pow(temp_arr(i,j,k), 4)) -
                             fac2 * (kp2_arr(i,j,k) - kp_arr(i,j,k)) * er_arr(i,j,k);

                Real frc = S_arr(i,j,k,URHO) * c_v_arr(i,j,k) + 1.0e-50_rt;
                Real etainv = frc / (alpha + frc);

                dcf_arr(i,j,k) = 2.e0_rt * etainv * (kp_arr(i,j,k) / kr_arr(i,j,k));
            });
        }
    }
}

void Radiation::update_dcf(MultiFab& dcf, MultiFab& etainv, MultiFab& kp, MultiFab& kr,
                           const Geometry& geom)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(dcf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        auto dcf_arr = dcf[mfi].array();
        auto etainv_arr = etainv[mfi].array();
        auto kp_arr = kp[mfi].array();
        auto kr_arr = kr[mfi].array();

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {
            dcf_arr(i,j,k) = 2.e0_rt * etainv_arr(i,j,k) * (kp_arr(i,j,k) / kr_arr(i,j,k));
        });
    }

    dcf.FillBoundary(geom.periodicity());
}

void Radiation::set_current_group(int igroup)
{
  current_group_number = igroup;
}


void Radiation::filter_prim(int level, MultiFab& State)
{
  Castro *castro = dynamic_cast<Castro*>(&parent->getLevel(level));
  const Geometry& geom = parent->Geom(level);
  auto geomdata = geom.data();

  int ngrow = filter_prim_T;
  int ncomp = State.nComp();
  Real time = castro->get_state_data(Rad_Type).curTime();

  FillPatchIterator fpi(*castro,State,ngrow,time,State_Type,0,ncomp);
  MultiFab& S_fp = fpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(State,true); mfi.isValid(); ++mfi)
  {
      const Box& bx = mfi.tilebox();

      auto S_fp_arr = S_fp[mfi].array();
      auto State_arr = State[mfi].array();

      int T = filter_prim_T;
      int S = filter_prim_S;

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
      {
          filt_prim(i, j, k,
                    S_fp_arr, State_arr,
                    T, S,
                    geomdata, time);
      });
  }
}
