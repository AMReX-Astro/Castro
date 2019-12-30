
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_ParmParse.H>
#include "Radiation.H"
#include "RadSolve.H"

#include "Castro_F.H"

#include "RadTests.H"

#include "RAD_F.H"
#include "AMReX_PROB_AMR_F.H"

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <sstream>

using namespace amrex;

// Radiation test problems and static parameters.  Some of these are
// initialized with inputs values in the function read_static_params.
// This is called from Castro::read_params in Castro.cpp.  The reason
// for initializing these early is that some of them are used in
// Castro::variableSetUp and affect the construction of the state.

int RadTests::do_thermal_wave_cgs = 0;
int RadTests::do_rad_sphere = 0;

Radiation::Solver_Type Radiation::SolverType = Radiation::InvalidSolver;
Radiation::Radiation_Type Radiation::radiation_type = Radiation::Unknown;

Real Radiation::radtoE = 0.;
//Real Radiation::radtoJ = 0.;
Real Radiation::Etorad = 0.;
Real Radiation::radfluxtoF = 0.;

int Radiation::do_multigroup = 0;
int Radiation::nGroups = NGROUPS;
#ifdef NEUTRINO
int Radiation::nNeutrinoSpecies = N_NEUTRINO_SPECIES;
#else
int Radiation::nNeutrinoSpecies = 0;
#endif
Vector<int> Radiation::nNeutrinoGroups(0);
int Radiation::plot_neutrino_group_energies_per_MeV = 1;
int Radiation::plot_neutrino_group_energies_total   = 0;
int Radiation::accelerate = 1;
int Radiation::rad_hydro_combined = 0;
int Radiation::comoving = 1;
int Radiation::Er_Lorentz_term = 1;
int Radiation::fspace_advection_type = 2;
int Radiation::do_inelastic_scattering = 0;
int Radiation::use_analytic_solution = 0;
int Radiation::plot_lambda   = 0;
int Radiation::plot_kappa_p  = 0;
int Radiation::plot_kappa_r  = 0;
int Radiation::plot_lab_Er   = 0;
int Radiation::plot_lab_flux = 0;
int Radiation::plot_com_flux = 0;
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

Real Radiation:: flatten_pp_threshold = -1.0;
int Radiation::pure_hydro = 0;
int Radiation ::do_real_eos = 1;

Real Radiation::const_c_v = -1.0;
Real Radiation::const_kappa_p = -1.0;
Real Radiation::const_kappa_r = -1.0;
Real Radiation::const_scattering = 0.0;

Real Radiation::c_v_exp_m = 0.0;
Real Radiation::c_v_exp_n = 0.0;

Real Radiation::kappa_p_exp_m = 0.0;
Real Radiation::kappa_p_exp_n = 0.0;
Real Radiation::kappa_p_exp_p = 0.0;

Real Radiation::kappa_r_exp_m = 0.0;
Real Radiation::kappa_r_exp_n = 0.0;
Real Radiation::kappa_r_exp_p = 0.0;

Real Radiation::scattering_exp_m = 0.0;
Real Radiation::scattering_exp_n = 0.0;
Real Radiation::scattering_exp_p = 0.0;

Real Radiation::prop_temp_floor = 0.0;


// static initialization, must be called before Castro::variableSetUp

void Radiation::read_static_params()
{
  ParmParse pp("radiation");

  pp.query("do_thermal_wave_cgs",  RadTests::do_thermal_wave_cgs);
  pp.query("do_rad_sphere",    RadTests::do_rad_sphere);

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

  if (SolverType == MGFLDSolver) {
    comoving = 1;
    fspace_advection_type = 2;
    pp.query("fspace_advection_type", fspace_advection_type);
  }
  else if (SolverType == SGFLDSolver) {
    comoving = 1;
    fspace_advection_type = 1;
    pp.query("comoving", comoving);
    if (comoving) {
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

  if (rad_hydro_combined) {
    if (Castro::riemann_solver == 1) {
      amrex::Error("The Colella and Glaz Riemann solver cannot be used with rad_hydro_combined.");
    }
  }

  if (Radiation::SolverType == Radiation::MGFLDSolver) {

#ifdef NEUTRINO
    radiation_type = Neutrino;
    std::string rad_type_s("Unknown");
    pp.query("radiation_type", rad_type_s);
    if (rad_type_s.compare("PHOTON") == 0 ||
	rad_type_s.compare("Photon") == 0 ||
	rad_type_s.compare("photon") == 0 ) {
      radiation_type = Photon;
    }

    if (radiation_type == Neutrino) {

      Radiation::nNeutrinoGroups.resize(Radiation::nNeutrinoSpecies);

      Radiation::nGroups = 0;
#if (N_NEUTRINO_SPECIES >= 1)
      Radiation::nGroups += N_NEUTRINO_GROUPS_1;
      Radiation::nNeutrinoGroups[0] = N_NEUTRINO_GROUPS_1;
#endif
#if (N_NEUTRINO_SPECIES >= 2)
      Radiation::nGroups += N_NEUTRINO_GROUPS_2;
      Radiation::nNeutrinoGroups[1] = N_NEUTRINO_GROUPS_2;
#endif
#if (N_NEUTRINO_SPECIES == 3)
      Radiation::nGroups += N_NEUTRINO_GROUPS_3;
      Radiation::nNeutrinoGroups[2] = N_NEUTRINO_GROUPS_3;
#endif

      pp.query("plot_neutrino_group_energies_per_MeV",
	       Radiation::plot_neutrino_group_energies_per_MeV);
      pp.query("plot_neutrino_group_energies_total",
	       Radiation::plot_neutrino_group_energies_total);
    }
    else {

      // if we are here, we are wanting to do photon radiation with
      // Neutrino compiled -- we no longer support this
      amrex::Error("We do not support photon radiation when compiled for neutrinos");
    }
#else
    // photon radiation
    radiation_type = Photon;
    Radiation::nGroups = NGROUPS;

    // sanity check -- in the old style, we allowed the number of groups to be
    // set at compile time, so ensure that, if the user set this, it matches what
    // we expect
    int test_groups = -1;
    pp.query("nGroups", test_groups);

    if (test_groups > 0 && test_groups != Radiation::nGroups) {
      amrex::Error("you set the number of groups at runtime, but this does not match the compiled value");
    }
#endif
  }
  else if (Radiation::SolverType != Radiation::SingleGroupSolver &&
	   Radiation::SolverType != Radiation::SGFLDSolver) {
      amrex::Error("Unknown Radiation::SolverType");
  }


#ifndef NEUTRINO
  if (Radiation::nGroups > 1)
      pp.query("do_inelastic_scattering", do_inelastic_scattering);
#endif

  pp.query("use_analytic_solution", use_analytic_solution);

  pp.query("plot_lambda", plot_lambda);
  pp.query("plot_kappa_p", plot_kappa_p);
  pp.query("plot_kappa_r", plot_kappa_r);
  if (comoving) pp.query("plot_lab_Er", plot_lab_Er);
  pp.query("plot_lab_flux", plot_lab_flux);
  pp.query("plot_com_flux", plot_com_flux);

  // set up the extra plot variables
  {
      if (use_analytic_solution) {
	  // This is for the special Thermal Wave test
	  plotvar_names.push_back("Texact");
	  plotvar_names.push_back("Terror");
      }
      if (plot_lambda) {
	  icomp_lambda = plotvar_names.size();

	  int limiter = 2;
	  pp.query("limiter", limiter);

	  if (!do_multigroup || limiter == 0) {
	      plotvar_names.push_back("lambda");
	  } else if (nNeutrinoSpecies == 0) {
	      for (int g=0; g<nGroups; ++g) {
		  std::ostringstream ss;
		  ss << "lambda" << g;
		  plotvar_names.push_back(ss.str());
	      }
	  } else {
	      for (int s = 0; s < nNeutrinoSpecies; ++s) {
		  for (int g=0; g<Radiation::nNeutrinoGroups[s]; ++g) {
		      std::ostringstream ss;
		      ss << "lambdas" << s << "g" << g;
		      plotvar_names.push_back(ss.str());
		  }
	      }
	  }
      }
      if (plot_kappa_p) {
	  icomp_kp = plotvar_names.size();
	  if (!do_multigroup) {
	      plotvar_names.push_back("kappa_P");
	  } else if (nNeutrinoSpecies == 0) {
	      for (int g=0; g<nGroups; ++g) {
		  std::ostringstream ss;
		  ss << "kappa_P" << g;
		  plotvar_names.push_back(ss.str());
	      }
	  } else {
	      for (int s = 0; s < nNeutrinoSpecies; ++s) {
		  for (int g=0; g<Radiation::nNeutrinoGroups[s]; ++g) {
		      std::ostringstream ss;
		      ss << "kappa_Ps" << s << "g" << g;
		      plotvar_names.push_back(ss.str());
		  }
	      }
	  }
      }
      if (plot_kappa_r) {
	  icomp_kr = plotvar_names.size();
	  if (!do_multigroup) {
	      plotvar_names.push_back("kappa_R");
	  } else if (nNeutrinoSpecies == 0) {
	      for (int g=0; g<nGroups; ++g) {
		  std::ostringstream ss;
		  ss << "kappa_R" << g;
		  plotvar_names.push_back(ss.str());
	      }
	  } else {
	      for (int s = 0; s < nNeutrinoSpecies; ++s) {
		  for (int g=0; g<Radiation::nNeutrinoGroups[s]; ++g) {
		      std::ostringstream ss;
		      ss << "kappa_Rs" << s << "g" << g;
		      plotvar_names.push_back(ss.str());
		  }
	      }
	  }
      }
      if (plot_lab_Er) {
	  icomp_lab_Er = plotvar_names.size();
	  if (!do_multigroup) {
	      plotvar_names.push_back("Erlab");
	  } else if (nNeutrinoSpecies == 0) {
	      for (int g=0; g<nGroups; ++g) {
		  std::ostringstream ss;
		  ss << "Erlab" << g;
		  plotvar_names.push_back(ss.str());
	      }
	  } else {
	      for (int s = 0; s < nNeutrinoSpecies; ++s) {
		  for (int g=0; g<nNeutrinoGroups[s]; ++g) {
		      std::ostringstream ss;
		      ss << "Erlab" << "s" << s << "g" << g;
		      plotvar_names.push_back(ss.str());
		  }
	      }
	  }
      }
      if (plot_lab_flux) {
	  icomp_lab_Fr = plotvar_names.size();
	  std::string frame = "lab";
	  Vector<std::string> dimname;
	  dimname.push_back("x");
	  dimname.push_back("y");
	  dimname.push_back("z");
	  if (!do_multigroup) {
	      for (int idim=0; idim<BL_SPACEDIM; ++idim) {
		  std::ostringstream ss;
		  ss << "Fr" << frame << dimname[idim];
		  plotvar_names.push_back(ss.str());
	      }
	  } else if (nNeutrinoSpecies == 0) {
	      for (int idim=0; idim<BL_SPACEDIM; ++idim) {
		  for (int g=0; g<nGroups; ++g) {
		      std::ostringstream ss;
		      ss << "Fr" << frame << g << dimname[idim];
		      plotvar_names.push_back(ss.str());
		  }
	      }
	  } else {
	      for (int idim=0; idim<BL_SPACEDIM; ++idim) {
		  for (int s = 0; s < nNeutrinoSpecies; ++s) {
		      for (int g=0; g<nNeutrinoGroups[s]; ++g) {
			  std::ostringstream ss;
			  ss << "Fr" << frame << "s" << s << "g" << g << dimname[idim];
			  plotvar_names.push_back(ss.str());
		      }
		  }
	      }
	  }
      }
      if (plot_com_flux) {
	  icomp_com_Fr = plotvar_names.size();
	  std::string frame = "com";
	  Vector<std::string> dimname;
	  dimname.push_back("x");
	  dimname.push_back("y");
	  dimname.push_back("z");
	  if (!do_multigroup) {
	      for (int idim=0; idim<BL_SPACEDIM; ++idim) {
		  std::ostringstream ss;
		  ss << "Fr" << frame << dimname[idim];
		  plotvar_names.push_back(ss.str());
	      }
	  } else if (nNeutrinoSpecies == 0) {
	      for (int idim=0; idim<BL_SPACEDIM; ++idim) {
		  for (int g=0; g<nGroups; ++g) {
		      std::ostringstream ss;
		      ss << "Fr" << frame << g << dimname[idim];
		      plotvar_names.push_back(ss.str());
		  }
	      }
	  } else {
	      for (int idim=0; idim<BL_SPACEDIM; ++idim) {
		  for (int s = 0; s < nNeutrinoSpecies; ++s) {
		      for (int g=0; g<nNeutrinoGroups[s]; ++g) {
			  std::ostringstream ss;
			  ss << "Fr" << frame << "s" << s << "g" << g << dimname[idim];
			  plotvar_names.push_back(ss.str());
		      }
		  }
	      }
	  }
      }
      nplotvar = plotvar_names.size();
  }

  std::cout << "radiation initialized, nGroups = " << Radiation::nGroups << std::endl;

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
    Real stefbol;

    int J_is_used = 0;

    ca_initradconstants(M_PI, clight, hPlanck, kBoltz, stefbol,
			Avogadro, convert_MeV_erg, J_is_used);

    aRad = 4.*stefbol/clight;

    c        = clight;
    sigma    = stefbol;

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
  limiter  = 2;              pp.query("limiter", limiter);
  if (SolverType == SGFLDSolver && limiter == 1) {
    amrex::Abort("SGFLDSolver does not support limiter = 1");
  }
  if (SolverType == MGFLDSolver && limiter == 1) {
    amrex::Abort("MGFLDSolver does not support limiter = 1");
  }

  closure = 3;
  pp.query("closure", closure);

  ca_initfluxlimiter(&limiter, &closure);

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
  surface_average = 2;       pp.query("surface_average", surface_average);

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
  dedY_fac = 1.0;
  pp.query("dedT_fac", dedT_fac);
  pp.query("dedY_fac", dedY_fac);

  inner_convergence_check = 2;
  pp.query("inner_convergence_check", inner_convergence_check);

  delta_e_rat_dt_tol = 100.0;
  pp.query("delta_e_rat_dt_tol", delta_e_rat_dt_tol);
  delta_T_rat_dt_tol = 100.0;
  pp.query("delta_T_rat_dt_tol", delta_T_rat_dt_tol);
  delta_Ye_dt_tol = 100.0;
  pp.query("delta_Ye_dt_tol", delta_Ye_dt_tol);

  underfac = 1.0;    pp.query("underfac", underfac);

  integrate_Planck = 1;
  pp.query("integrate_Planck", integrate_Planck);

  use_WiensLaw = 0;
  pp.query("use_WiensLaw", use_WiensLaw);
  Tf_Wien = -1.0;
  pp.query("Tf_Wien", Tf_Wien);

  verbose  = 0;      pp.query("v", verbose);  pp.query("verbose", verbose);


  use_opacity_table_module = 0;
  pp.query("use_opacity_table_module", use_opacity_table_module);

  do_kappa_stm_emission = 0;
  pp.query("do_kappa_stm_emission", do_kappa_stm_emission);

  pp.query("const_c_v",       const_c_v);
  pp.query("const_kappa_p",   const_kappa_p);
  pp.query("const_kappa_r",   const_kappa_r);
  pp.query("const_scattering", const_scattering);
  pp.query("c_v_exp_m",       c_v_exp_m);
  pp.query("c_v_exp_n",       c_v_exp_n);
  pp.query("kappa_p_exp_m",   kappa_p_exp_m);
  pp.query("kappa_p_exp_n",   kappa_p_exp_n);
  pp.query("kappa_r_exp_m",   kappa_r_exp_m);
  pp.query("kappa_r_exp_n",   kappa_r_exp_n);
  pp.query("scattering_exp_m", scattering_exp_m);
  pp.query("scattering_exp_n", scattering_exp_n);
  pp.query("prop_temp_floor", prop_temp_floor);
  if (nGroups > 1) {
    pp.query("kappa_p_exp_p",   kappa_p_exp_p);
    pp.query("kappa_r_exp_p",   kappa_r_exp_p);
    pp.query("scattering_exp_p", scattering_exp_p);
  }

  use_dkdT = (kappa_p_exp_n == 0.0) ? 0 : 1;
  pp.query("use_dkdT", use_dkdT);

  if (Radiation::nNeutrinoSpecies == 0 ||
      Radiation::nNeutrinoGroups[0] == 0) {
    // photon problem
    if (SolverType == SGFLDSolver || SolverType == MGFLDSolver) {
      kappa_r_floor = 0.0;
    }
    else {
      kappa_r_floor = 1.e-4;
    }
  }
  else {
    // neutrino problem
    kappa_r_floor = 0.0;
  }
  pp.query("kappa_r_floor", kappa_r_floor);

  temp_floor    = 1.e-10;  pp.query("temp_floor", temp_floor);

  if (verbose >= 1 && ParallelDescriptor::IOProcessor()) {
    std::cout << "const_c_v     = " << const_c_v << std::endl;
    std::cout << "const_kappa_p = " << const_kappa_p << std::endl;
    std::cout << "const_kappa_r = " << const_kappa_r << std::endl;
    std::cout << "const_scattering = " << const_scattering << std::endl;
    std::cout << "c_v_exp_m     = " << c_v_exp_m << std::endl;
    std::cout << "c_v_exp_n     = " << c_v_exp_n << std::endl;
    std::cout << "kappa_p_exp_m = " << kappa_p_exp_m << std::endl;
    std::cout << "kappa_p_exp_n = " << kappa_p_exp_n << std::endl;
    std::cout << "kappa_p_exp_p = " << kappa_p_exp_p << std::endl;
    std::cout << "kappa_r_exp_m = " << kappa_r_exp_m << std::endl;
    std::cout << "kappa_r_exp_n = " << kappa_r_exp_n << std::endl;
    std::cout << "kappa_r_exp_p = " << kappa_r_exp_p << std::endl;
    std::cout << "scattering_exp_m =" << scattering_exp_m << std::endl;
    std::cout << "scattering_exp_n =" << scattering_exp_n << std::endl;
    std::cout << "scattering_exp_p =" << scattering_exp_p << std::endl;
    std::cout << "prop_temp_floor =" << prop_temp_floor << std::endl;
    std::cout << "kappa_r_floor = " << kappa_r_floor << std::endl;
    std::cout << "temp_floor    = " << temp_floor << std::endl;
  }

  if (SolverType == MGFLDSolver ) {
    do_real_eos = 1;
  }
  else {
    do_real_eos = (const_c_v < 0.) ? 1 : 0;
  }
  pp.query("do_real_eos", do_real_eos);

  if (! do_real_eos) {
    BL_ASSERT(const_c_v > 0.);
  }

  if (verbose > 2) {
    Vector<int> temp;
    if (pp.queryarr("spot",temp,0,BL_SPACEDIM)) {
      IntVect tempi(temp);
      spot = tempi;
    }
    if (ParallelDescriptor::IOProcessor()) std::cout << "Spot: " << spot << std::endl;
  }

  // This call stores the value of surface_average in the kavg
  // routine.  The first three arguments are irrelevant here.
  Real foo=0.0;
  FORT_KAVG(foo, foo, foo, surface_average);

  if (verbose > 0 && ParallelDescriptor::IOProcessor()) {
    std::cout << "Creating Radiation object" << std::endl;
  }
  if (verbose >= 1 && ParallelDescriptor::IOProcessor()) {
    std::cout << "processors = " << ParallelDescriptor::NProcs() << std::endl;
    std::cout << "do_sync            = " << do_sync << std::endl;
    std::cout << "use_analytic_solution = " << use_analytic_solution << std::endl;

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
    std::cout << "delta_Ye_dt_tol    = " << delta_Ye_dt_tol    << std::endl;
    std::cout << "limiter  = " << limiter << std::endl;
    std::cout << "closure  = " << closure << std::endl;
    std::cout << "update_limiter   = " << update_limiter << std::endl;
    std::cout << "update_planck    = " << update_planck << std::endl;
    std::cout << "update_rosseland = " << update_rosseland << std::endl;
    std::cout << "delta_temp = " << dT << std::endl;
    std::cout << "surface_average = " << surface_average << std::endl;
    std::cout << "underfac = " << underfac << std::endl;
    std::cout << "do_real_eos = " << do_real_eos << std::endl;
    std::cout << "do_multigroup = " << do_multigroup << std::endl;
    std::cout << "accelerate = " << accelerate << std::endl;
    std::cout << "verbose  = " << verbose << std::endl;
    if (RadTests::do_thermal_wave_cgs)
      std::cout << "do_thermal_wave_cgs = " << RadTests::do_thermal_wave_cgs << std::endl;
    if (RadTests::do_rad_sphere)
      std::cout << "do_rad_sphere = " << RadTests::do_rad_sphere << std::endl;
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
      std::cout << "comoving = " << comoving << std::endl;
    }
    if (SolverType == MGFLDSolver) {
      std::cout << "fspace_advection_type = " << fspace_advection_type << std::endl;
      std::cout << "do_inelastic_scattering = " << do_inelastic_scattering << std::endl;
    }
    if (SolverType == SGFLDSolver && comoving == 0) {
      std::cout << "Er_Lorentz_term = " << Er_Lorentz_term << std::endl;
    }
  }

  if (do_multigroup) {

    get_groups(verbose);

#ifdef NEUTRINO
    if (SolverType == MGFLDSolver && radiation_type == Neutrino) {
      // load opacities from (Burrows) table
      int iverb = (verbose >= 1 && ParallelDescriptor::IOProcessor());
      //      int iverb = 0;
      if (iverb) {
        std::cout << "reading opacity tables..." << std::endl;
      }
      FORT_INIT_OPACITY_TABLE(iverb);
    }
#endif
  }
  else {
    ca_initsinglegroup(nGroups);

    // xnu is a dummy for single group
    xnu.resize(2, 1.0);
    nugroup.resize(1, 1.0);
  }

  // current implementation of the Radiation boundary condition reads
  // incoming flux information in the RadBndry constructor.  we just
  // set the boundary condition type here:

  Vector<int> lo_bc(BL_SPACEDIM), hi_bc(BL_SPACEDIM);
  pp.getarr("lo_bc",lo_bc,0,BL_SPACEDIM);
  pp.getarr("hi_bc",hi_bc,0,BL_SPACEDIM);
  for (int i = 0; i < BL_SPACEDIM; i++) {
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
  delta_Ye_level.resize(   levels, 0.0);

  Density   = castro->Density;
  Xmom      = castro->Xmom;
  Eden      = castro->Eden;
  Eint      = castro->Eint;
  Temp      = castro->Temp;
  FirstSpec = castro->FirstSpec;
  FirstAux  = castro->FirstAux;
  NUM_STATE = castro->NUM_STATE;

  pp.query("flatten_pp_threshold", flatten_pp_threshold);
  pp.query("pure_hydro", pure_hydro);

  if (pure_hydro || limiter == 0) {
      if (verbose > 1 && ParallelDescriptor::IOProcessor()) {
	  std::cout << "turning off inelastic scattering when (pure_hydro || limiter == 0)" << std::endl;
      }
      do_inelastic_scattering = 0;
  }

  ca_init_radhydro_pars(fspace_advection_type, do_inelastic_scattering,
			comoving, flatten_pp_threshold);
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
    else if (aString.find("delta_Ye") == 0) {
      is >> delta_Ye_level[level];
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
    sprintf(buf, "delta_e_rat_level[%d]= ", level);
    std::string DeltaString = buf;
    os << DeltaString << delta_e_rat_level[level] << '\n';
    sprintf(buf, "delta_T_rat_level[%d]= ", level);
    DeltaString = buf;
    os << DeltaString << delta_T_rat_level[level] << '\n';
    sprintf(buf, "delta_Ye_level[%d]= ", level);
    DeltaString = buf;
    os << DeltaString << delta_Ye_level[level] << '\n';
    os.precision(oldprec);
  }

  // Path name construction stolen from AmrLevel::checkPoint

  sprintf(buf, "Level_%d", level);
  std::string Level = buf;

  //
  // Write name of conservation flux register to header.
  //
  sprintf(buf, "/RadFlux");

  std::string PathNameInHeader = Level;
  PathNameInHeader += buf;
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
    FullPathName += buf;
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


void Radiation::analytic_solution(int level)
{
  if (!use_analytic_solution)
    return;

  Castro *castro        = dynamic_cast<Castro*>(&parent->getLevel(level));
  const BoxArray& grids = castro->boxArray();
  const DistributionMapping& dmap = castro->DistributionMap();

  Real time = castro->get_state_data(Rad_Type).curTime();

  if (RadTests::do_thermal_wave_cgs) {

    ParmParse pp("radiation");

    Real rhocv, Q;
    pp.get("thermal_wave_rhocv", rhocv);
    pp.get("thermal_wave_Eexp", Q);
    Q /= rhocv;

    Real a = (16.0 * sigma) / (3.0 * const_kappa_r) / rhocv;
    Real p = kappa_r_exp_n + 3.0;

    Real pfac = exp(lgamma(2.5 + 1.0 / p) - lgamma(1.0 + 1.0 / p) - lgamma(1.5));

    Real xi0  = pow((3.0 * p + 2.0) / (pow(2.0, p - 1.0) * p *
                                      pow(M_PI, p)),
                   1.0 / (3.0 * p + 2.0)) *
           pow(pfac, p / (3.0 * p + 2.0));
    Real xf = xi0 * pow(a * pow(Q,p) * time, 1.0 / (3.0 * p + 2.0));

    Real Tbar = Q / pow(xf,3);
    // Factor of (4 pi / 3) ignored in Tbar because it is cancelled in the following
    // In ZelDovich & Raizer, the \sqrt in the equation following 10.44 should
    //    /.
    Real Tc = pow(xi0,3)*pow((p*xi0*xi0/(6.*p+4.)),1./p)*Tbar;

    //std::cout << xi0 << " " << pfac << std::endl;
    if (verbose > 1 && level == 0 && ParallelDescriptor::IOProcessor()) {
      std::cout << "Front velocity = " << xf / (time * (3.0 * p + 2.0))
	   << "   Front position = " << xf << std::endl;
    }

    MultiFab& S_new = castro->get_new_data(State_Type);
    MultiFab& T_new = *plotvar[level];

    MultiFab temp(grids,dmap,1,0);

    {
      MultiFab fkp(grids,dmap,1,0);

      get_frhoe(temp, S_new);
      get_planck_and_temp(fkp, temp, S_new);
    }

    const Geometry& geom = parent->Geom(level);
    const Real *dx = geom.CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter ti(T_new); ti.isValid(); ++ti) {
      int i = ti.index();
      const Box &reg = grids[i];

      RealBox gridloc = RealBox(grids[i], geom.CellSize(), geom.ProbLo());

      anatw2(T_new[ti].dataPtr(), ARLIM(reg.loVect()), ARLIM(reg.hiVect()),
	     temp[ti].dataPtr(),
	     p, xf, Tc, dx, gridloc.lo(), reg.loVect());
    }

    return;
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
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter exi(exch,true); exi.isValid(); ++exi) {
	const Box& reg = exi.tilebox();

	cexch(ARLIM(reg.loVect()), ARLIM(reg.hiVect()),
	      BL_TO_FORTRAN(exch[exi]),
	      BL_TO_FORTRAN(Er[exi]),
	      BL_TO_FORTRAN(fkp[exi]),
	      sigma, c);
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

                if (use_opacity_table_module) {
#pragma gpu box(bx) sync
                    ca_compute_planck(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                                      BL_TO_FORTRAN_ANYD(eta[mfi]),
                                      BL_TO_FORTRAN_ANYD(state[mfi]),
                                      dT);
                }
                else {
#pragma gpu box(bx) sync
                    fkpn(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                         BL_TO_FORTRAN_ANYD(eta[mfi]),
                         const_kappa_p, kappa_p_exp_m, kappa_p_exp_n,
                         kappa_p_exp_p, nugroup[igroup], prop_temp_floor,
                         BL_TO_FORTRAN_ANYD(temp[mfi]),
                         BL_TO_FORTRAN_ANYD(state[mfi]),
                         dT);
                }

	    }

	    c_v.resize(bx);
            Elixir c_v_elix = c_v.elixir();

	    get_c_v(c_v, temp[mfi], state[mfi], bx);

#pragma gpu box(bx) sync
	    ceta2(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                  BL_TO_FORTRAN_ANYD(eta[mfi]),
		  BL_TO_FORTRAN_ANYD(etainv[mfi]),
                  BL_TO_FORTRAN_N_ANYD(state[mfi], Density),
                  BL_TO_FORTRAN_ANYD(temp[mfi]),
                  BL_TO_FORTRAN_ANYD(c_v),
                  BL_TO_FORTRAN_ANYD(fkp[mfi]),
                  BL_TO_FORTRAN_N_ANYD(Er[mfi], igroup),
                  dT, delta_t, sigma, c,
                  underrel, lag_planck);
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

  relative = 0.0;
  absolute = 0.0;
  Real theta = 1.0;

#ifdef _OPENMP
#pragma omp parallel reduction(max:relative, absolute)
#endif
  for (MFIter mfi(eta, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
      ceup(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
           AMREX_MFITER_REDUCE_MAX(&relative),
           AMREX_MFITER_REDUCE_MAX(&absolute),
           BL_TO_FORTRAN_ANYD(frhoes[mfi]),
           BL_TO_FORTRAN_ANYD(frhoem[mfi]),
           BL_TO_FORTRAN_ANYD(eta[mfi]),
           BL_TO_FORTRAN_ANYD(etainv[mfi]),
           BL_TO_FORTRAN_ANYD(dflux_old[mfi]),
           BL_TO_FORTRAN_ANYD(dflux_new[mfi]),
           BL_TO_FORTRAN_ANYD(exch[mfi]),
           delta_t, theta);
  }

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

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
      Real relative_priv = 0.0;
      Real absolute_priv = 0.0;
      Real theta = 1.0;

      for (MFIter mfi(eta,true); mfi.isValid(); ++mfi) {
	  const Box &reg = mfi.tilebox();
	  ceupdterm(ARLIM(reg.loVect()), ARLIM(reg.hiVect()),
		    relative_priv, absolute_priv,
		    BL_TO_FORTRAN(frhoes[mfi]),
		    frhoem[mfi].dataPtr(), eta[mfi].dataPtr(), etainv[mfi].dataPtr(),
		    dflux_old[mfi].dataPtr(), dflux_new[mfi].dataPtr(),
		    exch[mfi].dataPtr(), Dterm[mfi].dataPtr(), delta_t, theta);
      }
#ifdef _OPENMP
#pragma omp critical (rad_ceupdterm)
#endif
      {
	  relative = std::max(relative, relative_priv);
	  absolute = std::max(absolute, absolute_priv);
      }
  }

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

  relative = 0.0;
  absolute = 0.0;

#ifdef _OPENMP
#pragma omp parallel
#endif
  {
      Real relative_priv = 0.0;
      Real absolute_priv = 0.0;
      Real theta = 1.0;
      FArrayBox c_v;

      for (MFIter mfi(eta,true); mfi.isValid(); ++mfi) {
	  const Box &reg = mfi.tilebox();

	  c_v.resize(reg);
	  get_c_v(c_v, temp[mfi], state[mfi], reg);

	  nceup(ARLIM(reg.loVect()), ARLIM(reg.hiVect()),
		relative_priv, absolute_priv,
		BL_TO_FORTRAN(frhoes[mfi]),
		frhoem[mfi].dataPtr(), eta[mfi].dataPtr(), etainv[mfi].dataPtr(),
		BL_TO_FORTRAN(Er_new[mfi]),
		dflux_old[mfi].dataPtr(), dflux_new[mfi].dataPtr(),
		temp[mfi].dataPtr(), fkp[mfi].dataPtr(), c_v.dataPtr(),
		BL_TO_FORTRAN(state[mfi]),
		sigma, c, delta_t, theta);
      }
#ifdef _OPENMP
#pragma omp critical (rad_ceupdterm)
#endif
      {
	  relative = std::max(relative, relative_priv);
	  absolute = std::max(absolute, absolute_priv);
      }
  }

  ParallelDescriptor::ReduceRealMax(relative);
  ParallelDescriptor::ReduceRealMax(absolute);
}

void Radiation::state_update(MultiFab& state,
                             MultiFab& frhoes,
			     MultiFab& temp)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter si(state,true); si.isValid(); ++si) {
	const Box& fbox = frhoes[si].box();
	const Box& reg = si.tilebox();

	cetot(ARLIM(reg.loVect()), ARLIM(reg.hiVect()),
	      BL_TO_FORTRAN(state[si]),
	      BL_TO_FORTRAN(frhoes[si]));

	if (do_real_eos == 0) {

	    temp[si].copy(frhoes[si],reg);

#pragma gpu box(reg) sync
	    ca_compute_temp_given_cv
		(AMREX_INT_ANYD(reg.loVect()), AMREX_INT_ANYD(reg.hiVect()),
		 BL_TO_FORTRAN_ANYD(temp[si]),
		 BL_TO_FORTRAN_ANYD(state[si]),
		 const_c_v, c_v_exp_m, c_v_exp_n);

	    state[si].copy(temp[si],reg,0,reg,Temp,1);
	}
    }
}


void Radiation::extrapolateBorders(MultiFab& f, int indx)
{
  BL_PROFILE("Radiation::extrapolateBorders");

  BL_ASSERT(f.nGrow() >= 1);

#ifdef _OPENMP
#pragma omp parallel
#endif
  for(MFIter mfi(f); mfi.isValid(); ++mfi) {
    int i = mfi.index();

    const Box& reg  = f.box(i);

    bextrp(BL_TO_FORTRAN_N(f[mfi],indx),
	   ARLIM(reg.loVect()), ARLIM(reg.hiVect()));
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
      sold_tmp.copy(S_old, Rad, 0, 1);
      sold_tmp.FillBoundary(geom.periodicity());
    }

    if (need_new_data) {
      snew_tmp.define(grids, dmap, 1, n_grow);
      snew_tmp.setVal(0.0); // need legal numbers for linComb below
      snew_tmp.copy(S_new, Rad, 0, 1);
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

    if (do_real_eos == 1) {

#pragma gpu box(reg) sync
      ca_compute_c_v
          (AMREX_INT_ANYD(reg.loVect()), AMREX_INT_ANYD(reg.hiVect()),
           BL_TO_FORTRAN_ANYD(c_v),
           BL_TO_FORTRAN_ANYD(temp),
           BL_TO_FORTRAN_ANYD(state));

    }
    else if (do_real_eos == 0) {

        if (c_v_exp_m == 0.0 && c_v_exp_n == 0.0) {

            Array4<Real> const c_v_arr = c_v.array();
            Real c_v_new = const_c_v;
            AMREX_PARALLEL_FOR_3D(reg, i, j, k, { c_v_arr(i,j,k) = c_v_new; });
            Gpu::synchronize();

	}
	else {

#pragma gpu box(reg) sync
	    gcv(AMREX_INT_ANYD(reg.loVect()), AMREX_INT_ANYD(reg.hiVect()),
		BL_TO_FORTRAN_ANYD(c_v),
		BL_TO_FORTRAN_ANYD(temp),
		const_c_v, c_v_exp_m, c_v_exp_n,
                prop_temp_floor,
		BL_TO_FORTRAN_ANYD(state));

	}
    }
    else {
	amrex::Error("ERROR Radiation::get_c_v  do_real_eos < 0");
    }
}

void Radiation::get_rosseland_from_temp(FArrayBox& kappa_r,
                                        FArrayBox& temp,
                                        FArrayBox& state,
                                        const Box& reg,
					int igroup)
{
  BL_ASSERT(kappa_r.nComp() == Radiation::nGroups);

  state.copy(temp,reg,0,reg,Temp,1);

  if (use_opacity_table_module) {
#pragma gpu box(reg) sync
      ca_compute_rosseland(AMREX_INT_ANYD(reg.loVect()), AMREX_INT_ANYD(reg.hiVect()),
                           BL_TO_FORTRAN_ANYD(kappa_r),
                           BL_TO_FORTRAN_ANYD(state));
  }
  else if (const_scattering > 0.0) {
#pragma gpu box(reg) sync
      rosse1s(AMREX_ARLIM_ANYD(reg.loVect()), AMREX_ARLIM_ANYD(reg.hiVect()),
	      const_kappa_r, kappa_r_exp_m, kappa_r_exp_n,
	      kappa_r_exp_p, const_scattering,
              scattering_exp_m, scattering_exp_n,
	      scattering_exp_p, nugroup[igroup],
	      prop_temp_floor, kappa_r_floor,
	      BL_TO_FORTRAN_ANYD(state),
              BL_TO_FORTRAN_N_ANYD(kappa_r, igroup));
  }
  else {
#pragma gpu box(reg) sync
    rosse1(AMREX_ARLIM_ANYD(reg.loVect()), AMREX_ARLIM_ANYD(reg.hiVect()),
	   const_kappa_r, kappa_r_exp_m, kappa_r_exp_n,
	   kappa_r_exp_p, nugroup[igroup],
	   prop_temp_floor, kappa_r_floor,
	   BL_TO_FORTRAN_ANYD(state),
           BL_TO_FORTRAN_N_ANYD(kappa_r, igroup));
  }
}

void Radiation::get_frhoe(MultiFab& frhoe,
                          MultiFab& state)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter si(state,TilingIfNotGPU()); si.isValid(); ++si) {
	const Box& reg = si.tilebox();

#pragma gpu box(reg)
	cfrhoe(AMREX_INT_ANYD(reg.loVect()), AMREX_INT_ANYD(reg.hiVect()),
	       BL_TO_FORTRAN_ANYD(frhoe[si]),
	       BL_TO_FORTRAN_ANYD(state[si]));
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

        if (do_real_eos > 0) {
#pragma gpu box(bx)
            ca_compute_temp_given_rhoe
                (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                 BL_TO_FORTRAN_ANYD(temp[mfi]),
                 BL_TO_FORTRAN_ANYD(state[mfi]),
                 0);
        }
        else if (do_real_eos == 0) {
#pragma gpu box(bx)
            gtemp(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                  BL_TO_FORTRAN_ANYD(temp[mfi]),
                  const_c_v, c_v_exp_m, c_v_exp_n,
                  BL_TO_FORTRAN_ANYD(state[mfi]),
                  0);
        }
        else {
            amrex::Error("ERROR Radiation::get_planck_and_temp  do_real_eos < 0");
        }

        if (use_opacity_table_module) {
#pragma gpu box(bx)
            ca_compute_planck(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                              BL_TO_FORTRAN_ANYD(fkp[mfi]),
                              BL_TO_FORTRAN_ANYD(state[mfi]),
                              0.0);
        }
        else {
#pragma gpu box(bx) sync
            fkpn(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                 BL_TO_FORTRAN_ANYD(fkp[mfi]),
                 const_kappa_p, kappa_p_exp_m, kappa_p_exp_n,
                 kappa_p_exp_p, nugroup[igroup], prop_temp_floor,
                 BL_TO_FORTRAN_ANYD(temp[mfi]),
                 BL_TO_FORTRAN_ANYD(state[mfi]),
                 0.0);
        }

#pragma gpu box(bx) sync
        nfloor(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
               BL_TO_FORTRAN_ANYD(temp[mfi]),
               temp_floor, temp[mfi].nComp());

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

#pragma gpu box(bx)
          cfrhoe(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                 BL_TO_FORTRAN_ANYD(frhoe),
                 BL_TO_FORTRAN_ANYD(state[mfi]));

          if (do_real_eos > 0) {
#pragma gpu box(bx)
              ca_compute_temp_given_rhoe
                  (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                   BL_TO_FORTRAN_ANYD(frhoe),
                   BL_TO_FORTRAN_ANYD(state[mfi]),
                   1);
          }
          else if (do_real_eos == 0) {
#pragma gpu box(bx)
              gtemp(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                    BL_TO_FORTRAN_ANYD(frhoe),
                    const_c_v, c_v_exp_m, c_v_exp_n,
                    BL_TO_FORTRAN_ANYD(state[mfi]),
                    1);
          }
          else {
              amrex::Abort("ERROR Radiation::get_rosseland do_real_eos < 0");
          }

          if (use_opacity_table_module) {
#pragma gpu box(bx)
              ca_compute_rosseland(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                                   BL_TO_FORTRAN_ANYD(kappa_r[mfi]),
                                   BL_TO_FORTRAN_ANYD(state[mfi]));
          }
          else if (const_scattering > 0.0) {
#pragma gpu box(bx)
              rosse1s(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                      const_kappa_r, kappa_r_exp_m, kappa_r_exp_n,
                      const_scattering, kappa_r_exp_p,
                      scattering_exp_m, scattering_exp_n,
                      scattering_exp_p, nugroup[igroup],
                      prop_temp_floor, kappa_r_floor,
                      BL_TO_FORTRAN_ANYD(state[mfi]),
                      BL_TO_FORTRAN_N_ANYD(kappa_r[mfi], igroup));
          }
          else {
#pragma gpu box(bx)
              rosse1(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                     const_kappa_r, kappa_r_exp_m, kappa_r_exp_n,
                     kappa_r_exp_p, nugroup[igroup],
                     prop_temp_floor, kappa_r_floor,
                     BL_TO_FORTRAN_ANYD(state[mfi]),
                     BL_TO_FORTRAN_N_ANYD(kappa_r[mfi], igroup));
          }

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

  BL_ASSERT(kappa_r.nGrow() == 1);
  BL_ASSERT(temp.nGrow()    == 0);

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter si(state); si.isValid(); ++si) {
      get_rosseland_from_temp(kappa_r[si], temp[si], state[si], si.tilebox(), igroup);
  }

  kappa_r.FillBoundary(geom.periodicity());
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

        for (int dir = 0; dir < BL_SPACEDIM; dir++) {
          const Orientation lo_face(dir,Orientation::low);
          const Orientation hi_face(dir,Orientation::high);

          for (FabSetIter fsi(ref_sync_flux[lo_face]);
               fsi.isValid(); ++fsi) {

            rfface(BL_TO_FORTRAN(ref_sync_flux[lo_face][fsi]),
		   BL_TO_FORTRAN_N(crse_sync_flux[lo_face][fsi], indx),
		   dir, ref_rat.getVect());
          }

          for (FabSetIter fsi(ref_sync_flux[hi_face]);
               fsi.isValid(); ++fsi) {

            rfface(BL_TO_FORTRAN(ref_sync_flux[hi_face][fsi]),
		   BL_TO_FORTRAN_N(crse_sync_flux[hi_face][fsi], indx),
		   dir, ref_rat.getVect());
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
                               Array<MultiFab, BL_SPACEDIM>& R,
                               MultiFab& kappa_r, int kcomp,
                               MultiFab& Er, int igroup,
                               int limiter, int nGrow_Er, int Rcomp)
{
  BL_PROFILE("Radiation::scaledGradient");
  BL_ASSERT(kappa_r.nGrow() == 1);

  int Ercomp = igroup;

  MultiFab Erbtmp;
  if (nGrow_Er == 0) { // default value
    if (limiter > 0) {
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
    Ercomp = 0;
  }

  MultiFab& Erborder = (nGrow_Er==0) ? Erbtmp : Er;

  const Real* dx = parent->Geom(level).CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (int idim = 0; idim < BL_SPACEDIM; ++idim) {

      for (MFIter mfi(R[idim], TilingIfNotGPU()); mfi.isValid(); ++mfi) {

          const Box& nbx = mfi.tilebox();  // note that R is edge based
          const Box& bx = amrex::enclosedCells(nbx);

          Array4<Real> const R_arr = R[idim].array(mfi);

          if (limiter == 0) {

              int comp = Rcomp;
              AMREX_PARALLEL_FOR_3D(nbx, i, j, k, { R_arr(i,j,k,comp) = 0.0; });

          }
          else {

              int include_cross_terms = 0;

              if (limiter == 1) {
                  include_cross_terms = 0;
              } else if (limiter == 2) {
                  include_cross_terms = 1;
              } else {
                  amrex::Abort("Unknown limiter");
              }

#pragma gpu box(bx)
              scgrd(AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
                    BL_TO_FORTRAN_N_ANYD(R[idim][mfi], Rcomp),
                    idim, AMREX_REAL_ANYD(dx),
                    BL_TO_FORTRAN_N_ANYD(kappa_r[mfi], kcomp),
                    BL_TO_FORTRAN_ANYD(Erborder[mfi]),
                    include_cross_terms);

          }
      }
  }

}

// On input, lambda should contain scaled gradient.
// On output this will be overwritten with the flux limiter.

void Radiation::fluxLimiter(int level,
                            Array<MultiFab, BL_SPACEDIM>& lambda,
                            int limiter, int lamcomp)
{
    BL_PROFILE("Radiation:fluxLimiter");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
        for (MFIter mfi(lambda[idim], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
            flxlim(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                   BL_TO_FORTRAN_N_ANYD(lambda[idim][mfi], lamcomp),
                   limiter);
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

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	FArrayBox temp, c_v, kp, kp2;
	for (MFIter mfi(kappa_r,true); mfi.isValid(); ++mfi)
	{
	    const Box& reg = mfi.growntilebox();

	    temp.resize(reg);

#pragma gpu box(reg) sync
            cfrhoe(AMREX_INT_ANYD(reg.loVect()), AMREX_INT_ANYD(reg.hiVect()),
                   BL_TO_FORTRAN_ANYD(temp),
                   BL_TO_FORTRAN_ANYD(S[mfi]));

	    if (do_real_eos > 0) {
#pragma gpu box(reg) sync
	      ca_compute_temp_given_rhoe
                  (AMREX_INT_ANYD(reg.loVect()), AMREX_INT_ANYD(reg.hiVect()),
                   BL_TO_FORTRAN_ANYD(temp),
                   BL_TO_FORTRAN_ANYD(S[mfi]),
                   0);
	    }
	    else if (do_real_eos == 0) {
#pragma gpu box(reg) sync
		gtemp(AMREX_INT_ANYD(reg.loVect()), AMREX_INT_ANYD(reg.hiVect()),
		      BL_TO_FORTRAN_ANYD(temp),
		      const_c_v, c_v_exp_m, c_v_exp_n,
		      BL_TO_FORTRAN_ANYD(S[mfi]),
                      0);
	    }

	    c_v.resize(reg);
	    get_c_v(c_v, temp, S[mfi], reg);

	    S[mfi].copy(temp,reg,0,reg,Temp,1);

	    // compute rosseland
	    if (use_opacity_table_module) {
#pragma gpu box(reg) sync
                ca_compute_rosseland(AMREX_INT_ANYD(reg.loVect()), AMREX_INT_ANYD(reg.hiVect()),
                                     BL_TO_FORTRAN_ANYD(kappa_r[mfi]),
                                     BL_TO_FORTRAN_ANYD(S[mfi]));
	    }
	    else if (const_scattering > 0.0) {
#pragma gpu box(reg) sync
	      rosse1s(AMREX_ARLIM_ANYD(reg.loVect()), AMREX_ARLIM_ANYD(reg.hiVect()),
		      const_kappa_r, kappa_r_exp_m, kappa_r_exp_n,
		      kappa_r_exp_p, const_scattering,
                      scattering_exp_m, scattering_exp_n,
		      scattering_exp_p, nugroup[igroup],
		      prop_temp_floor, kappa_r_floor,
		      BL_TO_FORTRAN_ANYD(S[mfi]),
                      BL_TO_FORTRAN_N_ANYD(kappa_r[mfi], igroup));
	    }
	    else {
#pragma gpu box(reg) sync
	      rosse1(AMREX_ARLIM_ANYD(reg.loVect()), AMREX_ARLIM_ANYD(reg.hiVect()),
		     const_kappa_r, kappa_r_exp_m, kappa_r_exp_n,
		     kappa_r_exp_p, nugroup[igroup],
		     prop_temp_floor, kappa_r_floor,
		     BL_TO_FORTRAN_ANYD(S[mfi]),
                     BL_TO_FORTRAN_N_ANYD(kappa_r[mfi], igroup));
	    }

	    if (use_opacity_table_module) {
		std::cout << "SGFLDSolver Solver does not support both use_opacity_table_module=1 "
		     << "and Er_Lorentz_term=1; try comoving=1 or Er_Lorentz_term=0 "
		     << "when use_opacity_table_module=1"<<std::endl;
		amrex::Abort("SGFLDSolver: try comoving=1 or Er_Lorentz_term=0");
	    }


	    kp.resize(reg);
#pragma gpu box(reg) sync
	    fkpn(AMREX_INT_ANYD(reg.loVect()), AMREX_INT_ANYD(reg.hiVect()),
		 BL_TO_FORTRAN_ANYD(kp),
		 const_kappa_p, kappa_p_exp_m, kappa_p_exp_n,
		 kappa_p_exp_p, nugroup[igroup], prop_temp_floor,
		 BL_TO_FORTRAN_ANYD(temp),
		 BL_TO_FORTRAN_ANYD(S[mfi]),
                 0.0);

	    kp2.resize(reg);
	    temp.plus(dT, 0, 1);

#pragma gpu box(reg) sync
	    fkpn(AMREX_INT_ANYD(reg.loVect()), AMREX_INT_ANYD(reg.hiVect()),
		 BL_TO_FORTRAN_ANYD(kp2),
		 const_kappa_p, kappa_p_exp_m, kappa_p_exp_n,
		 kappa_p_exp_p, nugroup[igroup], prop_temp_floor,
		 BL_TO_FORTRAN_ANYD(temp),
		 BL_TO_FORTRAN_ANYD(S[mfi]),
                 0.0);

	    temp.plus(-dT, 0, 1);

	    ca_get_v_dcf(reg.loVect(), reg.hiVect(),
			 BL_TO_FORTRAN(Er[mfi]),
			 BL_TO_FORTRAN(S[mfi]),
			 BL_TO_FORTRAN(temp),
			 BL_TO_FORTRAN(c_v),
			 BL_TO_FORTRAN(kappa_r[mfi]),
			 BL_TO_FORTRAN(kp),
			 BL_TO_FORTRAN(kp2),
			 &dT, &delta_t, &sigma, &c,
			 BL_TO_FORTRAN(v[mfi]),
			 BL_TO_FORTRAN(dcf[mfi]));
	}
    }
}

void Radiation::update_dcf(MultiFab& dcf, MultiFab& etainv, MultiFab& kp, MultiFab& kr,
			   const Geometry& geom)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(dcf,true); mfi.isValid(); ++mfi) {
	const Box& bx = mfi.growntilebox();
	ca_update_dcf(bx.loVect(), bx.hiVect(),
		      BL_TO_FORTRAN(dcf[mfi]),
		      BL_TO_FORTRAN(etainv[mfi]),
		      BL_TO_FORTRAN(kp[mfi]),
		      BL_TO_FORTRAN(kr[mfi]));
    }

    dcf.FillBoundary(geom.periodicity());
}

void Radiation::EstTimeStep(Real & estdt, int level)
{
  Castro *castro        = dynamic_cast<Castro*>(&parent->getLevel(level));
  const Geometry& geom = parent->Geom(level);
  Real time = castro->get_state_data(Rad_Type).curTime();

  if (nNeutrinoSpecies > 0 &&
      nNeutrinoGroups[0] > 0) {
//    Real derat = deltaEnergyRatMax(level);
    Real dTrat = deltaTRatMax(level);
    Real dye   = deltaYeMax(level);

    //    Real fac = std::min(deltaEnergyTol() / (derat + 1.e-20),
    //		deltaTTol()      / (dTrat + 1.e-20));
    //    fac = std::min(fac, deltaYeTol()     / (dye   + 1.e-20));
    Real fac = deltaTTol()  / (dTrat + 1.e-20);
    fac = std::min(fac, deltaYeTol()     / (dye   + 1.e-20));

    Real estdt_rad = parent->dtLevel(level);
    estdt_rad *= fac;

    if (verbose && ParallelDescriptor::IOProcessor())
      std::cout << "radiation timestep at level " << level
		<< ":  estdt_rad = " << estdt_rad << std::endl;

    estdt = std::min(estdt, estdt_rad);
  }

  if (RadTests::do_thermal_wave_cgs) {
    ParmParse pp("radiation");
    Real cfl = -1.0;
    pp.query("thermal_wave_cfl", cfl);
    if (cfl > 0.0) {
      Real rhocv, Q;
      pp.get("thermal_wave_rhocv", rhocv);
      pp.get("thermal_wave_Eexp", Q);
      Q /= rhocv;

      Real a = (16.0 * sigma) / (3.0 * const_kappa_r) / rhocv;
      Real p = kappa_r_exp_n + 3.0;

      Real pfac = exp(lgamma(2.5 + 1.0 / p) - lgamma(1.0 + 1.0 / p) - lgamma(1.5));

      Real xi0  = pow((3.0 * p + 2.0) / (pow(2.0, p - 1.0) * p *
					 pow(M_PI, p)),
		      1.0 / (3.0 * p + 2.0)) *
	pow(pfac, p / (3.0 * p + 2.0));
      Real xf = xi0 * pow(a * pow(Q,p) * time, 1.0 / (3.0 * p + 2.0));
      Real vf = xf / (time * (3.0 * p + 2.0));

      const Real *dx = geom.CellSize();
      Real estdt_rad = cfl * dx[0] / vf;

      if (verbose && ParallelDescriptor::IOProcessor())
	std::cout << "radiation timestep at level " << level
		  << ":  estdt_rad = " << estdt_rad << std::endl;

      estdt = std::min(estdt, estdt_rad);
    }
  }
}

void Radiation::set_current_group(int igroup)
{
#ifdef NEUTRINO
  if (radiation_type == Photon) {
    current_group_number = igroup;
    current_group_name = "Photon";
  }
  else if (igroup < nNeutrinoGroups[0]) {
    current_group_number = igroup;
    current_group_name = "Electron";
  }
  else if (nNeutrinoSpecies >= 2 && igroup < nNeutrinoGroups[0]+nNeutrinoGroups[1]) {
    current_group_number = igroup - nNeutrinoGroups[0];
    current_group_name = "Anti Electron";
  }
  else if (nNeutrinoSpecies == 3 && igroup >= nNeutrinoGroups[0]+nNeutrinoGroups[1]) {
    current_group_number = igroup - (nNeutrinoGroups[0]+nNeutrinoGroups[1]);
    current_group_name = "Muon";
  }
  else {
    amrex::Abort("Something is wrong!  Maybe there are more than three neutrino flavors.");
  }
#else
  current_group_number = igroup;
#endif
}


void Radiation::filter_prim(int level, MultiFab& State)
{
  Castro *castro = dynamic_cast<Castro*>(&parent->getLevel(level));
  const BoxArray& grids = castro->boxArray();
  const DistributionMapping& dmap = castro->DistributionMap();
  const Geometry& geom = parent->Geom(level);

  const int*  domain_lo = geom.Domain().loVect();
  const int*  domain_hi = geom.Domain().hiVect();
  const Real* dx        = geom.CellSize();
  const Real* prob_lo   = geom.ProbLo();

  int ngrow = filter_prim_T;
  int ncomp = State.nComp();
  Real time = castro->get_state_data(Rad_Type).curTime();

  MultiFab mask(grids,dmap,1,ngrow);
  mask.setVal(-1.0,ngrow);
  mask.setVal( 0.0,0);
  mask.FillBoundary(geom.periodicity());

  if (level < parent->finestLevel())
  {
      BoxArray baf = parent->boxArray(level+1);
      baf.coarsen(parent->refRatio(level));

#ifdef _OPENMP
#pragma omp parallel
#endif
      for (MFIter mfi(mask); mfi.isValid(); ++mfi)
      {
	  FArrayBox& mask_fab = mask[mfi];
	  const Box& mask_box = mask_fab.box();

	  const std::vector< std::pair<int,Box> >& isects = baf.intersections(mask_box);

	  for (int ii = 0; ii < isects.size(); ii++) {
	      mask_fab.setVal(1.0, isects[ii].second, 0);
	  }
      }
  }

  FillPatchIterator fpi(*castro,State,ngrow,time,State_Type,0,ncomp);
  MultiFab& S_fp = fpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(State,true); mfi.isValid(); ++mfi)
  {
      const Box& bx = mfi.tilebox();

      const RealBox& gridloc = RealBox(bx, dx, prob_lo);
      const Real* xlo = gridloc.lo();

      ca_filt_prim(bx.loVect(), bx.hiVect(),
		   BL_TO_FORTRAN( S_fp[mfi]),
		   BL_TO_FORTRAN(State[mfi]),
		   BL_TO_FORTRAN( mask[mfi]),
		   &filter_prim_T, &filter_prim_S,
		   domain_lo, domain_hi,
		   dx, xlo, prob_lo,
		   &time, &level);
  }
}
