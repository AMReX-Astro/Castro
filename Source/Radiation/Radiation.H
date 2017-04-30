#ifndef _Radiation_H_
#define _Radiation_H_

#include <AMReX_AmrLevel.H>

#include "Castro.H"
#include "RadTypes.H"
#include "RadBndry.H"
#include "MGRadBndry.H"
#include "RadSolve.H"
#include <AMReX_FluxRegister.H>
#include <AMReX_Tuple.H>

class Radiation {

public:

  enum Solver_Type { InvalidSolver = -1,
                     SingleGroupSolver = 0,
                     // ProtoMultiGroupSolver = 1,  // deleted
                     //  NeutrinoSolver = 2,  // deleted
                     //  SGFLDHybrid = 3,  // deleted
		     //  SGFLDImplicit = 4, // deleted
                     SGFLDSolver = 5,
                     MGFLDSolver = 6};

  static Solver_Type SolverType;

  enum Radiation_Type {Unknown = 0, Neutrino = 1, Photon = 2};
  static Radiation_Type radiation_type;

  static amrex::Real radtoE, Etorad, radfluxtoF;

  static int do_multigroup;
  static int nGroups;
  static int nNeutrinoSpecies;
  static amrex::Array<int> nNeutrinoGroups;
  static int plot_neutrino_group_energies_per_MeV;
  static int plot_neutrino_group_energies_total;

  static int rad_hydro_combined;
  static int comoving;
  static int Er_Lorentz_term; 
  static int fspace_advection_type;
  static int do_inelastic_scattering;    // Compton scattering for photon only

  static int use_analytic_solution; // compare against an analytic solution

  static int plot_lambda, plot_kappa_p, plot_kappa_r;
  static int plot_lab_Er, plot_lab_flux, plot_com_flux;
  static int icomp_lambda, icomp_kp, icomp_kr;
  static int icomp_lab_Er, icomp_lab_Fr, icomp_com_Fr; 
  static int nplotvar;   // extra plot variables
  static amrex::Array<std::string> plotvar_names;

  static int filter_lambda_T, filter_lambda_S;
  static int filter_prim_int, filter_prim_T, filter_prim_S;

  static int accelerate;        // controls multigroup convergence acceleration

  static void read_static_params();

  // Physical constants for use in multigroup, made static and public
  // so that other modules such as RadSolve can get to them easily.
  // Values are set in the Radiation constructor.

  static amrex::Real convert_MeV_erg, clight, hPlanck, kBoltz, Avogadro;

  // Values used internally for speed of light, Stefan-Boltzmann constant,
  // radiation constant.  (For contrived test problems, these may differ
  // from the proper physical values, so be careful if you use them!)

  static amrex::Real c, sigma, aRad;

  static int current_group_number;
  static std::string current_group_name;

  static amrex::Real flatten_pp_threshold;
  static int pure_hydro;

  static int do_real_eos;  // if 0,  c_v = const_c_v * (rho**m) * (temp**(-n))


  static amrex::Real const_c_v, const_kappa_p, const_kappa_r, const_scattering;
  static amrex::Real c_v_exp_m, c_v_exp_n;
  static amrex::Real kappa_p_exp_m, kappa_p_exp_n;
  static amrex::Real kappa_r_exp_m, kappa_r_exp_n;
  static amrex::Real scattering_exp_m, scattering_exp_n;
  static amrex::Real kappa_p_exp_p, kappa_r_exp_p, scattering_exp_p;
  static amrex::Real prop_temp_floor;

  amrex::Array<std::unique_ptr<amrex::MultiFab> > plotvar;

  Radiation(amrex::Amr* Parent, class Castro* castro, int restart = 0);
  ~Radiation() { }

  void regrid(int level, const amrex::BoxArray& grids,
	      const amrex::DistributionMapping& dmap);

  void close(int level);

  void restart(int level, const amrex::BoxArray& grids,
	       const amrex::DistributionMapping& dmap,
               const std::string& dir, std::istream& is);

  void checkPoint(int level,
                  const std::string& dir,
                  std::ostream&  os,
                  amrex::VisMF::How     how);

  // access to group information:

  amrex::Real group_center(int i) {
    return nugroup[i];
  }

  amrex::Real group_weight(int i) {
    return dnugroup[i];
  }

  // diagnostics for adjusting timestep:

  amrex::Real deltaEnergyRatMax(int lev) {
    return delta_e_rat_level[lev];
  }

  amrex::Real deltaTRatMax(int lev) {
    return delta_T_rat_level[lev];
  }

  amrex::Real deltaYeMax(int lev) {
    return delta_Ye_level[lev];
  }

  amrex::Real deltaEnergyTol() {
    return delta_e_rat_dt_tol;
  }

  amrex::Real deltaTTol() {
    return delta_T_rat_dt_tol;
  }

  amrex::Real deltaYeTol() {
    return delta_Ye_dt_tol;
  }

  // access to conservation flux register for energy diagnostics:

  amrex::FluxRegister* consRegister(int level) {
      return flux_cons_old[level].get();
  }

  void analytic_solution(int level);

  void post_init(int level);

  void pre_timestep(int level);

  void init_flux(int level, int ncycle);

  void deferred_sync_setup(int level);

  void deferred_sync(int level, amrex::MultiFab& rhs, int indx);

  void single_group_update(int level, int iteration, int ncycle);

  void MGFLD_implicit_update(int level, int iteration, int ncycle);

  void inelastic_scattering(int level);

  void write_groups(ostream& os);

  void reflux(int level);

  void compute_limiter(int level, const amrex::BoxArray& grids,
		       const amrex::MultiFab &Sborder, 
		       const amrex::MultiFab &Erborder,
		       amrex::MultiFab &lamborder);

  void estimate_gamrPr(const amrex::FArrayBox& state, const amrex::FArrayBox& Er, 
		       amrex::FArrayBox& gPr, const amrex::Real* dx, const amrex::Box& box);

  void EstTimeStep(amrex::Real& estdt, int level);

  void filter_prim(int level, amrex::MultiFab& State);

protected:

  // internal utility functions:

  void compute_exchange(amrex::MultiFab& exch, amrex::MultiFab& Er, amrex::MultiFab& fkp,
			int igroup = -1);

  void compute_eta(amrex::MultiFab& eta, amrex::MultiFab& etainv, amrex::MultiFab& state,
		   amrex::MultiFab& temp, amrex::MultiFab& fkp, amrex::MultiFab& Er,
		   amrex::Real delta_t, amrex::Real c, amrex::Real underrel, int lag_planck,
		   int igroup = 0);

  void internal_energy_update(amrex::Real& relative, amrex::Real& absolute,
			      amrex::MultiFab& frhoes, amrex::MultiFab& frhoem,
			      amrex::MultiFab& eta, amrex::MultiFab& etainv,
			      amrex::MultiFab& dflux_old, amrex::MultiFab& dflux_new,
                              amrex::MultiFab& exchange, amrex::Real delta_t);

  void internal_energy_update(amrex::Real& relative, amrex::Real& absolute,
			      amrex::MultiFab& frhoes, amrex::MultiFab& frhoem,
			      amrex::MultiFab& eta, amrex::MultiFab& etainv,
			      amrex::MultiFab& dflux_old, amrex::MultiFab& dflux_new,
			      amrex::MultiFab& exchange, amrex::MultiFab& Dterm,
			      amrex::Real delta_t);

  void nonconservative_energy_update(amrex::Real& relative, amrex::Real& absolute,
				     amrex::MultiFab& frhoes, amrex::MultiFab& frhoem,
				     amrex::MultiFab& eta, amrex::MultiFab& etainv,
				     amrex::MultiFab& Er_new,
				     amrex::MultiFab& dflux_old, amrex::MultiFab& dflux_new,
				     amrex::MultiFab& temp, amrex::MultiFab& fkp,
				     amrex::MultiFab& state, amrex::Real delta_t);

  void state_update(amrex::MultiFab& state,  amrex::MultiFab& frhoes, amrex::MultiFab& temp);

  void getBndryData(RadBndry& bd, amrex::MultiFab& Er, amrex::Real time, int level);

  // multigroup version
  void getBndryDataMG(MGRadBndry& mgbd, amrex::MultiFab& Er, amrex::Real time, int level);
  void getBndryDataMG_ga(MGRadBndry& mgbd, amrex::MultiFab& Er, int level);

  void filBndry(amrex::BndryRegister& bdry, int level, amrex::Real time);

  // Flux limiter functions, potentially for use by all update modules

  // Computes the scaled gradient for use in flux limiters & Eddington factors

  void scaledGradient(int level,
                      amrex::Tuple<amrex::MultiFab, BL_SPACEDIM>& R,
                      amrex::MultiFab& kappa_r, int kcomp,
                      amrex::MultiFab& Er, int igroup,
                      int limiter, int nGrow_Er=0, int Rcomp=0);

  // On input, lambda should contain scaled gradient.
  // On output this will be overwritten with the flux limiter.

  void fluxLimiter(int level,
                   amrex::Tuple<amrex::MultiFab, BL_SPACEDIM>& lambda,
                   int limiter, int lamcomp=0);

  // Fab versions of conversion functions.  All except frhoe use eos data.

  void get_frhoe(amrex::FArrayBox& rhoe, amrex::FArrayBox& state, const amrex::Box& reg);

  void get_c_v(amrex::FArrayBox& c_v, amrex::FArrayBox& temp,
               amrex::FArrayBox& state, const amrex::Box& reg);

  // temp contains frhoe on input:

  void get_planck_and_temp(amrex::FArrayBox& fkp, amrex::FArrayBox& temp,
			   amrex::FArrayBox& state, const amrex::Box& reg,
			   int igroup = 0, amrex::Real delta_t = 0.0);

  void get_rosseland_and_temp(amrex::FArrayBox& kappa_r, amrex::FArrayBox& temp,
			      amrex::FArrayBox& state, const amrex::Box& reg,
			      int igroup = 0);

  // temp contains temp on input:
  void get_planck_from_temp(amrex::FArrayBox& fkp, amrex::FArrayBox& temp,
                            amrex::FArrayBox& state, const amrex::Box& reg,
			    int igroup = 0);

  void get_rosseland_from_temp(amrex::FArrayBox& kappa_r, amrex::FArrayBox& temp,
			       amrex::FArrayBox& state, const amrex::Box& reg,
			       int igroup = 0);

  // MultiFab versions of conversion functions call Fab versions for eos data:

  void get_frhoe(amrex::MultiFab& frhoe, amrex::MultiFab& state);

  // Planck mean is only needed in the grid interiors, but must be
  // computed separately for each fluid (fkp.nVar() == N fluids).
  // temp contains frhoe on input:

  void get_planck_and_temp(amrex::MultiFab& fkp, amrex::MultiFab& temp, amrex::MultiFab& state,
			   int igroup = 0, amrex::Real delta_t = 0.0);

  // Rosseland mean is needed on an extended region (nGrow == 1).
  // Fills patches with state information on extended region:

  void get_rosseland(amrex::MultiFab& kappa_r, amrex::AmrLevel* hclaw, int igroup = 0);

  // For the interior opacities, temperature is assumed to be known (input):

  // Updates Rosseland mean in interior of level, leaves ghost values
  // at physical or coarse-fine interfaces untouched:

  void update_rosseland_from_temp(amrex::MultiFab& kappa_r,
				  amrex::MultiFab& temp, amrex::MultiFab& state,
				  const amrex::Geometry& geom,
				  int igroup = 0);

  // RadHydro
  void get_rosseland_v_dcf(amrex::MultiFab& kappa_r, amrex::MultiFab& v, amrex::MultiFab& dcf, 
			   amrex::Real delta_t, amrex::Real c,
			   amrex::AmrLevel* hclaw,
			   int igroup = 0);  

  // multigroup methods
  void get_groups(int verbose);

  void extrapolateBorders(amrex::MultiFab& f, int indx);

  void SGFLD_compute_rosseland(amrex::MultiFab& kappa_r, const amrex::MultiFab& state);
  void SGFLD_compute_rosseland(amrex::FArrayBox& kappa_r, const amrex::FArrayBox& state);
  void update_dcf(amrex::MultiFab& dcf, amrex::MultiFab& etainv, amrex::MultiFab& kp, amrex::MultiFab& kr,
		  const amrex::Geometry& geom);

  // <MGFLD>

  int use_dkdT;

  // <Shestakov-Bolstad>
  int use_WiensLaw;
  amrex::Real Tf_Wien;
  // </ Shestakov-Bolstad>
  int integrate_Planck;

  void check_convergence_er(amrex::Real& relative_in, amrex::Real& absolute_in, amrex::Real& error_er,
  			    const amrex::MultiFab& Er_new, const amrex::MultiFab& Er_pi, 
			    const amrex::MultiFab& kappa_p, 
			    const amrex::MultiFab& etaTz, const amrex::MultiFab& etaYz, 
			    const amrex::MultiFab& theTz, const amrex::MultiFab& theYz,
			    const amrex::MultiFab& temp_new, const amrex::MultiFab& Ye_new,
			    const amrex::BoxArray& grids, amrex::Real delta_t);
  void check_convergence_matt(const amrex::MultiFab& rhoe_new, const amrex::MultiFab& rhoe_star, 
			      const amrex::MultiFab& rhoe_step, const amrex::MultiFab& Er_new, 
			      const amrex::MultiFab& temp_new, const amrex::MultiFab& temp_star, 
			      const amrex::MultiFab& rhoYe_new, const amrex::MultiFab& rhoYe_star, 
			      const amrex::MultiFab& rhoYe_step, const amrex::MultiFab& rho,
			      const amrex::MultiFab& kappa_p, const amrex::MultiFab& jg,
			      const amrex::MultiFab& dedT, const amrex::MultiFab& dedY,
			      amrex::Real& rel_rhoe, amrex::Real& abs_rhoe, 
			      amrex::Real& rel_FT,   amrex::Real& abs_FT, 
			      amrex::Real& rel_T,    amrex::Real& abs_T, 
			      amrex::Real& rel_FY,   amrex::Real& abs_FY, 
			      amrex::Real& rel_Ye,   amrex::Real& abs_Ye,
			      const amrex::BoxArray& grids, amrex::Real delta_t);
  void compute_coupling(amrex::MultiFab& coupT, amrex::MultiFab& coupY, 
			const amrex::MultiFab& kappa_p, const amrex::MultiFab& Er_pi,
			const amrex::MultiFab& jg);
  void compute_eta_theta(amrex::MultiFab& etaT, amrex::MultiFab& etaTz, amrex::MultiFab& etaY, amrex::MultiFab& etaYz, 
			 amrex::MultiFab& eta1,
			 amrex::MultiFab& thetaT, amrex::MultiFab& thetaTz, amrex::MultiFab& thetaY, amrex::MultiFab& thetaYz, 
			 amrex::MultiFab& theta1,
			 amrex::MultiFab& djdT, amrex::MultiFab& djdY, 
			 const amrex::MultiFab& dkdT, const amrex::MultiFab& dkdY, 
			 const amrex::MultiFab& dedT, const amrex::MultiFab& dedY, 
			 const amrex::MultiFab& Er_star, const amrex::MultiFab& rho, 
			 const amrex::BoxArray& grids, amrex::Real delta_t, amrex::Real ptc_tau);
  void eos_opacity_emissivity(const amrex::MultiFab& S_new, 
			      const amrex::MultiFab& temp_new, const amrex::MultiFab& Ye_new,
			      const amrex::MultiFab& temp_star, const amrex::MultiFab& Ye_star,
			      amrex::MultiFab& kappa_p, amrex::MultiFab& kappa_r, amrex::MultiFab& jg, 
			      amrex::MultiFab& djdT, amrex::MultiFab& djdY, 
			      amrex::MultiFab& dkdT, amrex::MultiFab& dkdY, 
			      amrex::MultiFab& dedT, amrex::MultiFab& dedY, 
			      int level, const amrex::BoxArray& grids, int it, int ngrow); 
  void gray_accel(amrex::MultiFab& Er_new, amrex::MultiFab& Er_pi, 
		  amrex::MultiFab& kappa_p, amrex::MultiFab& kappa_r,
		  amrex::MultiFab& etaT, amrex::MultiFab& etaY, amrex::MultiFab& eta1,
		  amrex::MultiFab& thetaT, amrex::MultiFab& thetaY,
		  amrex::MultiFab& mugT, amrex::MultiFab& mugY, 
		  amrex::Tuple<amrex::MultiFab, BL_SPACEDIM>& lambda,
		  RadSolve& solver, MGRadBndry& mgbd, 
		  const amrex::BoxArray& grids, int level, amrex::Real time, amrex::Real delta_t, amrex::Real ptc_tau);
  void local_accel(amrex::MultiFab& Er_new, const amrex::MultiFab& Er_pi, 
		   const amrex::MultiFab& kappa_p, 
		   const amrex::MultiFab& etaT, const amrex::MultiFab& etaY, 
		   const amrex::MultiFab& thetaT, const amrex::MultiFab& thetaY, 
		   const amrex::MultiFab& mugT, const amrex::MultiFab& mugY, 
		   amrex::Real delta_t, amrex::Real ptc_tau);
  void state_energy_update(amrex::MultiFab& state, const amrex::MultiFab& rhoe, const amrex::MultiFab& Ye,
			   const amrex::MultiFab& temp, const amrex::BoxArray& grids, 
			   amrex::Real& derat, amrex::Real& dT, amrex::Real&dye, int level);
  void update_matter(amrex::MultiFab& rhoe_new, amrex::MultiFab& temp_new, 
		     amrex::MultiFab& rhoYe_new, amrex::MultiFab& Ye_new, 
		     const amrex::MultiFab& Er_new, const amrex::MultiFab& Er_pi,
		     const amrex::MultiFab& rhoe_star,	const amrex::MultiFab& rhoYe_star, 
		     const amrex::MultiFab& rhoe_step, const amrex::MultiFab& rhoYe_step,
		     const amrex::MultiFab& etaT, const amrex::MultiFab& etaTz, 
		     const amrex::MultiFab& etaY, const amrex::MultiFab& etaYz, 
		     const amrex::MultiFab& eta1,
		     const amrex::MultiFab& thetaT, const amrex::MultiFab& thetaTz, 
		     const amrex::MultiFab& thetaY, const amrex::MultiFab& thetaYz, 
		     const amrex::MultiFab& theta1,
		     const amrex::MultiFab& coupT, const amrex::MultiFab& coupY, 
		     const amrex::MultiFab& kappa_p, const amrex::MultiFab& jg, 
		     const amrex::MultiFab& mugT, const amrex::MultiFab& mugY, 
		     const amrex::MultiFab& S_new, 
		     int level, amrex::Real delta_t, 
		     amrex::Real ptc_tau, int it, bool conservative_update);
  void bisect_matter(amrex::MultiFab& rhoe_new, amrex::MultiFab& temp_new, 
		     amrex::MultiFab& rhoYe_new, amrex::MultiFab& Ye_new, 
		     const amrex::MultiFab& rhoe_star, const amrex::MultiFab& temp_star, 
		     const amrex::MultiFab& rhoYe_star, const amrex::MultiFab& Ye_star, 
		     const amrex::MultiFab& S_new, const amrex::BoxArray& grids, int level);

  // for the hyperbolic solver
  void MGFLD_compute_rosseland(amrex::MultiFab& kappa_r, const amrex::MultiFab& state);
  void MGFLD_compute_rosseland(amrex::FArrayBox& kappa_r, const amrex::FArrayBox& state);

  void MGFLD_compute_scattering(amrex::FArrayBox& kappa_s, const amrex::FArrayBox& state);

  void rhstoEr(amrex::MultiFab& rhs, amrex::Real dt, int level);

  // </ MGFLD>

  void save_lambda_in_plotvar(int level, const amrex::Tuple<amrex::MultiFab,BL_SPACEDIM>& lambda);
  void save_lab_Er_in_plotvar(int level, const amrex::MultiFab& Snew,
			      const amrex::MultiFab& Ecom, const amrex::MultiFab& F, int iflx);
  void save_lab_flux_in_plotvar(int level, const amrex::MultiFab& Snew, 
				const amrex::Tuple<amrex::MultiFab,BL_SPACEDIM>& lambda,
				const amrex::MultiFab& Er, const amrex::MultiFab& F, int iflx);
  void save_com_flux_in_plotvar(int level, const amrex::MultiFab& Snew, 
				const amrex::Tuple<amrex::MultiFab,BL_SPACEDIM>& lambda,
				const amrex::MultiFab& Er, const amrex::MultiFab& F, int iflx);

protected:

  amrex::Amr* parent;
  amrex::BCRec rad_bc;          // types defined in LO_BCTYPES.H, not BCTYPES.H
  amrex::Real reltol, abstol;   // tolerances for implicit update loop
  int maxiter;           // iteration limit for implicit update loop
  int miniter;
  int convergence_check_type; // 0 (default): both 2 and 3, 
                              // 1: rhoe,  
                              // 2: residue of Eq. rhoe,  
                              // 3: T & Ye
  amrex::Real relInTol, absInTol; // tolerance for inner iternation of J equation
  int maxInIter;           // iteration limit for inner iteration of J equation
  int minInIter;
  int skipAccelAllowed;   // Skip acceleration if it doesn't help

  int matter_update_type; // 0: conservative  1: non-conservative  2: C and NC interwoven
                          // The last outer iteration is always conservative.

  int n_bisect;  // Bisection after n_bisect iterations
  amrex::Real dedT_fac, dedY_fac; // Make dedT and dedY larger for safety in Newton iteration

  int inner_convergence_check; 

  amrex::Real delta_e_rat_dt_tol; // energy change tolerance for adjusting timestep
  amrex::Real delta_T_rat_dt_tol; // T change tolerance for adjusting timestep
  amrex::Real delta_Ye_dt_tol;    // Ye change tolerance for adjusting timestep

  int limiter;           // 0: no limiter, 2: Lev-Pom limiter
                         // 12: Bruenn, 22: square root, 32: Minerbo
  int closure;           // 0: f = lambda, 1: f = 1/3, 
                         // 2: f = 1-2*lambda, 3: f = lambda+(lambda*R)^2
                         // 4: f = 1/3 + 2/3*(lambda*R)^2
  int update_planck;     // after this number of iterations, lag planck
  int update_rosseland;  // after this number of iterations, lag rosseland
  int update_opacity;
  int update_limiter;    // after this number of iterations, lag limiter
  int inner_update_limiter; // This is for MGFLD solver. 
                            // Stop updating limiter after ? inner iterations
                            // 0 means lagging by one outer iteration
  amrex::Real dT;               // temperature step for derivative estimate
  int surface_average;   // 0 = arithmetic, 1 = harmonic, 2 = surface formula

  amrex::Real underfac;         // factor controlling progressive underrelaxation

  int do_sync;               // perform sync (if false zero out sync source)

  amrex::Real kappa_r_floor, temp_floor;

  int use_opacity_table_module;  // Use opacity_table_module?

  int do_kappa_stm_emission;

  amrex::Array<amrex::Real> delta_e_rat_level, delta_T_rat_level, delta_Ye_level;

  // arrays of flux registers 
  // flux_cons[level] is a amrex::FluxRegister, which contains fluxes for all groups
  amrex::Array <std::unique_ptr<amrex::FluxRegister> > flux_cons;
  amrex::Array <std::unique_ptr<amrex::FluxRegister> > flux_cons_old;
  amrex::Array <std::unique_ptr<amrex::FluxRegister> > flux_trial;

  // for deferred sync
  amrex::Array<amrex::Real> delta_t_old;

  // divergence of flux
  amrex::Array<std::unique_ptr<amrex::MultiFab> > dflux;

  amrex::Array<amrex::Real> xnu, nugroup, dnugroup;
  std::string group_units;
  amrex::Real group_print_factor;

  amrex::IntVect spot;

  int Density, Xmom, Eden, Eint, Temp, FirstSpec, FirstAux, NUM_STATE;

  void set_current_group(int igroup);

public:
  int verbose;
};

#endif
