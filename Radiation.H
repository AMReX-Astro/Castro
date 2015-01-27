#ifndef _Radiation_H_
#define _Radiation_H_

#include <AmrLevel.H>

#include "LHH.H"

#include "Castro.H"
#include "RadTypes.H"
#include "RadBndry.H"
#include "MGRadBndry.H"
#include "RadSolve.H"
#include <FluxRegister.H>
#include <Tuple.H>

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

  static Real radtoE, Etorad, radfluxtoF;

  static int do_multigroup;
  static int nGroups;
  static int nNeutrinoSpecies;
  static Array<int> nNeutrinoGroups;
  static int plot_neutrino_group_energies_per_MeV;
  static int plot_neutrino_group_energies_total;

  static int rad_hydro_combined;
  static int comoving;
  static int Er_Lorentz_term; 
  static int fspace_advection_type;

  static int Test_Type_lambda, Test_Type_Flux, Test_Type_Flux_lab;

  static int filter_lambda_T, filter_lambda_S;
  static int filter_prim_int, filter_prim_T, filter_prim_S;

  static int accelerate;        // controls multigroup convergence acceleration

  static void read_static_params();

  // Physical constants for use in multigroup, made static and public
  // so that other modules such as RadSolve can get to them easily.
  // Values are set in the Radiation constructor.

  static Real convert_MeV_erg, clight, hPlanck, kBoltz, Avogadro;

  // Values used internally for speed of light, Stefan-Boltzmann constant,
  // radiation constant.  (For contrived test problems, these may differ
  // from the proper physical values, so be careful if you use them!)

  static Real c, sigma, aRad;

  static int current_group_number;
  static std::string current_group_name;

  static Real flatten_pp_threshold;
  static int first_order_hydro;
  static int pure_hydro;

  int do_real_eos;  // if 0,  c_v = const_c_v * (rho**m) * (temp**(-n))

  Radiation(Amr* Parent, class Castro* castro, int restart = 0);
  ~Radiation() { }

  void regrid(int level, const BoxArray& grids);

  void close(int level);

  void restart(int level,
               const std::string& dir,
               std::istream& is);

  void checkPoint(int level,
                  const std::string& dir,
                  std::ostream&  os,
                  VisMF::How     how);

  // access to group information:

  Real group_center(int i) {
    return nugroup[i];
  }

  Real group_weight(int i) {
    return dnugroup[i];
  }

  // diagnostics for adjusting timestep:

  Real deltaEnergyRatMax(int lev) {
    return delta_e_rat_level[lev];
  }

  Real deltaTRatMax(int lev) {
    return delta_T_rat_level[lev];
  }

  Real deltaYeMax(int lev) {
    return delta_Ye_level[lev];
  }

  Real deltaEnergyTol() {
    return delta_e_rat_dt_tol;
  }

  Real deltaTTol() {
    return delta_T_rat_dt_tol;
  }

  Real deltaYeTol() {
    return delta_Ye_dt_tol;
  }

  // access to conservation flux register for energy diagnostics:

  FluxRegister* consRegister(int level) {
    return (flux_cons_old.defined(level) ? &flux_cons_old[level] : 0);
  }

  void analytic_solution(int level);

  void post_init(int level);

  void pre_timestep(int level);

  void init_flux(int level, int ncycle);

  void deferred_sync_setup(int level);

  void deferred_sync(int level, MultiFab& rhs, int indx);

  void single_group_update(int level, int iteration, int ncycle);

  void MGFLD_implicit_update(int level, int iteration, int ncycle);

  void write_groups(ostream& os);

  void reflux(int level);

  void compute_limiter(int level, const BoxArray& grids,
		       const MultiFab &Sborder, 
		       const MultiFab &Erborder,
		       MultiFab &lmaborder);

  void estimate_gamrPr(const FArrayBox& state, const FArrayBox& Er, 
		       FArrayBox& gPr, const Real* dx, const Box& box);

  void EstTimeStep(Real& estdt, int level);

  void computeTemp(MultiFab& State, int resetEint);

  void filter_prim(int level, MultiFab& State);

protected:

  // internal utility functions:

  void clear_internal_borders(FluxRegister& fr);

  void compute_exchange(MultiFab& exch, MultiFab& Er, MultiFab& fkp,
			int igroup = -1);

  void compute_eta(MultiFab& eta, MultiFab& etainv, MultiFab& state,
		   MultiFab& temp, MultiFab& fkp, MultiFab& Er,
		   Real delta_t, Real c, Real underrel, int lag_planck,
		   int igroup = 0);

  void internal_energy_update(Real& relative, Real& absolute,
			      MultiFab& frhoes, MultiFab& frhoem,
			      MultiFab& eta, MultiFab& etainv,
			      MultiFab& dflux_old, MultiFab& dflux_new,
                              MultiFab& exchange, Real delta_t);

  void internal_energy_update(Real& relative, Real& absolute,
			      MultiFab& frhoes, MultiFab& frhoem,
			      MultiFab& eta, MultiFab& etainv,
			      MultiFab& dflux_old, MultiFab& dflux_new,
			      MultiFab& exchange, MultiFab& Dterm,
			      Real delta_t);

  void nonconservative_energy_update(Real& relative, Real& absolute,
				     MultiFab& frhoes, MultiFab& frhoem,
				     MultiFab& eta, MultiFab& etainv,
				     MultiFab& Er_new,
				     MultiFab& dflux_old, MultiFab& dflux_new,
				     MultiFab& temp, MultiFab& fkp,
				     MultiFab& state, Real delta_t);

  void state_update(MultiFab& state,  MultiFab& frhoes, MultiFab& temp);

  void getBndryData(RadBndry& bd, MultiFab& Er, Real time, int level);

  // multigroup version
  void getBndryDataMG(MGRadBndry& mgbd, MultiFab& Er, Real time, int level);
  void getBndryDataMG_ga(MGRadBndry& mgbd, MultiFab& Er, int level);

  void filBndry(BndryRegister& bdry, int level, Real time);

  // Flux limiter functions, potentially for use by all update modules

  // Computes the scaled gradient for use in flux limiters & Eddington factors

  void scaledGradient(int level,
                      Tuple<MultiFab, BL_SPACEDIM>& R,
                      MultiFab& kappa_r, int kcomp,
                      MultiFab& Er, int igroup,
                      int limiter, int nGrow_Er=0, int Rcomp=0);

  // On input, lambda should contain scaled gradient.
  // On output this will be overwritten with the flux limiter.

  void fluxLimiter(int level,
                   Tuple<MultiFab, BL_SPACEDIM>& lambda,
                   int limiter, int lamcomp=0);

  // Fab versions of conversion functions.  All except frhoe use eos data.

  void get_frhoe(Fab& rhoe, Fab& state, const Box& reg);

  void get_c_v(Fab& c_v, Fab& temp, Fab& state, const Box& reg);

  // temp contains frhoe on input:

  void get_planck_and_temp(Fab& fkp, Fab& temp,
			   Fab& state, const Box& reg,
			   int igroup = 0, Real delta_t = 0.0);

  void get_rosseland_and_temp(Fab& kappa_r, Fab& temp,
			      Fab& state, const Box& reg,
			      int igroup = 0);

  // temp contains temp on input:
  void get_planck_from_temp(Fab& fkp, Fab& temp,
			    Fab& state, const Box& reg,
			    int igroup = 0);

  void get_rosseland_from_temp(Fab& kappa_r, Fab& temp,
			       Fab& state, const Box& reg,
			       int igroup = 0);

  // MultiFab versions of conversion functions call Fab versions for eos data:

  void get_frhoe(MultiFab& frhoe, MultiFab& state);

  // Planck mean is only needed in the grid interiors, but must be
  // computed separately for each fluid (fkp.nVar() == N fluids).
  // temp contains frhoe on input:

  void get_planck_and_temp(MultiFab& fkp, MultiFab& temp, MultiFab& state,
			   int igroup = 0, Real delta_t = 0.0);

  // Rosseland mean is needed on an extended region (nGrow == 1).
  // Fills patches with state information on extended region:

  void get_rosseland(MultiFab& kappa_r, AmrLevel* hclaw, int igroup = 0);

  // For the interior opacities, temperature is assumed to be known (input):

  // Updates Rosseland mean in interior of level, leaves ghost values
  // at physical or coarse-fine interfaces untouched:

  void update_rosseland_from_temp(MultiFab& kappa_r,
				  MultiFab& temp, MultiFab& state,
				  const Geometry& geom,
				  int igroup = 0);

  // RadHydro
  void get_rosseland_v_dcf(MultiFab& kappa_r, MultiFab& v, MultiFab& dcf, 
			   Real delta_t, Real c,
			   AmrLevel* hclaw,
			   int igroup = 0);  

  // multigroup methods
  void get_groups(int verbose);

  void extrapolateBorders(MultiFab& f, int indx);

  void SGFLD_compute_rosseland(MultiFab& kappa_r, const MultiFab& state);
  void SGFLD_compute_rosseland(FArrayBox& kappa_r, const FArrayBox& state);
  void update_dcf(MultiFab& dcf, MultiFab& etainv, MultiFab& kp, MultiFab& kr,
		  const Geometry& geom);

  // <MGFLD>

  int use_dkdT;

  // <Shestakov-Bolstad>
  int use_WiensLaw;
  Real Tf_Wien;
  // </ Shestakov-Bolstad>
  int integrate_Planck;

  void check_convergence_er(Real& relative_in, Real& absolute_in, Real& error_er,
  			    const MultiFab& Er_new, const MultiFab& Er_pi, 
			    const MultiFab& kappa_p, 
			    const MultiFab& etaTz, const MultiFab& etaYz, 
			    const MultiFab& theTz, const MultiFab& theYz,
			    const MultiFab& temp_new, const MultiFab& Ye_new,
			    const BoxArray& grids, Real delta_t);
  void check_convergence_matt(const MultiFab& rhoe_new, const MultiFab& rhoe_star, 
			      const MultiFab& rhoe_step, const MultiFab& Er_new, 
			      const MultiFab& temp_new, const MultiFab& temp_star, 
			      const MultiFab& rhoYe_new, const MultiFab& rhoYe_star, 
			      const MultiFab& rhoYe_step, const MultiFab& rho,
			      const MultiFab& kappa_p, const MultiFab& jg,
			      const MultiFab& dedT, const MultiFab& dedY,
			      Real& rel_rhoe, Real& abs_rhoe, 
			      Real& rel_FT,   Real& abs_FT, 
			      Real& rel_T,    Real& abs_T, 
			      Real& rel_FY,   Real& abs_FY, 
			      Real& rel_Ye,   Real& abs_Ye,
			      const BoxArray& grids, Real delta_t);
  void compute_coupling(MultiFab& coupT, MultiFab& coupY, 
			const MultiFab& kappa_p, const MultiFab& Er_pi,
			const MultiFab& jg);
  void compute_eta_theta(MultiFab& etaT, MultiFab& etaTz, MultiFab& etaY, MultiFab& etaYz, 
			 MultiFab& eta1,
			 MultiFab& thetaT, MultiFab& thetaTz, MultiFab& thetaY, MultiFab& thetaYz, 
			 MultiFab& theta1,
			 MultiFab& djdT, MultiFab& djdY, 
			 const MultiFab& dkdT, const MultiFab& dkdY, 
			 const MultiFab& dedT, const MultiFab& dedY, 
			 const MultiFab& Er_star, const MultiFab& rho, 
			 const BoxArray& grids, Real delta_t, Real ptc_tau);
  void eos_opacity_emissivity(const MultiFab& S_new, 
			      const MultiFab& temp_new, const MultiFab& Ye_new,
			      const MultiFab& temp_star, const MultiFab& Ye_star,
			      MultiFab& kappa_p, MultiFab& kappa_r, MultiFab& jg, 
			      MultiFab& djdT, MultiFab& djdY, 
			      MultiFab& dkdT, MultiFab& dkdY, 
			      MultiFab& dedT, MultiFab& dedY, 
			      int level, const BoxArray& grids, int it, int ngrow); 
  void gray_accel(MultiFab& Er_new, MultiFab& Er_pi, 
		  MultiFab& kappa_p, MultiFab& kappa_r,
		  MultiFab& etaT, MultiFab& etaY, MultiFab& eta1,
		  MultiFab& thetaT, MultiFab& thetaY,
		  MultiFab& mugT, MultiFab& mugY, 
		  Tuple<MultiFab, BL_SPACEDIM>& lambda,
		  RadSolve& solver, MGRadBndry& mgbd, 
		  const BoxArray& grids, int level, Real time, Real delta_t, Real ptc_tau);
  void local_accel(MultiFab& Er_new, const MultiFab& Er_pi, 
		   const MultiFab& kappa_p, 
		   const MultiFab& etaT, const MultiFab& etaY, 
		   const MultiFab& thetaT, const MultiFab& thetaY, 
		   const MultiFab& mugT, const MultiFab& mugY, 
		   const BoxArray& grids, Real delta_t, Real ptc_tau);
  void state_energy_update(MultiFab& state, const MultiFab& rhoe, const MultiFab& Ye,
			   const MultiFab& temp, const BoxArray& grids, 
			   Real& derat, Real& dT, Real&dye, int level);
  void update_matter(MultiFab& rhoe_new, MultiFab& temp_new, 
		     MultiFab& rhoYe_new, MultiFab& Ye_new, 
		     const MultiFab& Er_new, const MultiFab& Er_pi,
		     const MultiFab& rhoe_star,	const MultiFab& rhoYe_star, 
		     const MultiFab& rhoe_step, const MultiFab& rhoYe_step,
		     const MultiFab& etaT, const MultiFab& etaTz, 
		     const MultiFab& etaY, const MultiFab& etaYz, 
		     const MultiFab& eta1,
		     const MultiFab& thetaT, const MultiFab& thetaTz, 
		     const MultiFab& thetaY, const MultiFab& thetaYz, 
		     const MultiFab& theta1,
		     const MultiFab& coupT, const MultiFab& coupY, 
		     const MultiFab& kappa_p, const MultiFab& jg, 
		     const MultiFab& mugT, const MultiFab& mugY, 
		     const MultiFab& S_new, 
		     int level, const BoxArray& grids, Real delta_t, 
		     Real ptc_tau, int it, bool conservative_update);
  void bisect_matter(MultiFab& rhoe_new, MultiFab& temp_new, 
		     MultiFab& rhoYe_new, MultiFab& Ye_new, 
		     const MultiFab& rhoe_star, const MultiFab& temp_star, 
		     const MultiFab& rhoYe_star, const MultiFab& Ye_star, 
		     const MultiFab& S_new, const BoxArray& grids, int level);

  // for the hyperbolic solver
  void MGFLD_compute_rosseland(MultiFab& kappa_r, const MultiFab& state);
  void MGFLD_compute_rosseland(FArrayBox& kappa_r, const FArrayBox& state);

  void rhstoEr(MultiFab& rhs, const BoxArray& grids, Real dt, int level);

  // </ MGFLD>

protected:

  Amr* parent;
  BCRec rad_bc;          // types defined in LO_BCTYPES.H, not BCTYPES.H
  Real reltol, abstol;   // tolerances for implicit update loop
  int maxiter;           // iteration limit for implicit update loop
  int miniter;
  int convergence_check_type; // 0 (default): both 2 and 3, 
                              // 1: rhoe,  
                              // 2: residue of Eq. rhoe,  
                              // 3: T & Ye
  Real relInTol, absInTol; // tolerance for inner iternation of J equation
  int maxInIter;           // iteration limit for inner iteration of J equation
  int minInIter;
  int skipAccelAllowed;   // Skip acceleration if it doesn't help

  int matter_update_type; // 0: conservative  1: non-conservative  2: C and NC interwoven
                          // The last outer iteration is always conservative.

  int n_bisect;  // Bisection after n_bisect iterations
  Real dedT_fac, dedY_fac; // Make dedT and dedY larger for safety in Newton iteration

  int inner_convergence_check; 

  Real delta_e_rat_dt_tol; // energy change tolerance for adjusting timestep
  Real delta_T_rat_dt_tol; // T change tolerance for adjusting timestep
  Real delta_Ye_dt_tol;    // Ye change tolerance for adjusting timestep

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
  Real dT;               // temperature step for derivative estimate
  int surface_average;   // 0 = arithmetic, 1 = harmonic, 2 = surface formula

  Real underfac;         // factor controlling progressive underrelaxation

  int do_sync;               // perform sync (if false zero out sync source)

  int use_analytic_solution; // compare against an analytic solution

  Real kappa_r_floor, temp_floor;

  int use_opacity_table_module;  // Use opacity_table_module?

  int do_kappa_stm_emission;
  Array<Real> const_c_v, const_kappa_p, const_kappa_r, const_scattering;
  Array<Real> c_v_exp_m, c_v_exp_n;
  Array<Real> kappa_p_exp_m, kappa_p_exp_n;
  Array<Real> kappa_r_exp_m, kappa_r_exp_n;
  Array<Real> scattering_exp_m, scattering_exp_n;
  Array<Real> kappa_p_exp_p, kappa_r_exp_p, scattering_exp_p;
  Array<Real> prop_temp_floor;

  Array<Real> delta_e_rat_level, delta_T_rat_level, delta_Ye_level;

  // arrays of flux registers 
  // flux_cons[level] is a FluxRegister, which contains fluxes for all groups
  PArray <FluxRegister> flux_cons;
  PArray <FluxRegister> flux_cons_old;
  PArray <FluxRegister> flux_trial;

  // for deferred sync
  Array<Real> delta_t_old;

  // divergence of flux
  PArray<MultiFab> dflux;

  Array<Real> xnu, nugroup, dnugroup;
  std::string group_units;
  Real group_print_factor;

  IntVect spot;

  int Density, Xmom, Eden, Eint, Temp, FirstSpec, FirstAux, NUM_STATE;

  void set_current_group(int igroup);

public:
  int verbose;
};

#endif
