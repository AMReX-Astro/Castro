#ifndef _Gravity_H_
#define _Gravity_H_

#include <AmrLevel.H>
#include <FluxRegister.H>

#ifdef PARTICLES
#include <Particles.H>
#endif

class Gravity {

public:

  Gravity(Amr* Parent, int _finest_level, BCRec* _phys_bc, int _Density);
  virtual ~Gravity ();
  void read_params ();

  void set_numpts_in_gravity(int numpts);

  void install_level (int                   level,
                      AmrLevel*             level_data,
                      MultiFab&             volume,
                      MultiFab*             area);

  std::string get_gravity_type ();
  Real get_const_grav ();
  int NoSync();
  int NoComposite();
  int test_results_of_solves ();
  
  void set_mass_offset(Real time, bool multi_level=true);

  PArray<MultiFab>& get_grad_phi_prev(int level);
  PArray<MultiFab>& get_grad_phi_curr(int level);

  MultiFab* get_grad_phi_prev_comp(int level,int comp);

  void plus_grad_phi_curr(int level, PArray<MultiFab>& addend);

  void swapTimeLevels (int level);

  void solve_for_phi (int               level,
		      MultiFab&         phi,
                      PArray<MultiFab>& grad_phi,
		      int               is_new);


  void solve_for_delta_phi (int                        crse_level, 
                            int                        fine_level,
                            MultiFab&                  CrseRhs,
                            PArray<MultiFab>&          delta_phi,
                            PArray<PArray<MultiFab> >& grad_delta_phi);

  void gravity_sync (int crse_level, int fine_level, int iteration, int ncycle,
                     const MultiFab& drho_and_drhoU, const MultiFab& dphi,
                     PArray<MultiFab>& grad_delta_phi_cc);

  void multilevel_solve_for_new_phi (int level, int finest_level,
                                     int use_previous_phi_as_guess = 0);
  void actual_multilevel_solve      (int level, int finest_level, 
                                     Array<PArray<MultiFab> >& grad_phi,
				     int is_new,
                                     int use_previous_phi_as_guess = 0);

  Real computeAvg    (int level, MultiFab* mf, bool mask=true);

  void create_comp_minus_level_grad_phi(int level, MultiFab& cmlp,
                                        PArray<MultiFab>& cmlgp);

  void GetCrseGradPhi(int level, 
                      PArray<MultiFab>& grad_phi_crse,
                      Real time);

  void GetCrsePhi(int level, 
                  MultiFab& phi_crse,
                  Real time);

  void zeroPhiFluxReg (int level);

  void get_old_grav_vector(int level, MultiFab& grav_vector, Real time);
  void get_new_grav_vector(int level, MultiFab& grav_vector, Real time);

  void test_level_grad_phi_prev(int level);
  void test_level_grad_phi_curr(int level);
  void test_composite_phi(int level);

  static void avgDown(MultiFab& crse, const MultiFab& fine, const IntVect& ratio);

  void average_fine_ec_onto_crse_ec(int level, int is_new);

  void add_to_fluxes(int level, int iteration, int ncycle);

  void reflux_phi(int level, MultiFab& dphi);

#if (BL_SPACEDIM == 1)
  void make_one_d_grav(int level, Real time, MultiFab& grav_vector);
#else
  void make_radial_gravity(int level, Real time, Array<Real>& radial_grav);
  void interpolate_monopole_grav(int level, 
                                 Array<Real>& radial_grav, MultiFab& grav_vector);
  void make_prescribed_grav(int level, Real time, MultiFab& grav);
#endif

  void make_radial_phi(int level, MultiFab& Rhs, MultiFab& phi, int fill_interior);

#if (BL_SPACEDIM == 3)
  void fill_multipole_BCs(int level, MultiFab& Rhs, MultiFab& phi);
  void fill_direct_sum_BCs(int level, MultiFab& Rhs, MultiFab& phi);
#endif

  void make_mg_bc();

protected:
  //
  // Pointers to amr,amrlevel.
  //
  Amr*             parent;
  PArray<AmrLevel> LevelData;
  //
  // Pointers to grad_phi at previous and current time
  //
  Array< PArray<MultiFab> > grad_phi_curr;
  Array< PArray<MultiFab> > grad_phi_prev;

  PArray<FluxRegister> phi_flux_reg;
  //
  // BoxArray at each level
  //
  Array<BoxArray> grids;
  //
  // Resnorm at each level
  //
  Array<Real> level_solver_resnorm;
  //
  // Volume and area fractions.
  //
  PArray<MultiFab> volume;
  Array<MultiFab*> area;

  const BoxArray& getEdgeBoxArray (int lev, int dir) const { return area[lev][dir].boxArray(); }

  int Density;
  int finest_level;
  int finest_level_allocated;

  int mg_bc[2*BL_SPACEDIM];

  BCRec*       phys_bc;

  int   numpts_at_level;

  static int   verbose;
  static int   test_solves;
  static int   no_sync;
  static int   no_composite;
  static int   drdxfac;
  static int   lnum;
  static int   direct_sum_bcs;
  static Real  const_grav;
  static Real  mass_offset;
  static Real  sl_tol;
  static Real  ml_tol;
  static Real  delta_tol;
  static Array< Array<Real> > radial_grav_old;
  static Array< Array<Real> > radial_grav_new;
  static Array< Array<Real> > radial_mass;
  static Array< Array<Real> > radial_vol;
#ifdef GR_GRAV
  static Array< Array<Real> > radial_pres;
#endif
  static std::string gravity_type;
  static int   stencil_type;

  static Real max_radius_all_in_domain;

  void fill_ec_grow (int level, PArray<MultiFab>& ecF, const PArray<MultiFab>& ecC) const;

#if (BL_SPACEDIM < 3)
  void applyMetricTerms(int level,MultiFab& Rhs, PArray<MultiFab>& coeffs);
  void unweight_cc(int level,MultiFab& cc);
  void unweight_edges(int level, PArray<MultiFab>& edges);
#endif

#ifdef POINTMASS
    void add_pointmass_to_gravity (int level, MultiFab& grav_vector, Real point_mass);
#endif

#ifdef PARTICLES
    void AddParticlesToRhs(int level, MultiFab& Rhs, int ngrow);
    void AddParticlesToRhs(int base_level, int finest_level, PArray<MultiFab>& Rhs_particles);
#endif

private:

    Real solve_phi_with_fmg (int crse_level, int fine_level,
			     PArray<MultiFab> & phi,
			     PArray<MultiFab> & rhs,
			     Array< PArray<MultiFab> >& grad_phi,
			     PArray<MultiFab>& res,
			     Real time);

    void get_rhs (int crse_level, int nlevs, PArray<MultiFab> & rhs_p, int is_new);

    void sanity_check (int level);
};
#endif

