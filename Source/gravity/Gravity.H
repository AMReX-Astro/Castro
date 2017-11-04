#ifndef _Gravity_H_
#define _Gravity_H_

#include <AMReX_AmrLevel.H>

#ifdef CASTRO_MLMG
#include <AMReX_MLLinOp.H>
#endif

class Gravity {

public:

  Gravity(amrex::Amr* Parent, int _finest_level, amrex::BCRec* _phys_bc, int _Density);
  ~Gravity ();
  void read_params ();

  void set_numpts_in_gravity(int numpts);

  void install_level (int                   level,
                      amrex::AmrLevel*             level_data,
                      amrex::MultiFab&             volume,
                      amrex::MultiFab*             area);

  std::string get_gravity_type ();
  int NoSync();
  int NoComposite();
  int DoCompositeCorrection();
  int test_results_of_solves ();
  
  void set_mass_offset(amrex::Real time, bool multi_level=true);

  amrex::Vector<std::unique_ptr<amrex::MultiFab> >& get_grad_phi_prev(int level);
  amrex::Vector<std::unique_ptr<amrex::MultiFab> >& get_grad_phi_curr(int level);

  amrex::MultiFab* get_grad_phi_prev_comp(int level,int comp);

  void plus_grad_phi_curr(int level, amrex::Vector<std::unique_ptr<amrex::MultiFab> >& addend);

  void swapTimeLevels (int level);

  void update_max_rhs();

  void solve_for_phi (int               level,
		      amrex::MultiFab&         phi,
                      const amrex::Vector<amrex::MultiFab*>& grad_phi,
		      int               is_new);


  void solve_for_delta_phi (int                        crse_level, 
                            int                        fine_level,
                            const amrex::Vector<amrex::MultiFab*>&          rhs,
                            const amrex::Vector<amrex::MultiFab*>&          delta_phi,
                            const amrex::Vector<amrex::Vector<amrex::MultiFab*> >&  grad_delta_phi);

  void gravity_sync (int crse_level, int fine_level,
		     const amrex::Vector<amrex::MultiFab*>& drho, const amrex::Vector<amrex::MultiFab*>& dphi);

  void multilevel_solve_for_new_phi (int level, int finest_level,
                                     int use_previous_phi_as_guess = 0);
  void actual_multilevel_solve      (int level, int finest_level, 
                                     const amrex::Vector<amrex::Vector<amrex::MultiFab*> >& grad_phi,
				     int is_new,
                                     int use_previous_phi_as_guess = 0);

  amrex::Real computeAvg    (int level, amrex::MultiFab* mf, bool mask=true);

  void create_comp_minus_level_grad_phi(int level,
					amrex::MultiFab& comp_phi,
					const amrex::Vector<amrex::MultiFab*>& comp_gphi,
					amrex::MultiFab& cml_phi,
					amrex::Vector<std::unique_ptr<amrex::MultiFab> >& cml_gphi);

  void GetCrsePhi(int level, 
                  amrex::MultiFab& phi_crse,
                  amrex::Real time);

  void get_old_grav_vector(int level, amrex::MultiFab& grav_vector, amrex::Real time);
  void get_new_grav_vector(int level, amrex::MultiFab& grav_vector, amrex::Real time);

  void test_level_grad_phi_prev(int level);
  void test_level_grad_phi_curr(int level);
  void test_composite_phi(int level);

  void average_fine_ec_onto_crse_ec(int level, int is_new);

  void make_radial_gravity(int level, amrex::Real time, amrex::Vector<amrex::Real>& radial_grav);
  void interpolate_monopole_grav(int level, amrex::Vector<amrex::Real>& radial_grav, amrex::MultiFab& grav_vector);

  void make_prescribed_grav(int level, amrex::Real time, amrex::MultiFab& grav, amrex::MultiFab& phi);
  void make_radial_phi(int level, const amrex::MultiFab& Rhs, amrex::MultiFab& phi, int fill_interior);

#if (BL_SPACEDIM > 1)
  void fill_multipole_BCs(int crse_level, int fine_level, const amrex::Vector<amrex::MultiFab*>& Rhs, amrex::MultiFab& phi);
  void init_multipole_grav();
#endif
#if (BL_SPACEDIM == 3)
  void fill_direct_sum_BCs(int crse_level, int fine_level, const amrex::Vector<amrex::MultiFab*>& Rhs, amrex::MultiFab& phi);
#endif

  void make_mg_bc();

protected:
  //
  // Pointers to amr,amrlevel.
  //
  amrex::Amr*             parent;
  amrex::Vector<amrex::AmrLevel*> LevelData;
  //
  // Pointers to grad_phi at previous and current time
  //
  amrex::Vector< amrex::Vector<std::unique_ptr<amrex::MultiFab> > > grad_phi_curr;
  amrex::Vector< amrex::Vector<std::unique_ptr<amrex::MultiFab> > > grad_phi_prev;

  //
  // BoxArray at each level
  //
  const amrex::Vector<amrex::BoxArray>& grids;
  //
  const amrex::Vector<amrex::DistributionMapping>& dmap;
  //
  // Absolute tolerance on each level
  //
  amrex::Vector<amrex::Real> abs_tol;
  //
  // Relative tolerance on each level
  //
  amrex::Vector<amrex::Real> rel_tol;
  //
  // Resnorm at each level
  //
  amrex::Vector<amrex::Real> level_solver_resnorm;
  //
  // Maximum value of the RHS (used for obtaining absolute tolerances)
  //
  amrex::Real max_rhs;
  //
  // Volume and area fractions.
  //
  amrex::Vector<amrex::MultiFab*> volume;
  amrex::Vector<amrex::MultiFab*> area;

  int Density;
  int finest_level;
  int finest_level_allocated;

  int mg_bc[2*BL_SPACEDIM];

  amrex::BCRec*       phys_bc;

#ifdef CASTRO_MLMG
  std::array<amrex::MLLinOp::BCType,AMREX_SPACEDIM> mlmg_lobc;
  std::array<amrex::MLLinOp::BCType,AMREX_SPACEDIM> mlmg_hibc;
#endif

  int   numpts_at_level;

  static int   test_solves;
  static amrex::Real  mass_offset;
  static amrex::Vector< amrex::Vector<amrex::Real> > radial_grav_old;
  static amrex::Vector< amrex::Vector<amrex::Real> > radial_grav_new;
  static amrex::Vector< amrex::Vector<amrex::Real> > radial_mass;
  static amrex::Vector< amrex::Vector<amrex::Real> > radial_vol;
#ifdef GR_GRAV
  static amrex::Vector< amrex::Vector<amrex::Real> > radial_pres;
#endif
  static int   stencil_type;

  static amrex::Real max_radius_all_in_domain;

#include "gravity_params.H"

#if (BL_SPACEDIM < 3)
  void applyMetricTerms(int level,amrex::MultiFab& Rhs, const amrex::Vector<amrex::MultiFab*>& coeffs);
  void unweight_cc(int level,amrex::MultiFab& cc);
  void unweight_edges(int level, const amrex::Vector<amrex::MultiFab*>& edges);
#endif

#ifdef POINTMASS
    void add_pointmass_to_gravity (int level, amrex::MultiFab& phi, amrex::MultiFab& grav_vector, amrex::Real point_mass);
#endif

private:

    amrex::Real solve_phi_with_fmg (int crse_level, int fine_level,
			     const amrex::Vector<amrex::MultiFab*>& phi,
			     const amrex::Vector<amrex::MultiFab*>& rhs,
			     const amrex::Vector<amrex::Vector<amrex::MultiFab*> >& grad_phi,
			     const amrex::Vector<amrex::MultiFab*>& res,
			     amrex::Real time);

    void solve_for_delta_phi_with_fmg (int crse_level, int fine_level,
                                       const amrex::Vector<amrex::MultiFab*>& rhs,
                                       const amrex::Vector<amrex::MultiFab*>& delta_phi,
                                       const amrex::Vector<amrex::Vector<amrex::MultiFab*> >&  grad_delta_phi);

    amrex::Vector<std::unique_ptr<amrex::MultiFab> > get_rhs (int crse_level, int nlevs, int is_new);

    void sanity_check (int level);

#ifdef CASTRO_MLMG

    amrex::Real actual_solve_with_mlmg (int crse_level, int fine_level,
                                        const amrex::Vector<amrex::MultiFab*>& phi,
                                        const amrex::Vector<const amrex::MultiFab*>& rhs,
                                        const amrex::Vector<std::array<amrex::MultiFab*,AMREX_SPACEDIM> >& grad_phi,
                                        const amrex::Vector<amrex::MultiFab*>& res,
                                        const amrex::MultiFab* const crse_bcdata,
                                        amrex::Real rel_eps, amrex::Real abs_eps);

    amrex::Real solve_phi_with_mlmg (int crse_level, int fine_level,
                                     const amrex::Vector<amrex::MultiFab*>& phi,
                                     const amrex::Vector<amrex::MultiFab*>& rhs,
                                     const amrex::Vector<amrex::Vector<amrex::MultiFab*> >& grad_phi,
                                     const amrex::Vector<amrex::MultiFab*>& res,
                                     amrex::Real time);

    void solve_for_delta_phi_with_mlmg (int crse_level, int fine_level,
                                        const amrex::Vector<amrex::MultiFab*>& rhs,
                                        const amrex::Vector<amrex::MultiFab*>& delta_phi,
                                        const amrex::Vector<amrex::Vector<amrex::MultiFab*> >&  grad_delta_phi);
#endif
};

class GradPhiPhysBCFunct
    : public amrex::PhysBCFunctBase
{
public:
    GradPhiPhysBCFunct ();
    virtual void FillBoundary (amrex::MultiFab& mf, int dcomp, int ncomp, amrex::Real time);
};
#endif
