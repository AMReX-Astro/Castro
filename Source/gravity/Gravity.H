#ifndef _Gravity_H_
#define _Gravity_H_

#include <AMReX_AmrLevel.H>
#include <AMReX_MLLinOp.H>

///
/// @class Gravity 
/// @brief
///
class Gravity {

public:

///
/// @param Parent
/// @param _finest_level
/// @param _phys_bc
/// @param _Density
///
  Gravity(amrex::Amr* Parent, int _finest_level, amrex::BCRec* _phys_bc, int _Density);
  ~Gravity ();
  void read_params ();


///
/// @param jobInfoFile
///
  void output_job_info_params(std::ostream& jobInfoFile);


///
/// @param numpts
///
  void set_numpts_in_gravity(int numpts);

  void install_level (int                   level,
                      amrex::AmrLevel*             level_data,
                      amrex::MultiFab&             volume,
                      amrex::MultiFab*             area);

  std::string get_gravity_type ();
  int get_max_solve_level();
  int NoSync();
  int NoComposite();
  int DoCompositeCorrection();
  int test_results_of_solves ();


///
/// @param time
/// @param multi_level
///
  void set_mass_offset(amrex::Real time, bool multi_level=true);

  amrex::Vector<std::unique_ptr<amrex::MultiFab> >& get_grad_phi_prev(int level);
  amrex::Vector<std::unique_ptr<amrex::MultiFab> >& get_grad_phi_curr(int level);


///
/// @param level
/// @param comp
///
  amrex::MultiFab* get_grad_phi_prev_comp(int level,int comp);


///
/// @param level
/// @param addend
///
  void plus_grad_phi_curr(int level, amrex::Vector<std::unique_ptr<amrex::MultiFab> >& addend);


///
/// @param level
///
  void swapTimeLevels (int level);

///
/// Calculate the maximum value of the RHS over all levels.
/// This should only be called at a synchronization point where
/// all Castro levels have valid new time data at the same simulation time.
/// The RHS we will use is the density multiplied by 4*pi*G and also
/// multiplied by the metric terms, just as it would be in a real solve.
///
  void update_max_rhs();

///
/// @param level
/// @param phi
/// @param grad_phi
/// @param is_new
///
  void solve_for_phi (int               level,
		      amrex::MultiFab&         phi,
                      const amrex::Vector<amrex::MultiFab*>& grad_phi,
		      int               is_new);

///
/// @param crse_level
/// @param fine_level
/// @param rhs
/// @param delta_phi
/// @param grad_delta_phi
///
  void solve_for_delta_phi (int                        crse_level,
                            int                        fine_level,
                            const amrex::Vector<amrex::MultiFab*>&          rhs,
                            const amrex::Vector<amrex::MultiFab*>&          delta_phi,
                            const amrex::Vector<amrex::Vector<amrex::MultiFab*> >&  grad_delta_phi);

///
/// @param crse_level
/// @param fine_level
/// @param drho
/// @param dphi
///
  void gravity_sync (int crse_level, int fine_level,
		     const amrex::Vector<amrex::MultiFab*>& drho, const amrex::Vector<amrex::MultiFab*>& dphi);


///
/// @param level
/// @param finest_level
/// @param use_previous_phi_as_guess
///
  void multilevel_solve_for_new_phi (int level, int finest_level,
                                     int use_previous_phi_as_guess = 0);

///
/// @param level
/// @param finest_level
/// @param grad_phi
/// @param is_new
/// @param use_previous_phi_as_guess
///
  void actual_multilevel_solve      (int level, int finest_level,
                                     const amrex::Vector<amrex::Vector<amrex::MultiFab*> >& grad_phi,
				     int is_new,
                                     int use_previous_phi_as_guess = 0);


///
/// @param level
/// @param mf
/// @param mask
///
  amrex::Real computeAvg    (int level, amrex::MultiFab* mf, bool mask=true);


///
/// @param level
/// @param comp_phi
/// @param comp_gphi
/// @param cml_phi
/// @param cml_gphi
///
  void create_comp_minus_level_grad_phi(int level,
					amrex::MultiFab& comp_phi,
					const amrex::Vector<amrex::MultiFab*>& comp_gphi,
					amrex::MultiFab& cml_phi,
					amrex::Vector<std::unique_ptr<amrex::MultiFab> >& cml_gphi);


///
/// @param level
/// @param phi_crse
/// @param time
///
  void GetCrsePhi(int level,
                  amrex::MultiFab& phi_crse,
                  amrex::Real time);


///
/// @param level
/// @param grav_vector
/// @param time
///
  void get_old_grav_vector(int level, amrex::MultiFab& grav_vector, amrex::Real time);

///
/// @param level
/// @param grav_vector
/// @param time
///
  void get_new_grav_vector(int level, amrex::MultiFab& grav_vector, amrex::Real time);


///
/// @param level
///
  void test_level_grad_phi_prev(int level);

///
/// @param level
///
  void test_level_grad_phi_curr(int level);

///
/// @param level
///
  void test_composite_phi(int level);


///
/// @param level
/// @param is_new
///
  void average_fine_ec_onto_crse_ec(int level, int is_new);


///
/// @param level
/// @param time
/// @param radial_grav
///
  void make_radial_gravity(int level, amrex::Real time, amrex::Vector<amrex::Real>& radial_grav);

///
/// @param level
/// @param radial_grav
/// @param grav_vector
///
  void interpolate_monopole_grav(int level, amrex::Vector<amrex::Real>& radial_grav, amrex::MultiFab& grav_vector);


///
/// @param level
/// @param time
/// @param grav
/// @param phi
///
  void make_prescribed_grav(int level, amrex::Real time, amrex::MultiFab& grav, amrex::MultiFab& phi);

///
/// @param level
/// @param Rhs
/// @param phi
/// @param fill_interior
///
  void make_radial_phi(int level, const amrex::MultiFab& Rhs, amrex::MultiFab& phi, int fill_interior);

#if (BL_SPACEDIM > 1)
///
/// @param crse_level
/// @param fine_level
/// @param Rhs
/// @param phi
///
  void fill_multipole_BCs(int crse_level, int fine_level, const amrex::Vector<amrex::MultiFab*>& Rhs, amrex::MultiFab& phi);
  void init_multipole_grav();
#endif
#if (BL_SPACEDIM == 3)

///
/// @param crse_level
/// @param fine_level
/// @param Rhs
/// @param phi
///
  void fill_direct_sum_BCs(int crse_level, int fine_level, const amrex::Vector<amrex::MultiFab*>& Rhs, amrex::MultiFab& phi);
#endif

  void make_mg_bc();

protected:

///
/// Pointers to amr,amrlevel.
///
  amrex::Amr*             parent;
  amrex::Vector<amrex::AmrLevel*> LevelData;

///
/// Pointers to grad_phi at previous and current time
///
  amrex::Vector< amrex::Vector<std::unique_ptr<amrex::MultiFab> > > grad_phi_curr;
  amrex::Vector< amrex::Vector<std::unique_ptr<amrex::MultiFab> > > grad_phi_prev;


///
/// BoxArray at each level
///
  const amrex::Vector<amrex::BoxArray>& grids;
  const amrex::Vector<amrex::DistributionMapping>& dmap;

///
/// Absolute tolerance on each level
///
  amrex::Vector<amrex::Real> abs_tol;

///
/// Relative tolerance on each level
///
  amrex::Vector<amrex::Real> rel_tol;

///
/// Resnorm at each level
///
  amrex::Vector<amrex::Real> level_solver_resnorm;

///
/// Maximum value of the RHS (used for obtaining absolute tolerances)
///
  amrex::Real max_rhs;

///
/// Volume and area fractions.
///
  amrex::Vector<amrex::MultiFab*> volume;
  amrex::Vector<amrex::MultiFab*> area;

  int Density;
  int finest_level;
  int finest_level_allocated;

  amrex::BCRec*       phys_bc;

  std::array<amrex::MLLinOp::BCType,AMREX_SPACEDIM> mlmg_lobc;
  std::array<amrex::MLLinOp::BCType,AMREX_SPACEDIM> mlmg_hibc;

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

///
/// @param level
/// @param Rhs
/// @param coeffs
///
  void applyMetricTerms(int level,amrex::MultiFab& Rhs, const amrex::Vector<amrex::MultiFab*>& coeffs);

///
/// @param level
/// @param cc
///
  void unweight_cc(int level,amrex::MultiFab& cc);

///
/// @param level
/// @param edges
///
  void unweight_edges(int level, const amrex::Vector<amrex::MultiFab*>& edges);
#endif

#ifdef POINTMASS

///
/// @param level
/// @param phi
/// @param grav_vector
/// @param point_mass
///
    void add_pointmass_to_gravity (int level, amrex::MultiFab& phi, amrex::MultiFab& grav_vector, amrex::Real point_mass);
#endif

private:

///
/// @param crse_level
/// @param nlevs
/// @param is_new
///
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > get_rhs (int crse_level, int nlevs, int is_new);


///
/// @param level
///
    void sanity_check (int level);

///
/// @param crse_level
/// @param fine_level
/// @param phi
/// @param rhs
/// @param amrex::Vector<std::array<amrex::MultiFab*
/// @param grad_phi
/// @param res
/// @param crse_bcdata
/// @param rel_eps
/// @param abs_eps
///
    amrex::Real actual_solve_with_mlmg (int crse_level, int fine_level,
                                        const amrex::Vector<amrex::MultiFab*>& phi,
                                        const amrex::Vector<const amrex::MultiFab*>& rhs,
                                        const amrex::Vector<std::array<amrex::MultiFab*,AMREX_SPACEDIM> >& grad_phi,
                                        const amrex::Vector<amrex::MultiFab*>& res,
                                        const amrex::MultiFab* const crse_bcdata,
                                        amrex::Real rel_eps, amrex::Real abs_eps);


///
/// @param crse_level
/// @param fine_level
/// @param phi
/// @param rhs
/// @param grad_phi
/// @param res
/// @param time
///

    amrex::Real solve_phi_with_mlmg (int crse_level, int fine_level,
                                     const amrex::Vector<amrex::MultiFab*>& phi,
                                     const amrex::Vector<amrex::MultiFab*>& rhs,
                                     const amrex::Vector<amrex::Vector<amrex::MultiFab*> >& grad_phi,
                                     const amrex::Vector<amrex::MultiFab*>& res,
                                     amrex::Real time);

};

///
/// @class GradPhiPhysBCFunct
/// @brief
///
class GradPhiPhysBCFunct
    : public amrex::PhysBCFunctBase
{
public:
    GradPhiPhysBCFunct ();

///
/// @param mf
/// @param dcomp
/// @param ncomp
/// @param time
/// @param bccomp
///
    virtual void FillBoundary (amrex::MultiFab& mf, int dcomp, int ncomp, amrex::Real time, int bccomp);
};
#endif
