#ifndef GRAVITY_H
#define GRAVITY_H

#include <AMReX_AmrLevel.H>
#include <AMReX_MLLinOp.H>

#include <gravity_params.H>

// This vector can be accessed on the GPU.
using RealVector = amrex::Gpu::ManagedVector<amrex::Real>;

///
/// Multipole gravity data
///
namespace multipole {
    const int lnum_max = 30;

    extern AMREX_GPU_MANAGED amrex::Real volumeFactor;
    extern AMREX_GPU_MANAGED amrex::Real parityFactor;

    extern AMREX_GPU_MANAGED amrex::Real rmax;

    extern AMREX_GPU_MANAGED amrex::Array1D<bool, 0, 2> doSymmetricAddLo;
    extern AMREX_GPU_MANAGED amrex::Array1D<bool, 0, 2> doSymmetricAddHi;
    extern AMREX_GPU_MANAGED bool doSymmetricAdd;

    extern AMREX_GPU_MANAGED amrex::Array1D<bool, 0, 2> doReflectionLo;
    extern AMREX_GPU_MANAGED amrex::Array1D<bool, 0, 2> doReflectionHi;

    extern AMREX_GPU_MANAGED amrex::Array2D<amrex::Real, 0, lnum_max, 0, lnum_max> factArray;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 0, lnum_max> parity_q0;
    extern AMREX_GPU_MANAGED amrex::Array2D<amrex::Real, 0, lnum_max, 0, lnum_max> parity_qC_qS;
}

///
/// @class Gravity
/// @brief
///
class Gravity {

public:

///
/// Constructor
///
/// @param Parent           Parent ``Amr`` object
/// @param _finest_level    finest level
/// @param _phys_bc         ``BCRec`` object, physical boundary conditions
/// @param _density         index of density component
///
  Gravity(amrex::Amr* Parent, int _finest_level, amrex::BCRec* _phys_bc, int _density);

///
/// Deconstructor
///
  ~Gravity ();

///
/// Remove copy/move constructors/assignment operators.
///
  Gravity (const Gravity&) = delete;
  Gravity (Gravity&&) = delete;
  Gravity& operator= (const Gravity&) = delete;
  Gravity& operator= (Gravity&&) = delete;

///
/// Read gravity-related parameters from parameter file
///
  void read_params ();


///
/// @param jobInfoFile      std::ostream object
///
  static void output_job_info_params(std::ostream& jobInfoFile);


///
/// Set ``Gravity`` object's ``numpts_at_level`` variable.
///
/// @param numpts       number of points
///
  void set_numpts_in_gravity(int numpts);

///
/// Setup gravity at level ``level``.
///
/// @param level        integer, level number
/// @param level_data   ``AmrLevel`` object containing state data at that level
/// @param volume       MultiFab, volume elements of cells
/// @param area         MultiFab, area elements of cells
///
  void install_level (int                   level,
                      amrex::AmrLevel*             level_data,
                      amrex::MultiFab&             volume,
                      amrex::MultiFab*             area);

///
/// Returns ``gravity_type``
///
  static std::string get_gravity_type ();

///
/// Returns ``max_solve_level``
///
  static int get_max_solve_level();

///
/// Returns ``no_sync``
///
  static int NoSync();

///
/// Returns ``do_composite_phi_correction``
///
  static int DoCompositeCorrection();

///
/// Returns ``test_solves``
///
  static int test_results_of_solves ();


///
/// Set the ``mass_offset``
///
/// @param time             Current time
/// @param multi_level      boolean, do we iterate over all levels or mask off fine grids?
///
  void set_mass_offset(amrex::Real time, bool multi_level=true) const;

///
/// Return ``grad_phi_prev`` at given level
///
/// @param level        level index
///
  amrex::Vector<std::unique_ptr<amrex::MultiFab> >& get_grad_phi_prev(int level);

///
/// Return ``grad_phi_curr`` at given level
///
/// @param level        level index
///
  amrex::Vector<std::unique_ptr<amrex::MultiFab> >& get_grad_phi_curr(int level);


///
/// Return given component of ``grad_phi_prev`` at given level
///
/// @param level        level index
/// @param comp         component index
///
  amrex::MultiFab* get_grad_phi_prev_comp(int level,int comp);


///
/// Return ``grad_phi_curr`` at given level plus the given vector
///
/// @param level        level index
/// @param addend       Vector of MultiFabs to add to grad phi
///
  void plus_grad_phi_curr(int level, amrex::Vector<std::unique_ptr<amrex::MultiFab> >& addend);


///
/// Swap ``grad_phi_prev`` with ``grad_phi_curr`` at given level at set new ``grad_phi_curr`` to 1.e50.
///
/// @param level        level index
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
/// Solve Poisson's equation to find the gravitational potential
///
/// @param level        level index
/// @param phi          MultiFab to store gravitational potential in
/// @param grad_phi     Vector of MultiFabs, \f$ \nabla \Phi \f$
/// @param is_new       do we use state data at current time (1) or old time (0)?
///
  void solve_for_phi (int               level,
                      amrex::MultiFab&         phi,
                      const amrex::Vector<amrex::MultiFab*>& grad_phi,
                      int               is_new);

///
/// Find delta phi
///
/// @param crse_level       index of coarse level
/// @param fine_level       index of fine level
/// @param rhs              Vector of MultiFabs with right hand side source terms
/// @param delta_phi        Vector of MultiFabs delta phi will be saved to
/// @param grad_delta_phi   Vector of MultiFabs, gradient of delta phi
///
  void solve_for_delta_phi (int                        crse_level,
                            int                        fine_level,
                            const amrex::Vector<amrex::MultiFab*>&          rhs,
                            const amrex::Vector<amrex::MultiFab*>&          delta_phi,
                            const amrex::Vector<amrex::Vector<amrex::MultiFab*> >&  grad_delta_phi);

///
/// Sync gravity across levels
///
/// @param crse_level       index of coarse level
/// @param fine_level       index of fine level
/// @param drho
/// @param dphi
///
  void gravity_sync (int crse_level, int fine_level,
                     const amrex::Vector<amrex::MultiFab*>& drho, const amrex::Vector<amrex::MultiFab*>& dphi);


///
/// Multilevel solve for new phi from base level to finest level
///
/// @param level                        Base level index
/// @param finest_level                 Fine level index
///
  void multilevel_solve_for_new_phi (int level, int finest_level);

///
/// Actually do the multilevel solve for new phi from base level to finest level
///
/// @param level            Base level index
/// @param finest_level     Fine level index
/// @param grad_phi         gradient of phi
/// @param is_new           Should we use the new state (1) or previous state (0)?
///
  void actual_multilevel_solve      (int level, int finest_level,
                                     const amrex::Vector<amrex::Vector<amrex::MultiFab*> >& grad_phi,
                                     int is_new);


///
/// Compute the difference between level and composite solves
///
/// @param level                       level index
/// @param comp_phi                    MultiFab containing computed phi
/// @param comp_gphi                   Vector of MultiFabs containing computed grad phi
/// @param comp_minus_level_phi        MultiFab, computed minus level phi
/// @param comp_minus_level_grad_phi   Vector of MultiFabs, computed minus level grad phi
///
  void create_comp_minus_level_grad_phi(int level,
                                        amrex::MultiFab& comp_phi,
                                        const amrex::Vector<amrex::MultiFab*>& comp_gphi,
                                        amrex::MultiFab& comp_minus_level_phi,
                                        amrex::Vector<std::unique_ptr<amrex::MultiFab> >& comp_minus_level_grad_phi);


///
/// Get coarse phi on level ``level``-1
///
/// @param level        level index of fine data
/// @param phi_crse     MultiFab to contain coarse phi
/// @param time         Current time
///
  void GetCrsePhi(int level,
                  amrex::MultiFab& phi_crse,
                  amrex::Real time);


///
/// Get old gravity vector
///
/// @param level        Level index
/// @param grav_vector  MultiFab containing gravity vector
/// @param time         Current time
///
  void get_old_grav_vector(int level, amrex::MultiFab& grav_vector, amrex::Real time);

///
/// Get new gravity vector
///
/// @param level        Level index
/// @param grav_vector  MultiFab containing gravity vector
/// @param time         Current time
///
  void get_new_grav_vector(int level, amrex::MultiFab& grav_vector, amrex::Real time);

///
/// Test whether using the edge-based gradients
/// to compute Div(Grad(Phi)) satisfies Lap(phi) = RHS
///
/// @param bx          box
/// @param rhs         right-hand-side
/// @param ecx         gradients wrt x
/// @param ecy         gradients wrt y
/// @param ecz         gradients wrt z
/// @param coord_type  coordinate system

  static void test_residual (const amrex::Box& bx,
                      amrex::Array4<amrex::Real> const& rhs,
                      amrex::Array4<amrex::Real> const& ecx,
#if AMREX_SPACEDIM >= 2
                      amrex::Array4<amrex::Real> const& ecy,
#endif
#if AMREX_SPACEDIM == 3
                      amrex::Array4<amrex::Real> const& ecz,
#endif
                      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
                      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> problo,
                      int coord_type);

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
/// @param level    Level index
/// @param is_new   Use new state data (1) or old state data (0)?
///
  void average_fine_ec_onto_crse_ec(int level, int is_new);


///
/// Make radial gravity
///
/// @param level        Level index
/// @param time         Current time
/// @param radial_grav  Vector containing gravity in radial direction
///
  void make_radial_gravity(int level, amrex::Real time, RealVector& radial_grav);

///
/// Interpolate monopole gravitational field in the radial direction
///
/// @param level        Level index
/// @param radial_grav  Radial gravity
/// @param grav_vector  Gravity vector
///
  void interpolate_monopole_grav(int level, RealVector& radial_grav, amrex::MultiFab& grav_vector) const;

///
/// Integrate radially outward to find radial mass distribution
///
/// @param bx           Box
/// @param u            Density (RHS of Poisson's equation)
/// @param radial_mass  Radially integrated mass
/// @param radial_vol   Radially integrated volume
/// @param radial_pres  Radially integrated pressure
/// @param n1d          Number of radial points in the domain
/// @param level        Level index
///
  void compute_radial_mass(const amrex::Box& bx,
                           amrex::Array4<amrex::Real const> const u,
                           RealVector& radial_mass,
                           RealVector& radial_vol,
#ifdef GR_GRAV
                           RealVector& radial_pres,
#endif
                           int n1d, int level) const;

///
/// Implement multipole boundary conditions
///
/// @param crse_level
/// @param fine_level
/// @param Rhs
/// @param phi
///
  void fill_multipole_BCs(int crse_level, int fine_level, const amrex::Vector<amrex::MultiFab*>& Rhs, amrex::MultiFab& phi);

///
/// Initialize multipole gravity
///
  void init_multipole_grav() const;

#if (AMREX_SPACEDIM == 3)

///
/// Compute and fill direct sum boundary conditions
///
/// @param crse_level   Index of coarse level
/// @param fine_level   Index of fine level
/// @param Rhs          Vector of MultiFabs, right hand side
/// @param phi          MultiFab, phi
///
  void fill_direct_sum_BCs(int crse_level, int fine_level, const amrex::Vector<amrex::MultiFab*>& Rhs, amrex::MultiFab& phi);
#endif

///
/// Make multigrid boundary conditions
///
  void make_mg_bc();

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
  int finest_level_allocated;

  amrex::BCRec*       phys_bc;

  std::array<amrex::MLLinOp::BCType,AMREX_SPACEDIM> mlmg_lobc;
  std::array<amrex::MLLinOp::BCType,AMREX_SPACEDIM> mlmg_hibc;

  int   numpts_at_level;

  static int   test_solves;
  static amrex::Real  mass_offset;
  amrex::Vector< RealVector > radial_grav_old;
  amrex::Vector< RealVector > radial_grav_new;
  amrex::Vector< RealVector > radial_mass;
  amrex::Vector< RealVector > radial_vol;
#ifdef GR_GRAV
  amrex::Vector< RealVector > radial_pres;
#endif
  static int   stencil_type;

  static amrex::Real max_radius_all_in_domain;

#if (AMREX_SPACEDIM < 3)

///
/// Modify Rhs and coeffs with the appropriate metric terms.
///
/// @param level    Level index
/// @param Rhs      MultiFab
/// @param coeffs   Vector of MultiFabs
///
  void applyMetricTerms(int level,amrex::MultiFab& Rhs, const amrex::Vector<amrex::MultiFab*>& coeffs) const;

///
///
///
/// @param level    Index of level
/// @param cc       Cell-centered data
///
  void unweight_cc(int level,amrex::MultiFab& cc) const;

///
///
///
/// @param level    Index of level
/// @param edges    Edge-based data
///
  void unweight_edges(int level, const amrex::Vector<amrex::MultiFab*>& edges) const;
#endif

///
/// @param level        Index of level
/// @param phi          Gravitational potential
/// @param grav_vector  Gravity vector
///
    void add_pointmass_to_gravity (int level, amrex::MultiFab& phi, amrex::MultiFab& grav_vector) const;

///
/// Get the rhs
///
/// @param crse_level   Index of coarse level
/// @param nlevs        Number of levels
/// @param is_new       Use new (1) or old (0) state data
///
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > get_rhs (int crse_level, int nlevs, int is_new);


///
/// This is a sanity check on whether we are trying to fill multipole boundary conditiosn
///  for grids at this level > 0 -- this case is not currently supported.
///  Here we shrink the domain at this level by 1 in any direction which is not symmetry or periodic,
///  then ask if the grids at this level are contained in the shrunken domain.  If not, then grids
///  at this level touch the domain boundary and we must abort.
///
/// @param level    Level index
///
    void sanity_check (int level);

///
/// Do multigrid solve
///
/// @param crse_level   Coarse level index
/// @param fine_level   Fine level index
/// @param phi          Gravitational potential
/// @param rhs          Right hand side
/// @param grad_phi     Grad phi
/// @param res
/// @param crse_bcdata
/// @param rel_eps      Relative tolerance
/// @param abs_eps      Absolute tolerance
///
    amrex::Real actual_solve_with_mlmg (int crse_level, int fine_level,
                                        const amrex::Vector<amrex::MultiFab*>& phi,
                                        const amrex::Vector<const amrex::MultiFab*>& rhs,
                                        const amrex::Vector<std::array<amrex::MultiFab*,AMREX_SPACEDIM> >& grad_phi,
                                        const amrex::Vector<amrex::MultiFab*>& res,
                                        const amrex::MultiFab* const crse_bcdata,
                                        amrex::Real rel_eps, amrex::Real abs_eps) const;


///
/// Do multigrid solve to find phi
///
/// @param crse_level   Coarse level index
/// @param fine_level   Fine level index
/// @param phi          Gravitational potential
/// @param rhs          Right hand side source term
/// @param grad_phi     Grad phi
/// @param res
/// @param time         Current time
///
    amrex::Real solve_phi_with_mlmg (int crse_level, int fine_level,
                                     const amrex::Vector<amrex::MultiFab*>& phi,
                                     const amrex::Vector<amrex::MultiFab*>& rhs,
                                     const amrex::Vector<amrex::Vector<amrex::MultiFab*> >& grad_phi,
                                     const amrex::Vector<amrex::MultiFab*>& res,
                                     amrex::Real time);

  static inline amrex::Real
  get_const_grav() {
    return gravity::const_grav;
  }


};

///
/// @class GradPhiPhysBCFunct
/// @brief A physical boundary condition function for grad phi
///
using GradPhiPhysBCFunct = amrex::PhysBCFunctNoOp;


#endif
