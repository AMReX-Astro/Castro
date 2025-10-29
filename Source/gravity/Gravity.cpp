#include <cmath>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_ParmParse.H>
#include <Gravity.H>
#include <Castro.H>

#include <AMReX_FillPatchUtil.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>

#include <castro_limits.H>

#include <fundamental_constants.H>

#include <Gravity_util.H>
#include <MGutils.H>

using namespace amrex;

#ifdef AMREX_DEBUG
int Gravity::test_solves  = 1;
#else
int Gravity::test_solves  = 0;
#endif
Real Gravity::mass_offset    =  0.0;

// ************************************************************************************** //

// Ggravity is defined as 4 * pi * G, where G is the gravitational constant.

// In CGS, this constant is currently
//      Gconst   =  6.67428e-8           cm^3/g/s^2 , which results in
//      Ggravity =  83.8503442814844e-8  cm^3/g/s^2

// ************************************************************************************** //

const Real Ggravity = 4.0 * M_PI * C::Gconst;

///
/// Multipole gravity data
///
AMREX_GPU_MANAGED Real multipole::volumeFactor;
AMREX_GPU_MANAGED Real multipole::parityFactor;

AMREX_GPU_MANAGED Real multipole::rmax;

AMREX_GPU_MANAGED Array1D<bool, 0, 2> multipole::doSymmetricAddLo;
AMREX_GPU_MANAGED Array1D<bool, 0, 2> multipole::doSymmetricAddHi;
AMREX_GPU_MANAGED bool multipole::doSymmetricAdd;

AMREX_GPU_MANAGED Array1D<bool, 0, 2> multipole::doReflectionLo;
AMREX_GPU_MANAGED Array1D<bool, 0, 2> multipole::doReflectionHi;

AMREX_GPU_MANAGED Array2D<Real, 0, multipole::lnum_max, 0, multipole::lnum_max> multipole::factArray;
AMREX_GPU_MANAGED Array1D<Real, 0, multipole::lnum_max> multipole::parity_q0;
AMREX_GPU_MANAGED Array2D<Real, 0, multipole::lnum_max, 0, multipole::lnum_max> multipole::parity_qC_qS;

Gravity::Gravity(Amr* Parent, int _finest_level, BCRec* _phys_bc, int _density)
  :
    parent(Parent),
    LevelData(MAX_LEV),
    grad_phi_curr(MAX_LEV),
    grad_phi_prev(MAX_LEV),
    grids(Parent->boxArray()),
    dmap(Parent->DistributionMap()),
    abs_tol(MAX_LEV),
    rel_tol(MAX_LEV),
    level_solver_resnorm(MAX_LEV),
    volume(MAX_LEV),
    area(MAX_LEV),
    phys_bc(_phys_bc)
{

     amrex::ignore_unused(_finest_level);

     AMREX_ALWAYS_ASSERT(parent->maxLevel() < MAX_LEV);

     Density = _density;
     read_params();
     finest_level_allocated = -1;

     radial_grav_old.resize(MAX_LEV);
     radial_grav_new.resize(MAX_LEV);
     radial_mass.resize(MAX_LEV);
     radial_vol.resize(MAX_LEV);
#ifdef GR_GRAV
     radial_pres.resize(MAX_LEV);
#endif

     if (gravity::gravity_type == "PoissonGrav") {
         make_mg_bc();
         init_multipole_grav();
     }
     max_rhs = 0.0;
     numpts_at_level = -1;
}

Gravity::~Gravity() = default;

void
Gravity::read_params ()
{
    static bool done = false;

    if (!done)
    {
        const Geometry& dgeom = DefaultGeometry();

        ParmParse pp("gravity");

        if ( (gravity::gravity_type != "ConstantGrav") &&
             (gravity::gravity_type != "PoissonGrav") &&
             (gravity::gravity_type != "MonopoleGrav") )
             {
                std::cout << "Sorry -- dont know this gravity type"  << std::endl;
                amrex::Abort("Options are ConstantGrav, PoissonGrav, or MonopoleGrav");
             }

#if (AMREX_SPACEDIM == 1)
        if (gravity::gravity_type == "PoissonGrav")
        {
          amrex::Abort(" gravity::gravity_type = PoissonGrav doesn't work well in 1-d -- please set gravity::gravity_type = MonopoleGrav");
        }
        else if (gravity::gravity_type == "MonopoleGrav" && !(dgeom.IsSPHERICAL()))
        {
          amrex::Abort("Only use MonopoleGrav in 1D spherical coordinates");
        }
        else if (gravity::gravity_type == "ConstantGrav" && dgeom.IsSPHERICAL())
        {
          amrex::Abort("Can't use constant gravity in 1D spherical coordinates");
        }

#elif (AMREX_SPACEDIM == 2)
        if (gravity::gravity_type == "MonopoleGrav" && dgeom.IsCartesian() )
        {
          amrex::Abort(" gravity::gravity_type = MonopoleGrav doesn't make sense in 2D Cartesian coordinates");
        }
#endif

        if (pp.contains("get_g_from_phi") && !gravity::get_g_from_phi && gravity::gravity_type == "PoissonGrav") {
            amrex::Print() << "Warning: gravity::gravity_type = PoissonGrav assumes get_g_from_phi is true" << std::endl;
        }

        int nlevs = parent->maxLevel() + 1;

        // Allow run-time input of solver tolerance. If the user
        // provides no value, set a reasonable default value on the
        // coarse level, and then increase it by ref_ratio**2 as the
        // levels get finer to account for the change in the absolute
        // scale of the Laplacian. If the user provides one value, use
        // that on the coarse level, and increase it the same way for
        // the fine levels. If the user provides more than one value,
        // we expect them to provide one for every level, and we do
        // not apply the ref_ratio effect.

        int n_abs_tol = pp.countval("abs_tol");

        if (n_abs_tol <= 1) {

            Real tol;

            if (n_abs_tol == 1) {

                pp.get("abs_tol", tol);

            } else {

                if (dgeom.IsCartesian()) {
                    tol = 1.e-11;
                } else {
                    tol = 1.e-10;
                }

            }

            abs_tol[0] = tol;

            // Account for the fact that on finer levels, the scale of the
            // Laplacian changes due to the zone size changing. We assume
            // dx == dy == dz, so it is fair to say that on each level the
            // tolerance should increase by the factor ref_ratio**2, since
            // in absolute terms the Laplacian increases by that ratio too.
            // The actual tolerance we'll send in is the effective tolerance
            // on the finest level that we solve for.

            for (int lev = 1; lev < nlevs; ++lev) {
                abs_tol[lev] = abs_tol[lev - 1] * std::pow(parent->refRatio(lev - 1)[0], 2);
            }

        } else if (n_abs_tol >= nlevs) {

            pp.getarr("abs_tol", abs_tol, 0, nlevs);

        } else {

            amrex::Abort("If you are providing multiple values for abs_tol, you must provide at least one value for every level up to amr.max_level.");

        }

        // For the relative tolerance, we can again accept a single
        // scalar (same for all levels) or one for all levels. The
        // default value is zero, so that we only use the absolute
        // tolerance.  The multigrid always chooses the looser of the
        // two criteria in determining whether the solve has
        // converged.

        // Note that the parameter rel_tol used to be known as ml_tol,
        // so if we detect that the user has set ml_tol but not
        // rel_tol, we'll accept that for specifying the relative
        // tolerance. ml_tol is now considered deprecated and will be
        // removed in a future release.

        std::string rel_tol_name = "rel_tol";

        if (pp.contains("ml_tol")) {

            amrex::Warning("The gravity parameter ml_tol has been renamed rel_tol. ml_tol is now deprecated.");

            if (!pp.contains("rel_tol")) {
                rel_tol_name = "ml_tol";
            }

        }

        int n_rel_tol = pp.countval(rel_tol_name);

        if (n_rel_tol <= 1) {

            Real tol;

            if (n_rel_tol == 1) {

                pp.get(rel_tol_name, tol);

            } else {

                tol = 0.0;

            }

            for (int lev = 0; lev < MAX_LEV; ++lev) {
                rel_tol[lev] = tol;
            }

        } else if (n_rel_tol >= nlevs) {

            pp.getarr(rel_tol_name, rel_tol, 0, nlevs);

        } else {

            amrex::Abort("If you are providing multiple values for rel_tol, you must provide at least one value for every level up to amr.max_level.");

        }

        done = true;
    }
}

void
Gravity::output_job_info_params(std::ostream& jobInfoFile)
{
#include <gravity_job_info_tests.H>
}


void
Gravity::set_numpts_in_gravity (int numpts)
{
  numpts_at_level = numpts;
}

void
Gravity::install_level (int                   level,
                        AmrLevel*             level_data,
                        MultiFab&             _volume,
                        MultiFab*             _area)
{
    if (gravity::verbose > 1) {
        amrex::Print() << "Installing Gravity level " << level << '\n';
    }

    LevelData[level] = level_data;

    volume[level] = &_volume;

    area[level] = _area;

    level_solver_resnorm[level] = 0.0;

    const Geometry& geom = level_data->Geom();

    if (gravity::gravity_type == "PoissonGrav") {

       const DistributionMapping& dm = level_data->DistributionMap();

       grad_phi_prev[level].resize(AMREX_SPACEDIM);
       for (int n=0; n<AMREX_SPACEDIM; ++n) {
           grad_phi_prev[level][n] = std::make_unique<MultiFab>(level_data->getEdgeBoxArray(n),dm,1,1);
       }
       grad_phi_curr[level].resize(AMREX_SPACEDIM);
       for (int n=0; n<AMREX_SPACEDIM; ++n) {
           grad_phi_curr[level][n] = std::make_unique<MultiFab>(level_data->getEdgeBoxArray(n),dm,1,1);
       }

    } else if (gravity::gravity_type == "MonopoleGrav") {

        if (!geom.isAllPeriodic())
        {
           int n1d = gravity::drdxfac*numpts_at_level;

           radial_grav_old[level].resize(n1d);
           radial_grav_new[level].resize(n1d);
           radial_mass[level].resize(n1d);
           radial_vol[level].resize(n1d);
#ifdef GR_GRAV
           radial_pres[level].resize(n1d);
#endif
        }

    }

    finest_level_allocated = level;
}

std::string Gravity::get_gravity_type()
{
  return gravity::gravity_type;
}

int Gravity::get_max_solve_level()
{
  return gravity::max_solve_level;
}

int Gravity::NoSync()
{
  return gravity::no_sync;
}

int Gravity::DoCompositeCorrection()
{
  return gravity::do_composite_phi_correction;
}

int Gravity::test_results_of_solves()
{
  return test_solves;
}

Vector<std::unique_ptr<MultiFab> >&
Gravity::get_grad_phi_prev(int level)
{
  return grad_phi_prev[level];
}

MultiFab*
Gravity::get_grad_phi_prev_comp(int level, int comp)
{
  return grad_phi_prev[level][comp].get();
}

Vector<std::unique_ptr<MultiFab> >&
Gravity::get_grad_phi_curr(int level)
{
  return grad_phi_curr[level];
}

void
Gravity::plus_grad_phi_curr(int level, Vector<std::unique_ptr<MultiFab> >& addend)
{
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
        grad_phi_curr[level][n]->plus(*addend[n],0,1,0);
    }
}

void
Gravity::swapTimeLevels (int level)
{
    BL_PROFILE("Gravity::swapTimeLevels()");

    if (gravity::gravity_type == "PoissonGrav") {
        for (int n=0; n < AMREX_SPACEDIM; n++) {
            std::swap(grad_phi_prev[level][n], grad_phi_curr[level][n]);
            grad_phi_curr[level][n]->setVal(1.e50);
        }
    }
}

void
Gravity::solve_for_phi (int               level,
                        MultiFab&         phi,
                        const Vector<MultiFab*>& grad_phi,
                        int               is_new)
{
    BL_PROFILE("Gravity::solve_for_phi()");

    if (gravity::verbose > 1) {
        amrex::Print() << " ... solve for phi at level " << level << std::endl;
    }

    const Real strt = ParallelDescriptor::second();

    if (is_new == 0) {
        sanity_check(level);
    }

    Real time;
    if (is_new == 1) {
      time = LevelData[level]->get_state_data(PhiGrav_Type).curTime();
    } else {
      time = LevelData[level]->get_state_data(PhiGrav_Type).prevTime();
    }

    // If we are below the max_solve_level, do the Poisson solve.
    // Otherwise, interpolate using a fillpatch from max_solve_level.

    if (level <= gravity::max_solve_level) {

        Vector<MultiFab*> phi_p(1, &phi);

        const auto& g_rhs = get_rhs(level, 1, is_new);

        Vector< Vector<MultiFab*> > grad_phi_p(1);
        grad_phi_p[0].resize(AMREX_SPACEDIM);
        for (int i = 0; i < AMREX_SPACEDIM ; i++) {
            grad_phi_p[0][i] = grad_phi[i];
        }

        Vector<MultiFab*> res_null;

        level_solver_resnorm[level] = solve_phi_with_mlmg(level, level,
                                                          phi_p,
                                                          amrex::GetVecOfPtrs(g_rhs),
                                                          grad_phi_p,
                                                          res_null,
                                                          time);

    }
    else {

        LevelData[level]->FillCoarsePatch(phi, 0, time, PhiGrav_Type, 0, 1, 1);

    }

    if (gravity::verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        amrex::Real end = ParallelDescriptor::second() - strt;
        amrex::Real llevel = level;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        amrex::Print() << "Gravity::solve_for_phi() time = " << end << " on level "
                       << llevel << std::endl << std::endl;
#ifdef BL_LAZY
        });
#endif
    }
}

void
Gravity::gravity_sync (int crse_level, int fine_level, const Vector<MultiFab*>& drho, const Vector<MultiFab*>& dphi)
{
    BL_PROFILE("Gravity::gravity_sync()");

    // There is no need to do a synchronization if
    // we didn't solve on the fine levels.

    if (fine_level > gravity::max_solve_level) {
        return;
    } else {
        fine_level = amrex::min(fine_level, gravity::max_solve_level);
    }

    BL_ASSERT(parent->finestLevel()>crse_level);
    if (gravity::verbose > 1 && ParallelDescriptor::IOProcessor()) {
          std::cout << " ... gravity_sync at crse_level " << crse_level << '\n';
          std::cout << " ...     up to finest_level     " << fine_level << '\n';
    }

    const Geometry& crse_geom = parent->Geom(crse_level);
    const Box& crse_domain = crse_geom.Domain();

    int nlevs = fine_level - crse_level + 1;

    // Construct delta(phi) and delta(grad_phi). delta(phi)
    // needs a ghost zone for holding the boundary condition
    // in the same way that phi does.

    Vector<std::unique_ptr<MultiFab> > delta_phi(nlevs);

    for (int lev = crse_level; lev <= fine_level; ++lev) {
        delta_phi[lev - crse_level] = std::make_unique<MultiFab>(grids[lev], dmap[lev], 1, 1);
        delta_phi[lev - crse_level]->setVal(0.0);
    }

    Vector< Vector<std::unique_ptr<MultiFab> > > ec_gdPhi(nlevs);

    for (int lev = crse_level; lev <= fine_level; ++lev) {
        ec_gdPhi[lev - crse_level].resize(AMREX_SPACEDIM);

        const DistributionMapping& dm = LevelData[lev]->DistributionMap();
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            ec_gdPhi[lev - crse_level][n] = std::make_unique<MultiFab>(LevelData[lev]->getEdgeBoxArray(n), dm, 1, 0);
            ec_gdPhi[lev - crse_level][n]->setVal(0.0);
        }
    }

    // Construct a container for the right-hand-side (4 * pi * G * drho + dphi).
    // dphi appears in the construction of the boundary conditions because it
    // indirectly represents a change in mass on the domain (the mass motion that
    // occurs on the fine grid, whose gravitational effects are now indirectly
    // being propagated to the coarse grid).

    // We will temporarily leave the RHS divided by (4 * pi * G) because that
    // is the form expected by the boundary condition routine.

    Vector<std::unique_ptr<MultiFab> > g_rhs(nlevs);

    for (int lev = crse_level; lev <= fine_level; ++lev) {
        g_rhs[lev - crse_level] = std::make_unique<MultiFab>(LevelData[lev]->boxArray(), LevelData[lev]->DistributionMap(), 1, 0);
        MultiFab::Copy(*g_rhs[lev - crse_level], *dphi[lev - crse_level], 0, 0, 1, 0);
        g_rhs[lev - crse_level]->mult(1.0 / Ggravity);
        MultiFab::Add(*g_rhs[lev - crse_level], *drho[lev - crse_level], 0, 0, 1, 0);
    }

    // Construct the boundary conditions for the Poisson solve.

    if (crse_level == 0 && !crse_geom.isAllPeriodic()) {

        if (gravity::verbose > 1) {
            amrex::Print() << " ... Making bc's for delta_phi at crse_level 0"  << std::endl;
        }

#if (AMREX_SPACEDIM == 3)
      if ( gravity::direct_sum_bcs )
          fill_direct_sum_BCs(crse_level,fine_level,amrex::GetVecOfPtrs(g_rhs),*delta_phi[crse_level]);
      else {
          fill_multipole_BCs(crse_level,fine_level,amrex::GetVecOfPtrs(g_rhs),*delta_phi[crse_level]);
      }
#elif (AMREX_SPACEDIM == 2)
      fill_multipole_BCs(crse_level,fine_level,amrex::GetVecOfPtrs(g_rhs),*delta_phi[crse_level]);
#else
      fill_multipole_BCs(crse_level,fine_level,amrex::GetVecOfPtrs(g_rhs),*delta_phi[crse_level]);
#endif

    }

    // Restore the factor of (4 * pi * G) for the Poisson solve.
    for (int lev = crse_level; lev <= fine_level; ++lev)
        g_rhs[lev - crse_level]->mult(Ggravity);

    // In the all-periodic case we enforce that the RHS sums to zero.
    // We only do this if we're periodic and the coarse level covers the whole domain.
    // In principle this could be true for level > 0, so we'll test on whether the number
    // of points on the level is equal to the number of points possible on the level.
    // Note that since we did the average-down, we can stick with the data on the coarse
    // level since the averaging down is conservative.

    if (crse_geom.isAllPeriodic() && (grids[crse_level].numPts() == crse_domain.numPts()))
    {

        // We assume that if we're fully periodic then we're going to be in Cartesian
        // coordinates, so to get the average value of the RHS we can divide the sum
        // of the RHS by the number of points. This correction should probably be
        // volume weighted if we somehow got here without being Cartesian.

        Real local_correction = g_rhs[0]->sum() / static_cast<Real>(grids[crse_level].numPts());

        if (gravity::verbose > 1) {
            amrex::Print() << "WARNING: Adjusting RHS in gravity_sync solve by " << local_correction << '\n';
        }

        for (int lev = fine_level; lev >= crse_level; --lev) {
            g_rhs[lev-crse_level]->plus(-local_correction, 0, 1, 0);
        }
    }

    // Do multi-level solve for delta_phi.

    solve_for_delta_phi(crse_level, fine_level,
                        amrex::GetVecOfPtrs(g_rhs),
                        amrex::GetVecOfPtrs(delta_phi),
                        amrex::GetVecOfVecOfPtrs(ec_gdPhi));

    // In the all-periodic case we enforce that delta_phi averages to zero.

    if (crse_geom.isAllPeriodic() && (grids[crse_level].numPts() == crse_domain.numPts()) ) {

        Real local_correction = delta_phi[0]->sum() / static_cast<Real>(grids[crse_level].numPts());

        for (int lev = crse_level; lev <= fine_level; ++lev) {
            delta_phi[lev - crse_level]->plus(-local_correction, 0, 1, 1);
        }

    }

    // Add delta_phi to phi_new, and grad(delta_phi) to grad(delta_phi_curr) on each level.
    // Update the cell-centered gravity too.

    for (int lev = crse_level; lev <= fine_level; lev++) {

        LevelData[lev]->get_new_data(PhiGrav_Type).plus(*delta_phi[lev - crse_level], 0, 1, 0);

        for (int n = 0; n < AMREX_SPACEDIM; n++) {
            grad_phi_curr[lev][n]->plus(*ec_gdPhi[lev - crse_level][n], 0, 1, 0);
        }

        get_new_grav_vector(lev, LevelData[lev]->get_new_data(Gravity_Type),
                            LevelData[lev]->get_state_data(State_Type).curTime());

    }

    int is_new = 1;

    for (int lev = fine_level-1; lev >= crse_level; --lev)
    {

        // Average phi_new from fine to coarse level

        const IntVect& ratio = parent->refRatio(lev);

        amrex::average_down(LevelData[lev+1]->get_new_data(PhiGrav_Type),
                             LevelData[lev  ]->get_new_data(PhiGrav_Type),
                             0, 1, ratio);

        // Average the edge-based grad_phi from finer to coarser level

        average_fine_ec_onto_crse_ec(lev, is_new);

        // Average down the gravitational acceleration too.

        amrex::average_down(LevelData[lev+1]->get_new_data(Gravity_Type),
                             LevelData[lev  ]->get_new_data(Gravity_Type),
                             0, 1, ratio);

    }

}

void
Gravity::GetCrsePhi(int level,
                    MultiFab& phi_crse,
                    Real      time      )
{
    BL_PROFILE("Gravity::GetCrsePhi()");

    BL_ASSERT(level!=0);

    phi_crse.clear();
    phi_crse.define(grids[level-1], dmap[level-1], 1, 1); // BUT NOTE we don't trust phi's ghost cells.

    const Real t_old = LevelData[level-1]->get_state_data(PhiGrav_Type).prevTime();
    const Real t_new = LevelData[level-1]->get_state_data(PhiGrav_Type).curTime();
    Real alpha = (time - t_old) / (t_new - t_old);
    Real omalpha = 1.0_rt - alpha;
    const Real threshold = 1.e-6_rt;

    MultiFab const& phi_new = LevelData[level-1]->get_new_data(PhiGrav_Type);

    if (std::abs(omalpha) < threshold) {
        MultiFab::Copy(phi_crse, phi_new, 0, 0, 1, 1);
    }
    else if (std::abs(alpha) < threshold) {
        // Note we only access the old time if it's actually needed, to guard against
        // scenarios where it may not be allocated yet, for example after a restart when
        // the old time was not dumped to the checkpoint.
        MultiFab const& phi_old = LevelData[level-1]->get_old_data(PhiGrav_Type);
        MultiFab::Copy(phi_crse, phi_old, 0, 0, 1, 1);
    }
    else {
        MultiFab const& phi_old = LevelData[level-1]->get_old_data(PhiGrav_Type);
        MultiFab::LinComb(phi_crse, alpha, phi_new, 0, omalpha, phi_old, 0, 0, 1, 1);
    }

    const Geometry& geom = parent->Geom(level-1);
    phi_crse.FillBoundary(geom.periodicity());
}

void
Gravity::multilevel_solve_for_new_phi (int level, int finest_level_in)
{
    BL_PROFILE("Gravity::multilevel_solve_for_new_phi()");

    if (gravity::verbose > 1) {
        amrex::Print() << "... multilevel solve for new phi at base level " << level << " to finest level " << finest_level_in << std::endl;
    }

    const Real strt = ParallelDescriptor::second();

    for (int lev = level; lev <= finest_level_in; lev++) {
       BL_ASSERT(grad_phi_curr[lev].size()==AMREX_SPACEDIM);
       for (int n=0; n<AMREX_SPACEDIM; ++n)
       {
           grad_phi_curr[lev][n] = std::make_unique<MultiFab>(LevelData[lev]->getEdgeBoxArray(n),
                                                              LevelData[lev]->DistributionMap(),1,1);
       }
    }

    int is_new = 1;
    actual_multilevel_solve(level, finest_level_in, amrex::GetVecOfVecOfPtrs(grad_phi_curr), is_new);

    if (gravity::verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        amrex::Print() << "Gravity::multilevel_solve_for_new_phi() time = " << end << std::endl << std::endl;
#ifdef BL_LAZY
        });
#endif
    }
}

void
Gravity::actual_multilevel_solve (int crse_level, int finest_level_in,
                                  const Vector<Vector<MultiFab*> >& grad_phi,
                                  int is_new)
{
    BL_PROFILE("Gravity::actual_multilevel_solve()");

    for (int ilev = crse_level; ilev <= finest_level_in ; ++ilev) {
        sanity_check(ilev);
    }

    int nlevels = finest_level_in - crse_level + 1;

    Vector<MultiFab*> phi_p(nlevels);
    for (int ilev = 0; ilev < nlevels; ilev++)
    {
        int amr_lev = ilev + crse_level;
        if (is_new == 1) {
            phi_p[ilev] = &LevelData[amr_lev]->get_new_data(PhiGrav_Type);
        } else {
            phi_p[ilev] = &LevelData[amr_lev]->get_old_data(PhiGrav_Type);
        }
    }

    const auto& g_rhs = get_rhs(crse_level, nlevels, is_new);

    Vector<Vector<MultiFab*> > grad_phi_p(nlevels);
    for (int ilev = 0; ilev < nlevels; ilev++)
    {
        int amr_lev = ilev + crse_level;
        grad_phi_p[ilev] = grad_phi[amr_lev];
    }

    Real time;
    if (is_new == 1) {
        time = LevelData[crse_level]->get_state_data(PhiGrav_Type).curTime();
    } else {
        time = LevelData[crse_level]->get_state_data(PhiGrav_Type).prevTime();
    }

    int fine_level = amrex::min(finest_level_in, gravity::max_solve_level);

    if (fine_level >= crse_level) {

        Vector<MultiFab*> res_null;
        solve_phi_with_mlmg(crse_level, fine_level,
                            phi_p, amrex::GetVecOfPtrs(g_rhs), grad_phi_p, res_null,
                            time);

        // Average phi from fine to coarse level
        for (int amr_lev = fine_level; amr_lev > crse_level; amr_lev--)
        {
            const IntVect& ratio = parent->refRatio(amr_lev-1);
            if (is_new == 1)
            {
                amrex::average_down(LevelData[amr_lev  ]->get_new_data(PhiGrav_Type),
                                    LevelData[amr_lev-1]->get_new_data(PhiGrav_Type),
                                    0, 1, ratio);
            }
            else if (is_new == 0)
            {
                amrex::average_down(LevelData[amr_lev  ]->get_old_data(PhiGrav_Type),
                                    LevelData[amr_lev-1]->get_old_data(PhiGrav_Type),
                                    0, 1, ratio);
            }
        }

        // Average grad_phi from fine to coarse level
        for (int amr_lev = fine_level; amr_lev > crse_level; amr_lev--) {
            average_fine_ec_onto_crse_ec(amr_lev-1,is_new);
        }

    }

    // For all levels on which we're not doing the solve, interpolate from
    // the coarsest level with correct data. Note that since FillCoarsePatch
    // fills from the coarse level just below it, we need to fill from the
    // lowest level upwards using successive interpolations.

    for (int amr_lev = gravity::max_solve_level+1; amr_lev <= finest_level_in; amr_lev++) {

        // Interpolate the potential.

        if (is_new == 1) {

            MultiFab& phi = LevelData[amr_lev]->get_new_data(PhiGrav_Type);

            LevelData[amr_lev]->FillCoarsePatch(phi,0,time,PhiGrav_Type,0,1,1);

        }
        else {

            MultiFab& phi = LevelData[amr_lev]->get_old_data(PhiGrav_Type);

            LevelData[amr_lev]->FillCoarsePatch(phi,0,time,PhiGrav_Type,0,1,1);

        }

        // Interpolate the grad_phi.

        // Instantiate a bare physical BC function for grad_phi. It doesn't do anything
        // since the fine levels for Poisson gravity do not touch the physical boundary.

        GradPhiPhysBCFunct gp_phys_bc;

        // We need to use a interpolater that works with data on faces.

        Interpolater* gp_interp = &face_linear_interp;

        // For the BCs, we will use the Gravity_Type BCs for convenience, but these will
        // not do anything because we do not fill on physical boundaries.

        const Vector<BCRec>& gp_bcs = LevelData[amr_lev]->get_desc_lst()[Gravity_Type].getBCs();

        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            amrex::InterpFromCoarseLevel(*grad_phi[amr_lev][n], time, *grad_phi[amr_lev-1][n],
                                         0, 0, 1,
                                         parent->Geom(amr_lev-1), parent->Geom(amr_lev),
                                         gp_phys_bc, 0, gp_phys_bc, 0, parent->refRatio(amr_lev-1),
                                         gp_interp, gp_bcs, 0);
        }

    }

}

void
Gravity::get_old_grav_vector(int level, MultiFab& grav_vector, Real time)
{
    BL_PROFILE("Gravity::get_old_grav_vector()");

    int ng = grav_vector.nGrow();

    // Fill data from the level below if we're not doing a solve on this level.

    if (level > gravity::max_solve_level) {

        LevelData[level]->FillCoarsePatch(grav_vector,0,time,Gravity_Type,0,3,ng);

        return;

    }

    // Note that grav_vector coming into this routine always has three components.
    // So we'll define a temporary MultiFab with AMREX_SPACEDIM dimensions.
    // Then at the end we'll copy in all AMREX_SPACEDIM dimensions from this into
    // the outgoing grav_vector, leaving any higher dimensions unchanged.

    MultiFab grav(grids[level], dmap[level], AMREX_SPACEDIM, ng);
    grav.setVal(0.0,ng);

    const Geometry& geom = parent->Geom(level);

    if (gravity::gravity_type == "ConstantGrav") {

        if (AMREX_SPACEDIM == 2 && geom.Coord() == 2) {
            // 2D spherical r-theta, we want g in the radial direction
            grav.setVal(gravity::const_grav, 0, 1, ng);
        } else {
            // Set to constant value in the AMREX_SPACEDIM direction and zero in all others.
            grav.setVal(gravity::const_grav, AMREX_SPACEDIM-1, 1, ng);
        }

    } else if (gravity::gravity_type == "MonopoleGrav") {

       const Real prev_time = LevelData[level]->get_state_data(State_Type).prevTime();
       make_radial_gravity(level,prev_time,radial_grav_old[level]);
       interpolate_monopole_grav(level,radial_grav_old[level],grav);

    } else if (gravity::gravity_type == "PoissonGrav") {

       amrex::average_face_to_cellcenter(grav, amrex::GetVecOfConstPtrs(grad_phi_prev[level]), geom);
       grav.mult(-1.0, ng); // g = - grad(phi)

    } else {
       amrex::Abort("Unknown gravity_type in get_old_grav_vector");
    }

    // Do the copy to the output vector.

    for (int dir = 0; dir < 3; dir++) {
        if (dir < AMREX_SPACEDIM) {
            MultiFab::Copy(grav_vector, grav, dir, dir, 1, ng);
        } else {
            grav_vector.setVal(0.,dir,1,ng);
        }
    }

#if (AMREX_SPACEDIM > 1)
    if (gravity::gravity_type != "ConstantGrav") {
        // Fill ghost cells
        AmrLevel* amrlev = &parent->getLevel(level) ;
        AmrLevel::FillPatch(*amrlev,grav_vector,ng,time,Gravity_Type,0,AMREX_SPACEDIM);
    }
#endif

    auto* cs = dynamic_cast<Castro*>(&parent->getLevel(level));
    if (cs->using_point_mass()) {
        MultiFab& phi = LevelData[level]->get_old_data(PhiGrav_Type);
        add_pointmass_to_gravity(level,phi,grav_vector);
    }
}

void
Gravity::get_new_grav_vector(int level, MultiFab& grav_vector, Real time)
{
    BL_PROFILE("Gravity::get_new_grav_vector()");

    int ng = grav_vector.nGrow();

    // Fill data from the level below if we're not doing a solve on this level.

    if (level > gravity::max_solve_level) {

        LevelData[level]->FillCoarsePatch(grav_vector,0,time,Gravity_Type,0,3,ng);

        return;

    }

    // Note that grav_vector coming into this routine always has three components.
    // So we'll define a temporary MultiFab with AMREX_SPACEDIM dimensions.
    // Then at the end we'll copy in all AMREX_SPACEDIM dimensions from this into
    // the outgoing grav_vector, leaving any higher dimensions unchanged.

    MultiFab grav(grids[level],dmap[level],AMREX_SPACEDIM,ng);
    grav.setVal(0.0,ng);
    const Geometry& geom = parent->Geom(level);

    if (gravity::gravity_type == "ConstantGrav") {

        if (AMREX_SPACEDIM == 2 && geom.Coord() == 2) {
            // 2D spherical r-theta, we want g in the radial direction
            grav.setVal(gravity::const_grav, 0, 1, ng);
        } else {
            // Set to constant value in the AMREX_SPACEDIM direction
            grav.setVal(gravity::const_grav, AMREX_SPACEDIM-1, 1, ng);
        }

    } else if (gravity::gravity_type == "MonopoleGrav") {

        // We always fill radial_grav_new (at every level)
        const Real cur_time = LevelData[level]->get_state_data(State_Type).curTime();
        make_radial_gravity(level,cur_time,radial_grav_new[level]);
        interpolate_monopole_grav(level,radial_grav_new[level],grav);

    } else if (gravity::gravity_type == "PoissonGrav") {

        amrex::average_face_to_cellcenter(grav, amrex::GetVecOfConstPtrs(grad_phi_curr[level]), geom);
        grav.mult(-1.0, ng); // g = - grad(phi)

    } else {
       amrex::Abort("Unknown gravity_type in get_new_grav_vector");
    }

    // Do the copy to the output vector.

    for (int dir = 0; dir < 3; dir++) {
        if (dir < AMREX_SPACEDIM) {
            MultiFab::Copy(grav_vector, grav, dir, dir, 1, ng);
        } else {
            grav_vector.setVal(0.,dir,1,ng);
        }
    }

#if (AMREX_SPACEDIM > 1)
    if (gravity::gravity_type != "ConstantGrav" && ng>0) {
        // Fill ghost cells
        AmrLevel* amrlev = &parent->getLevel(level) ;
        AmrLevel::FillPatch(*amrlev,grav_vector,ng,time,Gravity_Type,0,AMREX_SPACEDIM);
    }
#endif

    auto* cs = dynamic_cast<Castro*>(&parent->getLevel(level));
    if (cs->using_point_mass()) {
        MultiFab& phi = LevelData[level]->get_new_data(PhiGrav_Type);
        add_pointmass_to_gravity(level,phi,grav_vector);
    }
}

void
Gravity::test_residual (const Box& bx,
                        Array4<Real> const& g_rhs,
                        Array4<Real> const& ecx,
#if AMREX_SPACEDIM >= 2
                        Array4<Real> const& ecy,
#endif
#if AMREX_SPACEDIM == 3
                        Array4<Real> const& ecz,
#endif
                        GpuArray<Real, AMREX_SPACEDIM> dx,
                        GpuArray<Real, AMREX_SPACEDIM> problo,
                        int coord_type)
{
    // Test whether using the edge-based gradients
    // to compute Div(Grad(Phi)) satisfies Lap(phi) = RHS
    // Fill the RHS array with the residual

    AMREX_ALWAYS_ASSERT(coord_type >= 0 && coord_type <= 2);

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // Cartesian
        if (coord_type == 0) {

            Real lapphi = (ecx(i+1,j,k) - ecx(i,j,k)) / dx[0];
#if AMREX_SPACEDIM >= 2
            lapphi += (ecy(i,j+1,k) - ecy(i,j,k)) / dx[1];
#endif
#if AMREX_SPACEDIM == 3
            lapphi += (ecz(i,j,k+1) - ecz(i,j,k)) / dx[2];
#endif

            g_rhs(i,j,k) -= lapphi;

        // r-z
        } else if (coord_type == 1) {

            Real rlo  = problo[0] + static_cast<Real>(i) * dx[0];
            Real rhi  = rlo + dx[0];
            Real rcen = 0.5_rt * (rlo + rhi);

            Real lapphi = (rhi * ecx(i+1,j,k) - rlo * ecx(i,j,k)) / (rcen * dx[0]);
#if AMREX_SPACEDIM >= 2
            lapphi += (ecy(i,j+1,k) - ecy(i,j,k)) / dx[1];
#endif

            g_rhs(i,j,k) -= lapphi;

        // spherical
        } else if (coord_type == 2) {

            Real rlo  = problo[0] + static_cast<Real>(i) * dx[0];
            Real rhi  = rlo + dx[0];
            Real rcen = 0.5_rt * (rlo + rhi);

            Real lapphi = (rhi * rhi * ecx(i+1,j,k) - rlo * rlo * ecx(i,j,k)) / (rcen * rcen * dx[0]);
            g_rhs(i,j,k) -= lapphi;

        }
    });
}

void
Gravity::test_level_grad_phi_prev(int level)
{
    BL_PROFILE("Gravity::test_level_grad_phi_prev()");

    // Fill the RHS for the solve
    const MultiFab& S_old = LevelData[level]->get_old_data(State_Type);
    MultiFab Rhs(grids[level],dmap[level],1,0);
    MultiFab::Copy(Rhs,S_old, URHO,0,1,0);

    const Geometry& geom = parent->Geom(level);

    // This is a correction for fully periodic domains only
    if ( geom.isAllPeriodic() )
    {
       if (gravity::verbose > 1 && ParallelDescriptor::IOProcessor() && mass_offset != 0.0) {
          std::cout << " ... subtracting average density from RHS at level ... "
                    << level << " " << mass_offset << std::endl;
       }
       Rhs.plus(-mass_offset,0,1,0);
    }

    Rhs.mult(Ggravity);

    if (gravity::verbose > 1) {
       Real rhsnorm = Rhs.norm0();
       amrex::Print() << "... test_level_grad_phi_prev at level " << level << std::endl;
       amrex::Print() << "       norm of RHS             " << rhsnorm << std::endl;
    }

    auto dx     = parent->Geom(level).CellSizeArray();
    auto problo = parent->Geom(level).ProbLoArray();
    const int coord_type = geom.Coord();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(Rhs, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        test_residual(bx,
                      Rhs.array(mfi),
                      (*grad_phi_prev[level][0]).array(mfi),
#if AMREX_SPACEDIM >= 2
                      (*grad_phi_prev[level][1]).array(mfi),
#endif
#if AMREX_SPACEDIM == 3
                      (*grad_phi_prev[level][2]).array(mfi),
#endif
                      dx, problo, coord_type);
    }

    if (gravity::verbose > 1) {
       Real resnorm = Rhs.norm0();
       amrex::Print() << "       norm of residual        " << resnorm << std::endl;
    }
}

void
Gravity::test_level_grad_phi_curr(int level)
{
    BL_PROFILE("Gravity::test_level_grad_phi_curr()");

    // Fill the RHS for the solve
    const MultiFab& S_new = LevelData[level]->get_new_data(State_Type);
    MultiFab Rhs(grids[level],dmap[level],1,0);
    MultiFab::Copy(Rhs,S_new, URHO, 0,1,0);

    const Geometry& geom = parent->Geom(level);

    // This is a correction for fully periodic domains only
    if ( geom.isAllPeriodic() )
    {
       if (gravity::verbose > 1 && ParallelDescriptor::IOProcessor() && mass_offset != 0.0) {
           std::cout << " ... subtracting average density from RHS in solve ... " << mass_offset << std::endl;
       }
       Rhs.plus(-mass_offset,0,1,0);
    }

    Rhs.mult(Ggravity);

    if (gravity::verbose > 1) {
       Real rhsnorm = Rhs.norm0();
       if (ParallelDescriptor::IOProcessor()) {
          std::cout << "... test_level_grad_phi_curr at level " << level << std::endl;
          std::cout << "       norm of RHS             " << rhsnorm << std::endl;
        }
    }

    auto dx     = geom.CellSizeArray();
    auto problo = geom.ProbLoArray();
    const int coord_type = geom.Coord();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(Rhs, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        test_residual(bx,
                      Rhs.array(mfi),
                      (*grad_phi_curr[level][0]).array(mfi),
#if AMREX_SPACEDIM >= 2
                      (*grad_phi_curr[level][1]).array(mfi),
#endif
#if AMREX_SPACEDIM == 3
                      (*grad_phi_curr[level][2]).array(mfi),
#endif
                      dx, problo, coord_type);
    }

    if (gravity::verbose > 1) {
       Real resnorm = Rhs.norm0();
       amrex::Print() << "       norm of residual        " << resnorm << std::endl;
    }
}

void
Gravity::create_comp_minus_level_grad_phi(int level,
                                          MultiFab& comp_phi,
                                          const Vector<MultiFab*>& comp_gphi,
                                          MultiFab& comp_minus_level_phi,
                                          Vector<std::unique_ptr<MultiFab> >& comp_minus_level_grad_phi)
{
    BL_PROFILE("Gravity::create_comp_minus_level_grad_phi()");

    if (gravity::verbose > 1 && ParallelDescriptor::IOProcessor()) {
        std::cout << "\n";
        std::cout << "... compute difference between level and composite solves at level " << level << "\n";
        std::cout << "\n";
    }

    comp_minus_level_phi.define(LevelData[level]->boxArray(),
                                LevelData[level]->DistributionMap(),
                                1, 0);

    MultiFab::Copy(comp_minus_level_phi, comp_phi, 0, 0, 1, 0);
    comp_minus_level_phi.minus(parent->getLevel(level).get_old_data(PhiGrav_Type), 0, 1, 0);

    comp_minus_level_grad_phi.resize(AMREX_SPACEDIM);
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        comp_minus_level_grad_phi[n] = std::make_unique<MultiFab>(LevelData[level]->getEdgeBoxArray(n),
                                                                  LevelData[level]->DistributionMap(), 1, 0);
        MultiFab::Copy(*comp_minus_level_grad_phi[n], *comp_gphi[n], 0, 0, 1, 0);
        comp_minus_level_grad_phi[n]->minus(*grad_phi_prev[level][n], 0, 1, 0);
    }

}

void
Gravity::average_fine_ec_onto_crse_ec(int level, int is_new)
{
    BL_PROFILE("Gravity::average_fine_ec_onto_crse_ec()");

    // NOTE: this is called with level == the coarser of the two levels involved
    if (level == parent->finestLevel()) {
        return;
    }

    //
    // Coarsen() the fine stuff on processors owning the fine data.
    //
    BoxArray crse_gphi_fine_BA(grids[level+1].size());

    IntVect fine_ratio = parent->refRatio(level);

    for (int i = 0; i < crse_gphi_fine_BA.size(); ++i) {
        crse_gphi_fine_BA.set(i,amrex::coarsen(grids[level+1][i],fine_ratio));
    }

    Vector<std::unique_ptr<MultiFab> > crse_gphi_fine(AMREX_SPACEDIM);
    for (int n=0; n<AMREX_SPACEDIM; ++n)
    {
        BoxArray eba = crse_gphi_fine_BA;
        eba.surroundingNodes(n);
        crse_gphi_fine[n] = std::make_unique<MultiFab>(eba,dmap[level+1],1,0);
    }

    auto& grad_phi = (is_new) ? grad_phi_curr : grad_phi_prev;

    amrex::average_down_faces(amrex::GetVecOfConstPtrs(grad_phi[level+1]),
                               amrex::GetVecOfPtrs(crse_gphi_fine),
                               fine_ratio);

    const Geometry& cgeom = parent->Geom(level);

    for (int n = 0; n < AMREX_SPACEDIM; ++n)
    {
        grad_phi[level][n]->ParallelCopy(*crse_gphi_fine[n], cgeom.periodicity());
    }
}

void
Gravity::test_composite_phi (int crse_level)
{
    BL_PROFILE("Gravity::test_composite_phi()");

    if (gravity::verbose > 1 && ParallelDescriptor::IOProcessor()) {
        std::cout << "   " << '\n';
        std::cout << "... test_composite_phi at base level " << crse_level << '\n';
    }

    int finest_level_local = parent->finestLevel();
    int nlevels = finest_level_local - crse_level + 1;

    Vector<std::unique_ptr<MultiFab> > phi(nlevels);
    Vector<std::unique_ptr<MultiFab> > g_rhs(nlevels);
    Vector<std::unique_ptr<MultiFab> > res(nlevels);
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        int amr_lev = crse_level + ilev;

        phi[ilev] = std::make_unique<MultiFab>(grids[amr_lev],dmap[amr_lev],1,1);
        MultiFab::Copy(*phi[ilev],
                       LevelData[amr_lev]->get_new_data(PhiGrav_Type),
                       0,0,1,1);

        g_rhs[ilev] = std::make_unique<MultiFab>(grids[amr_lev],dmap[amr_lev],1,1);
        MultiFab::Copy(*g_rhs[ilev],
                       LevelData[amr_lev]->get_new_data(State_Type),
                       URHO, 0,1,0);

        res[ilev] = std::make_unique<MultiFab>(grids[amr_lev],dmap[amr_lev],1,0);
        res[ilev]->setVal(0.);
    }

    Real time = LevelData[crse_level]->get_state_data(PhiGrav_Type).curTime();

    Vector< Vector<MultiFab*> > grad_phi_null;
    solve_phi_with_mlmg(crse_level, finest_level_local,
                        amrex::GetVecOfPtrs(phi),
                        amrex::GetVecOfPtrs(g_rhs),
                        grad_phi_null,
                        amrex::GetVecOfPtrs(res),
                        time);

    // Average residual from fine to coarse level before printing the norm
    for (int amr_lev = finest_level_local-1; amr_lev >= 0; --amr_lev)
    {
        const IntVect& ratio = parent->refRatio(amr_lev);
        int ilev = amr_lev - crse_level;
        amrex::average_down(*res[ilev+1], *res[ilev],
                             0, 1, ratio);
    }

    for (int amr_lev = crse_level; amr_lev <= finest_level_local; ++amr_lev) {
        Real resnorm = res[amr_lev]->norm0();
        amrex::Print() << "      ... norm of composite residual at level "
                       << amr_lev << "  " << resnorm << '\n';
    }
    amrex::Print() << std::endl;
}

void
Gravity::interpolate_monopole_grav(int level, RealVector& radial_grav, MultiFab& grav_vector) const
{
    BL_PROFILE("Gravity::interpolate_monopole_grav()");

    int n1d = static_cast<int>(radial_grav.size());

    const Geometry& geom = parent->Geom(level);
    const auto dx = geom.CellSizeArray();
    const Real dr = dx[0] / static_cast<Real>(gravity::drdxfac);

    const auto problo = geom.ProbLoArray();
    const auto geomdata = geom.data();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(grav_vector, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox();

        auto grav = grav_vector.array(mfi);
        const Real* const radial_grav_ptr = radial_grav.dataPtr();

        // Note that we are interpolating onto the entire range of grav,
        // including the ghost cells.

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            GpuArray<Real, 3> loc;

            loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0];

#if AMREX_SPACEDIM >= 2
            loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1];
#else
            loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
            loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - problem::center[2];
#else
            loc[2] = 0.0_rt;
#endif

            Real r = distance(geomdata, loc);

            int index = static_cast<int>(r / dr);

            Real cen = (static_cast<Real>(index) + 0.5_rt) * dr;
            Real xi = r - cen;

            Real mag_grav;

            if (index == 0) {

                // Linear interpolation or extrapolation
                Real slope = (radial_grav_ptr[index+1] - radial_grav_ptr[index]) / dr;
                mag_grav = radial_grav_ptr[index] + slope * xi;

            } else if (index == n1d-1) {

                // Linear interpolation or extrapolation
                Real slope = (radial_grav_ptr[index] - radial_grav_ptr[index-1]) / dr;
                mag_grav = radial_grav_ptr[index] + slope * xi;

            } else if (index > n1d-1) {

#ifndef AMREX_USE_GPU
                if (level == 0) {
                    std::cout << "PUT_RADIAL_GRAV: INDEX TOO BIG " << index << " > " << n1d-1 << "\n";
                    std::cout << "AT (i,j,k) " << i << " " << j << " " << k << "\n";
                    std::cout << "R / DR IS " << r << " " << dr << "\n";
                    amrex::Abort("Error:: Gravity.cpp :: interpolate_monopole_grav");
                } else {
                    // NOTE: we don't do anything to this point if it's outside the
                    //       radial grid and level > 0
                }
#endif

            } else {

                // Quadratic interpolation
                Real ghi = radial_grav_ptr[index+1];
                Real gmd = radial_grav_ptr[index  ];
                Real glo = radial_grav_ptr[index-1];
                mag_grav = ( ghi -   2.0_rt * gmd + glo) * xi * xi / (2.0_rt * dr * dr) +
                           ( ghi                  - glo) * xi      / (2.0_rt * dr     ) +
                           (-ghi + 26.e0_rt * gmd - glo) / 24.e0_rt;

                Real minvar = amrex::min(gmd, amrex::min(glo, ghi));
                Real maxvar = amrex::max(gmd, amrex::max(glo, ghi));
                mag_grav = amrex::max(mag_grav, minvar);
                mag_grav = amrex::min(mag_grav, maxvar);

            }

            if (index <= n1d-1) {

                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    grav(i,j,k,n) = mag_grav * (loc[n] / r);
                }

            }

        });
    }
}

void
Gravity::compute_radial_mass(const Box& bx,
                             Array4<Real const> const u,
                             RealVector& radial_mass_local,
                             RealVector& radial_vol_local,
#ifdef GR_GRAV
                             RealVector& radial_pres_local,
#endif
                             int n1d, int level) const
{
    const Geometry& geom = parent->Geom(level);

    GpuArray<Real, 3> dx, problo;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        dx[i] = geom.CellSizeArray()[i];
        problo[i] = geom.ProbLoArray()[i];
    }
    for (int i = AMREX_SPACEDIM; i < 3; ++i) {
        dx[i] = 0.0_rt;
        problo[i] = 0.0_rt;
    }

    Real dr = dx[0] / static_cast<Real>(gravity::drdxfac);
    Real drinv = 1.0_rt / dr;

    const int coord_type = geom.Coord();
    const auto geomdata = geom.data();

    AMREX_ALWAYS_ASSERT(coord_type >= 0 && coord_type <= 2);

    Real octant_factor = 1.0_rt;

    if (coord_type == 0) {

        if ((std::abs(problem::center[0] - problo[0]) < 1.e-2_rt * dx[0]) &&
            (std::abs(problem::center[1] - problo[1]) < 1.e-2_rt * dx[1]) &&
            (std::abs(problem::center[2] - problo[2]) < 1.e-2_rt * dx[2])) {

            octant_factor = 8.0_rt;

        }

    } else if (coord_type == 1) {

        if (std::abs(problem::center[1] - problo[1]) < 1.e-2_rt * dx[1]) {

            octant_factor = 2.0_rt;

        }

    }

    Real fac = static_cast<Real>(gravity::drdxfac);

    Real dx_frac = dx[0] / fac;
    Real dy_frac = dx[1] / fac;
    Real dz_frac = dx[2] / fac;

    Real* const radial_mass_ptr = radial_mass_local.dataPtr();
    Real* const radial_vol_ptr = radial_vol_local.dataPtr();
#ifdef GR_GRAV
    Real* const radial_pres_ptr = radial_pres_local.dataPtr();
#endif

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0];
        Real lo_i = problo[0] + static_cast<Real>(i) * dx[0] - problem::center[0];

        loc[1]= problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1];
        Real lo_j = problo[1] + static_cast<Real>(j) * dx[1] - problem::center[1];

        loc[2]= problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - problem::center[2];
        Real lo_k = problo[2] + static_cast<Real>(k) * dx[2] - problem::center[2];

        Real r = distance(geomdata, loc);
        int index = static_cast<int>(r * drinv);

        // We may be coming in here with a masked out zone (in a zone on a coarse
        // level underlying a fine level). We don't want to be calling the EOS in
        // this case, so we'll skip these masked out zones (which will have rho
        // exactly equal to zero).

        if (u(i,j,k,URHO) == 0.0_rt) {
            return;
        }

#ifdef GR_GRAV
        Real rhoInv = 1.0_rt / u(i,j,k,URHO);

        eos_t eos_state;

        eos_state.rho = u(i,j,k,URHO);
        eos_state.e   = u(i,j,k,UEINT) * rhoInv;
        eos_state.T   = u(i,j,k,UTEMP);
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = u(i,j,k,UFS+n) * rhoInv;
        }
#if NAUX_NET > 0
        for (int n = 0; n < NumAux; ++n) {
            eos_state.aux[n] = u(i,j,k,UFX+n) * rhoInv;
        }
#endif

        // Compute pressure from the EOS
        eos(eos_input_re, eos_state);
#endif

        if (index > n1d - 1) {

#ifndef AMREX_USE_GPU
            if (level == 0) {
                std::cout << "   " << "\n";
                std::cout << ">>> Error: Gravity_nd::ca_compute_radial_mass " << i << " " << j << " " << k << "\n";
                std::cout << ">>> ... index too big: " << index << " > " << n1d-1 << "\n";
                std::cout << ">>> ... at (i,j,k)   : " << i << " " << j << " " << k << "\n";
                amrex::Abort("Error:: Gravity_nd.F90 :: ca_compute_radial_mass");
            }
#endif

        } else {

            for (int kk = 0; kk <= dg2 * (gravity::drdxfac - 1); ++kk) {
                Real zz   = lo_k + (static_cast<Real>(kk) + 0.5_rt) * dz_frac;
                Real zzsq = zz * zz;

                for (int jj = 0; jj <= dg1 * (gravity::drdxfac - 1); ++jj) {
                    Real yy   = lo_j + (static_cast<Real>(jj) + 0.5_rt) * dy_frac;
                    Real yysq = yy * yy;

                    for (int ii = 0; ii <= gravity::drdxfac - 1; ++ii) {
                        Real xx    = lo_i + (static_cast<Real>(ii) + 0.5_rt) * dx_frac;
                        Real xxsq  = xx * xx;

                        r     = std::sqrt(xxsq + yysq + zzsq);
                        index = static_cast<int>(r * drinv);

                        Real vol_frac{};

                        if (coord_type == 0) {

                            vol_frac = octant_factor * dx_frac * dy_frac * dz_frac;

                        } else if (coord_type == 1) {

                            vol_frac = 2.0_rt * M_PI * dx_frac * dy_frac * octant_factor * xx;

                        } else if (coord_type == 2) {

                            Real rlo = std::abs(lo_i + static_cast<Real>(ii  ) * dx_frac);
                            Real rhi = std::abs(lo_i + static_cast<Real>(ii+1) * dx_frac);
                            vol_frac = (4.0_rt / 3.0_rt) * M_PI * (rhi * rhi * rhi - rlo * rlo * rlo);

                        }

                        if (index <= n1d - 1) {
                            Gpu::Atomic::Add(&radial_mass_ptr[index], vol_frac * u(i,j,k,URHO));
                            Gpu::Atomic::Add(&radial_vol_ptr[index], vol_frac);
#ifdef GR_GRAV
                            Gpu::Atomic::Add(&radial_pres_ptr[index], vol_frac * eos_state.p);
#endif
                        }

                    }
                }
            }

        }

    });
}

void
Gravity::init_multipole_grav() const
{
    if (gravity::lnum < 0) {
        amrex::Abort("lnum negative");
    }

    if (gravity::lnum > multipole::lnum_max) {
        amrex::Abort("lnum greater than lnum_max");
    }

    int lo_bc[3];
    int hi_bc[3];

    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
      lo_bc[dir] = phys_bc->lo(dir);
      hi_bc[dir] = phys_bc->hi(dir);
    }
    for (int dir = AMREX_SPACEDIM; dir < 3; dir++)
    {
      lo_bc[dir] = -1;
      hi_bc[dir] = -1;
    }

    const auto problo = parent->Geom(0).ProbLoArray();
    const auto probhi = parent->Geom(0).ProbHiArray();

    // If any of the boundaries are symmetric, we need to account for the mass that is assumed
    // to lie on the opposite side of the symmetric axis. If the center in any direction
    // coincides with the boundary, then we can simply double the mass as a result of that reflection.
    // Otherwise, we need to do a more general solve. We include a logical that is set to true
    // if any boundary is symmetric, so that we can avoid unnecessary function calls.

    multipole::volumeFactor = 1.0_rt;
    multipole::parityFactor = 1.0_rt;

    multipole::doSymmetricAdd = false;

    for (int n = 0; n < 3; ++n) {
        multipole::doSymmetricAddLo(n) = false;
        multipole::doSymmetricAddHi(n) = false;

        multipole::doReflectionLo(n) = false;
        multipole::doReflectionHi(n) = false;
    }

    const Real edgeTolerance = 1.0e-2_rt;

    for (int b = 0; b < AMREX_SPACEDIM; ++b) {

        if ((lo_bc[b] == amrex::PhysBCType::symmetry) && (parent->Geom(0).Coord() == 0)) {
            if (std::abs(problem::center[b] - problo[b]) < edgeTolerance) {
                multipole::volumeFactor *= 2.0_rt;
                multipole::doReflectionLo(b) = true;
            }
            else {
                multipole::doSymmetricAddLo(b) = true;
                multipole::doSymmetricAdd      = true;
            }
        }

        if ((hi_bc[b] == amrex::PhysBCType::symmetry) && (parent->Geom(0).Coord() == 0)) {
            if (std::abs(problem::center[b] - probhi[b]) < edgeTolerance) {
                multipole::volumeFactor *= 2.0_rt;
                multipole::doReflectionHi(b) = true;
            }
            else {
                multipole::doSymmetricAddHi(b) = true;
                multipole::doSymmetricAdd      = true;
            }
        }

    }

    // Compute pre-factors now to save computation time, for qC and qS

    for (int l = 0; l <= multipole::lnum_max; ++l) {
        multipole::parity_q0(l) = 1.0_rt;
        for (int m = 0; m <= multipole::lnum_max; ++m) {
            multipole::factArray(l,m) = 0.0_rt;
            multipole::parity_qC_qS(l,m) = 1.0_rt;
        }
    }

    for (int l = 0; l <= multipole::lnum_max; ++l) {

        // The odd l Legendre polynomials are odd in their argument, so
        // a symmetric reflection about the z axis leads to a total cancellation.

        multipole::parity_q0(l) = 1.0_rt;

        if (l % 2 != 0) {
            if (AMREX_SPACEDIM == 3 && (multipole::doReflectionLo(2) || multipole::doReflectionHi(2))) {  // NOLINT(bugprone-branch-clone)
                multipole::parity_q0(l) = 0.0_rt;
            }
            else if (AMREX_SPACEDIM == 2 && parent->Geom(0).Coord() == 1) {
                multipole::parity_q0(l) = 0.0_rt;
            }
            else if (AMREX_SPACEDIM == 1 && parent->Geom(0).Coord() == 2) {
                multipole::parity_q0(l) = 0.0_rt;
            }
        }

        // In 1D spherical, every term above l == 0 cancels out since the integration
        // of the Legendre polynomial from 0 to pi is zero.

        if (l > 0) {
            if (AMREX_SPACEDIM == 1 && parent->Geom(0).Coord() == 2) {
                multipole::parity_q0(l) = 0.0_rt;
            }
        }

        for (int m = 1; m <= l; ++m) {

            // The parity properties of the associated Legendre polynomials are:
            // P_l^m (-x) = (-1)^(l+m) P_l^m (x)
            // Therefore, a complete cancellation occurs if l+m is odd and
            // we are reflecting about the z axis.

            // Additionally, the cosine and sine terms flip sign when reflected
            // about the x or y axis, so if we have a reflection about x or y
            // then the terms have a complete cancellation.

            if (AMREX_SPACEDIM == 3) {
                multipole::parity_qC_qS(l,m) = 1.0_rt;
            }
            else if (AMREX_SPACEDIM == 2 && parent->Geom(0).Coord() == 1) {  // NOLINT(bugprone-branch-clone)
                multipole::parity_qC_qS(l,m) = 0.0_rt;
            }
            else if (AMREX_SPACEDIM == 1 && parent->Geom(0).Coord() == 2) {
                multipole::parity_qC_qS(l,m) = 0.0_rt;
            }

            if ((l+m) % 2 != 0 && (multipole::doReflectionLo(2) || multipole::doReflectionHi(2))) {
                multipole::parity_qC_qS(l,m) = 0.0_rt;
            }

            if (multipole::doReflectionLo(0) || multipole::doReflectionLo(1) || multipole::doReflectionHi(0) || multipole::doReflectionHi(1)) {
                multipole::parity_qC_qS(l,m) = 0.0_rt;
            }

            multipole::factArray(l,m) = 2.0_rt * factorial(l-m) / factorial(l+m) * multipole::volumeFactor;

        }

    }

    // Now let's take care of a safety issue. The multipole calculation involves taking powers of r^l,
    // which can overflow the floating point exponent limit if lnum is very large. Therefore,
    // we will normalize all distances to the maximum possible physical distance from the center,
    // which is the diagonal from the center to the edge of the box. Then r^l will always be
    // less than or equal to one. For large enough lnum, this may still result in roundoff
    // errors that don't make your answer any more precise, but at least it avoids
    // possible NaN issues from having numbers that are too large for double precision.
    // We will put the rmax factor back in at the end of ca_put_multipole_phi.

    Real maxWidth = probhi[0] - problo[0];
    if (AMREX_SPACEDIM >= 2) {
        maxWidth = amrex::max(maxWidth, probhi[1] - problo[1]);
    }
    if (AMREX_SPACEDIM == 3) {
        maxWidth = amrex::max(maxWidth, probhi[2] - problo[2]);
    }

    multipole::rmax = 0.5_rt * maxWidth * std::sqrt(static_cast<Real>(AMREX_SPACEDIM));  // NOLINT(modernize-use-std-numbers)
}

void
Gravity::fill_multipole_BCs(int crse_level, int fine_level, const Vector<MultiFab*>& Rhs, MultiFab& phi)
{
    BL_PROFILE("Gravity::fill_multipole_BCs()");

    // Multipole BCs only make sense to construct if we are starting from the coarse level.

    BL_ASSERT(crse_level == 0);

    BL_ASSERT(gravity::lnum >= 0);

    const Real strt = ParallelDescriptor::second();

#if (AMREX_SPACEDIM == 3)
    const int npts = numpts_at_level;
#else
    const int npts = 1;
#endif

    // Storage arrays for the multipole moments.
    // We will initialize them to zero, and then
    // sum up the results over grids.
    // Note that since Boxes are defined with
    // AMREX_SPACEDIM dimensions, we cannot presently
    // use this array to fill the interior of the
    // domain in 2D, since we can only have one
    // radial index for calculating the multipole moments.

    Box boxq0( IntVect(AMREX_D_DECL(0, 0, 0)), IntVect(AMREX_D_DECL(gravity::lnum, 0,    npts-1)) );
    Box boxqC( IntVect(AMREX_D_DECL(0, 0, 0)), IntVect(AMREX_D_DECL(gravity::lnum, gravity::lnum, npts-1)) );
    Box boxqS( IntVect(AMREX_D_DECL(0, 0, 0)), IntVect(AMREX_D_DECL(gravity::lnum, gravity::lnum, npts-1)) );

    FArrayBox qL0(boxq0);
    FArrayBox qLC(boxqC);
    FArrayBox qLS(boxqS);

    FArrayBox qU0(boxq0);
    FArrayBox qUC(boxqC);
    FArrayBox qUS(boxqS);

    qL0.setVal<RunOn::Device>(0.0);
    qLC.setVal<RunOn::Device>(0.0);
    qLS.setVal<RunOn::Device>(0.0);
    qU0.setVal<RunOn::Device>(0.0);
    qUC.setVal<RunOn::Device>(0.0);
    qUS.setVal<RunOn::Device>(0.0);

    // This section needs to be generalized for computing
    // full multipole gravity, not just BCs. At present this
    // does nothing.

#if (AMREX_SPACEDIM == 3)
    int boundary_only = 1;
#else
    const int boundary_only = 1;
#endif

    // Use all available data in constructing the boundary conditions,
    // unless the user has indicated that a maximum level at which
    // to stop using the more accurate data.

    for (int lev = crse_level; lev <= fine_level; ++lev) {

        // Create a local copy of the RHS so that we can mask it.

        MultiFab source(Rhs[lev - crse_level]->boxArray(),
                        Rhs[lev - crse_level]->DistributionMap(), 1, 0);

        MultiFab::Copy(source, *Rhs[lev - crse_level], 0, 0, 1, 0);

        if (lev < fine_level) {
        auto *castro_level = dynamic_cast<Castro*>(&(parent->getLevel(lev+1)));
        if (castro_level != nullptr) {
        const MultiFab& mask = castro_level->build_fine_mask();
        MultiFab::Multiply(source, mask, 0, 0, 1, 0);
        } else {
                amrex::Abort("unable to access mask");
            }
        }

        // Loop through the grids and compute the individual contributions
        // to the various moments. The multipole moment constructor
        // is coded to only add to the moment arrays, so it is safe
        // to directly hand the arrays to them.

        const auto dx = parent->Geom(lev).CellSizeArray();
        const auto problo = parent->Geom(lev).ProbLoArray();
        const auto probhi = parent->Geom(lev).ProbHiArray();
        int coord_type = parent->Geom(lev).Coord();

#ifdef _OPENMP
        int nthreads = omp_get_max_threads();
        Vector<std::unique_ptr<FArrayBox> > priv_qL0(nthreads);
        Vector<std::unique_ptr<FArrayBox> > priv_qLC(nthreads);
        Vector<std::unique_ptr<FArrayBox> > priv_qLS(nthreads);
        Vector<std::unique_ptr<FArrayBox> > priv_qU0(nthreads);
        Vector<std::unique_ptr<FArrayBox> > priv_qUC(nthreads);
        Vector<std::unique_ptr<FArrayBox> > priv_qUS(nthreads);
        for (int i=0; i<nthreads; i++) {
            priv_qL0[i].reset(new FArrayBox(boxq0));
            priv_qLC[i].reset(new FArrayBox(boxqC));
            priv_qLS[i].reset(new FArrayBox(boxqS));
            priv_qU0[i].reset(new FArrayBox(boxq0));
            priv_qUC[i].reset(new FArrayBox(boxqC));
            priv_qUS[i].reset(new FArrayBox(boxqS));
        }
#pragma omp parallel
#endif
        {
#ifdef _OPENMP
            int tid = omp_get_thread_num();
            priv_qL0[tid]->setVal<RunOn::Device>(0.0);
            priv_qLC[tid]->setVal<RunOn::Device>(0.0);
            priv_qLS[tid]->setVal<RunOn::Device>(0.0);
            priv_qU0[tid]->setVal<RunOn::Device>(0.0);
            priv_qUC[tid]->setVal<RunOn::Device>(0.0);
            priv_qUS[tid]->setVal<RunOn::Device>(0.0);
#endif
            for (MFIter mfi(source, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();

#ifdef _OPENMP
                auto qL0_arr = priv_qL0[tid]->array();
                auto qLC_arr = priv_qLC[tid]->array();
                auto qLS_arr = priv_qLS[tid]->array();
                auto qU0_arr = priv_qU0[tid]->array();
                auto qUC_arr = priv_qUC[tid]->array();
                auto qUS_arr = priv_qUS[tid]->array();
#else
                auto qL0_arr = qL0.array();
                auto qLC_arr = qLC.array();
                auto qLS_arr = qLS.array();
                auto qU0_arr = qU0.array();
                auto qUC_arr = qUC.array();
                auto qUS_arr = qUS.array();
#endif

                auto rho = source[mfi].array();
                auto vol = (*volume[lev])[mfi].array();

                amrex::ParallelFor(amrex::Gpu::KernelInfo().setReduction(true), bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::Gpu::Handler const& handler) noexcept
                {
                    // If we're using this to construct boundary values, then only fill
                    // the outermost bin.

                    int nlo = 0;
                    if (boundary_only == 1) {
                        nlo = npts-1;
                    }

                    // Note that we don't currently support dx != dy != dz, so this is acceptable.

                    Real drInv = multipole::rmax / dx[0];

                    Real rmax_cubed_inv = 1.0_rt / (multipole::rmax * multipole::rmax * multipole::rmax);

                    Real x = (problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0]) / multipole::rmax;

#if AMREX_SPACEDIM >= 2
                    Real y = (problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1]) / multipole::rmax;
#else
                    Real y = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
                    Real z = (problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - problem::center[2]) / multipole::rmax;
#else
                    Real z = 0.0_rt;
#endif

                    Real r = std::sqrt(x * x + y * y + z * z);

                    Real cosTheta{}, phiAngle{};
                    int index{};

                    if (AMREX_SPACEDIM == 3) {
                        index = static_cast<int>(r * drInv);
                        cosTheta = z / r;
                        phiAngle = std::atan2(y, x);
                    }
                    else if (AMREX_SPACEDIM == 2 && coord_type == 1) {
                        index = nlo; // We only do the boundary potential in 2D.
                        cosTheta = y / r;
                        phiAngle = z;
                    }
                    else if (AMREX_SPACEDIM == 1 && coord_type == 2) {
                        index = nlo; // We only do the boundary potential in 1D.
                        cosTheta = 1.0_rt;
                        phiAngle = 0.0_rt;
                    }

                    // Now, compute the multipole moments.

                    multipole_add(cosTheta, phiAngle, r, rho(i,j,k), vol(i,j,k) * rmax_cubed_inv,
                                  qL0_arr, qLC_arr, qLS_arr, qU0_arr, qUC_arr, qUS_arr,
                                  npts, nlo, index, handler, true);

                    // Now add in contributions if we have any symmetric boundaries in 3D.
                    // The symmetric boundary in 2D axisymmetric is handled separately.

                    if (multipole::doSymmetricAdd) {

                        multipole_symmetric_add(x, y, z, problo, probhi,
                                                rho(i,j,k), vol(i,j,k) * rmax_cubed_inv,
                                                qL0_arr, qLC_arr, qLS_arr, qU0_arr, qUC_arr, qUS_arr,
                                                npts, nlo, index, handler);

                    }
                });
            }

#ifdef _OPENMP
            int np0 = boxq0.numPts();
            int npC = boxqC.numPts();
            int npS = boxqS.numPts();
            Real* pL0 = qL0.dataPtr();
            Real* pLC = qLC.dataPtr();
            Real* pLS = qLS.dataPtr();
            Real* pU0 = qU0.dataPtr();
            Real* pUC = qUC.dataPtr();
            Real* pUS = qUS.dataPtr();
#pragma omp barrier
#pragma omp for nowait
            for (int i=0; i<np0; ++i)
            {
                for (int it=0; it<nthreads; it++) {
                    const Real* pp = priv_qL0[it]->dataPtr();
                    pL0[i] += pp[i];
                }
            }
#pragma omp for nowait
            for (int i=0; i<npC; ++i)
            {
                for (int it=0; it<nthreads; it++) {
                    const Real* pp = priv_qLC[it]->dataPtr();
                    pLC[i] += pp[i];
                }
            }
#pragma omp for nowait
            for (int i=0; i<npS; ++i)
            {
                for (int it=0; it<nthreads; it++) {
                    const Real* pp = priv_qLS[it]->dataPtr();
                    pLS[i] += pp[i];
                }
            }
#pragma omp for nowait
            for (int i=0; i<np0; ++i)
            {
                for (int it=0; it<nthreads; it++) {
                  const Real* pp = priv_qU0[it]->dataPtr();
                  pU0[i] += pp[i];
                }
            }
#pragma omp for nowait
            for (int i=0; i<npC; ++i)
            {
                for (int it=0; it<nthreads; it++) {
                  const Real* pp = priv_qUC[it]->dataPtr();
                  pUC[i] += pp[i];
                }
            }
#pragma omp for nowait
            for (int i=0; i<npS; ++i)
            {
                for (int it=0; it<nthreads; it++) {
                    const Real* pp = priv_qUS[it]->dataPtr();
                    pUS[i] += pp[i];
                }
            }
#endif

        } // end OpenMP parallel loop

    } // end loop over levels

    // Now, do a global reduce over all processes.

    Real* qL0_ptr = qL0.dataPtr();
    Real* qLC_ptr = qLC.dataPtr();
    Real* qLS_ptr = qLS.dataPtr();

    // Create temporary pinned host containers in case we need them.

    FArrayBox qL0_host(The_Pinned_Arena());
    FArrayBox qLC_host(The_Pinned_Arena());
    FArrayBox qLS_host(The_Pinned_Arena());

    if (!ParallelDescriptor::UseGpuAwareMpi()) {
        if (The_Arena() == The_Managed_Arena()) {
            qL0.prefetchToHost();
            qLC.prefetchToHost();
            qLS.prefetchToHost();
        }
        else if (The_Arena() == The_Device_Arena()) {
            qL0_host.resize(boxq0);
            qLC_host.resize(boxqC);
            qLS_host.resize(boxqS);

            qL0_host.copy<RunOn::Device>(qL0, boxq0);
            qLC_host.copy<RunOn::Device>(qLC, boxqC);
            qLS_host.copy<RunOn::Device>(qLS, boxqS);

            qL0_ptr = qL0_host.dataPtr();
            qLC_ptr = qLC_host.dataPtr();
            qLS_ptr = qLS_host.dataPtr();
        }
    }

    Gpu::synchronize();

    ParallelDescriptor::ReduceRealSum(qL0_ptr, static_cast<int>(boxq0.numPts()));
    ParallelDescriptor::ReduceRealSum(qLC_ptr, static_cast<int>(boxqC.numPts()));
    ParallelDescriptor::ReduceRealSum(qLS_ptr, static_cast<int>(boxqS.numPts()));

    if (!ParallelDescriptor::UseGpuAwareMpi()) {
        if (The_Arena() == The_Managed_Arena()) {
            qL0.prefetchToDevice();
            qLC.prefetchToDevice();
            qLS.prefetchToDevice();
        }
        else if (The_Arena() == The_Device_Arena()) {
            qL0.copy<RunOn::Device>(qL0_host, boxq0);
            qLC.copy<RunOn::Device>(qLC_host, boxqC);
            qLS.copy<RunOn::Device>(qLS_host, boxqS);
        }
    }

    if (boundary_only != 1) {

        Real* qU0_ptr = qU0.dataPtr();
        Real* qUC_ptr = qUC.dataPtr();
        Real* qUS_ptr = qUS.dataPtr();

        FArrayBox qU0_host(The_Pinned_Arena());
        FArrayBox qUC_host(The_Pinned_Arena());
        FArrayBox qUS_host(The_Pinned_Arena());

        if (!ParallelDescriptor::UseGpuAwareMpi()) {
            if (The_Arena() == The_Managed_Arena()) {
                qU0.prefetchToHost();
                qUC.prefetchToHost();
                qUS.prefetchToHost();
            }
            else if (The_Arena() == The_Device_Arena()) {
                qU0_host.resize(boxq0);
                qUC_host.resize(boxqC);
                qUS_host.resize(boxqS);

                qU0_host.copy<RunOn::Device>(qU0, boxq0);
                qUC_host.copy<RunOn::Device>(qUC, boxqC);
                qUS_host.copy<RunOn::Device>(qUS, boxqS);

                qU0_ptr = qU0_host.dataPtr();
                qUC_ptr = qUC_host.dataPtr();
                qUS_ptr = qUS_host.dataPtr();
            }
        }

        Gpu::synchronize();

        ParallelDescriptor::ReduceRealSum(qU0_ptr, static_cast<int>(boxq0.numPts()));
        ParallelDescriptor::ReduceRealSum(qUC_ptr, static_cast<int>(boxqC.numPts()));
        ParallelDescriptor::ReduceRealSum(qUS_ptr, static_cast<int>(boxqS.numPts()));

        if (!ParallelDescriptor::UseGpuAwareMpi()) {
            if (The_Arena() == The_Managed_Arena()) {
                qU0.prefetchToDevice();
                qUC.prefetchToDevice();
                qUS.prefetchToDevice();
            }
            else if (The_Arena() == The_Device_Arena()) {
                qU0.copy<RunOn::Device>(qU0_host, boxq0);
                qUC.copy<RunOn::Device>(qUC_host, boxqC);
                qUS.copy<RunOn::Device>(qUS_host, boxqS);
            }
        }

    }

    // Finally, construct the boundary conditions using the
    // complete multipole moments, for all points on the
    // boundary that are held on this process.

    const Box& domain = parent->Geom(crse_level).Domain();
    const auto dx = parent->Geom(crse_level).CellSizeArray();
    const auto problo = parent->Geom(crse_level).ProbLoArray();
    int coord_type = parent->Geom(crse_level).Coord();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox();

        auto qL0_arr = qL0.array();
        auto qLC_arr = qLC.array();
        auto qLS_arr = qLS.array();
        auto phi_arr = phi[mfi].array();

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int* domlo = domain.loVect();
            const int* domhi = domain.hiVect();

            // If we're using this to construct boundary values, then only use
            // the outermost bin.

            int nlo = 0;
            if (boundary_only == 1) {
                nlo = npts-1;
            }

            Real rmax_cubed = multipole::rmax * multipole::rmax * multipole::rmax;

            Real x;
            if (i > domhi[0]) {
                x = problo[0] + (static_cast<Real>(i  )         ) * dx[0] - problem::center[0];
            }
            else if (i < domlo[0]) {
                x = problo[0] + (static_cast<Real>(i+1)         ) * dx[0] - problem::center[0];
            }
            else {
                x = problo[0] + (static_cast<Real>(i  ) + 0.5_rt) * dx[0] - problem::center[0];
            }

            x = x / multipole::rmax;

#if AMREX_SPACEDIM >= 2
            Real y;
            if (j > domhi[1]) {
                y = problo[1] + (static_cast<Real>(j  )         ) * dx[1] - problem::center[1];
            }
            else if (j < domlo[1]) {
                y = problo[1] + (static_cast<Real>(j+1)         ) * dx[1] - problem::center[1];
            }
            else {
                y = problo[1] + (static_cast<Real>(j  ) + 0.5_rt) * dx[1] - problem::center[1];
            }
#else
            Real y = 0.0_rt;
#endif

            y = y / multipole::rmax;

#if AMREX_SPACEDIM == 3
            Real z;
            if (k > domhi[2]) {
                z = problo[2] + (static_cast<Real>(k  )         ) * dx[2] - problem::center[2];
            }
            else if (k < domlo[2]) {
                z = problo[2] + (static_cast<Real>(k+1)         ) * dx[2] - problem::center[2];
            }
            else {
                z = problo[2] + (static_cast<Real>(k  ) + 0.5_rt) * dx[2] - problem::center[2];
            }
#else
            Real z = 0.0;
#endif

            z = z / multipole::rmax;

            // Only adjust ghost zones here

            if (i < domlo[0] || i > domhi[0]
#if AMREX_SPACEDIM >= 2
                || j < domlo[1] || j > domhi[1]
#endif
#if AMREX_SPACEDIM >= 3
                || k < domlo[2] || k > domhi[2]
#endif
                ) {

                // There are some cases where r == 0. This might occur, for example,
                // when we have symmetric BCs and our corner is at one edge.
                // In this case, we'll set phi to zero for safety, to avoid NaN issues.
                // These cells should not be accessed anyway during the gravity solve.

                Real r = std::sqrt(x * x + y * y + z * z);

                if (r < 1.0e-12_rt) {
                    phi_arr(i,j,k) = 0.0_rt;
                    return;
                }

                Real cosTheta{}, phiAngle{};
                if (AMREX_SPACEDIM == 3) {
                    cosTheta = z / r;
                    phiAngle = std::atan2(y, x);
                }
                else if (AMREX_SPACEDIM == 2 && coord_type == 1) {
                    cosTheta = y / r;
                    phiAngle = 0.0_rt;
                }

                phi_arr(i,j,k) = 0.0_rt;

                // Compute the potentials on the ghost cells.

                Real legPolyL, legPolyL1, legPolyL2;
                Real assocLegPolyLM, assocLegPolyLM1, assocLegPolyLM2;

                for (int n = nlo; n <= npts - 1; ++n) {

                    for (int l = 0; l <= gravity::lnum; ++l) {

                        calcLegPolyL(l, legPolyL, legPolyL1, legPolyL2, cosTheta);

                        Real r_U = std::pow(r, -l-1);

                        // Make sure we undo the volume scaling here.

                        phi_arr(i,j,k) += qL0_arr(l,0,n) * legPolyL * r_U * rmax_cubed;

                    }

                    for (int m = 1; m <= gravity::lnum; ++m) {
                        for (int l = 1; l <= gravity::lnum; ++l) {

                            if (m > l) {
                                continue;
                            }

                            calcAssocLegPolyLM(l, m, assocLegPolyLM, assocLegPolyLM1, assocLegPolyLM2, cosTheta);

                            Real r_U = std::pow(r, -l-1);

                            // Make sure we undo the volume scaling here.

                            phi_arr(i,j,k) += (qLC_arr(l,m,n) * std::cos(m * phiAngle) + qLS_arr(l,m,n) * std::sin(m * phiAngle)) *
                                              assocLegPolyLM * r_U * rmax_cubed;

                        }
                    }

                }

                phi_arr(i,j,k) = -C::Gconst * phi_arr(i,j,k) / multipole::rmax;
            }
        });
    }

    if (gravity::verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        amrex::Print() << "Gravity::fill_multipole_BCs() time = " << end << std::endl << std::endl;
#ifdef BL_LAZY
        });
#endif
    }

}

#if (AMREX_SPACEDIM == 3)
void
Gravity::fill_direct_sum_BCs(int crse_level, int fine_level, const Vector<MultiFab*>& Rhs, MultiFab& phi)
{
    BL_PROFILE("Gravity::fill_direct_sum_BCs()");

    BL_ASSERT(crse_level==0);

    const Real strt = ParallelDescriptor::second();

    const Geometry& crse_geom = parent->Geom(crse_level);

    // Storage arrays for the BCs.

    const int* domlo = crse_geom.Domain().loVect();
    const int* domhi = crse_geom.Domain().hiVect();

    const int loVectXY[3] = {domlo[0]-1, domlo[1]-1, 0         };
    const int hiVectXY[3] = {domhi[0]+1, domhi[1]+1, 0         };

    const int loVectXZ[3] = {domlo[0]-1, 0         , domlo[2]-1};
    const int hiVectXZ[3] = {domhi[0]+1, 0         , domhi[2]+1};

    const int loVectYZ[3] = {0         , domlo[1]-1, domlo[2]-1};
    const int hiVectYZ[3] = {0         , domhi[1]+1, domhi[1]+1};

    const int bc_lo[3] = {domlo[0]-1, domlo[1]-1, domlo[2]-1};
    const int bc_hi[3] = {domhi[0]+1, domhi[1]+1, domhi[2]+1};

    GpuArray<Real, 3> bc_dx;
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        bc_dx[n] = crse_geom.CellSizeArray()[n];
    }
    for (int n = AMREX_SPACEDIM; n < 3; ++n) {
        bc_dx[n] = 0.0_rt;
    }

    GpuArray<Real, 3> problo;
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        problo[n] = crse_geom.ProbLoArray()[n];
    }
    for (int n = AMREX_SPACEDIM; n < 3; ++n) {
        problo[n] = 0.0_rt;
    }

    GpuArray<Real, 3> probhi;
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        probhi[n] = crse_geom.ProbHiArray()[n];
    }
    for (int n = AMREX_SPACEDIM; n < 3; ++n) {
        probhi[n] = 0.0_rt;
    }

    IntVect smallEndXY( loVectXY );
    IntVect bigEndXY  ( hiVectXY );
    IntVect smallEndXZ( loVectXZ );
    IntVect bigEndXZ  ( hiVectXZ );
    IntVect smallEndYZ( loVectYZ );
    IntVect bigEndYZ  ( hiVectYZ );

    Box boxXY(smallEndXY, bigEndXY);
    Box boxXZ(smallEndXZ, bigEndXZ);
    Box boxYZ(smallEndYZ, bigEndYZ);

    const long nPtsXY = boxXY.numPts();
    const long nPtsXZ = boxXZ.numPts();
    const long nPtsYZ = boxYZ.numPts();

    FArrayBox bcXYLo(boxXY);
    FArrayBox bcXYHi(boxXY);
    FArrayBox bcXZLo(boxXZ);
    FArrayBox bcXZHi(boxXZ);
    FArrayBox bcYZLo(boxYZ);
    FArrayBox bcYZHi(boxYZ);

    bcXYLo.setVal<RunOn::Device>(0.0);
    bcXYHi.setVal<RunOn::Device>(0.0);
    bcXZLo.setVal<RunOn::Device>(0.0);
    bcXZHi.setVal<RunOn::Device>(0.0);
    bcYZLo.setVal<RunOn::Device>(0.0);
    bcYZHi.setVal<RunOn::Device>(0.0);

    // Loop through the grids and compute the individual contributions
    // to the BCs. The BC constructor is coded to only add to the
    // BCs, so it is safe to directly hand the arrays to them.

    int physbc_lo[3];
    int physbc_hi[3];

    for (int dir = 0; dir < 3; dir++) {
        physbc_lo[dir] = phys_bc->lo(dir);
        physbc_hi[dir] = phys_bc->hi(dir);
    }

    for (int lev = crse_level; lev <= fine_level; ++lev) {

        // Create a local copy of the RHS so that we can mask it.

        MultiFab source(Rhs[lev - crse_level]->boxArray(),
                        Rhs[lev - crse_level]->DistributionMap(),
                        1, 0);

        MultiFab::Copy(source, *Rhs[lev - crse_level], 0, 0, 1, 0);

        if (lev < fine_level) {
            const MultiFab& mask = dynamic_cast<Castro*>(&(parent->getLevel(lev+1)))->build_fine_mask();
            MultiFab::Multiply(source, mask, 0, 0, 1, 0);
        }

        const auto dx = parent->Geom(lev).CellSizeArray();

#ifdef _OPENMP
        int nthreads = omp_get_max_threads();
        Vector<std::unique_ptr<FArrayBox> > priv_bcXYLo(nthreads);
        Vector<std::unique_ptr<FArrayBox> > priv_bcXYHi(nthreads);
        Vector<std::unique_ptr<FArrayBox> > priv_bcXZLo(nthreads);
        Vector<std::unique_ptr<FArrayBox> > priv_bcXZHi(nthreads);
        Vector<std::unique_ptr<FArrayBox> > priv_bcYZLo(nthreads);
        Vector<std::unique_ptr<FArrayBox> > priv_bcYZHi(nthreads);
        for (int i=0; i<nthreads; i++) {
            priv_bcXYLo[i].reset(new FArrayBox(boxXY));
            priv_bcXYHi[i].reset(new FArrayBox(boxXY));
            priv_bcXZLo[i].reset(new FArrayBox(boxXZ));
            priv_bcXZHi[i].reset(new FArrayBox(boxXZ));
            priv_bcYZLo[i].reset(new FArrayBox(boxYZ));
            priv_bcYZHi[i].reset(new FArrayBox(boxYZ));
        }
#pragma omp parallel
#endif
        {
#ifdef _OPENMP
            int tid = omp_get_thread_num();
            priv_bcXYLo[tid]->setVal<RunOn::Gpu>(0.0);
            priv_bcXYHi[tid]->setVal<RunOn::Gpu>(0.0);
            priv_bcXZLo[tid]->setVal<RunOn::Gpu>(0.0);
            priv_bcXZHi[tid]->setVal<RunOn::Gpu>(0.0);
            priv_bcYZLo[tid]->setVal<RunOn::Gpu>(0.0);
            priv_bcYZHi[tid]->setVal<RunOn::Gpu>(0.0);
#endif
            for (MFIter mfi(source, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box bx = mfi.tilebox();

                const auto rho = source[mfi].array();
                const auto vol = (*volume[lev])[mfi].array();

                // Determine if we need to add contributions from any symmetric boundaries.

                GpuArray<bool, 3> doSymmetricAddLo {false};
                GpuArray<bool, 3> doSymmetricAddHi {false};
                bool doSymmetricAdd {false};

                for (int b = 0; b < 3; ++b) {
                    if (physbc_lo[b] == amrex::PhysBCType::symmetry) {
                        doSymmetricAddLo[b] = true;
                        doSymmetricAdd      = true;
                    }

                    if (physbc_hi[b] == amrex::PhysBCType::symmetry) {
                        doSymmetricAddHi[b] = true;
                        doSymmetricAdd      = true;
                    }
                }

#ifdef _OPENMP
                auto bcXYLo_arr = priv_bcXYLo[tid]->array();
                auto bcXYHi_arr = priv_bcXYHi[tid]->array();
                auto bcXZLo_arr = priv_bcXZLo[tid]->array();
                auto bcXZHi_arr = priv_bcXZHi[tid]->array();
                auto bcYZLo_arr = priv_bcYZLo[tid]->array();
                auto bcYZHi_arr = priv_bcYZHi[tid]->array();
#else
                auto bcXYLo_arr = bcXYLo.array();
                auto bcXYHi_arr = bcXYHi.array();
                auto bcXZLo_arr = bcXZLo.array();
                auto bcXZHi_arr = bcXZHi.array();
                auto bcYZLo_arr = bcYZLo.array();
                auto bcYZHi_arr = bcYZHi.array();
#endif

                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    GpuArray<Real, 3> loc, locb;
                    loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0];

#if AMREX_SPACEDIM >= 2
                    loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1];
#else
                    loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
                    loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2];
#else
                    loc[2] = 0.0_rt;
#endif

                    // Do xy interfaces first. Note that the boundary conditions
                    // on phi are expected to live directly on the interface.
                    // We also have to handle the domain corners correctly. We are
                    // assuming that bc_lo = domlo - 1 and bc_hi = domhi + 1, where
                    // domlo and domhi are the coarse domain extent.

                    for (int m = bc_lo[1]; m <= bc_hi[1]; ++m) {
                        if (m == bc_lo[1]) {
                            locb[1] = problo[1];
                        }
                        else if (m == bc_hi[1]) {
                            locb[1] = probhi[1];
                        }
                        else {
                            locb[1] = problo[1] + (static_cast<Real>(m) + 0.5_rt) * bc_dx[1];
                        }
                        Real dy2 = (loc[1] - locb[1]) * (loc[1] - locb[1]);

                        for (int l = bc_lo[0]; l <= bc_hi[0]; ++l) {
                            if (l == bc_lo[0]) {
                                locb[0] = problo[0];
                            }
                            else if (l == bc_hi[0]) {
                                locb[0] = probhi[1];
                            }
                            else {
                                locb[0] = problo[0] + (static_cast<Real>(l) + 0.5_rt) * bc_dx[0];
                            }
                            Real dx2 = (loc[0] - locb[0]) * (loc[0] - locb[0]);

                            locb[2] = problo[2];
                            Real dz2 = (loc[2] - locb[2]) * (loc[2] - locb[2]);

                            Real r = std::sqrt(dx2 + dy2 + dz2);

                            Real dbc = -C::Gconst * rho(i,j,k) * vol(i,j,k) / r;

                            // Now, add any contributions from mass that is hidden behind
                            // a symmetric boundary.

                            if (doSymmetricAdd) {

                                dbc += direct_sum_symmetric_add(loc, locb, problo, probhi,
                                                                rho(i,j,k), vol(i,j,k),
                                                                doSymmetricAddLo, doSymmetricAddHi);

                            }

                            Gpu::Atomic::Add(&bcXYLo_arr(l,m,0), dbc);

                            locb[2] = probhi[2];
                            dz2 = (loc[2] - locb[2]) * (loc[2] - locb[2]);

                            r = std::sqrt(dx2 + dy2 + dz2);

                            dbc = -C::Gconst * rho(i,j,k) * vol(i,j,k) / r;

                            if (doSymmetricAdd) {

                                dbc += direct_sum_symmetric_add(loc, locb, problo, probhi,
                                                                rho(i,j,k), vol(i,j,k),
                                                                doSymmetricAddLo, doSymmetricAddHi);

                            }

                            Gpu::Atomic::Add(&bcXYHi_arr(l,m,0), dbc);

                        }

                    }

                    // Now do xz interfaces.

                    for (int n = bc_lo[2]; n <= bc_hi[2]; ++n) {
                        if (n == bc_lo[2]) {
                            locb[2] = problo[2];
                        }
                        else if (n == bc_hi[2]) {
                            locb[2] = probhi[2];
                        }
                        else {
                            locb[2] = problo[2] + (static_cast<Real>(n) + 0.5_rt) * bc_dx[2];
                        }
                        Real dz2 = (loc[2] - locb[2]) * (loc[2] - locb[2]);

                        for (int l = bc_lo[0]; l <= bc_hi[0]; ++l) {
                            if (l == bc_lo[0]) {
                                locb[0] = problo[0];
                            }
                            else if (l == bc_hi[0]) {
                                locb[0] = probhi[0];
                            }
                            else {
                                locb[0] = problo[0] + (static_cast<Real>(l) + 0.5_rt) * bc_dx[0];
                            }
                            Real dx2 = (loc[0] - locb[0]) * (loc[0] - locb[0]);

                            locb[1] = problo[1];
                            Real dy2 = (loc[1] - locb[1]) * (loc[1] - locb[1]);

                            Real r = std::sqrt(dx2 + dy2 + dz2);

                            Real dbc = -C::Gconst * rho(i,j,k) * vol(i,j,k) / r;

                            if (doSymmetricAdd) {

                                dbc += direct_sum_symmetric_add(loc, locb, problo, probhi,
                                                                rho(i,j,k), vol(i,j,k),
                                                                doSymmetricAddLo, doSymmetricAddHi);

                            }

                            Gpu::Atomic::Add(&bcXZLo_arr(l,0,n), dbc);

                            locb[1] = probhi[1];
                            dy2 = (loc[1] - locb[1]) * (loc[1] - locb[1]);

                            r = std::sqrt(dx2 + dy2 + dz2);

                            dbc = -C::Gconst * rho(i,j,k) * vol(i,j,k) / r;

                            if (doSymmetricAdd) {

                                dbc += direct_sum_symmetric_add(loc, locb, problo, probhi,
                                                                rho(i,j,k), vol(i,j,k),
                                                                doSymmetricAddLo, doSymmetricAddHi);

                            }

                            Gpu::Atomic::Add(&bcXZHi_arr(l,0,n), dbc);

                        }

                    }

                    // Finally, do yz interfaces.

                    for (int n = bc_lo[2]; n <= bc_hi[2]; ++n) {
                        if (n == bc_lo[2]) {
                            locb[2] = problo[2];
                        }
                        else if (n == bc_hi[2]) {
                            locb[2] = probhi[2];
                        }
                        else {
                            locb[2] = problo[2] + (static_cast<Real>(n) + 0.5_rt) * bc_dx[2];
                        }
                        Real dz2 = (loc[2] - locb[2]) * (loc[2] - locb[2]);

                        for (int m = bc_lo[1]; m <= bc_hi[1]; ++m) {
                            if (m == bc_lo[1]) {
                                locb[1] = problo[1];
                            }
                            else if (m == bc_hi[1]) {
                                locb[1] = probhi[1];
                            }
                            else {
                                locb[1] = problo[1] + (static_cast<Real>(m) + 0.5_rt) * bc_dx[1];
                            }
                            Real dy2 = (loc[1] - locb[1]) * (loc[1] - locb[1]);

                            locb[0] = problo[0];
                            Real dx2 = (loc[0] - locb[0]) * (loc[0] - locb[0]);

                            Real r = std::sqrt(dx2 + dy2 + dz2);

                            Real dbc = -C::Gconst * rho(i,j,k) * vol(i,j,k) / r;

                            if (doSymmetricAdd) {

                                dbc += direct_sum_symmetric_add(loc, locb, problo, probhi,
                                                                rho(i,j,k), vol(i,j,k),
                                                                doSymmetricAddLo, doSymmetricAddHi);

                            }

                            Gpu::Atomic::Add(&bcYZLo_arr(0,m,n), dbc);

                            locb[0] = probhi[0];
                            dx2 = (loc[0] - locb[0]) * (loc[0] - locb[0]);

                            r = std::sqrt(dx2 + dy2 + dz2);

                            dbc = -C::Gconst * rho(i,j,k) * vol(i,j,k) / r;

                            if (doSymmetricAdd) {

                                dbc += direct_sum_symmetric_add(loc, locb, problo, probhi,
                                                                rho(i,j,k), vol(i,j,k),
                                                                doSymmetricAddLo, doSymmetricAddHi);

                            }

                            Gpu::Atomic::Add(&bcYZHi_arr(0,m,n), dbc);

                        }

                    }

                });

            }

#ifdef _OPENMP
            Real* pXYLo = bcXYLo.dataPtr();
            Real* pXYHi = bcXYHi.dataPtr();
            Real* pXZLo = bcXZLo.dataPtr();
            Real* pXZHi = bcXZHi.dataPtr();
            Real* pYZLo = bcYZLo.dataPtr();
            Real* pYZHi = bcYZHi.dataPtr();
#pragma omp barrier
#pragma omp for nowait
            for (int i=0; i<nPtsXY; i++) {
                for (int it=0; it<nthreads; it++) {
                    const Real* pl = priv_bcXYLo[it]->dataPtr();
                    const Real* ph = priv_bcXYHi[it]->dataPtr();
                    pXYLo[i] += pl[i];
                    pXYHi[i] += ph[i];
                }
            }
#pragma omp for nowait
            for (int i=0; i<nPtsXZ; i++) {
                for (int it=0; it<nthreads; it++) {
                    const Real* pl = priv_bcXZLo[it]->dataPtr();
                    const Real* ph = priv_bcXZHi[it]->dataPtr();
                    pXZLo[i] += pl[i];
                    pXZHi[i] += ph[i];
                }
            }
#pragma omp for nowait
            for (int i=0; i<nPtsYZ; i++) {
                for (int it=0; it<nthreads; it++) {
                    const Real* pl = priv_bcYZLo[it]->dataPtr();
                    const Real* ph = priv_bcYZHi[it]->dataPtr();
                    pYZLo[i] += pl[i];
                    pYZHi[i] += ph[i];
                }
            }
#endif
        }

    } // end loop over levels

    // because the number of elements in mpi_reduce is int
    BL_ASSERT(nPtsXY <= std::numeric_limits<int>::max());
    BL_ASSERT(nPtsXZ <= std::numeric_limits<int>::max());
    BL_ASSERT(nPtsYZ <= std::numeric_limits<int>::max());

    ParallelDescriptor::ReduceRealSum(bcXYLo.dataPtr(), static_cast<int>(nPtsXY));
    ParallelDescriptor::ReduceRealSum(bcXYHi.dataPtr(), static_cast<int>(nPtsXY));
    ParallelDescriptor::ReduceRealSum(bcXZLo.dataPtr(), static_cast<int>(nPtsXZ));
    ParallelDescriptor::ReduceRealSum(bcXZHi.dataPtr(), static_cast<int>(nPtsXZ));
    ParallelDescriptor::ReduceRealSum(bcYZLo.dataPtr(), static_cast<int>(nPtsYZ));
    ParallelDescriptor::ReduceRealSum(bcYZHi.dataPtr(), static_cast<int>(nPtsYZ));

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx= mfi.growntilebox();

        auto p = phi[mfi].array();

        auto bcXYLo_arr = bcXYLo.array();
        auto bcXYHi_arr = bcXYHi.array();
        auto bcXZLo_arr = bcXZLo.array();
        auto bcXZHi_arr = bcXZHi.array();
        auto bcYZLo_arr = bcYZLo.array();
        auto bcYZHi_arr = bcYZHi.array();

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (i == bc_lo[0]) {
                p(i,j,k) = bcYZLo_arr(0,j,k);
            }

            if (i == bc_hi[0]) {
                p(i,j,k) = bcYZHi_arr(0,j,k);
            }

            if (j == bc_lo[1]) {
                p(i,j,k) = bcXZLo_arr(i,0,k);
            }

            if (j == bc_hi[1]) {
                p(i,j,k) = bcXZHi_arr(i,0,k);
            }

            if (k == bc_lo[2]) {
                p(i,j,k) = bcXYLo_arr(i,j,0);
            }

            if (k == bc_hi[2]) {
                p(i,j,k) = bcXYHi_arr(i,j,0);
            }
        });
    }

    if (gravity::verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        amrex::Print() << "Gravity::fill_direct_sum_BCs() time = " << end << std::endl << std::endl;
#ifdef BL_LAZY
        });
#endif
    }

}
#endif

#if (AMREX_SPACEDIM < 3)
void
Gravity::applyMetricTerms(int level, MultiFab& Rhs, const Vector<MultiFab*>& coeffs) const
{
    BL_PROFILE("Gravity::applyMetricTerms()");

    auto dx = parent->Geom(level).CellSizeArray();
    int coord_type = parent->Geom(level).Coord();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(Rhs, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(1);
        const Box& xbx = mfi.nodaltilebox(0);
#if AMREX_SPACEDIM >= 2
        const Box& ybx = mfi.nodaltilebox(1);
#endif

        // Modify Rhs and coeffs with the appropriate metric terms.
        apply_metric(bx,
                     Rhs.array(mfi), (Rhs[mfi]).box(),
                     (*coeffs[0]).array(mfi), xbx,
#if AMREX_SPACEDIM >= 2
                     (*coeffs[1]).array(mfi), ybx,
#endif
                     dx, coord_type);
    }
}

#endif

void
Gravity::make_mg_bc ()
{
    const Geometry& geom = parent->Geom(0);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (geom.isPeriodic(idim)) {
            mlmg_lobc[idim] = MLLinOp::BCType::Periodic;
            mlmg_hibc[idim] = MLLinOp::BCType::Periodic;
        } else {
            if (phys_bc->lo(idim) == amrex::PhysBCType::symmetry) {
                mlmg_lobc[idim] = MLLinOp::BCType::Neumann;
            } else {
                mlmg_lobc[idim] = MLLinOp::BCType::Dirichlet;
            }
            if (phys_bc->hi(idim) == amrex::PhysBCType::symmetry) {
                mlmg_hibc[idim] = MLLinOp::BCType::Neumann;
            } else {
                mlmg_hibc[idim] = MLLinOp::BCType::Dirichlet;
            }
        }
    }

    // Set Neumann bc at r=0.
    if (geom.IsSPHERICAL() || geom.IsRZ() ) {
        mlmg_lobc[0] = MLLinOp::BCType::Neumann;
    }
}

void
Gravity::set_mass_offset (Real time, bool multi_level) const
{
    BL_PROFILE("Gravity::set_mass_offset()");

    const Geometry& geom = parent->Geom(0);

    if (!geom.isAllPeriodic()) {
        mass_offset = 0.0;
    }
    else {
        Real old_mass_offset = mass_offset;
        mass_offset = 0.0;

        if (multi_level)
        {
            for (int lev = 0; lev <= parent->finestLevel(); lev++) {
                auto* cs = dynamic_cast<Castro*>(&parent->getLevel(lev));
                if (cs != nullptr) {
                    mass_offset += cs->volWgtSum("density", time);
                } else {
                    amrex::Abort("unable to access volWgtSum");
                }
            }
        }
        else
        {
            auto* cs = dynamic_cast<Castro*>(&parent->getLevel(0));
            if (cs != nullptr) {
                mass_offset = cs->volWgtSum("density", time, false, false);  // do not mask off fine grids
            } else {
                amrex::Abort("unable to access volWgtSum");
            }
        }

        mass_offset = mass_offset / geom.ProbSize();
        if (gravity::verbose > 1) {
            amrex::Print() << "Defining average density to be " << mass_offset << std::endl;
        }

        Real diff = std::abs(mass_offset - old_mass_offset);
        Real eps = 1.e-10 * std::abs(old_mass_offset);
        if (diff > eps && old_mass_offset > 0)
        {
            if (ParallelDescriptor::IOProcessor())
            {
                std::cout << " ... new vs old mass_offset " << mass_offset << " " << old_mass_offset
                          << " ... diff is " << diff <<  std::endl;
                std::cout << " ... Gravity::set_mass_offset -- total mass has changed!" << std::endl;;
            }
        }
    }
}

void
Gravity::add_pointmass_to_gravity (int level, MultiFab& phi, MultiFab& grav_vector) const
{
    BL_PROFILE("Gravity::add_pointmass_to_gravity()");

    const auto dx     = parent->Geom(level).CellSizeArray();
    const auto problo = parent->Geom(level).ProbLoArray();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(grav_vector, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox();

        Array4<Real> const grav_arr = grav_vector.array(mfi);
        Array4<Real> const phi_arr = phi.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Compute radial gravity due to a point mass at center[:].

            Real x = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0];
#if AMREX_SPACEDIM >= 2
            Real y = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1];
#else
            Real y = 0.0;
#endif
#if AMREX_SPACEDIM == 3
            Real z = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - problem::center[2];
#else
            Real z = 0.0;
#endif


            if(castro::point_mass_offset_is_true == 1)
            {
                Real star_radius = castro::point_mass_location_offset;

                if(AMREX_SPACEDIM == 1)
                {
                    x += star_radius;
                }
                else if(AMREX_SPACEDIM ==2)
                {
                    y += star_radius;
                }
                else if(AMREX_SPACEDIM == 3)
                {
                    z += star_radius;
                }

            }

            Real rsq = x * x + y * y + z * z;
            Real radial_force = -C::Gconst * castro::point_mass / rsq;

            Real rinv = 1.e0_rt / std::sqrt(rsq);

            // Note that grav may have more ghost zones than
            // phi, so we need to check that we're doing
            // valid indexing here.

            if (phi_arr.contains(i,j,k)) {
                phi_arr(i,j,k) -= C::Gconst * castro::point_mass * rinv;
            }

            grav_arr(i,j,k,0) += radial_force * (x * rinv);
            grav_arr(i,j,k,1) += radial_force * (y * rinv);
            grav_arr(i,j,k,2) += radial_force * (z * rinv);
        });
    }
}

void
Gravity::make_radial_gravity(int level, Real time, RealVector& radial_grav)
{
    BL_PROFILE("Gravity::make_radial_gravity()");

    const Real strt = ParallelDescriptor::second();

    // This is just here in case we need to debug ...
    int do_diag = 0;

    Real sum_over_levels = 0.;

    for (int lev = 0; lev <= level; lev++)
    {
        const Real t_old = LevelData[lev]->get_state_data(State_Type).prevTime();
        const Real t_new = LevelData[lev]->get_state_data(State_Type).curTime();
        const Real eps   = (t_new - t_old) * 1.e-6;

        // Create MultiFab with NUM_STATE components and no ghost cells
        MultiFab S(grids[lev],dmap[lev],NUM_STATE,0);

        if ( eps == 0.0 ) {  // NOLINT(bugprone-branch-clone,-warnings-as-errors)
            // Old and new time are identical; this should only happen if
            // dt is smaller than roundoff compared to the current time,
            // in which case we're probably in trouble anyway,
            // but we will still handle it gracefully here.
            MultiFab::Copy(S, LevelData[lev]->get_new_data(State_Type), 0, 0, NUM_STATE, 0);
        }
        else if ( std::abs(time-t_old) < eps)
        {
            MultiFab::Copy(S, LevelData[lev]->get_old_data(State_Type), 0, 0, NUM_STATE, 0);
        }
        else if ( std::abs(time-t_new) < eps)
        {
            MultiFab::Copy(S, LevelData[lev]->get_new_data(State_Type), 0, 0, NUM_STATE, 0);
        }
        else if (time > t_old && time < t_new)
        {
            Real alpha   = (time - t_old)/(t_new - t_old);
            Real omalpha = 1.0 - alpha;

            MultiFab::Copy(S, LevelData[lev]->get_old_data(State_Type), 0, 0, NUM_STATE, 0);
            S.mult(omalpha);

            MultiFab S_new(grids[lev],dmap[lev],NUM_STATE,0);
            MultiFab::Copy(S_new, LevelData[lev]->get_new_data(State_Type), 0, 0, NUM_STATE, 0);
            S_new.mult(alpha);

            S.plus(S_new,0,NUM_STATE,0);
        }
        else
        {
            std::cout << " Level / Time in make_radial_gravity is: " << lev << " " << time  << std::endl;
            std::cout << " but old / new time      are: " << t_old << " " << t_new << std::endl;
            amrex::Abort("Problem in Gravity::make_radial_gravity");
        }

        if (lev < level)
        {
            auto* fine_level = dynamic_cast<Castro*>(&(parent->getLevel(lev+1)));
        if (fine_level != nullptr) {
        const MultiFab& mask = fine_level->build_fine_mask();
        for (int n = 0; n < NUM_STATE; ++n) {
            MultiFab::Multiply(S, mask, 0, n, 1, 0);
        }
        } else {
                amrex::Abort("unable to create mask");
            }
        }

        int n1d = static_cast<int>(radial_mass[lev].size());

#ifdef GR_GRAV
        Real* const lev_pres = radial_pres[lev].dataPtr();
#endif
        Real* const lev_vol = radial_vol[lev].dataPtr();
        Real* const lev_mass = radial_mass[lev].dataPtr();

        amrex::ParallelFor(n1d,
        [=] AMREX_GPU_DEVICE (int i) noexcept
        {
#ifdef GR_GRAV
            lev_pres[i] = 0.;
#endif
            lev_vol[i] = 0.;
            lev_mass[i] = 0.;
        });

#ifdef _OPENMP
        int nthreads = omp_get_max_threads();
#ifdef GR_GRAV
        Vector< RealVector > priv_radial_pres(nthreads);
#endif
        Vector< RealVector > priv_radial_mass(nthreads);
        Vector< RealVector > priv_radial_vol (nthreads);
        for (int i=0; i<nthreads; i++) {
#ifdef GR_GRAV
            priv_radial_pres[i].resize(n1d,0.0);
#endif
            priv_radial_mass[i].resize(n1d,0.0);
            priv_radial_vol [i].resize(n1d,0.0);
        }
#pragma omp parallel
#endif
        {
#ifdef _OPENMP
            int tid = omp_get_thread_num();
#endif
            for (MFIter mfi(S, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                FArrayBox& fab = S[mfi];

                compute_radial_mass(bx,
                                    fab.array(),
#ifdef _OPENMP
                                    priv_radial_mass[tid],
                                    priv_radial_vol[tid],
#ifdef GR_GRAV
                                    priv_radial_pres[tid],
#endif
#else
                                    radial_mass[lev],
                                    radial_vol[lev],
#ifdef GR_GRAV
                                    radial_pres[lev],
#endif
#endif
                                    n1d, lev);
            }

#ifdef _OPENMP
#pragma omp barrier
#pragma omp for
            for (int i=0; i<n1d; i++) {
                for (int it=0; it<nthreads; it++) {
#ifdef GR_GRAV
                    radial_pres[lev][i] += priv_radial_pres[it][i];
#endif
                    radial_mass[lev][i] += priv_radial_mass[it][i];
                    radial_vol [lev][i] += priv_radial_vol [it][i];
                }
            }
#endif
        }

        if (!ParallelDescriptor::UseGpuAwareMpi()) {
            Gpu::prefetchToHost(radial_mass[lev].begin(), radial_mass[lev].end());
            Gpu::prefetchToHost(radial_vol[lev].begin(), radial_vol[lev].end());
#ifdef GR_GRAV
            Gpu::prefetchToHost(radial_pres[lev].begin(), radial_pres[lev].end());
#endif
        }

        ParallelDescriptor::ReduceRealSum(radial_mass[lev].dataPtr() ,n1d);
        ParallelDescriptor::ReduceRealSum(radial_vol[lev].dataPtr()  ,n1d);
#ifdef GR_GRAV
        ParallelDescriptor::ReduceRealSum(radial_pres[lev].dataPtr()  ,n1d);
#endif

        if (!ParallelDescriptor::UseGpuAwareMpi()) {
            Gpu::prefetchToDevice(radial_mass[lev].begin(), radial_mass[lev].end());
            Gpu::prefetchToDevice(radial_vol[lev].begin(), radial_vol[lev].end());
#ifdef GR_GRAV
            Gpu::prefetchToDevice(radial_pres[lev].begin(), radial_pres[lev].end());
#endif
        }

        if (do_diag > 0)
        {
            ReduceOps<ReduceOpSum> reduce_op;
            ReduceData<Real> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;

            reduce_op.eval(n1d, reduce_data,
            [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
            {
                return {lev_mass[i]};
            });

            ReduceTuple hv = reduce_data.value();
            Real sum = amrex::get<0>(hv);

            sum_over_levels += sum;
        }
    }

    if (do_diag > 0) {
        amrex::Print() << "Gravity::make_radial_gravity: Sum of mass over all levels " << sum_over_levels << std::endl;
    }

    int n1d = static_cast<int>(radial_mass[level].size());
    RealVector radial_mass_summed(n1d,0);

    const Real* const level_mass = radial_mass[level].dataPtr();
    Real* const mass_summed = radial_mass_summed.dataPtr();

    // First add the contribution from this level
    amrex::ParallelFor(n1d,
    [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        mass_summed[i] = level_mass[i];
    });

    // Now add the contribution from coarser levels
    if (level > 0)
    {
        int ratio = parent->refRatio(level-1)[0];
        for (int lev = level-1; lev >= 0; lev--)
        {
            if (lev < level-1) {
                ratio *= parent->refRatio(lev)[0];
            }

            Real* const lev_mass = radial_mass[lev].dataPtr();

            amrex::ParallelFor(n1d/ratio,
            [=] AMREX_GPU_DEVICE (int i) noexcept
            {
                for (int n = 0; n < ratio; n++)
                {
                    mass_summed[ratio*i+n] += 1./double(ratio) * lev_mass[i];
                }
            });
        }
    }

    if (do_diag > 0 && ParallelDescriptor::IOProcessor())
    {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        reduce_op.eval(n1d, reduce_data,
        [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
        {
            return {mass_summed[i]};
        });

        ReduceTuple hv = reduce_data.value();
        Real sum_added = amrex::get<0>(hv);

        std::cout << "Gravity::make_radial_gravity: Sum of combined mass " << sum_added << std::endl;
    }

    const Geometry& geom = parent->Geom(level);
    const Real* dx = geom.CellSize();
    Real dr        = dx[0] / static_cast<Real>(gravity::drdxfac);

    RealVector radial_vol_summed(n1d,0);
    RealVector radial_den_summed(n1d,0);

    Real* const vol_summed = radial_vol_summed.dataPtr();
    Real* const den_summed = radial_den_summed.dataPtr();

    const Real* const level_vol = radial_vol[level].dataPtr();

    // First add the contribution from this level
    amrex::ParallelFor(n1d,
    [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        vol_summed[i] = level_vol[i];
    });

    // Now add the contribution from coarser levels
    if (level > 0)
    {
        int ratio = parent->refRatio(level-1)[0];
        for (int lev = level-1; lev >= 0; lev--)
        {
            if (lev < level-1) {
                ratio *= parent->refRatio(lev)[0];
            }

            const Real* lev_vol = radial_vol[lev].dataPtr();

            amrex::ParallelFor(n1d/ratio,
            [=] AMREX_GPU_DEVICE (int i) noexcept
            {
                for (int n = 0; n < ratio; n++)
                {
                    vol_summed[ratio*i+n]  += 1./double(ratio) * lev_vol[i];
                }
            });
        }
    }

    amrex::ParallelFor(n1d,
    [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        den_summed[i] = mass_summed[i];
        if (vol_summed[i] > 0.) {
            den_summed[i] /= vol_summed[i];
        }
    });

#ifdef GR_GRAV
    RealVector radial_pres_summed(n1d,0);

    Real* const pres_summed = radial_pres_summed.dataPtr();
    Real* const level_pres = radial_pres[level].dataPtr();

    // First add the contribution from this level
    amrex::ParallelFor(n1d,
    [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        pres_summed[i] = level_pres[i];
    });

    // Now add the contribution from coarser levels
    if (level > 0)
    {
        int ratio = parent->refRatio(level-1)[0];
        for (int lev = level-1; lev >= 0; lev--)
        {
            if (lev < level-1) ratio *= parent->refRatio(lev)[0];

            const Real* lev_pres = radial_pres[lev].dataPtr();

            amrex::ParallelFor(n1d/ratio,
            [=] AMREX_GPU_DEVICE (int i) noexcept
            {
                for (int n = 0; n < ratio; n++) {
                    pres_summed[ratio*i+n] += 1./double(ratio) * lev_pres[i];
                }
            });
        }
    }

    amrex::ParallelFor(n1d,
    [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        if (vol_summed[i] > 0.) {
            pres_summed[i] /= vol_summed[i];
        }
    });
#endif

    // Integrate radially outward to define the gravity

    const Real* const mass = radial_mass_summed.dataPtr();
#ifdef GR_GRAV
    Real* const den = radial_den_summed.dataPtr();
    Real* const pres = radial_pres_summed.dataPtr();
#endif
    Real* const grav = radial_grav.dataPtr();

    // At a given radius r corresponding to an index i, the gravity is
    // g(r) = -G * M(r) / r**2. The enclosed mass can be computed via
    // an inclusive prefix sum, which yields the same result one would
    // obtain with a serial calculation from 0 to i but with the ability
    // to be parallelized. The first lambda returns the mass at radial
    // zone index i (which includes the upper shell of zone i-1 and the
    // lower shell of zone i) and the second lambda accepts the current
    // enclosed mass as an argument and computes the gravity at the
    // corresponding radius.

    Scan::PrefixSum<Real> (n1d,
        [=] AMREX_GPU_DEVICE (int i) -> Real
        {
            Real dM = 0.0;

            if (i > 0) {
                // The mass at (i-1) is distributed into an upper and lower shell; the
                // contribution to the mass at zone center i is from the upper shell.
                Real rlo = (static_cast<Real>(i-1)         ) * dr;
                Real rc  = (static_cast<Real>(i-1) + 0.5_rt) * dr;
                Real rhi = (static_cast<Real>(i-1) + 1.0_rt) * dr;

                Real vol_shell = (4.0_rt / 3.0_rt * M_PI) * dr * (rhi * rhi * rhi - rc  * rc  * rc);
                Real vol_zone  = (4.0_rt / 3.0_rt * M_PI) * dr * (rhi * rhi * rhi - rlo * rlo * rlo);
                dM = dM + (vol_shell / vol_zone) * mass[i-1];
            }

            Real rlo = (static_cast<Real>(i)         ) * dr;
            Real rc  = (static_cast<Real>(i) + 0.5_rt) * dr;
            Real rhi = (static_cast<Real>(i) + 1.0_rt) * dr;

            // The mass at (i) is distributed into an upper and lower shell; the
            // contribution to the mass at zone center i is from the lower shell.
            Real vol_shell = (4.0_rt / 3.0_rt * M_PI) * dr * (rc  * rc  * rc  - rlo * rlo * rlo);
            Real vol_zone  = (4.0_rt / 3.0_rt * M_PI) * dr * (rhi * rhi * rhi - rlo * rlo * rlo);
            dM = dM + (vol_shell / vol_zone) * mass[i];

            return dM;
        },
        [=] AMREX_GPU_DEVICE (int i, Real const& mass_encl_local)
        {
            Real rc = (static_cast<Real>(i) + 0.5_rt) * dr;

            grav[i] = -C::Gconst * mass_encl_local / (rc * rc);

#ifdef GR_GRAV
            // Tolman-Oppenheimer-Volkoff (TOV) post-Newtonian correction

            if (den[i] > 0.0_rt) {
                Real ga = (1.0_rt + pres[i] / (den[i] * C::c_light * C::c_light));
                Real gb = (1.0_rt + (4.0_rt * M_PI) * rc * rc * rc * pres[i] / (mass_encl_local * C::c_light * C::c_light));
                Real gc = 1.0_rt / (1.0_rt - 2.0_rt * C::Gconst * mass_encl_local / (rc * C::c_light * C::c_light));

                grav[i] = grav[i] * ga * gb * gc;
            }
#endif
        },
        Scan::Type::inclusive);

    if (gravity::verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        amrex::Print() << "Gravity::make_radial_gravity() time = " << end << std::endl << std::endl;
#ifdef BL_LAZY
        });
#endif
    }
}

Vector<std::unique_ptr<MultiFab> >
Gravity::get_rhs (int crse_level, int nlevs, int is_new)
{
    Vector<std::unique_ptr<MultiFab> > g_rhs(nlevs);

    for (int ilev = 0; ilev < nlevs; ++ilev)
    {
        int amr_lev = ilev + crse_level;
        g_rhs[ilev] = std::make_unique<MultiFab>(grids[amr_lev],dmap[amr_lev],1,0);
        MultiFab& state = (is_new == 1) ?
            LevelData[amr_lev]->get_new_data(State_Type) :
            LevelData[amr_lev]->get_old_data(State_Type);
        MultiFab::Copy(*g_rhs[ilev], state, URHO, 0,1,0);
    }
    return g_rhs;
}

void
Gravity::sanity_check (int level)
{
    // This is a sanity check on whether we are trying to fill multipole boundary conditiosn
    //  for grids at this level > 0 -- this case is not currently supported.
    //  Here we shrink the domain at this level by 1 in any direction which is not symmetry or periodic,
    //  then ask if the grids at this level are contained in the shrunken domain.  If not, then grids
    //  at this level touch the domain boundary and we must abort.

    const Geometry& geom = parent->Geom(level);

    if (level > 0  && !geom.isAllPeriodic())
    {
        Box shrunk_domain(geom.Domain());
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
        {
            if (!geom.isPeriodic(dir))
            {
                if (phys_bc->lo(dir) != amrex::PhysBCType::symmetry) {
                    shrunk_domain.growLo(dir,-1);
                }
                if (phys_bc->hi(dir) != amrex::PhysBCType::symmetry) {
                    shrunk_domain.growHi(dir,-1);
                }
            }
        }
        if (!shrunk_domain.contains(grids[level].minimalBox())) {
            amrex::Error("Oops -- don't know how to set boundary conditions for grids at this level that touch the domain boundary!");
        }
    }
}

void
Gravity::update_max_rhs()
{
    BL_PROFILE("Gravity::update_max_rhs()");

    // Calculate the maximum value of the RHS over all levels.
    // This should only be called at a synchronization point where
    // all Castro levels have valid new time data at the same simulation time.
    // The RHS we will use is the density multiplied by 4*pi*G and also
    // multiplied by the metric terms, just as it would be in a real solve.

    int crse_level = 0;
    int nlevs = parent->finestLevel() + 1;
    int is_new = 1;

    const auto& g_rhs = get_rhs(crse_level, nlevs, is_new);

    const Geometry& geom0 = parent->Geom(0);

#if (AMREX_SPACEDIM == 3)
    if ( geom0.isAllPeriodic() )
    {
        for (int lev = 0; lev < nlevs; ++lev) {
            g_rhs[lev]->plus(-mass_offset,0,1,0);
        }
    }
#endif

    for (int lev = 0; lev < nlevs; ++lev)
    {
        g_rhs[lev]->mult(Ggravity);
    }

#if (AMREX_SPACEDIM < 3)
    if (geom0.IsSPHERICAL() || geom0.IsRZ() )
    {
        Vector<Vector<std::unique_ptr<MultiFab> > > coeffs(nlevs);

        for (int lev = 0; lev < nlevs; ++lev) {

            // We need to include this bit about the coefficients because
            // it's required by applyMetricTerms.

            coeffs[lev].resize(AMREX_SPACEDIM);

            for (int i = 0; i < AMREX_SPACEDIM ; i++) {
                coeffs[lev][i] = std::make_unique<MultiFab>(amrex::convert(grids[lev],
                                                                           IntVect::TheDimensionVector(i)),
                                                            dmap[lev], 1, 0);

                coeffs[lev][i]->setVal(1.0);
            }

            applyMetricTerms(lev, *g_rhs[lev], amrex::GetVecOfPtrs(coeffs[lev]));
        }
    }
#endif

    max_rhs = 0.0;

    for (int lev = 0; lev < nlevs; ++lev)
        max_rhs = std::max(max_rhs, g_rhs[lev]->max(0));

}

Real
Gravity::solve_phi_with_mlmg (int crse_level, int fine_level,
                              const Vector<MultiFab*>& phi,
                              const Vector<MultiFab*>& g_rhs,
                              const Vector<Vector<MultiFab*> >& grad_phi,
                              const Vector<MultiFab*>& res,
                              Real time)
{
    BL_PROFILE("Gravity::solve_phi_with_mlmg()");

    int nlevs = fine_level-crse_level+1;

    if (crse_level == 0 && !(parent->Geom(0).isAllPeriodic()))
    {
        if (gravity::verbose > 1) {
            amrex::Print() << " ... Making bc's for phi at level 0\n";
        }

#if (AMREX_SPACEDIM == 3)
        if ( gravity::direct_sum_bcs ) {
            fill_direct_sum_BCs(crse_level, fine_level, g_rhs, *phi[0]);
        } else {
            fill_multipole_BCs(crse_level, fine_level, g_rhs, *phi[0]);
        }
#elif (AMREX_SPACEDIM == 2)
        fill_multipole_BCs(crse_level, fine_level, g_rhs, *phi[0]);
#else
        fill_multipole_BCs(crse_level, fine_level, g_rhs, *phi[0]);
#endif
    }

    for (int ilev = 0; ilev < nlevs; ++ilev)
    {
        g_rhs[ilev]->mult(Ggravity);
    }

    MultiFab CPhi;
    const MultiFab* crse_bcdata = nullptr;
    if (crse_level > 0)
    {
        GetCrsePhi(crse_level, CPhi, time);
        crse_bcdata = &CPhi;
    }

    Real rel_eps = rel_tol[fine_level];

    // The absolute tolerance is determined by the error tolerance
    // chosen by the user (tol) multiplied by the maximum value of
    // the RHS (4 * pi * G * rho). If we're doing periodic BCs, we
    // subtract off the mass_offset corresponding to the average
    // density on the domain. This will automatically be zero for
    // non-periodic BCs. And this also accounts for the metric
    // terms that are applied in non-Cartesian coordinates.

    Real abs_eps = abs_tol[fine_level] * max_rhs;

    Vector<const MultiFab*> crhs{g_rhs.begin(), g_rhs.end()};
    Vector<std::array<MultiFab*,AMREX_SPACEDIM> > gp;
    for (const auto& x : grad_phi) {
        gp.push_back({AMREX_D_DECL(x[0],x[1],x[2])});
    }

    return actual_solve_with_mlmg(crse_level, fine_level, phi, crhs, gp, res,
                                  crse_bcdata, rel_eps, abs_eps);
}

void
Gravity::solve_for_delta_phi(int crse_level, int fine_level,
                             const Vector<MultiFab*>& g_rhs,
                             const Vector<MultiFab*>& delta_phi,
                             const Vector<Vector<MultiFab*> >&  grad_delta_phi)
{
    BL_PROFILE("Gravity::solve_for_delta_phi");

    BL_ASSERT(grad_delta_phi.size() == fine_level - crse_level + 1);
    BL_ASSERT(delta_phi.size() == fine_level - crse_level + 1);

    if (gravity::verbose > 1 && ParallelDescriptor::IOProcessor()) {
      std::cout << "... solving for delta_phi at crse_level = " << crse_level << std::endl;
      std::cout << "...                    up to fine_level = " << fine_level << std::endl;
    }

    Vector<const MultiFab*> crhs{g_rhs.begin(), g_rhs.end()};
    Vector<std::array<MultiFab*,AMREX_SPACEDIM> > gp;
    for (const auto& x : grad_delta_phi) {
        gp.push_back({AMREX_D_DECL(x[0],x[1],x[2])});
    }

    Real rel_eps = 0.0;
    Real abs_eps = *(std::max_element(level_solver_resnorm.begin() + crse_level,
                                      level_solver_resnorm.begin() + fine_level+1));

    actual_solve_with_mlmg(crse_level, fine_level, delta_phi, crhs, gp, {},
                           nullptr, rel_eps, abs_eps);
}

Real
Gravity::actual_solve_with_mlmg (int crse_level, int fine_level,
                                 const amrex::Vector<amrex::MultiFab*>& phi,
                                 const amrex::Vector<const amrex::MultiFab*>& g_rhs,
                                 const amrex::Vector<std::array<amrex::MultiFab*,AMREX_SPACEDIM> >& grad_phi,
                                 const amrex::Vector<amrex::MultiFab*>& res,
                                 const amrex::MultiFab* const crse_bcdata,
                                 amrex::Real rel_eps, amrex::Real abs_eps) const
{
    BL_PROFILE("Gravity::actual_solve_with_mlmg()");

    Real final_resnorm = -1.0;

    int nlevs = fine_level-crse_level+1;

    Vector<Geometry> gmv;
    Vector<BoxArray> bav;
    Vector<DistributionMapping> dmv;
    for (int ilev = 0; ilev < nlevs; ++ilev)
    {
        gmv.push_back(parent->Geom(ilev+crse_level));
        bav.push_back(g_rhs[ilev]->boxArray());
        dmv.push_back(g_rhs[ilev]->DistributionMap());
    }

    LPInfo info;
    info.setAgglomeration(gravity::mlmg_agglomeration);
    info.setConsolidation(gravity::mlmg_consolidation);

    MLPoisson mlpoisson(gmv, bav, dmv, info);

    // BC
    mlpoisson.setDomainBC(mlmg_lobc, mlmg_hibc);
    if (mlpoisson.needsCoarseDataForBC())
    {
        mlpoisson.setCoarseFineBC(crse_bcdata, parent->refRatio(crse_level-1)[0]);
    }

    for (int ilev = 0; ilev < nlevs; ++ilev)
    {
        mlpoisson.setLevelBC(ilev, phi[ilev]);
    }

    MLMG mlmg(mlpoisson);
    mlmg.setVerbose(gravity::verbose - 1); // With normal verbosity we don't want MLMG information
    if (crse_level == 0) {
        mlmg.setMaxFmgIter(gravity::mlmg_max_fmg_iter);
    } else {
        mlmg.setMaxFmgIter(0); // Vcycle
    }

    AMREX_ALWAYS_ASSERT( !grad_phi.empty() or !res.empty() );
    AMREX_ALWAYS_ASSERT(  grad_phi.empty() or  res.empty() );

    if (!grad_phi.empty())
    {
        if (!gmv[0].isAllPeriodic()) {
            mlmg.setConvergenceNormType(MLMGNormType::bnorm);
        }

        mlmg.setNSolve(gravity::mlmg_nsolve);
        final_resnorm = mlmg.solve(phi, g_rhs, rel_eps, abs_eps);

        mlmg.getGradSolution(grad_phi);
    }
    else if (!res.empty())
    {
        mlmg.compResidual(res, phi, g_rhs);
    }

    return final_resnorm;
}
