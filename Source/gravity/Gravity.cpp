#include <cmath>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_ParmParse.H>
#include "Gravity.H"
#include "Castro.H"
#include <Gravity_F.H>
#include <Castro_F.H>

#include <AMReX_FillPatchUtil.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>

#define MAX_LEV 30

#include "gravity_defaults.H"
#include "fundamental_constants.H"

using namespace amrex;

#ifdef AMREX_DEBUG
int Gravity::test_solves  = 1;
#else
int Gravity::test_solves  = 0;
#endif
Real Gravity::max_radius_all_in_domain =  0.0;
Real Gravity::mass_offset    =  0.0;

// ************************************************************************************** //

// Ggravity is defined as 4 * pi * G, where G is the gravitational constant.

// In CGS, this constant is currently
//      Gconst   =  6.67428e-8           cm^3/g/s^2 , which results in
//      Ggravity =  83.8503442814844e-8  cm^3/g/s^2

// ************************************************************************************** //

static Real Ggravity = 0.;

Gravity::Gravity(Amr* Parent, int _finest_level, BCRec* _phys_bc, int _Density)
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
     Density = _Density;
     read_params();
     finest_level_allocated = -1;

     radial_grav_old.resize(MAX_LEV);
     radial_grav_new.resize(MAX_LEV);
     radial_mass.resize(MAX_LEV);
     radial_vol.resize(MAX_LEV);
#ifdef GR_GRAV
     radial_pres.resize(MAX_LEV);
#endif

     if (gravity_type == "PoissonGrav") make_mg_bc();
#if (BL_SPACEDIM > 1)
     if (gravity_type == "PoissonGrav") init_multipole_grav();
#endif
     max_rhs = 0.0;
}

Gravity::~Gravity() {}

void
Gravity::read_params ()
{
    static bool done = false;

    if (!done)
    {
        const Geometry& dgeom = DefaultGeometry();

        ParmParse pp("gravity");

#include "gravity_queries.H"

        if ( (gravity_type != "ConstantGrav") &&
             (gravity_type != "PoissonGrav") &&
             (gravity_type != "MonopoleGrav") &&
             (gravity_type != "PrescribedGrav") )
             {
                std::cout << "Sorry -- dont know this gravity type"  << std::endl;
                amrex::Abort("Options are ConstantGrav, PoissonGrav, MonopoleGrav, or PrescribedGrav");
             }

        if (  gravity_type == "ConstantGrav")
        {
          if ( dgeom.IsSPHERICAL() )
              amrex::Abort("Can't use constant direction gravity with non-Cartesian coordinates");
        }

#if (BL_SPACEDIM == 1)
        if (gravity_type == "PoissonGrav")
        {
          amrex::Abort(" gravity_type = PoissonGrav doesn't work well in 1-d -- please set gravity_type = MonopoleGrav");
        }
        else if (gravity_type == "MonopoleGrav" && !(dgeom.IsSPHERICAL()))
        {
          amrex::Abort("Only use MonopoleGrav in 1D spherical coordinates");
        }
        else if (gravity_type == "ConstantGrav" && dgeom.IsSPHERICAL())
        {
          amrex::Abort("Can't use constant gravity in 1D spherical coordinates");
        }

#elif (BL_SPACEDIM == 2)
        if (gravity_type == "MonopoleGrav" && dgeom.IsCartesian() )
        {
          amrex::Abort(" gravity_type = MonopoleGrav doesn't make sense in 2D Cartesian coordinates");
        }
#endif

        if (pp.contains("get_g_from_phi") && !get_g_from_phi && gravity_type == "PoissonGrav")
          if (ParallelDescriptor::IOProcessor())
            std::cout << "Warning: gravity_type = PoissonGrav assumes get_g_from_phi is true" << std::endl;

        int nlevs = parent->maxLevel() + 1;

        // Allow run-time input of solver tolerance. If the user provides no value, set a reasonable default
        // value on the coarse level, and then increase it by ref_ratio**2 as the levels get finer to account
        // for the change in the absolute scale of the Laplacian. If the user provides one value, use that
        // on the coarse level, and increase it the same way for the fine levels. If the user provides more than
        // one value, we expect them to provide one for every level, and we do not apply the ref_ratio effect.

        int n_abs_tol = pp.countval("abs_tol");

        if (n_abs_tol <= 1) {

            Real tol;

            if (n_abs_tol == 1) {

                pp.get("abs_tol", tol);

            } else {

                if (dgeom.IsCartesian())
                    tol = 1.e-11;
                else
                    tol = 1.e-10;

            }

            abs_tol[0] = tol;

            // Account for the fact that on finer levels, the scale of the
            // Laplacian changes due to the zone size changing. We assume
            // dx == dy == dz, so it is fair to say that on each level the
            // tolerance should increase by the factor ref_ratio**2, since
            // in absolute terms the Laplacian increases by that ratio too.
            // The actual tolerance we'll send in is the effective tolerance
            // on the finest level that we solve for.

            for (int lev = 1; lev < nlevs; ++lev)
                abs_tol[lev] = abs_tol[lev - 1] * std::pow(parent->refRatio(lev - 1)[0], 2);

        } else if (n_abs_tol >= nlevs) {

            pp.getarr("abs_tol", abs_tol, 0, nlevs);

        } else {

            amrex::Abort("If you are providing multiple values for abs_tol, you must provide at least one value for every level up to amr.max_level.");

        }

        // For the relative tolerance, we can again accept a single scalar (same for all levels)
        // or one for all levels. The default value is zero, so that we only use the absolute tolerance.
        // The multigrid always chooses the looser of the two criteria in determining whether the solve
        // has converged.

        // Note that the parameter rel_tol used to be known as ml_tol, so if we detect that the user has
        // set ml_tol but not rel_tol, we'll accept that for specifying the relative tolerance. ml_tol
        // is now considered deprecated and will be removed in a future release.

        std::string rel_tol_name = "rel_tol";

        if (pp.contains("ml_tol")) {

            amrex::Warning("The gravity parameter ml_tol has been renamed rel_tol. ml_tol is now deprecated.");

            if (!pp.contains("rel_tol"))
                rel_tol_name = "ml_tol";

        }

        int n_rel_tol = pp.countval(rel_tol_name.c_str());

        if (n_rel_tol <= 1) {

            Real tol;

            if (n_rel_tol == 1) {

                pp.get(rel_tol_name.c_str(), tol);

            } else {

                tol = 0.0;

            }

            for (int lev = 0; lev < MAX_LEV; ++lev)
                rel_tol[lev] = tol;

        } else if (n_rel_tol >= nlevs) {

            pp.getarr(rel_tol_name.c_str(), rel_tol, 0, nlevs);

        } else {

            amrex::Abort("If you are providing multiple values for rel_tol, you must provide at least one value for every level up to amr.max_level.");

        }

        // Warn user about obsolete tolerance parameters.

        if (pp.contains("delta_tol"))
            amrex::Warning("The gravity parameter delta_tol is no longer used.");

        if (pp.contains("sl_tol"))
            amrex::Warning("The gravity parameter sl_tol is no longer used.");
        Ggravity = 4.0 * M_PI * Gconst;
        if (verbose > 1 && ParallelDescriptor::IOProcessor())
        {
           std::cout << "Getting Gconst from constants: " << Gconst << std::endl;
           std::cout << "Using " << Ggravity << " for 4 pi G in Gravity.cpp " << std::endl;
        }

        done = true;
    }
}

void
Gravity::output_job_info_params(std::ostream& jobInfoFile)
{
#include "gravity_job_info_tests.H"
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
    if (verbose > 1 && ParallelDescriptor::IOProcessor())
        std::cout << "Installing Gravity level " << level << '\n';

    LevelData[level] = level_data;

    volume[level] = &_volume;

    area[level] = _area;

    level_solver_resnorm[level] = 0.0;

    const Geometry& geom = level_data->Geom();

    if (gravity_type == "PoissonGrav") {

       const DistributionMapping& dm = level_data->DistributionMap();

       grad_phi_prev[level].resize(BL_SPACEDIM);
       for (int n=0; n<BL_SPACEDIM; ++n)
           grad_phi_prev[level][n].reset(new MultiFab(level_data->getEdgeBoxArray(n),dm,1,1));

       grad_phi_curr[level].resize(BL_SPACEDIM);
       for (int n=0; n<BL_SPACEDIM; ++n)
           grad_phi_curr[level][n].reset(new MultiFab(level_data->getEdgeBoxArray(n),dm,1,1));

    } else if (gravity_type == "MonopoleGrav") {

        if (!geom.isAllPeriodic())
        {
           int n1d = drdxfac*numpts_at_level;

           radial_grav_old[level].resize(n1d);
           radial_grav_new[level].resize(n1d);
           radial_mass[level].resize(n1d);
           radial_vol[level].resize(n1d);
#ifdef GR_GRAV
           radial_pres[level].resize(n1d);
#endif
        }

    }

    // Compute the maximum radius at which all the mass at that radius is in the domain,
    //   assuming that the "hi" side of the domain is away from the center.
#if (BL_SPACEDIM > 1)
    if (level == 0)
    {
        Real center[3];
        ca_get_center(center);
        Real x = geom.ProbHi(0) - center[0];
        Real y = geom.ProbHi(1) - center[1];
        max_radius_all_in_domain = std::min(x,y);
#if (BL_SPACEDIM == 3)
        Real z = geom.ProbHi(2) - center[2];
        max_radius_all_in_domain = std::min(max_radius_all_in_domain,z);
#endif
        if (verbose > 1 && ParallelDescriptor::IOProcessor())
            std::cout << "Maximum radius for which the mass is contained in the domain: "
                      << max_radius_all_in_domain << std::endl;
    }
#endif

    finest_level_allocated = level;
}

std::string Gravity::get_gravity_type()
{
  return gravity_type;
}

int Gravity::get_max_solve_level()
{
  return max_solve_level;
}

int Gravity::NoSync()
{
  return no_sync;
}

int Gravity::NoComposite()
{
  return no_composite;
}

int Gravity::DoCompositeCorrection()
{
  return do_composite_phi_correction;
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
  for (int n = 0; n < BL_SPACEDIM; n++)
    grad_phi_curr[level][n]->plus(*addend[n],0,1,0);
}

void
Gravity::swapTimeLevels (int level)
{
    BL_PROFILE("Gravity::swapTimeLevels()");
    
    if (gravity_type == "PoissonGrav") {
        for (int n=0; n < BL_SPACEDIM; n++) {
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

    if (verbose > 1 && ParallelDescriptor::IOProcessor())
        std::cout << " ... solve for phi at level " << level << std::endl;

    const Real strt = ParallelDescriptor::second();

    if (is_new == 0) sanity_check(level);

    Real time;
    if (is_new == 1) {
      time = LevelData[level]->get_state_data(PhiGrav_Type).curTime();
    } else {
      time = LevelData[level]->get_state_data(PhiGrav_Type).prevTime();
    }

    // If we are below the max_solve_level, do the Poisson solve.
    // Otherwise, interpolate using a fillpatch from max_solve_level.

    if (level <= max_solve_level) {

        Vector<MultiFab*> phi_p(1, &phi);

        const auto& rhs = get_rhs(level, 1, is_new);

        Vector< Vector<MultiFab*> > grad_phi_p(1);
        grad_phi_p[0].resize(BL_SPACEDIM);
        for (int i = 0; i < BL_SPACEDIM ; i++) {
            grad_phi_p[0][i] = grad_phi[i];
        }

        Vector<MultiFab*> res_null;

        level_solver_resnorm[level] = solve_phi_with_mlmg(level, level,
                                                          phi_p,
                                                          amrex::GetVecOfPtrs(rhs),
                                                          grad_phi_p,
                                                          res_null,
                                                          time);

    }
    else {

        LevelData[level]->FillCoarsePatch(phi, 0, time, PhiGrav_Type, 0, 1, 1);

    }

    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "Gravity::solve_for_phi() time = " << end << std::endl << std::endl;
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

    if (fine_level > max_solve_level) {
        return;
    } else {
        fine_level = std::min(fine_level, max_solve_level);
    }

    BL_ASSERT(parent->finestLevel()>crse_level);
    if (verbose > 1 && ParallelDescriptor::IOProcessor()) {
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
        delta_phi[lev - crse_level].reset(new MultiFab(grids[lev], dmap[lev], 1, 1));
        delta_phi[lev - crse_level]->setVal(0.0);
    }

    Vector< Vector<std::unique_ptr<MultiFab> > > ec_gdPhi(nlevs);

    for (int lev = crse_level; lev <= fine_level; ++lev) {
        ec_gdPhi[lev - crse_level].resize(BL_SPACEDIM);

        const DistributionMapping& dm = LevelData[lev]->DistributionMap();
        for (int n = 0; n < BL_SPACEDIM; ++n) {
            ec_gdPhi[lev - crse_level][n].reset(new MultiFab(LevelData[lev]->getEdgeBoxArray(n), dm, 1, 0));
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

    Vector<std::unique_ptr<MultiFab> > rhs(nlevs);

    for (int lev = crse_level; lev <= fine_level; ++lev) {
        rhs[lev - crse_level].reset( new MultiFab(LevelData[lev]->boxArray(), LevelData[lev]->DistributionMap(), 1, 0));
        MultiFab::Copy(*rhs[lev - crse_level], *dphi[lev - crse_level], 0, 0, 1, 0);
        rhs[lev - crse_level]->mult(1.0 / Ggravity);
        MultiFab::Add(*rhs[lev - crse_level], *drho[lev - crse_level], 0, 0, 1, 0);
    }

    // Construct the boundary conditions for the Poisson solve.

    if (crse_level == 0 && !crse_geom.isAllPeriodic()) {

        if (verbose > 1 && ParallelDescriptor::IOProcessor())
         std::cout << " ... Making bc's for delta_phi at crse_level 0"  << std::endl;

#if (BL_SPACEDIM == 3)
      if ( direct_sum_bcs )
          fill_direct_sum_BCs(crse_level,fine_level,amrex::GetVecOfPtrs(rhs),*delta_phi[crse_level]);
      else {
          if (lnum >= 0) {
              fill_multipole_BCs(crse_level,fine_level,amrex::GetVecOfPtrs(rhs),*delta_phi[crse_level]);
          } else {
              int fill_interior = 0;
              make_radial_phi(crse_level,*rhs[0],*delta_phi[crse_level],fill_interior);
          }
      }
#elif (BL_SPACEDIM == 2)
      if (lnum >= 0) {
          fill_multipole_BCs(crse_level,fine_level,amrex::GetVecOfPtrs(rhs),*delta_phi[crse_level]);
      } else {
          int fill_interior = 0;
          make_radial_phi(crse_level,*rhs[0],*delta_phi[crse_level],fill_interior);
      }
#else
      int fill_interior = 0;
      make_radial_phi(crse_level,*rhs[0],*delta_phi[crse_level],fill_interior);
#endif

    }

    // Restore the factor of (4 * pi * G) for the Poisson solve.
    for (int lev = crse_level; lev <= fine_level; ++lev)
        rhs[lev - crse_level]->mult(Ggravity);

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

        Real local_correction = rhs[0]->sum() / grids[crse_level].numPts();

        if (verbose > 1 && ParallelDescriptor::IOProcessor())
            std::cout << "WARNING: Adjusting RHS in gravity_sync solve by " << local_correction << '\n';

        for (int lev = fine_level; lev >= crse_level; --lev)
            rhs[lev-crse_level]->plus(-local_correction, 0, 1, 0);

    }

    // Do multi-level solve for delta_phi.

    solve_for_delta_phi(crse_level, fine_level,
                        amrex::GetVecOfPtrs(rhs),
                        amrex::GetVecOfPtrs(delta_phi),
                        amrex::GetVecOfVecOfPtrs(ec_gdPhi));

    // In the all-periodic case we enforce that delta_phi averages to zero.

    if (crse_geom.isAllPeriodic() && (grids[crse_level].numPts() == crse_domain.numPts()) ) {

        Real local_correction = delta_phi[0]->sum() / grids[crse_level].numPts();

        for (int lev = crse_level; lev <= fine_level; ++lev)
            delta_phi[lev - crse_level]->plus(-local_correction, 0, 1, 1);

    }

    // Add delta_phi to phi_new, and grad(delta_phi) to grad(delta_phi_curr) on each level.
    // Update the cell-centered gravity too.

    for (int lev = crse_level; lev <= fine_level; lev++) {

        LevelData[lev]->get_new_data(PhiGrav_Type).plus(*delta_phi[lev - crse_level], 0, 1, 0);

        for (int n = 0; n < BL_SPACEDIM; n++)
            grad_phi_curr[lev][n]->plus(*ec_gdPhi[lev - crse_level][n], 0, 1, 0);

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

    const Real t_old = LevelData[level-1]->get_state_data(PhiGrav_Type).prevTime();
    const Real t_new = LevelData[level-1]->get_state_data(PhiGrav_Type).curTime();
    Real alpha = (time - t_old)/(t_new - t_old);
    Real omalpha = 1.0 - alpha;

    MultiFab const& phi_old = LevelData[level-1]->get_old_data(PhiGrav_Type);
    MultiFab const& phi_new = LevelData[level-1]->get_new_data(PhiGrav_Type);

    phi_crse.clear();
    phi_crse.define(grids[level-1], dmap[level-1], 1, 1); // BUT NOTE we don't trust phi's ghost cells.

    MultiFab::LinComb(phi_crse,
                      alpha  , phi_new, 0,
                      omalpha, phi_old, 0,
                      0, 1, 1);

    const Geometry& geom = parent->Geom(level-1);
    phi_crse.FillBoundary(geom.periodicity());
}

void
Gravity::multilevel_solve_for_new_phi (int level, int finest_level_in, int use_previous_phi_as_guess)
{
    BL_PROFILE("Gravity::multilevel_solve_for_new_phi()");

    if (verbose > 1 && ParallelDescriptor::IOProcessor())
      std::cout << "... multilevel solve for new phi at base level " << level << " to finest level " << finest_level_in << std::endl;

    for (int lev = level; lev <= finest_level_in; lev++) {
       BL_ASSERT(grad_phi_curr[lev].size()==BL_SPACEDIM);
       for (int n=0; n<BL_SPACEDIM; ++n)
       {
           grad_phi_curr[lev][n].reset(new MultiFab(LevelData[lev]->getEdgeBoxArray(n),
                                                    LevelData[lev]->DistributionMap(),1,1));
       }
    }

    int is_new = 1;
    actual_multilevel_solve(level,finest_level_in,amrex::GetVecOfVecOfPtrs(grad_phi_curr),
                            is_new,use_previous_phi_as_guess);
}

void
Gravity::actual_multilevel_solve (int crse_level, int finest_level_in,
                                  const Vector<Vector<MultiFab*> >& grad_phi,
                                  int is_new,
                                  int use_previous_phi_as_guess)
{
    BL_PROFILE("Gravity::actual_multilevel_solve()");

    for (int ilev = crse_level; ilev <= finest_level_in ; ++ilev)
        sanity_check(ilev);

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

        if (!use_previous_phi_as_guess)
            phi_p[ilev]->setVal(0.);
    }

    const auto& rhs = get_rhs(crse_level, nlevels, is_new);

    if (!use_previous_phi_as_guess && crse_level == 0 && !(parent->Geom(0).isAllPeriodic()))
    {
        make_radial_phi(0, *rhs[0], *phi_p[0], 1);
    }

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

    int fine_level = std::min(finest_level_in, max_solve_level);

    if (fine_level >= crse_level) {

        Vector<MultiFab*> res_null;
        solve_phi_with_mlmg(crse_level, fine_level,
                            phi_p, amrex::GetVecOfPtrs(rhs), grad_phi_p, res_null,
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
        for (int amr_lev = fine_level; amr_lev > crse_level; amr_lev--)
            average_fine_ec_onto_crse_ec(amr_lev-1,is_new);

    }

    // For all levels on which we're not doing the solve, interpolate from
    // the coarsest level with correct data. Note that since FillCoarsePatch
    // fills from the coarse level just below it, we need to fill from the
    // lowest level upwards using successive interpolations.

    for (int amr_lev = max_solve_level+1; amr_lev <= finest_level_in; amr_lev++) {

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

        // We need to use a nodal interpolater.

        Interpolater* gp_interp = &node_bilinear_interp;

        // For the BCs, we will use the Gravity_Type BCs for convenience, but these will
        // not do anything because we do not fill on physical boundaries.

        const Vector<BCRec>& gp_bcs = LevelData[amr_lev]->get_desc_lst()[Gravity_Type].getBCs();

        for (int n = 0; n < BL_SPACEDIM; ++n) {
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

    if (level > max_solve_level) {

        LevelData[level]->FillCoarsePatch(grav_vector,0,time,Gravity_Type,0,3,ng);

        return;

    }

    // Note that grav_vector coming into this routine always has three components.
    // So we'll define a temporary MultiFab with BL_SPACEDIM dimensions.
    // Then at the end we'll copy in all BL_SPACEDIM dimensions from this into
    // the outgoing grav_vector, leaving any higher dimensions unchanged.

    MultiFab grav(grids[level], dmap[level], BL_SPACEDIM, ng);
    grav.setVal(0.0,ng);

    if (gravity_type == "ConstantGrav") {

       // Set to constant value in the BL_SPACEDIM direction and zero in all others.

       grav.setVal(const_grav,BL_SPACEDIM-1,1,ng);

    } else if (gravity_type == "MonopoleGrav") {

       const Real prev_time = LevelData[level]->get_state_data(State_Type).prevTime();
       make_radial_gravity(level,prev_time,radial_grav_old[level]);
       interpolate_monopole_grav(level,radial_grav_old[level],grav);

    } else if (gravity_type == "PrescribedGrav") {

        MultiFab& phi = LevelData[level]->get_old_data(PhiGrav_Type);
      make_prescribed_grav(level,time,grav,phi);

    } else if (gravity_type == "PoissonGrav") {

       const Geometry& geom = parent->Geom(level);
       amrex::average_face_to_cellcenter(grav, amrex::GetVecOfConstPtrs(grad_phi_prev[level]), geom);
       grav.mult(-1.0, ng); // g = - grad(phi)

    } else {
       amrex::Abort("Unknown gravity_type in get_old_grav_vector");
    }

    // Do the copy to the output vector.

    for (int dir = 0; dir < 3; dir++) {
        if (dir < BL_SPACEDIM) {
            MultiFab::Copy(grav_vector, grav, dir, dir, 1, ng);
        } else {
            grav_vector.setVal(0.,dir,1,ng);
        }
    }

#if (BL_SPACEDIM > 1)
    if (gravity_type != "ConstantGrav") {
        // Fill ghost cells
        AmrLevel* amrlev = &parent->getLevel(level) ;
        AmrLevel::FillPatch(*amrlev,grav_vector,ng,time,Gravity_Type,0,BL_SPACEDIM);
    }
#endif

    Castro* cs = dynamic_cast<Castro*>(&parent->getLevel(level));
    if (cs->using_point_mass()) {
        Real point_mass = cs->get_point_mass();
        MultiFab& phi = LevelData[level]->get_old_data(PhiGrav_Type);
        add_pointmass_to_gravity(level,phi,grav_vector,point_mass);
    }
}

void
Gravity::get_new_grav_vector(int level, MultiFab& grav_vector, Real time)
{
    BL_PROFILE("Gravity::get_new_grav_vector()");

    int ng = grav_vector.nGrow();

    // Fill data from the level below if we're not doing a solve on this level.

    if (level > max_solve_level) {

        LevelData[level]->FillCoarsePatch(grav_vector,0,time,Gravity_Type,0,3,ng);

        return;

    }

    // Note that grav_vector coming into this routine always has three components.
    // So we'll define a temporary MultiFab with BL_SPACEDIM dimensions.
    // Then at the end we'll copy in all BL_SPACEDIM dimensions from this into
    // the outgoing grav_vector, leaving any higher dimensions unchanged.

    MultiFab grav(grids[level],dmap[level],BL_SPACEDIM,ng);
    grav.setVal(0.0,ng);

    if (gravity_type == "ConstantGrav") {

       // Set to constant value in the BL_SPACEDIM direction
       grav.setVal(const_grav,BL_SPACEDIM-1,1,ng);

    } else if (gravity_type == "MonopoleGrav") {

        // We always fill radial_grav_new (at every level)
        const Real cur_time = LevelData[level]->get_state_data(State_Type).curTime();
        make_radial_gravity(level,cur_time,radial_grav_new[level]);
        interpolate_monopole_grav(level,radial_grav_new[level],grav);

    } else if (gravity_type == "PrescribedGrav") {

    MultiFab& phi = LevelData[level]->get_new_data(PhiGrav_Type);
    make_prescribed_grav(level,time,grav,phi);

    } else if (gravity_type == "PoissonGrav") {

        const Geometry& geom = parent->Geom(level);
        amrex::average_face_to_cellcenter(grav, amrex::GetVecOfConstPtrs(grad_phi_curr[level]), geom);
        grav.mult(-1.0, ng); // g = - grad(phi)

    } else {
       amrex::Abort("Unknown gravity_type in get_new_grav_vector");
    }

    // Do the copy to the output vector.

    for (int dir = 0; dir < 3; dir++) {
        if (dir < BL_SPACEDIM) {
            MultiFab::Copy(grav_vector, grav, dir, dir, 1, ng);
        } else {
            grav_vector.setVal(0.,dir,1,ng);
        }
    }

#if (BL_SPACEDIM > 1)
    if (gravity_type != "ConstantGrav" && ng>0) {
        // Fill ghost cells
        AmrLevel* amrlev = &parent->getLevel(level) ;
        AmrLevel::FillPatch(*amrlev,grav_vector,ng,time,Gravity_Type,0,BL_SPACEDIM);
    }
#endif

    Castro* cs = dynamic_cast<Castro*>(&parent->getLevel(level));
    if (cs->using_point_mass()) {
        Real point_mass = cs->get_point_mass();
        MultiFab& phi = LevelData[level]->get_new_data(PhiGrav_Type);
        add_pointmass_to_gravity(level,phi,grav_vector,point_mass);
    }
}

void
Gravity::test_level_grad_phi_prev(int level)
{
    BL_PROFILE("Gravity::test_level_grad_phi_prev()");

    // Fill the RHS for the solve
    MultiFab& S_old = LevelData[level]->get_old_data(State_Type);
    MultiFab Rhs(grids[level],dmap[level],1,0);
    MultiFab::Copy(Rhs,S_old, URHO,0,1,0);

    const Geometry& geom = parent->Geom(level);

    // This is a correction for fully periodic domains only
    if ( geom.isAllPeriodic() )
    {
       if (verbose > 1 && ParallelDescriptor::IOProcessor() && mass_offset != 0.0)
          std::cout << " ... subtracting average density from RHS at level ... "
                    << level << " " << mass_offset << std::endl;
       Rhs.plus(-mass_offset,0,1,0);
    }

    Rhs.mult(Ggravity);

    if (verbose > 1) {
       Real rhsnorm = Rhs.norm0();
       if (ParallelDescriptor::IOProcessor()) {
          std::cout << "... test_level_grad_phi_prev at level " << level << std::endl;
          std::cout << "       norm of RHS             " << rhsnorm << std::endl;
       }
    }

    const Real* dx     = parent->Geom(level).CellSize();
    const Real* problo = parent->Geom(level).ProbLo();
    const int coord_type     = geom.Coord();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(Rhs, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        // Test whether using the edge-based gradients
        //   to compute Div(Grad(Phi)) satisfies Lap(phi) = RHS
        // Fill the RHS array with the residual
        ca_test_residual(bx.loVect(), bx.hiVect(),
                         BL_TO_FORTRAN(Rhs[mfi]),
                         D_DECL(BL_TO_FORTRAN((*grad_phi_prev[level][0])[mfi]),
                                BL_TO_FORTRAN((*grad_phi_prev[level][1])[mfi]),
                                BL_TO_FORTRAN((*grad_phi_prev[level][2])[mfi])),
                         dx,problo,&coord_type);
    }
    if (verbose > 1) {
       Real resnorm = Rhs.norm0();
//     Real gppxnorm = grad_phi_prev[level][0]->norm0();
#if (BL_SPACEDIM > 1)
//     Real gppynorm = grad_phi_prev[level][1]->norm0();
#endif
#if (BL_SPACEDIM > 2)
//     Real gppznorm = grad_phi_prev[level][2]->norm0();
#endif
      if (ParallelDescriptor::IOProcessor())
        std::cout << "       norm of residual        " << resnorm << std::endl;
//      std::cout << "       norm of grad_phi_prev_x " << gppxnorm << std::endl;
#if (BL_SPACEDIM > 1)
//      std::cout << "       norm of grad_phi_prev_y " << gppynorm << std::endl;
#endif
#if (BL_SPACEDIM > 2)
//      std::cout << "       norm of grad_phi_prev_z " << gppznorm << std::endl;
#endif
    }
}

void
Gravity::test_level_grad_phi_curr(int level)
{
    BL_PROFILE("Gravity::test_level_grad_phi_curr()");

    // Fill the RHS for the solve
    MultiFab& S_new = LevelData[level]->get_new_data(State_Type);
    MultiFab Rhs(grids[level],dmap[level],1,0);
    MultiFab::Copy(Rhs,S_new, URHO, 0,1,0);

    const Geometry& geom = parent->Geom(level);

    // This is a correction for fully periodic domains only
    if ( geom.isAllPeriodic() )
    {
       if (verbose > 1 && ParallelDescriptor::IOProcessor() && mass_offset != 0.0)
          std::cout << " ... subtracting average density from RHS in solve ... " << mass_offset << std::endl;
       Rhs.plus(-mass_offset,0,1,0);
    }

    Rhs.mult(Ggravity);

    if (verbose > 1) {
       Real rhsnorm = Rhs.norm0();
       if (ParallelDescriptor::IOProcessor()) {
          std::cout << "... test_level_grad_phi_curr at level " << level << std::endl;
          std::cout << "       norm of RHS             " << rhsnorm << std::endl;
        }
    }

    const Real*     dx   = geom.CellSize();
    const Real* problo   = geom.ProbLo();
    const int coord_type = geom.Coord();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(Rhs, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        // Test whether using the edge-based gradients
        //   to compute Div(Grad(Phi)) satisfies Lap(phi) = RHS
        // Fill the RHS array with the residual
        ca_test_residual(bx.loVect(), bx.hiVect(),
                         BL_TO_FORTRAN(Rhs[mfi]),
                         D_DECL(BL_TO_FORTRAN((*grad_phi_curr[level][0])[mfi]),
                                BL_TO_FORTRAN((*grad_phi_curr[level][1])[mfi]),
                                BL_TO_FORTRAN((*grad_phi_curr[level][2])[mfi])),
                         dx,problo,&coord_type);
    }
    if (verbose > 1) {
       Real resnorm = Rhs.norm0();
//     Real gppxnorm = grad_phi_curr[level][0]->norm0();
#if (BL_SPACEDIM > 1)
//     Real gppynorm = grad_phi_curr[level][1]->norm0();
#endif
#if (BL_SPACEDIM > 2)
//     Real gppznorm = grad_phi_curr[level][2]->norm0();
#endif
       if (ParallelDescriptor::IOProcessor())
          std::cout << "       norm of residual        " << resnorm << std::endl;
//        std::cout << "       norm of grad_phi_curr_x " << gppxnorm << std::endl;
#if (BL_SPACEDIM > 1)
//        std::cout << "       norm of grad_phi_curr_y " << gppynorm << std::endl;
#endif
#if (BL_SPACEDIM > 2)
//        std::cout << "       norm of grad_phi_curr_z " << gppznorm << std::endl;
#endif
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

    if (verbose > 1 && ParallelDescriptor::IOProcessor()) {
        std::cout << "\n";
        std::cout << "... compute difference between level and composite solves at level " << level << "\n";
        std::cout << "\n";
    }

    comp_minus_level_phi.define(LevelData[level]->boxArray(),
                                LevelData[level]->DistributionMap(),
                                1, 0);

    comp_minus_level_phi.copy(comp_phi, 0, 0, 1);
    comp_minus_level_phi.minus(parent->getLevel(level).get_old_data(PhiGrav_Type), 0, 1, 0);

    comp_minus_level_grad_phi.resize(BL_SPACEDIM);
    for (int n = 0; n < BL_SPACEDIM; ++n) {
        comp_minus_level_grad_phi[n].reset(new MultiFab(LevelData[level]->getEdgeBoxArray(n),
                                                        LevelData[level]->DistributionMap(), 1, 0));
        comp_minus_level_grad_phi[n]->copy(*comp_gphi[n], 0, 0, 1);
        comp_minus_level_grad_phi[n]->minus(*grad_phi_prev[level][n], 0, 1, 0);
    }

}

void
Gravity::average_fine_ec_onto_crse_ec(int level, int is_new)
{
    BL_PROFILE("Gravity::average_fine_ec_onto_crse_ec()");

    // NOTE: this is called with level == the coarser of the two levels involved
    if (level == parent->finestLevel()) return;

    //
    // Coarsen() the fine stuff on processors owning the fine data.
    //
    BoxArray crse_gphi_fine_BA(grids[level+1].size());

    IntVect fine_ratio = parent->refRatio(level);

    for (int i = 0; i < crse_gphi_fine_BA.size(); ++i)
        crse_gphi_fine_BA.set(i,amrex::coarsen(grids[level+1][i],fine_ratio));

    Vector<std::unique_ptr<MultiFab> > crse_gphi_fine(BL_SPACEDIM);
    for (int n=0; n<BL_SPACEDIM; ++n)
    {
        BoxArray eba = crse_gphi_fine_BA;
        eba.surroundingNodes(n);
        crse_gphi_fine[n].reset(new MultiFab(eba,dmap[level+1],1,0));
    }

    auto& grad_phi = (is_new) ? grad_phi_curr : grad_phi_prev;

    amrex::average_down_faces(amrex::GetVecOfConstPtrs(grad_phi[level+1]),
                               amrex::GetVecOfPtrs(crse_gphi_fine),
                               fine_ratio);

    const Geometry& cgeom = parent->Geom(level);

    for (int n = 0; n < BL_SPACEDIM; ++n)
    {
        grad_phi[level][n]->copy(*crse_gphi_fine[n], cgeom.periodicity());
    }
}

void
Gravity::test_composite_phi (int crse_level)
{
    BL_PROFILE("Gravity::test_composite_phi()");

    if (verbose > 1 && ParallelDescriptor::IOProcessor()) {
        std::cout << "   " << '\n';
        std::cout << "... test_composite_phi at base level " << crse_level << '\n';
    }

    int finest_level_local = parent->finestLevel();
    int nlevels = finest_level_local - crse_level + 1;

    Vector<std::unique_ptr<MultiFab> > phi(nlevels);
    Vector<std::unique_ptr<MultiFab> > rhs(nlevels);
    Vector<std::unique_ptr<MultiFab> > res(nlevels);
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        int amr_lev = crse_level + ilev;

        phi[ilev].reset(new MultiFab(grids[amr_lev],dmap[amr_lev],1,1));
        MultiFab::Copy(*phi[ilev],
                       LevelData[amr_lev]->get_new_data(PhiGrav_Type),
                       0,0,1,1);

        rhs[ilev].reset(new MultiFab(grids[amr_lev],dmap[amr_lev],1,1));
        MultiFab::Copy(*rhs[ilev],
                       LevelData[amr_lev]->get_new_data(State_Type),
                       URHO, 0,1,0);

        res[ilev].reset(new MultiFab(grids[amr_lev],dmap[amr_lev],1,0));
        res[ilev]->setVal(0.);
    }

    Real time = LevelData[crse_level]->get_state_data(PhiGrav_Type).curTime();

    Vector< Vector<MultiFab*> > grad_phi_null;
    solve_phi_with_mlmg(crse_level, finest_level_local,
                        amrex::GetVecOfPtrs(phi),
                        amrex::GetVecOfPtrs(rhs),
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
        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "      ... norm of composite residual at level "
                      << amr_lev << "  " << resnorm << '\n';
        }
    }
    if (ParallelDescriptor::IOProcessor()) std::cout << std::endl;
}

void
Gravity::make_prescribed_grav(int level, Real time, MultiFab& grav_vector, MultiFab& phi)
{
    BL_PROFILE("Gravity::make_prescribed_grav()");
    
    const Real strt = ParallelDescriptor::second();

    const Geometry& geom = parent->Geom(level);
    const Real* dx   = geom.CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.growntilebox();
       ca_prescribe_phi(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                        BL_TO_FORTRAN_ANYD(phi[mfi]),dx);
    }

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(grav_vector, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.growntilebox();
       ca_prescribe_grav(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                         BL_TO_FORTRAN_ANYD(grav_vector[mfi]),dx);
    }

    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "Gravity::make_prescribed_grav() time = " << end << std::endl << std::endl;
#ifdef BL_LAZY
        });
#endif
    }
}

void
Gravity::interpolate_monopole_grav(int level, RealVector& radial_grav, MultiFab& grav_vector)
{
    BL_PROFILE("Gravity::interpolate_monopole_grav()");
    
    int n1d = radial_grav.size();

    const Geometry& geom = parent->Geom(level);
    const Real* dx = geom.CellSize();
    const Real dr        = dx[0] / double(drdxfac);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(grav_vector, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.growntilebox();
       ca_put_radial_grav(bx.loVect(),bx.hiVect(),dx,&dr,
                          BL_TO_FORTRAN(grav_vector[mfi]),
                          radial_grav.dataPtr(),geom.ProbLo(),
                          &n1d,&level);
    }
}

void
Gravity::make_radial_phi(int level, const MultiFab& Rhs, MultiFab& phi, int fill_interior)
{
    BL_PROFILE("Gravity::make_radial_phi()");

    BL_ASSERT(level==0);

    const Real strt = ParallelDescriptor::second();

    int n1d = drdxfac*numpts_at_level;

    RealVector radial_mass(n1d,0.0);
    RealVector radial_vol(n1d,0.0);
    RealVector radial_phi(n1d,0.0);
    RealVector radial_grav(n1d,0.0);

    const Geometry& geom = parent->Geom(level);
    const Real* dx   = geom.CellSize();
    Real dr = dx[0] / double(drdxfac);

    // Define total mass in each shell
    // Note that RHS = density (we have not yet multiplied by G)

#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
    Vector< RealVector > priv_radial_mass(nthreads);
    Vector< RealVector > priv_radial_vol (nthreads);
    for (int i=0; i<nthreads; i++) {
        priv_radial_mass[i].resize(n1d,0.0);
        priv_radial_vol [i].resize(n1d,0.0);
    }
#pragma omp parallel
#endif
    {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#endif
        for (MFIter mfi(Rhs, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
            ca_compute_radial_mass(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                                   AMREX_REAL_ANYD(dx), dr,
                                   BL_TO_FORTRAN_ANYD(Rhs[mfi]),
#ifdef _OPENMP
                                   priv_radial_mass[tid].dataPtr(),
                                   priv_radial_vol[tid].dataPtr(),
#else
                                   radial_mass.dataPtr(),
                                   radial_vol.dataPtr(),
#endif
                                   AMREX_REAL_ANYD(geom.ProbLo()),
                                   n1d, drdxfac, level);
        }

#ifdef _OPENMP
#pragma omp barrier
#pragma omp for
        for (int i=0; i<n1d; i++) {
            for (int it=0; it<nthreads; it++) {
                radial_mass[i] += priv_radial_mass[it][i];
                radial_vol [i] += priv_radial_vol [it][i];
            }
        }
#endif
    }

    ParallelDescriptor::ReduceRealSum(radial_mass.dataPtr(),n1d);

    RealVector radial_den(n1d, 0.0);

    for (int i = 0; i < n1d; ++i)
    {
        radial_den[i] = radial_mass[i];
        if (radial_vol[i] > 0.0) radial_den[i] /= radial_vol[i];
    }

    // Integrate radially outward to define the gravity
    ca_integrate_grav(radial_mass.dataPtr(), radial_den.dataPtr(),
                      radial_grav.dataPtr(), &max_radius_all_in_domain, &dr, &n1d);

    // Integrate radially inward to define the potential
    ca_integrate_phi(radial_mass.dataPtr(),radial_grav.dataPtr(),
                     radial_phi.dataPtr(),&dr,&n1d);

    Box domain(parent->Geom(level).Domain());
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox();

#pragma gpu box(bx)
        ca_put_radial_phi(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                          AMREX_INT_ANYD(domain.loVect()), AMREX_INT_ANYD(domain.hiVect()),
                          AMREX_REAL_ANYD(dx), dr,
                          BL_TO_FORTRAN_ANYD(phi[mfi]),
                          radial_phi.dataPtr(),
                          AMREX_REAL_ANYD(geom.ProbLo()),
                          n1d, fill_interior);
    }

    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "Gravity::make_radial_phi() time = " << end << std::endl << std::endl;
#ifdef BL_LAZY
        });
#endif
    }

}



#if (BL_SPACEDIM > 1)
void
Gravity::init_multipole_grav()
{

    int lo_bc[3];
    int hi_bc[3];

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      lo_bc[dir] = phys_bc->lo(dir);
      hi_bc[dir] = phys_bc->hi(dir);
    }
    for (int dir = BL_SPACEDIM; dir < 3; dir++)
    {
      lo_bc[dir] = -1;
      hi_bc[dir] = -1;
    }

    if (lnum >= 0) {
        init_multipole_gravity(&lnum, lo_bc, hi_bc);
    }
}

void
Gravity::fill_multipole_BCs(int crse_level, int fine_level, const Vector<MultiFab*>& Rhs, MultiFab& phi)
{
    BL_PROFILE("Gravity::fill_multipole_BCs()");

    // Multipole BCs only make sense to construct if we are starting from the coarse level.

    BL_ASSERT(crse_level == 0);

    BL_ASSERT(lnum >= 0);

    const Real strt = ParallelDescriptor::second();

#if (BL_SPACEDIM == 3)
    const int npts = numpts_at_level;
#else
    const int npts = 1;
#endif

    // Storage arrays for the multipole moments.
    // We will initialize them to zero, and then
    // sum up the results over grids.
    // Note that since Boxes are defined with
    // BL_SPACEDIM dimensions, we cannot presently
    // use this array to fill the interior of the
    // domain in 2D, since we can only have one
    // radial index for calculating the multipole moments.

    Box boxq0( IntVect(D_DECL(0, 0, 0)), IntVect(D_DECL(lnum, 0,    npts-1)) );
    Box boxqC( IntVect(D_DECL(0, 0, 0)), IntVect(D_DECL(lnum, lnum, npts-1)) );
    Box boxqS( IntVect(D_DECL(0, 0, 0)), IntVect(D_DECL(lnum, lnum, npts-1)) );

    FArrayBox qL0(boxq0);
    FArrayBox qLC(boxqC);
    FArrayBox qLS(boxqS);

    FArrayBox qU0(boxq0);
    FArrayBox qUC(boxqC);
    FArrayBox qUS(boxqS);

    Array4<Real> const& qL0_arr = qL0.array();
    Array4<Real> const& qU0_arr = qU0.array();

    amrex::ParallelFor(boxq0,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        qL0_arr(i,j,k) = 0.0;
        qU0_arr(i,j,k) = 0.0;
    });

    Array4<Real> const& qLC_arr = qLC.array();
    Array4<Real> const& qUC_arr = qUC.array();

    amrex::ParallelFor(boxqC,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        qLC_arr(i,j,k) = 0.0;
        qUC_arr(i,j,k) = 0.0;
    });

    Array4<Real> const& qLS_arr = qLS.array();
    Array4<Real> const& qUS_arr = qUS.array();

    amrex::ParallelFor(boxqS,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        qLS_arr(i,j,k) = 0.0;
        qUS_arr(i,j,k) = 0.0;
    });

    // This section needs to be generalized for computing
    // full multipole gravity, not just BCs. At present this
    // does nothing.

#if (BL_SPACEDIM == 3)
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
            const MultiFab& mask = dynamic_cast<Castro*>(&(parent->getLevel(lev+1)))->build_fine_mask();
            MultiFab::Multiply(source, mask, 0, 0, 1, 0);
        }

        // Loop through the grids and compute the individual contributions
        // to the various moments. The multipole moment constructor
        // is coded to only add to the moment arrays, so it is safe
        // to directly hand the arrays to them.

        const Box& domain = parent->Geom(lev).Domain();
        const Real* dx = parent->Geom(lev).CellSize();

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
            priv_qL0[tid]->setVal(0.0);
            priv_qLC[tid]->setVal(0.0);
            priv_qLS[tid]->setVal(0.0);
            priv_qU0[tid]->setVal(0.0);
            priv_qUC[tid]->setVal(0.0);
            priv_qUS[tid]->setVal(0.0);
#endif
            for (MFIter mfi(source, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
                ca_compute_multipole_moments(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                                             AMREX_INT_ANYD(domain.loVect()), AMREX_INT_ANYD(domain.hiVect()),
                                             AMREX_REAL_ANYD(dx), BL_TO_FORTRAN_ANYD(source[mfi]),
                                             BL_TO_FORTRAN_ANYD((*volume[lev])[mfi]),
                                             lnum,
#ifdef _OPENMP
                                             priv_qL0[tid]->dataPtr(),
                                             priv_qLC[tid]->dataPtr(),priv_qLS[tid]->dataPtr(),
                                             priv_qU0[tid]->dataPtr(),
                                             priv_qUC[tid]->dataPtr(),priv_qUS[tid]->dataPtr(),
#else
                                             qL0.dataPtr(), qLC.dataPtr(), qLS.dataPtr(),
                                             qU0.dataPtr(), qUC.dataPtr(), qUS.dataPtr(),
#endif
                                             npts, boundary_only);
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

    if (!ParallelDescriptor::UseGpuAwareMpi()) {
        qL0.prefetchToHost();
        qLC.prefetchToHost();
        qLS.prefetchToHost();
    }

    ParallelDescriptor::ReduceRealSum(qL0.dataPtr(),boxq0.numPts());
    ParallelDescriptor::ReduceRealSum(qLC.dataPtr(),boxqC.numPts());
    ParallelDescriptor::ReduceRealSum(qLS.dataPtr(),boxqS.numPts());

    if (!ParallelDescriptor::UseGpuAwareMpi()) {
        qL0.prefetchToDevice();
        qLC.prefetchToDevice();
        qLS.prefetchToDevice();
    }

    if (boundary_only != 1) {

      if (!ParallelDescriptor::UseGpuAwareMpi()) {
          qU0.prefetchToHost();
          qUC.prefetchToHost();
          qUS.prefetchToHost();
      }

      ParallelDescriptor::ReduceRealSum(qU0.dataPtr(),boxq0.numPts());
      ParallelDescriptor::ReduceRealSum(qUC.dataPtr(),boxqC.numPts());
      ParallelDescriptor::ReduceRealSum(qUS.dataPtr(),boxqS.numPts());

      if (!ParallelDescriptor::UseGpuAwareMpi()) {
          qU0.prefetchToDevice();
          qUC.prefetchToDevice();
          qUS.prefetchToDevice();
      }

    }

    // Finally, construct the boundary conditions using the
    // complete multipole moments, for all points on the
    // boundary that are held on this process.

    const Box& domain = parent->Geom(crse_level).Domain();
    const Real* dx = parent->Geom(crse_level).CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox();

#pragma gpu box(bx)
        ca_put_multipole_phi(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                             AMREX_INT_ANYD(domain.loVect()), AMREX_INT_ANYD(domain.hiVect()),
                             AMREX_REAL_ANYD(dx), BL_TO_FORTRAN_ANYD(phi[mfi]),
                             lnum,
                             qL0.dataPtr(), qLC.dataPtr(), qLS.dataPtr(),
                             qU0.dataPtr(), qUC.dataPtr(), qUS.dataPtr(),
                             npts, boundary_only);
    }

    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "Gravity::fill_multipole_BCs() time = " << end << std::endl << std::endl;
#ifdef BL_LAZY
        });
#endif
    }

}
#endif

#if (BL_SPACEDIM == 3)
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

    const Real* bc_dx = crse_geom.CellSize();

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

    Array4<Real> const& bcXYLo_arr = bcXYLo.array();
    Array4<Real> const& bcXYHi_arr = bcXYHi.array();

    amrex::ParallelFor(boxXY,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        bcXYLo_arr(i,j,k) = 0.0;
        bcXYHi_arr(i,j,k) = 0.0;
    });

    Array4<Real> const& bcXZLo_arr = bcXZLo.array();
    Array4<Real> const& bcXZHi_arr = bcXZHi.array();

    amrex::ParallelFor(boxXZ,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        bcXZLo_arr(i,j,k) = 0.0;
        bcXZHi_arr(i,j,k) = 0.0;
    });

    Array4<Real> const& bcYZLo_arr = bcYZLo.array();
    Array4<Real> const& bcYZHi_arr = bcYZHi.array();

    amrex::ParallelFor(boxYZ,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        bcYZLo_arr(i,j,k) = 0.0;
        bcYZHi_arr(i,j,k) = 0.0;
    });

    // Loop through the grids and compute the individual contributions
    // to the BCs. The BC constructor is coded to only add to the
    // BCs, so it is safe to directly hand the arrays to them.

    int physbc_lo[3];
    int physbc_hi[3];

    for (int dir = 0; dir < 3; dir++)
    {
      physbc_lo[dir] = phys_bc->lo(dir);
      physbc_lo[dir] = phys_bc->hi(dir);
    }

    int symmetry_type = Symmetry;

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

        const Real* dx = parent->Geom(lev).CellSize();

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
            priv_bcXYLo[tid]->setVal(0.0);
            priv_bcXYHi[tid]->setVal(0.0);
            priv_bcXZLo[tid]->setVal(0.0);
            priv_bcXZHi[tid]->setVal(0.0);
            priv_bcYZLo[tid]->setVal(0.0);
            priv_bcYZHi[tid]->setVal(0.0);
#endif
            for (MFIter mfi(source, TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box bx = mfi.tilebox();

                const FArrayBox& r = source[mfi];
                const FArrayBox& v = (*volume[lev])[mfi];

#pragma gpu box(bx)
                ca_compute_direct_sum_bc(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()), AMREX_REAL_ANYD(dx),
                                         symmetry_type, AMREX_INT_ANYD(physbc_lo), AMREX_INT_ANYD(physbc_hi),
                                         BL_TO_FORTRAN_ANYD(r),
                                         BL_TO_FORTRAN_ANYD(v),
                                         AMREX_REAL_ANYD(crse_geom.ProbLo()), AMREX_REAL_ANYD(crse_geom.ProbHi()),
#ifdef _OPENMP
                                         priv_bcXYLo[tid]->dataPtr(),
                                         priv_bcXYHi[tid]->dataPtr(),
                                         priv_bcXZLo[tid]->dataPtr(),
                                         priv_bcXZHi[tid]->dataPtr(),
                                         priv_bcYZLo[tid]->dataPtr(),
                                         priv_bcYZHi[tid]->dataPtr(),
#else
                                         bcXYLo.dataPtr(), bcXYHi.dataPtr(),
                                         bcXZLo.dataPtr(), bcXZHi.dataPtr(),
                                         bcYZLo.dataPtr(), bcYZHi.dataPtr(),
#endif
                                         AMREX_INT_ANYD(bc_lo), AMREX_INT_ANYD(bc_hi), AMREX_REAL_ANYD(bc_dx));
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

    // because the number of elments in mpi_reduce is int
    BL_ASSERT(nPtsXY <= std::numeric_limits<int>::max());
    BL_ASSERT(nPtsXZ <= std::numeric_limits<int>::max());
    BL_ASSERT(nPtsYZ <= std::numeric_limits<int>::max());

    ParallelDescriptor::ReduceRealSum(bcXYLo.dataPtr(), nPtsXY);
    ParallelDescriptor::ReduceRealSum(bcXYHi.dataPtr(), nPtsXY);
    ParallelDescriptor::ReduceRealSum(bcXZLo.dataPtr(), nPtsXZ);
    ParallelDescriptor::ReduceRealSum(bcXZHi.dataPtr(), nPtsXZ);
    ParallelDescriptor::ReduceRealSum(bcYZLo.dataPtr(), nPtsYZ);
    ParallelDescriptor::ReduceRealSum(bcYZHi.dataPtr(), nPtsYZ);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx= mfi.growntilebox();

        FArrayBox& p = phi[mfi];

#pragma gpu box(bx)
        ca_put_direct_sum_bc(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                             BL_TO_FORTRAN_ANYD(p),
                             bcXYLo.dataPtr(), bcXYHi.dataPtr(),
                             bcXZLo.dataPtr(), bcXZHi.dataPtr(),
                             bcYZLo.dataPtr(), bcYZHi.dataPtr(),
                             AMREX_INT_ANYD(bc_lo), AMREX_INT_ANYD(bc_hi));
    }

    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "Gravity::fill_direct_sum_BCs() time = " << end << std::endl << std::endl;
#ifdef BL_LAZY
        });
#endif
    }

}
#endif

#if (BL_SPACEDIM < 3)
void
Gravity::applyMetricTerms(int level, MultiFab& Rhs, const Vector<MultiFab*>& coeffs)
{
    BL_PROFILE("Gravity::applyMetricTerms()");
    
    const Real* dx = parent->Geom(level).CellSize();
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
#pragma gpu box(bx)
        ca_apply_metric(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                        AMREX_INT_ANYD(xbx.loVect()), AMREX_INT_ANYD(xbx.hiVect()),
#if AMREX_SPACEDIM >= 2
                        AMREX_INT_ANYD(ybx.loVect()), AMREX_INT_ANYD(ybx.hiVect()),
#endif
                        BL_TO_FORTRAN_ANYD(Rhs[mfi]),
                        BL_TO_FORTRAN_ANYD((*coeffs[0])[mfi]),
#if AMREX_SPACEDIM >= 2
                        BL_TO_FORTRAN_ANYD((*coeffs[1])[mfi]),
#endif
                        AMREX_REAL_ANYD(dx), coord_type);
    }
}

void
Gravity::unweight_cc(int level, MultiFab& cc)
{
    BL_PROFILE("Gravity::unweight_cc()");
    
    const Real* dx = parent->Geom(level).CellSize();
    const int coord_type = parent->Geom(level).Coord();
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(cc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
        ca_unweight_cc(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                       BL_TO_FORTRAN_ANYD(cc[mfi]),
                       AMREX_REAL_ANYD(dx), coord_type);
    }
}

void
Gravity::unweight_edges(int level, const Vector<MultiFab*>& edges)
{
    BL_PROFILE("Gravity::unweight_edges()");
    
    const Real* dx = parent->Geom(level).CellSize();
    const int coord_type = parent->Geom(level).Coord();
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (int idir=0; idir<BL_SPACEDIM; ++idir) {
        for (MFIter mfi(*edges[idir], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
            ca_unweight_edges(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                              BL_TO_FORTRAN_ANYD((*edges[idir])[mfi]),
                              AMREX_REAL_ANYD(dx),
                              coord_type, idir);
        }
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
            if (phys_bc->lo(idim) == Symmetry) {
                mlmg_lobc[idim] = MLLinOp::BCType::Neumann;
            } else {
                mlmg_lobc[idim] = MLLinOp::BCType::Dirichlet;
            }
            if (phys_bc->hi(idim) == Symmetry) {
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
Gravity::set_mass_offset (Real time, bool multi_level)
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
                Castro* cs = dynamic_cast<Castro*>(&parent->getLevel(lev));
                mass_offset += cs->volWgtSum("density", time);
            }
        }
        else
        {
            Castro* cs = dynamic_cast<Castro*>(&parent->getLevel(0));
            mass_offset = cs->volWgtSum("density", time, false, false);  // do not mask off fine grids
        }

        mass_offset = mass_offset / geom.ProbSize();
        if (verbose > 1 && ParallelDescriptor::IOProcessor())
            std::cout << "Defining average density to be " << mass_offset << std::endl;

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
Gravity::add_pointmass_to_gravity (int level, MultiFab& phi, MultiFab& grav_vector, Real point_mass)
{
    BL_PROFILE("Gravity::add_pointmass_to_gravity()");
    
    const Real* dx     = parent->Geom(level).CellSize();
    const Real* problo = parent->Geom(level).ProbLo();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(grav_vector, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox();

#pragma gpu box(bx)
        pm_add_to_grav(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                       point_mass, BL_TO_FORTRAN_ANYD(phi[mfi]),
                       BL_TO_FORTRAN_ANYD(grav_vector[mfi]),
                       AMREX_REAL_ANYD(problo), AMREX_REAL_ANYD(dx));
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

        if ( eps == 0.0 )
        {
            // Old and new time are identical; this should only happen if
            // dt is smaller than roundoff compared to the current time,
            // in which case we're probably in trouble anyway,
            // but we will still handle it gracefully here.
            S.copy(LevelData[lev]->get_new_data(State_Type),0,0,NUM_STATE);
        }
        else if ( std::abs(time-t_old) < eps)
        {
            S.copy(LevelData[lev]->get_old_data(State_Type),0,0,NUM_STATE);
        }
        else if ( std::abs(time-t_new) < eps)
        {
            S.copy(LevelData[lev]->get_new_data(State_Type),0,0,NUM_STATE);
        }
        else if (time > t_old && time < t_new)
        {
            Real alpha   = (time - t_old)/(t_new - t_old);
            Real omalpha = 1.0 - alpha;

            S.copy(LevelData[lev]->get_old_data(State_Type),0,0,NUM_STATE);
            S.mult(omalpha);

            MultiFab S_new(grids[lev],dmap[lev],NUM_STATE,0);
            S_new.copy(LevelData[lev]->get_new_data(State_Type),0,0,NUM_STATE);
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
            Castro* fine_level = dynamic_cast<Castro*>(&(parent->getLevel(lev+1)));
            const MultiFab& mask = fine_level->build_fine_mask();
            for (int n = 0; n < NUM_STATE; ++n)
                MultiFab::Multiply(S, mask, 0, n, 1, 0);
        }

        int n1d = radial_mass[lev].size();

#ifdef GR_GRAV
        for (int i = 0; i < n1d; i++) radial_pres[lev][i] = 0.;
#endif
        for (int i = 0; i < n1d; i++) radial_vol[lev][i] = 0.;
        for (int i = 0; i < n1d; i++) radial_mass[lev][i] = 0.;

        const Geometry& geom = parent->Geom(lev);
        const Real* dx   = geom.CellSize();
        Real dr = dx[0] / double(drdxfac);

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

#pragma gpu box(bx)
                ca_compute_radial_mass(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                                       AMREX_REAL_ANYD(dx), dr,
                                       BL_TO_FORTRAN_ANYD(fab),
#ifdef _OPENMP
                                       priv_radial_mass[tid].dataPtr(),
                                       priv_radial_vol[tid].dataPtr(),
#else
                                       radial_mass[lev].dataPtr(),
                                       radial_vol[lev].dataPtr(),
#endif
                                       AMREX_REAL_ANYD(geom.ProbLo()),
                                       n1d, drdxfac, lev);

#ifdef GR_GRAV
                ca_compute_avgpres(bx.loVect(), bx.hiVect(), dx, &dr,
                                   BL_TO_FORTRAN(fab),
#ifdef _OPENMP
                                   priv_radial_pres[tid].dataPtr(),
#else
                                   radial_pres[lev].dataPtr(),
#endif
                                   geom.ProbLo(),&n1d,&drdxfac,&lev);
#endif
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

        ParallelDescriptor::ReduceRealSum(radial_mass[lev].dataPtr() ,n1d);
        ParallelDescriptor::ReduceRealSum(radial_vol[lev].dataPtr()  ,n1d);
#ifdef GR_GRAV
        ParallelDescriptor::ReduceRealSum(radial_pres[lev].dataPtr()  ,n1d);
#endif

        if (do_diag > 0)
        {
            Real sum = 0.;
            for (int i = 0; i < n1d; i++) sum += radial_mass[lev][i];
            sum_over_levels += sum;
        }
    }

    if (do_diag > 0 && ParallelDescriptor::IOProcessor())
        std::cout << "Gravity::make_radial_gravity: Sum of mass over all levels " << sum_over_levels << std::endl;

    int n1d = radial_mass[level].size();
    RealVector radial_mass_summed(n1d,0);

    // First add the contribution from this level
    for (int i = 0; i < n1d; i++)
    {
        radial_mass_summed[i] = radial_mass[level][i];
    }

    // Now add the contribution from coarser levels
    if (level > 0)
    {
        int ratio = parent->refRatio(level-1)[0];
        for (int lev = level-1; lev >= 0; lev--)
        {
            if (lev < level-1) ratio *= parent->refRatio(lev)[0];
            for (int i = 0; i < n1d/ratio; i++)
            {
                for (int n = 0; n < ratio; n++)
                {
                   radial_mass_summed[ratio*i+n] += 1./double(ratio) * radial_mass[lev][i];
                }
            }
        }
    }

    if (do_diag > 0 && ParallelDescriptor::IOProcessor())
    {
        Real sum_added = 0.;
        for (int i = 0; i < n1d; i++) sum_added += radial_mass_summed[i];
        std::cout << "Gravity::make_radial_gravity: Sum of combined mass " << sum_added << std::endl;
    }

    const Geometry& geom = parent->Geom(level);
    const Real* dx = geom.CellSize();
    Real dr        = dx[0] / double(drdxfac);

    // ***************************************************************** //
    // Compute the average density to use at the radius above
    //   max_radius_all_in_domain so we effectively count mass outside
    //   the domain.
    // ***************************************************************** //

    RealVector radial_vol_summed(n1d,0);
    RealVector radial_den_summed(n1d,0);

    // First add the contribution from this level
    for (int i = 0; i < n1d; i++)
         radial_vol_summed[i] =  radial_vol[level][i];

    // Now add the contribution from coarser levels
    if (level > 0)
    {
        int ratio = parent->refRatio(level-1)[0];
        for (int lev = level-1; lev >= 0; lev--)
        {
            if (lev < level-1) ratio *= parent->refRatio(lev)[0];
            for (int i = 0; i < n1d/ratio; i++)
            {
                for (int n = 0; n < ratio; n++)
                {
                   radial_vol_summed[ratio*i+n]  += 1./double(ratio) * radial_vol[lev][i];
                }
            }
        }
    }

    for (int i = 0; i < n1d; i++)
    {
        radial_den_summed[i] = radial_mass_summed[i];
        if (radial_vol_summed[i] > 0.) radial_den_summed[i]  /= radial_vol_summed[i];
    }

#ifdef GR_GRAV
    RealVector radial_pres_summed(n1d,0);

    // First add the contribution from this level
    for (int i = 0; i < n1d; i++)
        radial_pres_summed[i] = radial_pres[level][i];

    // Now add the contribution from coarser levels
    if (level > 0)
    {
        int ratio = parent->refRatio(level-1)[0];
        for (int lev = level-1; lev >= 0; lev--)
        {
            if (lev < level-1) ratio *= parent->refRatio(lev)[0];
            for (int i = 0; i < n1d/ratio; i++)
                for (int n = 0; n < ratio; n++)
                   radial_pres_summed[ratio*i+n] += 1./double(ratio) * radial_pres[lev][i];
        }
    }

    for (int i = 0; i < n1d; i++)
        if (radial_vol_summed[i] > 0.) radial_pres_summed[i] /= radial_vol_summed[i];

    // Integrate radially outward to define the gravity -- here we add the post-Newtonian correction
    ca_integrate_gr_grav(radial_den_summed.dataPtr(),radial_mass_summed.dataPtr(),
                         radial_pres_summed.dataPtr(),radial_grav.dataPtr(),&dr,&n1d);

#else
    // Integrate radially outward to define the gravity
    ca_integrate_grav(radial_mass_summed.dataPtr(),radial_den_summed.dataPtr(),
                      radial_grav.dataPtr(),&max_radius_all_in_domain,&dr,&n1d);
#endif

    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = ParallelDescriptor::second() - strt;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(end,IOProc);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "Gravity::make_radial_gravity() time = " << end << std::endl << std::endl;
#ifdef BL_LAZY
        });
#endif
    }
}

Vector<std::unique_ptr<MultiFab> >
Gravity::get_rhs (int crse_level, int nlevs, int is_new)
{
    Vector<std::unique_ptr<MultiFab> > rhs(nlevs);

    for (int ilev = 0; ilev < nlevs; ++ilev)
    {
        int amr_lev = ilev + crse_level;
        rhs[ilev].reset(new MultiFab(grids[amr_lev],dmap[amr_lev],1,0));
        MultiFab& state = (is_new == 1) ?
            LevelData[amr_lev]->get_new_data(State_Type) :
            LevelData[amr_lev]->get_old_data(State_Type);
        MultiFab::Copy(*rhs[ilev], state, URHO, 0,1,0);
    }
    return rhs;
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
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (!geom.isPeriodic(dir))
            {
                if (phys_bc->lo(dir) != Symmetry)
                    shrunk_domain.growLo(dir,-1);
                if (phys_bc->hi(dir) != Symmetry)
                    shrunk_domain.growHi(dir,-1);
            }
        }
        if (!shrunk_domain.contains(grids[level].minimalBox()))
            amrex::Error("Oops -- don't know how to set boundary conditions for grids at this level that touch the domain boundary!");
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

    const auto& rhs = get_rhs(crse_level, nlevs, is_new);

    const Geometry& geom0 = parent->Geom(0);

#if (BL_SPACEDIM == 3)
    if ( geom0.isAllPeriodic() )
    {
        for (int lev = 0; lev < nlevs; ++lev)
            rhs[lev]->plus(-mass_offset,0,1,0);
    }
#endif

    for (int lev = 0; lev < nlevs; ++lev)
    {
        rhs[lev]->mult(Ggravity);
    }

#if (BL_SPACEDIM < 3)
    if (geom0.IsSPHERICAL() || geom0.IsRZ() )
    {
        Vector<Vector<std::unique_ptr<MultiFab> > > coeffs(nlevs);

        for (int lev = 0; lev < nlevs; ++lev) {

            // We need to include this bit about the coefficients because
            // it's required by applyMetricTerms.

            coeffs[lev].resize(BL_SPACEDIM);

            for (int i = 0; i < BL_SPACEDIM ; i++) {
                coeffs[lev][i].reset(new MultiFab(amrex::convert(grids[lev],
                                                  IntVect::TheDimensionVector(i)),
                                                  dmap[lev], 1, 0));

                coeffs[lev][i]->setVal(1.0);
            }

            applyMetricTerms(lev, *rhs[lev], amrex::GetVecOfPtrs(coeffs[lev]));
        }
    }
#endif

    max_rhs = 0.0;

    for (int lev = 0; lev < nlevs; ++lev)
        max_rhs = std::max(max_rhs, rhs[lev]->max(0));

}

Real
Gravity::solve_phi_with_mlmg (int crse_level, int fine_level,
                              const Vector<MultiFab*>& phi,
                              const Vector<MultiFab*>& rhs,
                              const Vector<Vector<MultiFab*> >& grad_phi,
                              const Vector<MultiFab*>& res,
                              Real time)
{
    BL_PROFILE("Gravity::solve_phi_with_mlmg()");

    int nlevs = fine_level-crse_level+1;

    if (crse_level == 0 && !(parent->Geom(0).isAllPeriodic()))
    {
        if (verbose > 1) {
            amrex::Print() << " ... Making bc's for phi at level 0\n";
        }

#if (BL_SPACEDIM == 3)
        if ( direct_sum_bcs ) {
            fill_direct_sum_BCs(crse_level, fine_level, rhs, *phi[0]);
        } else {
            if (lnum >= 0) {
                fill_multipole_BCs(crse_level, fine_level, rhs, *phi[0]);
            } else {
                int fill_interior = 0;
                make_radial_phi(crse_level, *rhs[0], *phi[0], fill_interior);
            }
        }
#elif (BL_SPACEDIM == 2)
        if (lnum >= 0) {
            fill_multipole_BCs(crse_level, fine_level, rhs, *phi[0]);
        } else {
            int fill_interior = 0;
            make_radial_phi(crse_level, *rhs[0], *phi[0], fill_interior);
        }
#else
        int fill_interior = 0;
        make_radial_phi(crse_level, *rhs[0], *phi[0], fill_interior);
#endif
    }

    for (int ilev = 0; ilev < nlevs; ++ilev)
    {
        rhs[ilev]->mult(Ggravity);
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

    Vector<const MultiFab*> crhs{rhs.begin(), rhs.end()};
    Vector<std::array<MultiFab*,AMREX_SPACEDIM> > gp;
    for (const auto& x : grad_phi) {
        gp.push_back({AMREX_D_DECL(x[0],x[1],x[2])});
    }

    return actual_solve_with_mlmg(crse_level, fine_level, phi, crhs, gp, res,
                                  crse_bcdata, rel_eps, abs_eps);
}

void
Gravity::solve_for_delta_phi(int crse_level, int fine_level,
                             const Vector<MultiFab*>& rhs,
                             const Vector<MultiFab*>& delta_phi,
                             const Vector<Vector<MultiFab*> >&  grad_delta_phi)
{
    BL_PROFILE("Gravity::solve_for_delta_phi");

    BL_ASSERT(grad_delta_phi.size() == fine_level - crse_level + 1);
    BL_ASSERT(delta_phi.size() == fine_level - crse_level + 1);

    if (verbose > 1 && ParallelDescriptor::IOProcessor()) {
      std::cout << "... solving for delta_phi at crse_level = " << crse_level << std::endl;
      std::cout << "...                    up to fine_level = " << fine_level << std::endl;
    }

    Vector<const MultiFab*> crhs{rhs.begin(), rhs.end()};
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
                                 const amrex::Vector<const amrex::MultiFab*>& rhs,
                                 const amrex::Vector<std::array<amrex::MultiFab*,AMREX_SPACEDIM> >& grad_phi,
                                 const amrex::Vector<amrex::MultiFab*>& res,
                                 const amrex::MultiFab* const crse_bcdata,
                                 amrex::Real rel_eps, amrex::Real abs_eps)
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
        bav.push_back(rhs[ilev]->boxArray());
        dmv.push_back(rhs[ilev]->DistributionMap());
    }

    LPInfo info;
    info.setAgglomeration(mlmg_agglomeration);
    info.setConsolidation(mlmg_consolidation);

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
    mlmg.setVerbose(verbose - 1); // With normal verbosity we don't want MLMG information
    if (crse_level == 0) {
        mlmg.setMaxFmgIter(mlmg_max_fmg_iter);
    } else {
        mlmg.setMaxFmgIter(0); // Vcycle
    }

    AMREX_ALWAYS_ASSERT( !grad_phi.empty() or !res.empty() );
    AMREX_ALWAYS_ASSERT(  grad_phi.empty() or  res.empty() );

    if (!grad_phi.empty())
    {
        if (!gmv[0].isAllPeriodic()) mlmg.setAlwaysUseBNorm(true);

        mlmg.setNSolve(mlmg_nsolve);
        final_resnorm = mlmg.solve(phi, rhs, rel_eps, abs_eps);

        mlmg.getGradSolution(grad_phi);
    }
    else if (!res.empty())
    {
        mlmg.compResidual(res, phi, rhs);
    }

    return final_resnorm;
}
