#ifndef RUNTIME_PARAMETERS_H
#define RUNTIME_PARAMETERS_H

#include <AMReX_ParmParse.H>
#include <castro_limits.H>

#ifdef DIFFUSION
#include <diffusion_params.H>
#endif
#ifdef GRAVITY
#include <gravity_params.H>
#endif
#ifdef RADIATION
#include <radsolve_params.H>
#include <radiation_params.H>
#endif

///
/// Initialize all of the runtime parameters defined in _cpp_parameters
/// regardless of the namespace
///
AMREX_INLINE
void
initialize_cpp_runparams()
{

    {
        ParmParse pp("castro");
#include <castro_queries.H>
    }

#ifdef AMREX_PARTICLES
    {
        ParmParse pp("particles");
#include <particles_queries.H>
    }
#endif

#ifdef DIFFUSION
    {
        ParmParse pp("diffusion");
#include <diffusion_queries.H>
    }
#endif

#ifdef GRAVITY
    {
        ParmParse pp("gravity");
#include <gravity_queries.H>
    }
#endif

#ifdef RADIATION
    {
        ParmParse pp("radiation");
#include <radiation_queries.H>
    }
#endif

}


///
/// Check to ensure that any runtime parameters set in inputs or on
/// the command line are valid -- everything should have been accessed
/// by now.
///
AMREX_INLINE
void
validate_runparams()
{

    amrex::Vector<std::string> check_namespaces = {"castro", "problem"};
#ifdef AMREX_PARTICLES
    check_namespaces.push_back("particles");
#endif
#ifdef DIFFUSION
    check_namespaces.push_back("diffusion");
#endif
#ifdef GRAVITY
    check_namespaces.push_back("gravity");
#endif
#ifdef RADIATION
    //check_namespaces.push_back("radiation");  // this doesn't yet work because RadBndry isn't setup yet
    check_namespaces.push_back("radsolve");
#endif

    // now Microphysics

    check_namespaces.push_back("eos");
#ifdef REACTIONS
    check_namespaces.push_back("network");
    check_namespaces.push_back("integrator");
#endif
#ifdef RADIATION
    check_namespaces.push_back("opacity");
#endif
#ifdef DIFFUSION
    check_namespaces.push_back("conductivity");
#endif

    for (const auto& nm: check_namespaces)
    {
        // "castro"
        if (ParmParse::hasUnusedInputs(nm)) {
            amrex::Print() << "Warning: the following " + nm + ".* parameters are ignored\n";
            auto unused = ParmParse::getUnusedInputs(nm);
            for (const auto& p: unused) {
                amrex::Print() << p << "\n";
            }
            amrex::Print() << std::endl;
            if (abort_on_invalid_params) {
                amrex::Error("Error: invalid parameters");
            }
        }

    }

}

///
/// Parse the AMReX runtime parameters to compute dx on the finest level.
/// This is needed for problem_initialize.H, since we have not yet
/// created the level hierarchy at that point.
///
AMREX_INLINE
Vector<Real>
get_fine_grid_dx() {

    Vector<Real> prob_lo(AMREX_SPACEDIM);
    Vector<Real> prob_hi(AMREX_SPACEDIM);

    {
        ParmParse pp("geometry");

        pp.getarr("prob_lo", prob_lo, 0, AMREX_SPACEDIM);
        pp.getarr("prob_hi", prob_hi, 0, AMREX_SPACEDIM);
    }

    int max_level{-1};
    Vector<int> ratios;
    Vector<int> n_cell(AMREX_SPACEDIM);

    {
        ParmParse pp("amr");

        pp.get("max_level", max_level);

        const int got_int = pp.queryarr("ref_ratio", ratios);

        if (! got_int) {
            amrex::Error("cannot read ref_ratio");
        }

        pp.getarr("n_cell", n_cell, 0, AMREX_SPACEDIM);
    }

    if (ratios.size() < max_level) {
        amrex::Error("not enough ref_ratio defined");
    }

    // now compute the dx
    Vector<Real> dx(AMREX_SPACEDIM);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        dx[idim] = (prob_hi[idim] - prob_lo[idim]) /
            static_cast<Real>(n_cell[idim]);
    }

    for (int ilev = 0; ilev < max_level; ++ilev) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            dx[idim] /= static_cast<Real>(ratios[ilev]);
        }
    }

    return dx;
}
#endif
