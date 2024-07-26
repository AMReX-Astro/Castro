#ifndef RUNTIME_PARAMETERS_H
#define RUNTIME_PARAMETERS_H

#include <AMReX_Print.H>
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
        amrex::ParmParse pp("castro");
#include <castro_queries.H>
    }

#ifdef AMREX_PARTICLES
    {
        amrex::ParmParse pp("particles");
#include <particles_queries.H>
    }
#endif

#ifdef DIFFUSION
    {
        amrex::ParmParse pp("diffusion");
#include <diffusion_queries.H>
    }
#endif

#ifdef GRAVITY
    {
        amrex::ParmParse pp("gravity");
#include <gravity_queries.H>
    }
#endif

#ifdef RADIATION
    {
        amrex::ParmParse pp("radiation");
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
        if (amrex::ParmParse::hasUnusedInputs(nm)) {
            amrex::Print() << "Warning: the following " + nm + ".* parameters are ignored\n";
            auto unused = amrex::ParmParse::getUnusedInputs(nm);
            for (const auto& p: unused) {
                amrex::Print() << p << "\n";
            }
            amrex::Print() << std::endl;
            if (castro::abort_on_invalid_params) {
                amrex::Error("Error: invalid parameters");
            }
        }

    }

}

#endif
