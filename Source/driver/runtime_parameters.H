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
/// Initialize all of the runtime parameters defined in _cpp_parameters
/// regardless of the namespace
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


    // we don't validate radiation parameters yet, since a lot of the reads
    // are hardcoded into the radiation constructor

    for (auto nm: check_namespaces)
    {
        // "castro"
        if (ParmParse::hasUnusedInputs(nm)) {
            amrex::Print() << "Warning: the following " + nm + ".* parameters are ignored\n";
            auto unused = ParmParse::getUnusedInputs(nm); 
            for (auto p: unused) {
                amrex::Print() << p << "\n";
            }
            amrex::Print() << std::endl;
            if (abort_on_invalid_params) {
                amrex::Error("Error: invalid parameters");
            }
        }

    }

}

#endif
