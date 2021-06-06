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

    {
        // "castro"
        if (ParmParse::hasUnusedInputs("castro")) {
            amrex::Print() << "Warning: the following castro.* parameters are ignored\n";
            auto unused = ParmParse::getUnusedInputs("castro"); 
            for (auto p: unused) {
                amrex::Print() << p << "\n";
            }
            amrex::Print() << std::endl;
            if (abort_on_invalid_params) {
                amrex::Error("Error: invalid parameters");
            }
        }

    }

#ifdef AMREX_PARTICLES
    {
        // "particles"
        if (ParmParse::hasUnusedInputs("particles")) {
            amrex::Print() << "Warning: the following particles.* parameters are ignored\n";
            auto unused = ParmParse::getUnusedInputs("particles"); 
            for (auto p: unused) {
                amrex::Print() << p << "\n";
            }
            amrex::Print() << std::endl;
            if (abort_on_invalid_params) {
                amrex::Error("Error: invalid parameters");
            }
        }

    }
#endif

#ifdef DIFFUSION
    {
        // "diffusion"
        if (ParmParse::hasUnusedInputs("diffusion")) {
            amrex::Print() << "Warning: the following diffusion.* parameters are ignored\n";
            auto unused = ParmParse::getUnusedInputs("diffusion"); 
            for (auto p: unused) {
                amrex::Print() << p << "\n";
            }
            amrex::Print() << std::endl;
            if (abort_on_invalid_params) {
                amrex::Error("Error: invalid parameters");
            }
        }

    }
#endif

#ifdef GRAVITY
    {
        // "gravity"
        if (ParmParse::hasUnusedInputs("gravity")) {
            amrex::Print() << "Warning: the following gravity.* parameters are ignored\n";
            auto unused = ParmParse::getUnusedInputs("gravity"); 
            for (auto p: unused) {
                amrex::Print() << p << "\n";
            }
            amrex::Print() << std::endl;
            if (abort_on_invalid_params) {
                amrex::Error("Error: invalid parameters");
            }
        }

    }
#endif

    // we don't validate radiation parameters yet, since a lot of the reads
    // are hardcoded into the radiation constructor

}

#endif
