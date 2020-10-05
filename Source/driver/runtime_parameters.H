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
    {
        ParmParse pp("radsolve");
#include <radsolve_queries.H>
    }
#endif

}

#endif
