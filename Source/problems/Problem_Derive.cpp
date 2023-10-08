#include <AMReX_REAL.H>

#include <Derive.H>
#include <Castro.H>

using namespace amrex;

#ifdef __cplusplus
extern "C"
{
#endif

    // Note that in the following routines, we are NOT passing
    // several variables to Fortran that would be unused.

    // These routines are called in an MFIter loop, so we do not
    // need to explicitly synchronize after GPU kernels.

#ifdef __cplusplus
}
#endif
