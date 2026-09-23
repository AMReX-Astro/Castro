#ifndef MODEL_PARSER_DATA_H
#define MODEL_PARSER_DATA_H

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuContainers.H>

#include <network_properties.H>

using namespace amrex;

namespace model
{

    // integer keys for indexing the model_state array

    constexpr int nvars = 4 + NumSpec + NumAux;
    constexpr int idens = 0;
    constexpr int itemp = 1;
    constexpr int ipres = 2;
    constexpr int ivelr = 3;
    constexpr int ispec = 4;
#if NAUX_NET > 0
    constexpr int iaux = 4 + NumSpec;
#else
    constexpr int iaux = -1;
#endif

#if DIM_MODEL == 1
    extern AMREX_GPU_MANAGED int npts;
#elif DIM_MODEL == 2
    extern AMREX_GPU_MANAGED int npts_x;
    extern AMREX_GPU_MANAGED int npts_y;
#endif
    extern AMREX_GPU_MANAGED bool initialized;

    struct initial_model_t {
#if DIM_MODEL == 1
        amrex::Array2D<amrex::Real, 0, NPTS_MODEL-1, 0, nvars-1> state;
        amrex::Array1D<amrex::Real, 0, NPTS_MODEL-1> r;
#elif DIM_MODEL == 2
        amrex::Array3D<amrex::Real, 0, NPTS_X_MODEL-1,
                       0, NPTS_Y_MODEL-1, 0, nvars-1> state;
        amrex::Array1D<amrex::Real, 0, NPTS_X_MODEL-1> x;
        amrex::Array1D<amrex::Real, 0, NPTS_Y_MODEL-1> y;
#endif
    };

    // Tolerance used for getting the total star mass equal to the desired mass.
    // It can be reasonably small, since there will always be a central density
    // value that can give the desired mass on the grid.

    const amrex::Real mass_tol = 1.0e-6_rt;

    // Tolerance used when iterating over a zone to force it into HSE by adjusting
    // the current density (and possibly temperature). Should be very small (~1.e-10).

    const amrex::Real hse_tol = 1.0e-10_rt;

    extern AMREX_GPU_MANAGED amrex::Array1D<initial_model_t, 0, NUM_MODELS-1> profile;
}
#endif
