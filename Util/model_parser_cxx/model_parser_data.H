#ifndef MODEL_PARSER_DATA_H
#define MODEL_PARSER_DATA_H

namespace model
{

    // integer keys for indexing the model_state array

    constexpr int nvars = 3 + NumSpec;
    constexpr int idens = 0;
    constexpr int itemp = 1;
    constexpr int ipres = 2;
    constexpr int ispec = 3;

    extern AMREX_GPU_MANAGED int npts;
    extern AMREX_GPU_MANAGED bool initialized;

    struct initial_model_t {
        amrex::Array2D<amrex::Real, 0, NPTS_MODEL-1, 0, nvars-1> state;
        amrex::Array1D<amrex::Real, 0, NPTS_MODEL-1> r;
    };

    extern AMREX_GPU_MANAGED amrex::Array1D<initial_model_t, 0, NUM_MODELS-1> profile;
}
#endif
