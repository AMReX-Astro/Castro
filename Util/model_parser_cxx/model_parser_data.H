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

    extern amrex::Gpu::ManagedVector<amrex::Real> state;
    extern amrex::Gpu::ManagedVector<amrex::Real> r;
}
#endif
