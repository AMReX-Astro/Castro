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

    AMREX_GPU_MANAGED int npts;
    AMREX_GPU_MANAGED bool initialized;

    amrex::Gpu::ManagedVector<amrex::Real> state;
    amrex::Gpu::ManagedVector<amrex::Real> r;
}
#endif
