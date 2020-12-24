#include <AMReX.H>
#include <model_parser_data.H>

namespace model
{

    AMREX_GPU_MANAGED int npts;
    AMREX_GPU_MANAGED bool initialized;

    AMREX_GPU_MANAGED amrex::Array1D<initial_model_t, 0, NUM_MODELS-1> profile;

}
