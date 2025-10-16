#include <AMReX.H>
#include <model_parser_data.H>

namespace model
{

#if DIM_MODEL >= 1
    AMREX_GPU_MANAGED int npts_x;
#endif
#if DIM_MODEL >= 2
    AMREX_GPU_MANAGED int npts_y;
#endif
    AMREX_GPU_MANAGED bool initialized;

    AMREX_GPU_MANAGED amrex::Array1D<initial_model_t, 0, NUM_MODELS-1> profile;

}
