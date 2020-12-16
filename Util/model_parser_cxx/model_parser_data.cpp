#include <AMReX_Vector.H>
#include <AMReX_GpuContainers.H>
#include <model_parser.H>

namespace model
{

    AMREX_GPU_MANAGED int npts;
    AMREX_GPU_MANAGED bool initialized;

    AMREX_GPU_MANAGED amrex::Array1D<initial_model_t, 0, NUM_MODELS-1> profile;

}
