#include <AMReX_Vector.H>
#include <AMReX_GpuContainers.H>
#include <model_parser.H>

namespace model
{

    AMREX_GPU_MANAGED int npts;
    AMREX_GPU_MANAGED bool initialized;

    AMREX_GPU_MANAGED amrex::Array2D<Real, 0, NPTS_MODEL-1, 0, nvars-1> state;
    AMREX_GPU_MANAGED amrex::Array1D<Real, 0, NPTS_MODEL-1> r;

}
