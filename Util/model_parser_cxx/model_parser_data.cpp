#include <AMReX_Vector.H>
#include <AMReX_GpuContainers.H>
#include <model_parser.H>

namespace model
{

    AMREX_GPU_MANAGED int npts;
    AMREX_GPU_MANAGED bool initialized;

    amrex::Gpu::ManagedVector<amrex::Real> state;
    amrex::Gpu::ManagedVector<amrex::Real> r;

}
