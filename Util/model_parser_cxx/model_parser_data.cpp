#include <AMReX_Vector.H>
#include <AMReX_GpuContainers.H>
#include <model_parser.H>

namespace model
{

    extern AMREX_GPU_MANAGED int npts;
    extern AMREX_GPU_MANAGED bool initialized;

    extern amrex::Gpu::ManagedVector<amrex::Real> state;
    extern amrex::Gpu::ManagedVector<amrex::Real> r;

}
