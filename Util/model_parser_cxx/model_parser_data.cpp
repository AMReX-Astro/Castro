#include <AMReX_Vector.H>
#include "nse.H"

namespace model
{

    extern AMREX_GPU_MANAGED int npts;
    extern AMREX_GPU_MANAGED bool initialized;

    extern amrex::Gpu::ManagedVector<Real> state;
    extern amrex::Gpu::ManagedVector<Real> r;

}
