#include <AMReX_REAL.H>

#include <Derive.H>
#include <Castro.H>
#include <model_parser.H>

using namespace amrex;

using RealVector = amrex::Gpu::ManagedVector<amrex::Real>;

void ca_derpi(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
              const FArrayBox& datfab, const Geometry& geomdata,
              Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

    // derive the dynamic pressure

    const auto dx = geomdata.CellSizeArray();
    const auto problo = geomdata.ProbLoArray();

    auto const dat = datfab.array();
    auto const der = derfab.array();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

        Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt) - problem::center[0];

        Real y = 0.0;
#if AMREX_SPACEDIM >= 2
        y = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt) - problem::center[1];
#endif

        Real z = 0.0;
#if AMREX_SPACEDIM == 3
        z = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt) - problem::center[2];
#endif

#if AMREX_SPACEDIM == 1
        Real dist = x;
#elif AMREX_SPACEDIM == 2
        Real dist = std::sqrt(x*x + y*y);
#else
        Real dist = std::sqrt(x*x + y*y + z*z);
#endif

        Real rhoInv = 1.0_rt / dat(i,j,k,URHO);
        Real e = dat(i,j,k,UEINT) * rhoInv;
        Real temp = dat(i,j,k,UTEMP);

        eos_t eos_state;

        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = dat(i,j,k,UFS+n) / dat(i,j,k,URHO);
        }
#if NAUX_NET > 0
        for (int n = 0; n < NumAux; n++) {
            eos_state.aux[n] = dat(i,j,k,UFX+n) / dat(i,j,k,URHO);
        }
#endif

        // Protect against negative internal energy

        if (e <= 0.0_rt) {
            eos_state.rho = dat(i,j,k,URHO);
            eos_state.T = temp;

            eos(eos_input_rt, eos_state);

            der(i,j,k,0) = eos_state.p;

        } else {

            eos_state.rho = dat(i,j,k,URHO);
            eos_state.T = temp;
            eos_state.e = e;

            eos(eos_input_re, eos_state);

            der(i,j,k,0) = eos_state.p;
        }

        der(i,j,k,0) -= interpolate(dist, model::ipres);

    });
}


void ca_derpioverp0(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                    const FArrayBox& datfab, const Geometry& geomdata,
                    Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

    const auto dx = geomdata.CellSizeArray();
    const auto problo = geomdata.ProbLoArray();

    auto const dat = datfab.array();
    auto const der = derfab.array();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

        Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt) - problem::center[0];

        Real y = 0.0;
#if AMREX_SPACEDIM >= 2
        y = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt) - problem::center[1];
#endif

        Real z = 0.0;
#if AMREX_SPACEDIM == 3
        z = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt) - problem::center[2];
#endif

#if AMREX_SPACEDIM == 1
        Real dist = x;
#elif AMREX_SPACEDIM == 2
        Real dist = std::sqrt(x*x + y*y);
#else
        Real dist = std::sqrt(x*x + y*y + z*z);
#endif

        Real rhoInv = 1.0_rt / dat(i,j,k,URHO);
        Real e = dat(i,j,k,UEINT) * rhoInv;
        Real temp = dat(i,j,k,UTEMP);

        eos_t eos_state;

        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = dat(i,j,k,UFS+n) / dat(i,j,k,URHO);
        }
#if NAUX_NET > 0
        for (int n = 0; n < NumAux; n++) {
            eos_state.aux[n] = dat(i,j,k,UFX+n) / dat(i,j,k,URHO);
        }
#endif

        // Protect against negative internal energy

        if (e <= 0.0_rt) {
            eos_state.rho = dat(i,j,k,URHO);
            eos_state.T = temp;

            eos(eos_input_rt, eos_state);

            der(i,j,k,0) = eos_state.p;

        } else {

            eos_state.rho = dat(i,j,k,URHO);
            eos_state.T = temp;
            eos_state.e = e;

            eos(eos_input_re, eos_state);

            der(i,j,k,0) = eos_state.p;
        }

        Real p0 = interpolate(dist, model::ipres);
        der(i,j,k,0) = std::abs(der(i,j,k,0) - p0) / p0;

    });

}
