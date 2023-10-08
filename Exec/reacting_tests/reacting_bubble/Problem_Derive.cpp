#include <AMReX_REAL.H>

#include <Derive.H>
#include <Castro.H>
#include <model_parser.H>

using namespace amrex;

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

        Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

        Real y = 0.0;
#if AMREX_SPACEDIM >= 2
        y = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
#endif

        Real z = 0.0;
#if AMREX_SPACEDIM == 3
        z = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
#endif

#if AMREX_SPACEDIM == 2
        Real height = y;
#else
        Real height = z;
#endif

        Real pres = interpolate(height, model::ipres);

        eos_t eos_state;

        eos_state.rho = dat(i,j,k,URHO);
        eos_state.T = dat(i,j,k,UTEMP);
        eos_state.e = dat(i,j,k,UEINT) / dat(i,j,k,URHO);
        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = dat(i,j,k,UFS+n) / dat(i,j,k,URHO);
        }
#if NAUX_NET > 0
        for (int n = 0; n < NumAux; n++) {
            eos_state.aux[n] = dat(i,j,k,UFX+n) / dat(i,j,k,URHO);
        }
#endif

        eos(eos_input_rt, eos_state);
        der(i,j,k,0) = eos_state.p - pres;

    });

}


void ca_derpioverp0(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                    const FArrayBox& datfab, const Geometry& geomdata,
                    Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

    // derive the dynamic pressure / thermodynamic pressure

    const auto dx = geomdata.CellSizeArray();
    const auto problo = geomdata.ProbLoArray();

    auto const dat = datfab.array();
    auto const der = derfab.array();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

        Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

        Real y = 0.0;
#if AMREX_SPACEDIM >= 2
        y = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
#endif

        Real z = 0.0;
#if AMREX_SPACEDIM == 3
        z = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
#endif

#if AMREX_SPACEDIM == 2
        Real height = y;
#else
        Real height = z;
#endif

        Real pres = interpolate(height, model::ipres);

        eos_t eos_state;

        eos_state.rho = dat(i,j,k,URHO);
        eos_state.T = dat(i,j,k,UTEMP);
        eos_state.e = dat(i,j,k,UEINT) / dat(i,j,k,URHO);
        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = dat(i,j,k,UFS+n) / dat(i,j,k,URHO);
        }
#if NAUX_NET > 0
        for (int n = 0; n < NumAux; n++) {
            eos_state.aux[n] = dat(i,j,k,UFX+n) / dat(i,j,k,URHO);
        }
#endif

        eos(eos_input_rt, eos_state);
        der(i,j,k,0) = (eos_state.p - pres) / pres;

    });

}


void ca_derrhopert(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const FArrayBox& datfab, const Geometry& geomdata,
                   Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

    // derive the dynamic pressure / thermodynamic pressure

    const auto dx = geomdata.CellSizeArray();
    const auto problo = geomdata.ProbLoArray();

    auto const dat = datfab.array();
    auto const der = derfab.array();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

        Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

        Real y = 0.0;
#if AMREX_SPACEDIM >= 2
        y = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
#endif

        Real z = 0.0;
#if AMREX_SPACEDIM == 3
        z = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
#endif

#if AMREX_SPACEDIM == 2
        Real height = y;
#else
        Real height = z;
#endif

        Real dens = interpolate(height, model::idens);

        der(i,j,k,0) = dat(i,j,k,URHO) - dens;

    });

}


void ca_dertpert(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                 const FArrayBox& datfab, const Geometry& geomdata,
                 Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

    // derive the temperature perturbation

    const auto dx = geomdata.CellSizeArray();
    const auto problo = geomdata.ProbLoArray();

    auto const dat = datfab.array();
    auto const der = derfab.array();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

        Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);

        Real y = 0.0;
#if AMREX_SPACEDIM >= 2
        y = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
#endif

        Real z = 0.0;
#if AMREX_SPACEDIM == 3
        z = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
#endif

#if AMREX_SPACEDIM == 2
        Real height = y;
#else
        Real height = z;
#endif

        Real temp = interpolate(height, model::itemp);

        der(i,j,k,0) = dat(i,j,k,UTEMP) - temp;

    });

}

