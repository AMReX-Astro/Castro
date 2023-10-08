#include <AMReX_REAL.H>

#include <Derive.H>
#include <Castro.H>
#include <model_parser.H>

using namespace amrex;


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

