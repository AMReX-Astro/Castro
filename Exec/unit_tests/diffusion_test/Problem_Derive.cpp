#include "AMReX_REAL.H"

#include "Derive.H"
#include "Castro.H"
#include "Castro_F.H"
#include "prob_parameters.H"

using namespace amrex;


void deranalytic(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                 const FArrayBox& datfab, const Geometry& geomdata,
                 Real time, const int* /*bcrec*/, int /*level*/)
{

  // derive the dynamic pressure

  const auto dx = geomdata.CellSizeArray();
  const auto problo = geomdata.ProbLoArray();

  const int coord_type = geomdata.Coord();

  auto const dat = datfab.array();
  auto const der = derfab.array();

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    Real r[3] = {0.0};
    r[0] = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt);
#if AMREX_SPACEDIM >= 2
    r[1] = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt);
#endif
#if AMREX_SPACEDIM == 3
    r[2] = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt);
#endif

    Real exponent = 0.0;

    if (AMREX_SPACEDIM == 1 && coord_type == 2) {
      // Handle spherical coordinates
       exponent = 3.0_rt / 2.0_rt;

    } else if (AMREX_SPACEDIM == 2 && coord_type == 1) {
      // Handle cylindrical coordinates
      exponent = 3.0_rt / 2.0_rt;

    } else {
      exponent = AMREX_SPACEDIM / 2.0_rt;
    }

    Real dist2 = 0.0;
    for (int d = 0; d < AMREX_SPACEDIM; d++) {
      dist2 += (r[d] - problem::center[d]) * (r[d] - problem::center[d]);
    }

    der(i,j,k,0) = problem::T1 + (problem::T2 - problem::T1) *
        std::pow(problem::t_0 / (time + problem::t_0), exponent) *
        std::exp(-0.25_rt * dist2 / (problem::diff_coeff * (time + problem::t_0)));

  });

}
