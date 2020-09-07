#include "AMReX_REAL.H"

#include "Derive.H"
#include "Castro.H"
#include "Castro_F.H"
#include "prob_parameters.H"
#include "extern_parameters.H"
#include "gravity_params.H"

using namespace amrex;


void deranalytic(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                 const FArrayBox& datfab, const Geometry& geomdata,
                 Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

  // derive the dynamic pressure

  const auto dx = geomdata.CellSizeArray();
  const auto problo = geomdata.ProbLoArray();

  GpuArray<Real, 3> center;
  ca_get_center(center.begin());

  auto const dat = datfab.array();
  auto const der = derfab.array();


  // Compute pressure from the EOS

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    Real loc[3];

    loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_Rt) * dx[0] - center[0];
    loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_Rt) * dx[1] - center[1];
    loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_Rt) * dx[2] - center[2];

    Real r = std::sqrt(loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2]);

    // compute the analytic radiating sphere solution for each point
    for (int g = 0; g < NGROUPS; g++) {
      Real F = F_radsphere(r, time, nugroup[g]);
      Real E = planck(nugroup[g], T_0) +
        (R_sphere / r) * (planck(nugroup[g], T_sphere) - planck(nugroup[g], T_0)) * F;
      der(i,j,k,g) = E * dnugroup(g);
    }

  });

}
