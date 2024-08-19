#include <AMReX_REAL.H>

#include <Derive.H>
#include <Castro.H>
#include <Radiation.H>
#include <prob_parameters.H>
#include <extern_parameters.H>
#include <global.H>
#include <gravity_params.H>
#include <problem_util.H>

using namespace amrex;


void deranalytic(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                 const FArrayBox& datfab, const Geometry& geomdata,
                 Real time, const int* /*bcrec*/, int /*level*/)
{

  const auto dx = geomdata.CellSizeArray();
  const auto problo = geomdata.ProbLoArray();

  auto const dat = datfab.array();
  auto const der = derfab.array();

  GpuArray<Real, NGROUPS> nugroup = {0.0};
  GpuArray<Real, NGROUPS> dnugroup = {0.0};

  for (int i = 0; i < NGROUPS; ++i) {
      nugroup[i] = global::the_radiation_ptr->nugroup[i];
      dnugroup[i] = global::the_radiation_ptr->dnugroup[i];
  }

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
  {

    Real loc[3] = {0.0};

    loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0];
#if AMREX_SPACEDIM >= 2
    loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1];
#endif
#if AMREX_SPACEDIM == 3
    loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - problem::center[2];
#endif

    Real r = std::sqrt(loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2]);

    // compute the analytic radiating sphere solution for each point
    for (int g = 0; g < NGROUPS; g++) {
      Real F = F_radsphere(r, time, nugroup[g]);
      Real E = planck(nugroup[g], problem::T_0) +
        (R_sphere / r) * (planck(nugroup[g], T_sphere) - planck(nugroup[g], problem::T_0)) * F;
      der(i,j,k,g) = E * dnugroup[g];
    }

  });

}
