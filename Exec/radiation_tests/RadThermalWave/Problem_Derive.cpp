#include "AMReX_REAL.H"

#include "Derive.H"
#include "Problem_Derive_F.H"
#include "Castro.H"
#include "Castro_F.H"
#include "fundamental_constants.H"
#include "prob_parameters.H"

using namespace amrex;


void dertexact(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
               const FArrayBox& datfab, const Geometry& geomdata,
               Real time, const int* /*bcrec*/, int /*level*/)
{

  const auto dx = geomdata.CellSizeArray();
  const auto problo = geomdata.ProbLoArray();

  auto const dat = datfab.array();
  auto const der = derfab.array();

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    Real loc[3];

    loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0];
    loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1];
    loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2];

    Real r2 = loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2];

    Real Q = Eexp / rhocv;

    Real a = (16.0_rt * C::sigma_SB) / (3.0_rt * const_kappa_r) / rhocv;
    Real pe = kappa_r_exp_n + 3.0_rt;

    Real pfac = std::exp(std::lgamma(2.5_rt + 1.0_rt / pe) -
                         std::lgamma(1.0_rt + 1.0_rt / pe) -
                         std::lgamma(1.5_rt));

    Real num = 3.0_rt * pe + 2.0_rt;
    Real fac1 = std::pow(2.0_rt, pe - 1.0_rt) * pe * std::pow(M_PI, pe);
    Real ex1 = 1.0_rt / (3.0_rt * pe + 2.0_rt);
    Real fac2 = std::pow(pfac, pe / (3.0_rt * pe + 2.0_rt));

    Real xi0 = std::pow(num / fac1, ex1) * fac2;

    fac1 = a * std::pow(Q, pe) * max(time, 1.e-50_rt);
    ex1 = 1.0_rt / (3.0_rt * pe + 2.0_rt);

    Real xf = xi0 * std::pow(fac1, ex1);

    Real Tbar = Q / std::pow(xf, 3);

    Real Tc = std::pow(xi0, 3) * std::pow(pe * xi0 * xi0 / (6.0_rt * pe + 4.0_rt), 1.0_rt / pe) * Tbar;

    der(i,j,k,0) = Tc * std::pow(amrex::max((1.e0_rt - r2 / (xf * xf)), 0.0_rt), 1.0_rt / pe);

  });
}

void derterror(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
               const FArrayBox& datfab, const Geometry& geomdata,
               Real time, const int* /*bcrec*/, int /*level*/)

{

  const auto dx = geomdata.CellSizeArray();
  const auto problo = geomdata.ProbLoArray();

  GpuArray<Real, 3> center;
  ca_get_center(center.begin());

  auto const dat = datfab.array();
  auto const der = derfab.array();

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
  {

    Real loc[3];

    loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0];
    loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1];
    loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2];

    Real r2 = loc[0] * loc[0] + loc[1] * loc[1] + loc[2] * loc[2];

    Real Q = Eexp / rhocv;

    Real a = (16.0_rt * C::sigma_SB) / (3.0_rt * const_kappa_r) / rhocv;
    Real pe = kappa_r_exp_n + 3.0_rt;

    Real pfac = std::exp(std::lgamma(2.5_rt + 1.0_rt / pe) -
                         std::lgamma(1.0_rt + 1.0_rt / pe) -
                         std::lgamma(1.5_rt));

    Real num = 3.0_rt * pe + 2.0_rt;
    Real fac1 = std::pow(2.0_rt, pe - 1.0_rt) * pe * std::pow(M_PI, pe);
    Real ex1 = 1.0_rt / (3.0_rt * pe + 2.0_rt);
    Real fac2 = std::pow(pfac, pe / (3.0_rt * pe + 2.0_rt));

    Real xi0 = std::pow(num / fac1, ex1) * fac2;

    fac1 = a * std::pow(Q, pe) * max(time, 1.e-50_rt);
    ex1 = 1.0_rt / (3.0_rt * pe + 2.0_rt);

    Real xf = xi0 * std::pow(fac1, ex1);

    Real Tbar = Q / std::pow(xf, 3);

    Real Tc = std::pow(xi0, 3) * std::pow(pe * xi0 * xi0 / (6.0_rt * pe + 4.0_rt), 1.0_rt / pe) * Tbar;

    der(i,j,k,0) = dat(i,j,k,0) - Tc * std::pow(amrex::max((1.e0_rt - r2 / (xf * xf)), 0.0_rt), 1.0_rt / pe);

  });

}
