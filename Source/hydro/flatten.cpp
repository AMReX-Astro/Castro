#include <Castro.H>

#include <cmath>

#ifdef RADIATION
#include <Radiation.H>
#endif

using namespace amrex;

void
Castro::uflatten(const Box& bx,
                 Array4<Real const> const& q_arr,
                 Array4<Real> const& flatn, const int pres_comp) {

  constexpr Real small_pres_flatn = 1.e-200_rt;

  // Knobs for detection of strong shock
  constexpr Real shktst = 0.33_rt;
  constexpr Real zcut1 = 0.75_rt;
  constexpr Real zcut2 = 0.85_rt;
  constexpr Real dzcut = 1.0_rt / (zcut2-zcut1);

  amrex::ParallelFor(bx,
  [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
  {

    // x-direction flattening coef

    Real dp = q_arr(i+1,j,k,pres_comp) - q_arr(i-1,j,k,pres_comp);

    int ishft = dp > 0.0_rt ? 1 : -1;

    Real denom = std::max(small_pres_flatn, std::abs(q_arr(i+2,j,k,pres_comp) - q_arr(i-2,j,k,pres_comp)));
    Real zeta = std::abs(dp) / denom;
    Real z = std::clamp(dzcut * (zeta - zcut1), 0.0_rt, 1.0_rt);

    Real tst = 0.0_rt;
    if (q_arr(i-1,j,k,QU) - q_arr(i+1,j,k,QU) >= 0.0_rt) {
      tst = 1.0_rt;
    }

    Real tmp = std::min(q_arr(i+1,j,k,pres_comp), q_arr(i-1,j,k,pres_comp));

    Real chi = 0.0_rt;
    if (std::abs(dp) > shktst*tmp) {
      chi = tst;
    }


    dp = q_arr(i+1-ishft,j,k,pres_comp) - q_arr(i-1-ishft,j,k,pres_comp);

    denom = std::max(small_pres_flatn, std::abs(q_arr(i+2-ishft,j,k,pres_comp)-q_arr(i-2-ishft,j,k,pres_comp)));
    zeta = std::abs(dp) / denom;
    Real z2 = std::clamp(dzcut * (zeta - zcut1), 0.0_rt, 1.0_rt);

    tst = 0.0_rt;
    if (q_arr(i-1-ishft,j,k,QU) - q_arr(i+1-ishft,j,k,QU) >= 0.0_rt) {
      tst = 1.0_rt;
    }

    tmp = std::min(q_arr(i+1-ishft,j,k,pres_comp), q_arr(i-1-ishft,j,k,pres_comp));

    Real chi2 = 0.0_rt;
    if (std::abs(dp) > shktst*tmp) {
      chi2 = tst;
    }

    flatn(i,j,k) = 1.0_rt - std::max(chi2 * z2, chi * z);


#if AMREX_SPACEDIM >= 2
    // y-direction flattening coef

    dp = q_arr(i,j+1,k,pres_comp) - q_arr(i,j-1,k,pres_comp);

    ishft = dp > 0.0_rt ? 1 : -1;

    denom = std::max(small_pres_flatn, std::abs(q_arr(i,j+2,k,pres_comp) - q_arr(i,j-2,k,pres_comp)));
    zeta = std::abs(dp) / denom;
    z = std::clamp(dzcut * (zeta - zcut1), 0.0_rt, 1.0_rt);

    tst = 0.0_rt;
    if (q_arr(i,j-1,k,QV) - q_arr(i,j+1,k,QV) >= 0.0_rt) {
      tst = 1.0_rt;
    }

    tmp = std::min(q_arr(i,j+1,k,pres_comp), q_arr(i,j-1,k,pres_comp));

    chi = 0.0_rt;
    if (std::abs(dp) > shktst*tmp) {
      chi = tst;
    }


    dp = q_arr(i,j+1-ishft,k,pres_comp) - q_arr(i,j-1-ishft,k,pres_comp);

    denom = std::max(small_pres_flatn, std::abs(q_arr(i,j+2-ishft,k,pres_comp) - q_arr(i,j-2-ishft,k,pres_comp)));
    zeta = std::abs(dp) / denom;
    z2 = std::clamp(dzcut * (zeta - zcut1), 0.0_rt, 1.0_rt);

    tst = 0.0_rt;
    if (q_arr(i,j-1-ishft,k,QV) - q_arr(i,j+1-ishft,k,QV) >= 0.0_rt) {
      tst = 1.0_rt;
    }

    tmp = std::min(q_arr(i,j+1-ishft,k,pres_comp), q_arr(i,j-1-ishft,k,pres_comp));

    chi2 = 0.0_rt;
    if (std::abs(dp) > shktst*tmp) {
      chi2 = tst;
    }

    flatn(i,j,k) = std::min(flatn(i,j,k), 1.0_rt - std::max(chi2 * z2, chi * z));
#endif


#if AMREX_SPACEDIM == 3
    // z-direction flattening coef

    dp = q_arr(i,j,k+1,pres_comp) - q_arr(i,j,k-1,pres_comp);

    ishft = dp > 0.0_rt ? 1: -1;

    denom = std::max(small_pres_flatn, std::abs(q_arr(i,j,k+2,pres_comp) - q_arr(i,j,k-2,pres_comp)));
    zeta = std::abs(dp) / denom;
    z = std::clamp(dzcut * (zeta - zcut1), 0.0_rt, 1.0_rt);

    tst = 0.0_rt;
    if (q_arr(i,j,k-1,QW) - q_arr(i,j,k+1,QW) >= 0.0_rt) {
      tst = 1.0_rt;
    }

    tmp = std::min(q_arr(i,j,k+1,pres_comp), q_arr(i,j,k-1,pres_comp));

    chi = 0.0_rt;
    if (std::abs(dp) > shktst*tmp) {
      chi = tst;
    }


    dp = q_arr(i,j,k+1-ishft,pres_comp) - q_arr(i,j,k-1-ishft,pres_comp);

    denom = std::max(small_pres_flatn, std::abs(q_arr(i,j,k+2-ishft,pres_comp) - q_arr(i,j,k-2-ishft,pres_comp)));
    zeta = std::abs(dp) / denom;
    z2 = std::clamp(dzcut * (zeta - zcut1), 0.0_rt, 1.0_rt);

    tst = 0.0_rt;
    if (q_arr(i,j,k-1-ishft,QW) - q_arr(i,j,k+1-ishft,QW) >= 0.0_rt) {
      tst = 1.0_rt;
    }

    tmp = std::min(q_arr(i,j,k+1-ishft,pres_comp), q_arr(i,j,k-1-ishft,pres_comp));

    chi2 = 0.0_rt;
    if (std::abs(dp) > shktst*tmp) {
      chi2 = tst;
    }

    flatn(i,j,k) = std::min(flatn(i,j,k),  1.0_rt - std::max(chi2 * z2, chi * z));
#endif

  });

}
